import gzip, io, json, re, sys, tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from matplotlib.image import imread

import config as C

KINDS = [
    "filtered_feature_bc_matrix", "raw_feature_bc_matrix", "matrix", "barcodes", "features", "genes",
    "tissue_positions_list", "tissue_positions", "scalefactors_json",
    "tissue_hires_image", "tissue_lowres_image", "aligned_fiducials", "detected_tissue_image",
]
PAT = re.compile(r"^(GSM\d+)_(.+?)_(" + "|".join(sorted(KINDS, key=len, reverse=True)) + r")\.(.+)$")


def _open(p: Path, mode="rb"):
    return gzip.open(p, mode) if p.suffix == ".gz" else open(p, mode)


def inventory(folder: Path) -> dict:
    # unpack any tar archives into their own subfolder first
    for t in list(folder.glob("*.tar")) + list(folder.glob("*.tar.gz")):
        dest = folder / (t.name.split(".tar")[0] + "_untar")
        if not dest.exists():
            with tarfile.open(t) as tf:
                tf.extractall(dest, filter="data")
    samples = defaultdict(dict)
    unmatched = []
    for f in sorted(folder.rglob("*")):
        if not f.is_file():
            continue
        m = PAT.match(f.name)
        if m:
            gsm, name, kind, _ = m.groups()
            samples[(gsm, name)][kind] = f
        else:
            unmatched.append(f.name)
    if unmatched:
        print(f"[inventory] {len(unmatched)} files not recognised, e.g. {unmatched[:5]}")
    return samples


def parse_sample(name):
    m = re.match(C.SAMPLE_PATTERN, name)
    if not m:
        print(f"[parse] '{name}' does not match SAMPLE_PATTERN - patient/region set to the name")
        return name, "Unknown", ""
    patient = f"{m['letter'].upper()}{int(m['number']):03d}"
    return patient, C.TYPE_REGION.get(m["type"].upper(), "Unknown"), m["part"]


def read_counts(files: dict) -> ad.AnnData:
    if "filtered_feature_bc_matrix" in files:
        a = sc.read_10x_h5(files["filtered_feature_bc_matrix"])
    elif "matrix" in files:
        import scipy.io
        with _open(files["matrix"]) as fh:
            X = scipy.io.mmread(fh).T.tocsr()
        bc = pd.read_csv(files["barcodes"], header=None, sep="\t")[0].values
        feat_f = files.get("features", files.get("genes"))
        feat = pd.read_csv(feat_f, header=None, sep="\t")
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=bc),
                       var=pd.DataFrame({"gene_ids": feat[0].values}, index=feat[1].values))
    elif "raw_feature_bc_matrix" in files:
        a = sc.read_10x_h5(files["raw_feature_bc_matrix"])  # filtered to in-tissue below
    else:
        raise FileNotFoundError("no count matrix")
    a.var_names_make_unique()
    return a


def read_positions(files: dict) -> pd.DataFrame:
    f = files.get("tissue_positions", files.get("tissue_positions_list"))
    if f is None:
        raise FileNotFoundError("no tissue_positions file")
    with _open(f, "rt") as fh:
        first = fh.readline()
    header = 0 if first.lower().startswith("barcode") else None
    pos = pd.read_csv(f, header=header)
    pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                   "pxl_row_in_fullres", "pxl_col_in_fullres"]
    return pos.set_index("barcode")


def read_image(f):
    if f is None:
        return None
    if f.suffix == ".gz":
        with gzip.open(f) as fh:
            return imread(io.BytesIO(fh.read()), format="png")
    return imread(f)


def main():
    C.OUT_DIR.mkdir(exist_ok=True, parents=True)
    samples = inventory(C.VISIUM_DIR)
    if not samples:
        sys.exit(f"No GSM*_ files found in {C.VISIUM_DIR.resolve()}")

    rows, adatas = [], []
    for (gsm, name), files in samples.items():
        row = {"gsm": gsm, "sample": name, "files": ",".join(sorted(files))}
        if name in getattr(C, "EXCLUDE_SAMPLES", []):
            row["status"] = "EXCLUDED (config.EXCLUDE_SAMPLES)"
            rows.append(row)
            continue
        try:
            a = read_counts(files)
            pos = read_positions(files)
            with _open(files["scalefactors_json"], "rt") as fh:
                scalef = json.load(fh)
        except Exception as e:  # keep going, report at the end
            row["status"] = f"SKIPPED: {e}"
            rows.append(row)
            continue
        pos = pos[pos["in_tissue"] == 1]
        common = a.obs_names.intersection(pos.index)
        a = a[common].copy()
        a.obs = a.obs.join(pos)
        a.obsm["spatial"] = a.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
        a.obs["sample"] = name
        a.obs["gsm"] = gsm
        patient, region, part = parse_sample(name)
        a.obs["patient"] = patient
        a.obs["region_label"] = region
        a.obs["section_part"] = part
        a.obs_names = [f"{name}_{b}" for b in a.obs_names]
        a.uns["spatial"] = {name: {
            "images": {k: v for k, v in {
                "hires": read_image(files.get("tissue_hires_image")),
                "lowres": read_image(files.get("tissue_lowres_image"))}.items() if v is not None},
            "scalefactors": scalef,
            "metadata": {"gsm": gsm},
        }}
        a.var["SYMBOL"] = a.var_names
        a.var["MT"] = a.var_names.str.upper().str.startswith("MT-")
        row.update(status="ok", spots_in_tissue=a.n_obs,
                   median_umi=float(np.median(np.asarray(a.X.sum(1)).ravel())),
                   median_genes=float(np.median(np.asarray((a.X > 0).sum(1)).ravel())),
                   pct_mt=float(100 * a[:, a.var["MT"]].X.sum() / a.X.sum()),
                   patient=patient, region_label=region, section_part=part)
        rows.append(row)
        adatas.append(a)

    inv = pd.DataFrame(rows)
    inv.to_csv(C.OUT_DIR / "visium_inventory.csv", index=False)
    print(inv.drop(columns="files").to_string(index=False))
    if not adatas:
        sys.exit("No sample could be loaded - see visium_inventory.csv")

    spatial = {k: v for a in adatas for k, v in a.uns["spatial"].items()}
    adata = ad.concat(adatas, join="outer", merge="first", fill_value=0)
    adata.uns["spatial"] = spatial
    adata.var["SYMBOL"] = adata.var_names
    adata.var["MT"] = adata.var_names.str.upper().str.startswith("MT-")
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.X = adata.X.astype(np.float32)
    adata.write_h5ad(C.OUT_DIR / "visium_raw.h5ad", compression="gzip")
    print(f"\nWrote {C.OUT_DIR/'visium_raw.h5ad'}: {adata.n_obs} spots x {adata.n_vars} genes, "
          f"{adata.obs['sample'].nunique()} sections")


if __name__ == "__main__":
    main()
