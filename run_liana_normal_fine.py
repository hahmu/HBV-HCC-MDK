import os, glob, shutil, time
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad
import scanpy as sc
import liana as li

BASE = os.environ.get("HBV_BASE", "./HBV_HCC")
LDIR = f"{BASE}/results/liana"
IN_DIR = f"{LDIR}/input"


SETTINGS = dict(
    label_col=os.environ.get("LABEL_COL", "liana_label"),   # subtype_fine + Malignant
    tissue_col=os.environ.get("TISSUE_COL", "tissue"),
    normalise="total_1e4_log1p",       
    expr_prop=0.10,                    
    min_cells=5,                       
    n_perms=1000,
    seed=1337,
)

DROP_MALIGNANT_IN_NORMAL = os.environ.get("DROP_MALIGNANT_IN_NORMAL", "1") == "1"
MAG, SPEC = 0.01, 0.05
EXCLUDE_GENES = {"B2M", "ALB", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3",
                 "APOE", "APOH", "HP", "HPX", "TF", "TTR", "FGA", "FGB", "FGG", "SERPINA1",
                 "AHSG", "ORM1", "ORM2", "RBP4", "AMBP", "CP", "VTN", "FTL", "FTH1"}


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def first(pattern):
    hits = sorted(glob.glob(f"{IN_DIR}/{pattern}"))
    if not hits:
        raise SystemExit(f"No file matching {IN_DIR}/{pattern}")
    return hits[0]


def read_names(path):
    sep = "," if path.endswith(".csv") else "\t"
    col = [ln.rstrip("\n").split(sep)[0].strip('"') for ln in open(path) if ln.strip()]
    if col[0].lower() in ("x", "gene", "genes", "cell", "cells", "barcode", "features"):
        col = col[1:]
    return col


def load_input():
    X = scipy.io.mmread(first("*.mtx")).tocsr()
    genes = read_names(first("gene*"))
    cells = read_names(first("cell*"))
    meta = pd.read_csv(first("meta*.csv"), index_col=0)
    if X.shape == (len(genes), len(cells)):
        X = X.T.tocsr()
    elif X.shape != (len(cells), len(genes)):
        raise SystemExit(f"matrix {X.shape} does not match {len(cells)} cells x {len(genes)} genes")
    meta = meta.reindex(cells)
    for c in (SETTINGS["label_col"], SETTINGS["tissue_col"]):
        if c not in meta.columns:
            raise SystemExit(f"column '{c}' not in meta.csv; columns are {list(meta.columns)} "
                             f"(set LABEL_COL / TISSUE_COL)")
    a = ad.AnnData(X=X.astype(np.float32), obs=meta, var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    log(f"loaded {a.n_obs} cells x {a.n_vars} genes; tissues: {a.obs[SETTINGS['tissue_col']].value_counts().to_dict()}")
    return a


def tissue_subset(a, tissue):
    t = a.obs[SETTINGS["tissue_col"]].astype(str)
    hit = [v for v in t.unique() if v.lower().startswith(tissue.lower()[:3])]   # Tumor/Tumour, Normal
    if len(hit) != 1:
        raise SystemExit(f"cannot match tissue '{tissue}' among {t.unique().tolist()}")
    s = a[t == hit[0]].copy()
    if tissue == "Normal" and DROP_MALIGNANT_IN_NORMAL:
        n = (s.obs[SETTINGS["label_col"]] == "Malignant").sum()
        s = s[s.obs[SETTINGS["label_col"]] != "Malignant"].copy()
        log(f"Normal: dropped {n} adjacent cells labelled Malignant")
    return s


def run(a, tissue):
    s = tissue_subset(a, tissue)
    if SETTINGS["normalise"] == "total_1e4_log1p":
        sc.pp.normalize_total(s, target_sum=1e4)
        sc.pp.log1p(s)
    s.obs[SETTINGS["label_col"]] = s.obs[SETTINGS["label_col"]].astype(str).astype("category")
    li.mt.rank_aggregate(s, groupby=SETTINGS["label_col"], resource_name=SETTINGS["resource_name"],
                         expr_prop=SETTINGS["expr_prop"], min_cells=SETTINGS["min_cells"],
                         n_perms=SETTINGS["n_perms"], seed=SETTINGS["seed"], use_raw=False,
                         verbose=True, key_added="liana_res")
    res = s.uns["liana_res"].copy()
    res["condition"] = tissue
    counts = s.obs[SETTINGS["label_col"]].value_counts().rename(tissue)
    return res, counts


def sig_filtered(d):
    bad = lambda c: any(g in EXCLUDE_GENES for g in str(c).split("_"))
    d = d[~d["ligand_complex"].map(bad) & ~d["receptor_complex"].map(bad)]
    return d[(d["magnitude_rank"] <= MAG) & (d["specificity_rank"] <= SPEC)]


def key(d):
    return set(zip(d["source"], d["target"], d["ligand_complex"], d["receptor_complex"]))


def main():
    a = load_input()
    report = [f"{k} = {v}" for k, v in SETTINGS.items()] + \
             [f"drop_malignant_in_normal = {DROP_MALIGNANT_IN_NORMAL}",
              f"significant = magnitude_rank <= {MAG} & specificity_rank <= {SPEC}, ambient genes removed",
              f"liana version = {li.__version__}"]
    counts = []

    if os.environ.get("RUN_BOTH") == "1":
        tum = f"{LDIR}/liana_Tumor.csv"
        if os.path.exists(tum) and not os.path.exists(f"{LDIR}/liana_Tumor_35labels.csv"):
            shutil.copy(tum, f"{LDIR}/liana_Tumor_35labels.csv")
            log("old liana_Tumor.csv kept as liana_Tumor_35labels.csv")
        new_t, ct = run(a, "Tumor")
        counts.append(ct)
        new_t.to_csv(tum, index=False)
        nt = sorted(set(new_t["source"]) | set(new_t["target"]))
        report += [f"TUMOUR (RUN_BOTH): {len(nt)} labels tested; {len(sig_filtered(new_t))} significant "
                   "non-ambient interactions", "TUMOUR labels tested: " + ", ".join(nt)]
        log(f"Tumour: {len(sig_filtered(new_t))} significant interactions")
    elif os.environ.get("SKIP_CHECK") != "1":
        old_t = pd.read_csv(f"{LDIR}/liana_Tumor.csv")
        new_t, ct = run(a, "Tumor")
        counts.append(ct)
        new_t.to_csv(f"{LDIR}/liana_Tumor_check.csv", index=False)
        so, sn = key(sig_filtered(old_t)), key(sig_filtered(new_t))
        jac = len(so & sn) / max(1, len(so | sn))
        m = old_t.merge(new_t, on=["source", "target", "ligand_complex", "receptor_complex"], suffixes=("_old", "_new"))
        rho = m["magnitude_rank_old"].corr(m["magnitude_rank_new"], method="spearman")
        labels_same = set(old_t["source"]) == set(new_t["source"])
        ok = labels_same and jac >= 0.95 and rho >= 0.99
        msg = (f"TUMOUR CHECK: old {len(so)} vs new {len(sn)} significant; Jaccard {jac:.3f}; "
               f"magnitude_rank Spearman {rho:.4f}; same labels {labels_same} -> {'PASS' if ok else 'FAIL'}")
        log(msg); report.append(msg)
        if not ok:
            report.append("Settings do not reproduce liana_Tumor.csv: Normal NOT run. Edit SETTINGS and rerun.")
            open(f"{LDIR}/liana_run_settings.txt", "w").write("\n".join(report) + "\n")
            raise SystemExit("Tumour check failed - see liana_run_settings.txt. Normal was not run.")

    old_n = f"{LDIR}/liana_Normal.csv"
    if os.path.exists(old_n) and not os.path.exists(f"{LDIR}/liana_Normal_OLDLABELS.csv"):
        shutil.copy(old_n, f"{LDIR}/liana_Normal_OLDLABELS.csv")
        log("old liana_Normal.csv kept as liana_Normal_OLDLABELS.csv")
    new_n, cn = run(a, "Normal")
    counts.append(cn)
    new_n.to_csv(old_n, index=False)
    tested = sorted(set(new_n["source"]) | set(new_n["target"]))
    n_sig = len(sig_filtered(new_n))
    report += [f"NORMAL: {len(tested)} labels tested; {n_sig} significant non-ambient interactions",
               "NORMAL labels tested: " + ", ".join(tested)]

    tab = pd.concat(counts, axis=1).fillna(0).astype(int)
    tab["tested_Normal"] = tab.index.isin(tested)
    tab.to_csv(f"{LDIR}/cells_per_label_by_tissue.csv")
    below = tab.index[(tab.get("Normal", 0) < SETTINGS["min_cells"])].tolist()
    report.append(f"Labels not tested in Normal (< {SETTINGS['min_cells']} non-tumour cells; "
                  f"adjacent Malignant excluded by design): {', '.join(below) or 'none'}")
    open(f"{LDIR}/liana_run_settings.txt", "w").write("\n".join(report) + "\n")
    log(f"Normal: {n_sig} significant interactions. Next: python scripts/fig4_liana_overview.py")


if __name__ == "__main__":
    main()
