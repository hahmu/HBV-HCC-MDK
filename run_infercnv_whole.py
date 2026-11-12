import os, sys, time, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- settings ----------------
BASE      = os.environ.get("HBV_BASE", "./HBV_HCC")
IN_DIR    = f"{BASE}/results/infercnv/input"
OUT_DIR   = f"{BASE}/results/infercnv/whole"
FIG_DIR   = f"{OUT_DIR}/figures"
GENE_POS  = f"{BASE}/reference/hg38_gencode_v27.txt"
CACHE     = f"{OUT_DIR}/whole_cnv.h5ad"

QUERY_TYPES = ["Hepatocyte", "Cholangiocyte", "Proliferating"]
REF_TYPES   = ["T_cell", "NK_cell", "B_cell", "Plasma_cell", "Myeloid", "Dendritic_cell",
               "Endothelial", "Fibroblast_HSC", "Mast_cell"]
WINDOW_SIZE     = 250
LEIDEN_RES      = 1.0
MIN_QUERY_FRAC  = 0.80     
MIN_PATIENT_FRAC = 0.80    
MIN_TUMOR_FRAC  = 0.60     
SCORE_FOLD      = 1.5      
N_JOBS = max(1, (os.cpu_count() or 2) - 2)
MASK_REGIONS = [("chr6",  28_500_000,  33_500_000),   
                ("chr14", 105_550_000, 106_900_000),  
                ("chr2",  88_800_000,  90_300_000),   
                ("chr22", 22_000_000,  22_950_000)]   

plt.rcParams["svg.fonttype"] = "none"
sc.settings.verbosity = 0
os.makedirs(FIG_DIR, exist_ok=True)

def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)

def savefig(name, svg=True):
    plt.savefig(f"{FIG_DIR}/{name}.png", dpi=200, bbox_inches="tight")
    if svg:
        plt.savefig(f"{FIG_DIR}/{name}.svg", dpi=200, bbox_inches="tight")
    plt.close("all")

def read_gene_pos(path):
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                f = line.split()
            if len(f) < 4:
                continue
            g, ch, s, e = (x.strip() for x in f[:4])
            if not ch.startswith("chr"):
                ch = "chr" + ch
            if s.isdigit() and e.isdigit():
                rows.append((g, ch, int(s), int(e)))
    gp = pd.DataFrame(rows, columns=["gene", "chromosome", "start", "end"])
    if len(gp) < 10000:
        sys.exit(f"ERROR: only {len(gp)} usable rows in {path}")
    return gp.drop_duplicates("gene").set_index("gene")


def build_and_run():
    log("Loading counts ...")
    X     = scipy.io.mmread(f"{IN_DIR}/counts.mtx").T.tocsr()
    genes = pd.read_csv(f"{IN_DIR}/genes.txt", header=None)[0].values
    cells = pd.read_csv(f"{IN_DIR}/cells.txt", header=None)[0].values
    meta  = pd.read_csv(f"{IN_DIR}/meta.csv", index_col=0).loc[cells]
    a = ad.AnnData(X=X.astype(np.float32), obs=meta, var=pd.DataFrame(index=genes))
    a.var_names_make_unique()

    ct = a.obs["celltype"].astype(str)
    a = a[ct.isin(QUERY_TYPES + REF_TYPES).values].copy()
    log(f"{a.n_obs} cells x {a.n_vars} genes")

    a.var = a.var.join(read_gene_pos(GENE_POS))
    keep = a.var["chromosome"].isin([f"chr{i}" for i in list(range(1, 23)) + ["X"]]).values
    a = a[:, keep].copy()
    a.var["start"] = a.var["start"].astype(int); a.var["end"] = a.var["end"].astype(int)
    masked = np.zeros(a.n_vars, dtype=bool)
    for ch, s, e in MASK_REGIONS:
        masked |= ((a.var["chromosome"] == ch) & (a.var["end"] >= s) & (a.var["start"] <= e)).values
    masked |= np.asarray(a.var_names.str.startswith(("HLA-", "IGH", "IGK", "IGL")))
    a = a[:, ~masked].copy()
    log(f"{a.n_vars} genes with position used ({masked.sum()} HLA/Ig genes masked)")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    is_ref = a.obs["celltype"].astype(str).isin(REF_TYPES).values
    a.obs["cnv_group"] = np.where(is_ref, "Ref_" + a.obs["dataset"].astype(str),
                                  a.obs["celltype"].astype(str))
    ref_cats = sorted(set(a.obs["cnv_group"][is_ref]))
    log(f"Query cells: {(~is_ref).sum()} | reference cells: {is_ref.sum()} in {ref_cats}")

    log("Running infercnv (whole dataset) ...")
    cnv.tl.infercnv(a, reference_key="cnv_group", reference_cat=ref_cats,
                    window_size=WINDOW_SIZE, n_jobs=N_JOBS)
    log("PCA / neighbours / leiden on CNV profiles ...")
    cnv.tl.pca(a)
    cnv.pp.neighbors(a)
    cnv.tl.leiden(a, resolution=LEIDEN_RES)
    cnv.tl.umap(a)
    cnv.tl.cnv_score(a)
    a.write_h5ad(CACHE)
    return a


def main():
    if os.path.exists(CACHE):
        log("Loading cached whole-dataset CNV result")
        a = ad.read_h5ad(CACHE)
    else:
        a = build_and_run()

    obs = a.obs
    is_ref = obs["cnv_group"].astype(str).str.startswith("Ref_").values
    M = a.obsm["X_cnv"]
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    obs["cnv_cell_score"] = np.abs(M).mean(axis=1)
    obs["is_query"] = ~is_ref

    # ---------- cluster table ----------
    rows = []
    for cl, d in obs.groupby("cnv_leiden", observed=True):
        q = d[d["is_query"]]
        top_pt = q["patient"].value_counts()
        rows.append(dict(
            cnv_leiden=cl, n_cells=len(d),
            frac_query=round(d["is_query"].mean(), 3),
            n_query=len(q),
            top_patient=top_pt.index[0] if len(q) else "",
            frac_top_patient=round(top_pt.iloc[0] / len(q), 3) if len(q) else 0,
            frac_tumor=round((q["tissue"] == "Tumor").mean(), 3) if len(q) else 0,
            main_celltype=q["celltype"].value_counts().index[0] if len(q) else
                          d["celltype"].value_counts().index[0],
            cnv_score=round(float(d["cnv_score"].iloc[0]), 5)))
    cl = pd.DataFrame(rows)
    ref_level = cl.loc[cl["frac_query"] <= 0.2, "cnv_score"].max()
    if not np.isfinite(ref_level):
        ref_level = obs.loc[is_ref, "cnv_cell_score"].median()
    thr = SCORE_FOLD * ref_level
    cl["malignant_cluster"] = ((cl["frac_query"] >= MIN_QUERY_FRAC) &
                               (cl["frac_top_patient"] >= MIN_PATIENT_FRAC) &
                               (cl["frac_tumor"] >= MIN_TUMOR_FRAC) &
                               (cl["cnv_score"] >= thr))
    cl = cl.sort_values(["malignant_cluster", "cnv_score"], ascending=False)
    cl.to_csv(f"{OUT_DIR}/cnv_cluster_summary.csv", index=False)
    log(f"Reference-cluster CNV score = {ref_level:.4f}; malignant threshold = {thr:.4f}")
    log("\n" + cl.to_string(index=False))

    
    mal = cl[cl["malignant_cluster"]]
    owner_of = dict(zip(mal["cnv_leiden"].astype(str), mal["top_patient"]))
    owners = sorted(set(owner_of.values()))
    leid = obs["cnv_leiden"].astype(str).values
    pat  = obs["patient"].astype(str).values
    tis  = obs["tissue"].astype(str).values
    q    = ~is_ref
    in_own_mal = np.array([owner_of.get(l) == p for l, p in zip(leid, pat)])

    
    null_mask = q & (tis == "Normal") & ~np.isin(leid, list(owner_of))
    
    ctv = obs["celltype"].astype(str).values
    M = M.copy()
    for t in ["Hepatocyte", "Cholangiocyte"]:
        nb = null_mask & (ctv == t)
        if nb.sum() >= 20:
            M[q & (ctv == t)] -= np.median(M[nb], axis=0)
            log(f"Baseline {t}: median of {nb.sum()} adjacent-normal cells subtracted")
    corr = np.full(len(obs), np.nan)
    call = np.where(is_ref, "reference", "non-malignant").astype(object)
    log(f"Owner patients (own malignant cluster): {owners}")
    for pt in owners:
        core = in_own_mal & (pat == pt)
        sig  = M[core].mean(axis=0)
        sel  = np.abs(sig) >= np.quantile(np.abs(sig), 0.80)
        S, s = M[:, sel], sig[sel]
        Sc, sc_ = S - S.mean(1, keepdims=True), s - s.mean()
        r = (Sc @ sc_) / (np.linalg.norm(Sc, axis=1) * np.linalg.norm(sc_) + 1e-12)
        nm = null_mask & (pat != pt)
        hi, lo = np.quantile(r[nm], 0.99), np.quantile(r[nm], 0.95)
        m = q & (pat == pt)
        corr[m] = r[m]
        call[m & (r > hi)] = "malignant"
        call[m & (r > lo) & (r <= hi)] = "unresolved"
        call[core] = "malignant"
        log(f"{pt}: signature corr threshold = {hi:.3f} (null n={nm.sum()}); "
            f"{(core).sum()} cells in own cluster, "
            f"{(m & ~core & (call == 'malignant')).sum()} rescued from other clusters")
    non_owner = q & ~np.isin(pat, owners)
    call[non_owner & (tis == "Tumor")] = "undetermined"
    call[non_owner & (tis == "Normal")] = "non-malignant"
    call = call.astype(str)
    obs["cnv_call"] = call
    obs["cnv_corr"] = corr

    # ---------- per-patient summary ----------
    summ = []
    for pt, d in obs[obs["is_query"]].groupby("patient", observed=True):
        t = d[d["tissue"] == "Tumor"]; n = d[d["tissue"] == "Normal"]
        summ.append(dict(patient=pt, n_tumor_query=len(t), n_normal_query=len(n),
                         n_malignant=int((d["cnv_call"] == "malignant").sum()),
                         pct_tumor_query_malignant=round(100 * (t["cnv_call"] == "malignant").mean(), 1) if len(t) else np.nan,
                         pct_normal_query_malignant=round(100 * (n["cnv_call"] == "malignant").mean(), 1) if len(n) else np.nan))
    summ = pd.DataFrame(summ)
    summ.to_csv(f"{OUT_DIR}/infercnv_summary_by_patient.csv", index=False)
    summ.to_csv(f"{BASE}/results/infercnv/infercnv_summary_by_patient.csv", index=False)

    out = obs.loc[:, ["patient", "tissue", "celltype", "cnv_group",
                                           "cnv_leiden", "cnv_score", "cnv_cell_score",
                                           "cnv_corr", "cnv_call"]].copy()
    out.index.name = "cell"
    out.to_csv(f"{OUT_DIR}/infercnv_calls.csv")
    out.to_csv(f"{BASE}/results/infercnv/infercnv_calls.csv")      # used by infercnv_3_import.R

    # ---------- figures ----------
    log("Drawing figures ...")
    a.obs["cnv_call"] = pd.Categorical(call, categories=["reference", "non-malignant", "undetermined", "unresolved", "malignant"])
    for g in ["patient", "tissue", "celltype", "cnv_call", "cnv_leiden"]:
        sc.pl.embedding(a, basis="cnv_umap", color=g, show=False, frameon=False,
                        legend_loc="on data" if g == "cnv_leiden" else "right margin",
                        legend_fontsize=6)
        savefig(f"CNV_UMAP_{g}")

    q = a[~is_ref].copy()
    q.obs["patient_call"] = (q.obs["patient"].astype(str) + "_" +
                             q.obs["cnv_call"].astype(str)).astype("category")
    cnv.pl.chromosome_heatmap(q, groupby="patient_call", show=False, figsize=(14, 12))
    savefig("Heatmap_query_by_patient_call", svg=False)
    cnv.pl.chromosome_heatmap(a, groupby="cnv_leiden", show=False, figsize=(14, 12))
    savefig("Heatmap_all_by_cnv_cluster", svg=False)

    log("Finished. Results in", OUT_DIR)
    log("\n" + summ.to_string(index=False))


if __name__ == "__main__":
    main()
