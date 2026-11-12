import os, glob, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import liana as li
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.environ.get("HBV_BASE", "./HBV_HCC")
DEG_DIR = f"{BASE}/results/DEG"
IN_DIR  = f"{BASE}/results/liana/input"
OUT_DIR = f"{BASE}/results/liana/de_informed"
FIG_DIR = f"{OUT_DIR}/figures"
os.makedirs(FIG_DIR, exist_ok=True)

GROUP     = "de_group"
EXPR_PROP = 0.1          # ligand / receptor detected in >= 10 % of the sender / receiver cells
STAT_KEYS = ["stat", "pvalue", "padj", "log2FoldChange"]
TOP_N     = 30
PADJ      = 0.05         # for the "top" tables: ligand AND receptor padj below this
EXCLUDE_GENES = {"B2M", "ALB", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3",
                 "APOE", "APOH", "HP", "HPX", "TF", "TTR", "FGA", "FGB", "FGG", "SERPINA1",
                 "AHSG", "ORM1", "ORM2", "RBP4", "AMBP", "CP", "VTN", "FTL", "FTH1"}
FONT = 8
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": FONT})


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def save(fig, name):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{FIG_DIR}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ inputs
def load_dea():
    tabs = []
    for f in sorted(glob.glob(f"{DEG_DIR}/DEG_TumorVsNormal_*.csv")):
        ct = os.path.basename(f)[len("DEG_TumorVsNormal_"):-4]
        if ct == "summary":
            continue
        tabs.append(pd.read_csv(f).assign(**{GROUP: ct}))
    f = f"{DEG_DIR}/DEG_Malignant_vs_NormalHepatocyte.csv"
    if os.path.exists(f):
        tabs.append(pd.read_csv(f).assign(**{GROUP: "Malignant"}))
    if not tabs:
        raise SystemExit(f"No DESeq2 tables in {DEG_DIR} - run DESeq2.R first")
    dea = pd.concat(tabs, ignore_index=True)
    dea = dea.dropna(subset=["stat"]).set_index("gene")
    log("DESeq2 groups:", dea[GROUP].value_counts().to_dict())
    return dea


def load_adata(groups):
    X = scipy.io.mmread(f"{IN_DIR}/counts.mtx").T.tocsr()
    genes = pd.read_csv(f"{IN_DIR}/genes.txt", header=None)[0].values
    cells = pd.read_csv(f"{IN_DIR}/cells.txt", header=None)[0].values
    meta = pd.read_csv(f"{IN_DIR}/meta.csv", index_col=0).loc[cells]
    a = ad.AnnData(X=X.astype(np.float32), obs=meta, var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    # major type, with inferCNV-malignant cells as their own group (matches DESeq2 part B)
    g = a.obs["celltype"].astype(str).where(a.obs["liana_label"].astype(str) != "Malignant", "Malignant")
    a.obs[GROUP] = g
    a = a[a.obs[GROUP].isin(groups)].copy()
    a.obs[GROUP] = a.obs[GROUP].astype("category")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    log(f"{a.n_obs} cells in {a.obs[GROUP].nunique()} groups")
    return a


# ------------------------------------------------------------------ fallback (same logic as df_to_lr)
def manual_df_to_lr(a, dea):
    """Used only if li.multi.df_to_lr is unavailable in the installed LIANA version."""
    res = li.resource.select_resource("consensus")
    groups = list(a.obs[GROUP].cat.categories)
    Xb = a.X > 0
    props = pd.DataFrame({g: np.asarray(Xb[(a.obs[GROUP] == g).values].mean(0)).ravel() for g in groups},
                         index=a.var_names)

    def side(col):
        rows = []
        for g in groups:
            d = dea[dea[GROUP] == g]
            d = d[~d.index.duplicated()]
            for cpx in res[col].unique():
                subs = cpx.split("_")
                if not all(s in d.index and s in props.index for s in subs):
                    continue
                k = min(subs, key=lambda s: abs(d.at[s, "stat"]))   # complex = its least-changed subunit
                r = d.loc[k]
                rows.append({"group": g, "complex": cpx, "stat": r["stat"], "pvalue": r["pvalue"],
                             "padj": r["padj"], "log2FoldChange": r["log2FoldChange"],
                             "props": min(props.at[s, g] for s in subs)})
        return pd.DataFrame(rows)

    L = side("ligand").add_prefix("ligand_").rename(columns={"ligand_group": "source",
                                                             "ligand_complex": "ligand_complex"})
    R = side("receptor").add_prefix("receptor_").rename(columns={"receptor_group": "target",
                                                                 "receptor_complex": "receptor_complex"})
    pairs = res.rename(columns={"ligand": "ligand_complex", "receptor": "receptor_complex"})
    out = pairs.merge(L, on="ligand_complex").merge(R, on="receptor_complex")
    out = out[(out["ligand_props"] >= EXPR_PROP) & (out["receptor_props"] >= EXPR_PROP)]
    out["interaction_stat"] = (out["ligand_stat"] + out["receptor_stat"]) / 2
    return out


def run_df_to_lr(a, dea):
    try:
        out = li.multi.df_to_lr(a, dea_df=dea, resource_name="consensus", expr_prop=EXPR_PROP,
                                groupby=GROUP, stat_keys=STAT_KEYS, use_raw=False,
                                complex_col="stat", verbose=True, return_all_lrs=False)
        log("li.multi.df_to_lr used")
    except Exception as e:                                            # noqa: BLE001
        log(f"li.multi.df_to_lr failed ({type(e).__name__}: {e}) - using the equivalent manual version")
        out = manual_df_to_lr(a, dea)
    if "interaction_stat" not in out:
        out["interaction_stat"] = (out["ligand_stat"] + out["receptor_stat"]) / 2
    return out


# ------------------------------------------------------------------ figure
def plot_top(df, title, name):
    if df.empty:
        log(f"{name}: nothing to plot"); return
    rows = [f"{l} → {r}   {s} → {t}" for l, r, s, t in
            zip(df.ligand_complex, df.receptor_complex, df.source, df.target)]
    st = df[["ligand_stat", "receptor_stat"]].to_numpy(float)
    pa = df[["ligand_padj", "receptor_padj"]].to_numpy(float)
    size = 12 + 18 * np.clip(-np.log10(np.clip(pa, 1e-300, 1)), 0, 10)
    v = np.nanmax(np.abs(st))
    fig, ax = plt.subplots(figsize=(6.4, 0.22 * len(df) + 1.2))
    yy, xx = np.meshgrid(np.arange(len(df)), [0, 1], indexing="ij")
    s_ = ax.scatter(xx.ravel(), yy.ravel(), c=st.ravel(), s=size.ravel(), cmap="RdBu_r",
                    vmin=-v, vmax=v, edgecolors="#555555", linewidths=0.4)
    ax.set_xticks([0, 1], ["ligand\n(in sender)", "receptor\n(in receiver)"])
    ax.set_yticks(range(len(df)), rows, fontsize=FONT - 1)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(len(df) - 0.5, -0.5)
    ax.grid(color="#eeeeee", lw=0.5); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(s_, ax=ax, shrink=0.4, pad=0.02)
    cb.set_label("DESeq2 Wald stat\nTumor vs Normal", fontsize=FONT - 1)
    ax.set_title(title + "\n(dot size = -log10 padj, capped at 10)", fontsize=FONT + 1, loc="left")
    save(fig, name)


def main():
    dea = load_dea()
    a = load_adata(dea[GROUP].unique())
    out = run_df_to_lr(a, dea)
    bad = lambda c: any(g in EXCLUDE_GENES for g in str(c).split("_"))
    out = out[~out["ligand_complex"].map(bad) & ~out["receptor_complex"].map(bad)]
    out = out.sort_values("interaction_stat", ascending=False)
    out.to_csv(f"{OUT_DIR}/lr_de_informed_all.csv", index=False)
    log(f"{len(out)} sender-receiver-LR rows scored")

    sig = out[(out["ligand_padj"] < PADJ) & (out["receptor_padj"] < PADJ)]
    up = sig[(sig["ligand_stat"] > 0) & (sig["receptor_stat"] > 0)].head(TOP_N)
    dn = sig[(sig["ligand_stat"] < 0) & (sig["receptor_stat"] < 0)].sort_values("interaction_stat").head(TOP_N)
    up.to_csv(f"{OUT_DIR}/lr_de_informed_top_up.csv", index=False)
    dn.to_csv(f"{OUT_DIR}/lr_de_informed_top_down.csv", index=False)
    log(f"ligand and receptor both padj < {PADJ}: {len(sig)} rows "
        f"({(sig.ligand_stat > 0).mul(sig.receptor_stat > 0).sum()} both up, "
        f"{(sig.ligand_stat < 0).mul(sig.receptor_stat < 0).sum()} both down)")
    plot_top(up, "Interactions up in tumour (ligand and receptor both up, padj < 0.05)", "LIANA_DE_top_up")
    plot_top(dn, "Interactions down in tumour (ligand and receptor both down, padj < 0.05)", "LIANA_DE_top_down")
    if len(up):
        log("\nTop up in tumour:\n" + up[["source", "target", "ligand_complex", "receptor_complex",
                                          "ligand_stat", "receptor_stat", "interaction_stat"]]
            .head(15).round(2).to_string(index=False))
    log("Finished. Results in", OUT_DIR)


if __name__ == "__main__":
    main()
