import os, sys, time, warnings, logging
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.getLogger("tensorflow").setLevel(logging.ERROR)
logging.getLogger("fontTools").setLevel(logging.WARNING)
import numpy as np
import pandas as pd

import tensorflow as tf
try:
    tf.config.set_visible_devices([], "GPU")
except Exception:
    pass
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.signal
if not hasattr(scipy.signal, "gaussian"):          # older arviz + scipy >= 1.13
    from scipy.signal.windows import gaussian as _gaussian
    scipy.signal.gaussian = _gaussian
import anndata as ad
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

BASE    = os.environ.get("HBV_BASE", "./HBV_HCC")
H5AD    = f"{BASE}/results/HBV_HCC_subclustered.h5ad"
OUT_DIR = f"{BASE}/results/abundance"
FIG_DIR = f"{BASE}/results/figures"
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(FIG_DIR, exist_ok=True)

RUNS        = ["major", "major_noEpi", "subtype", "subtype_noEpi",
               "stromal", "myeloid", "tnk"]   # or set env SCCODA_RUNS=stromal,myeloid,tnk
# Per-compartment runs 
COMPARTMENT_SUB = {
    "stromal": ["LSEC", "Venous EC", "Arterial EC", "Tumour EC", "Lymphatic EC", "Proliferating EC",
                "cyHSC", "Mesenchymal-like HSC", "myHSC", "vSMC"],
    "myeloid": ["Kupffer", "FOLR2+ TAM", "SPP1+ TAM", "Classical monocyte", "Non-classical monocyte",
                "cDC1", "cDC2", "mregDC", "pDC", "Proliferating myeloid"],
    "tnk":     ["CD4 T", "Tfh-like CD4 T", "Treg", "CD8 Tem", "CD8 Trm", "CD8 Temra", "MAIT", "gdT",
                "NK", "Liver-resident NK", "ILC", "Proliferating T/NK"],
}
MIN_COMP_CELLS = 20   
PAIRED      = True
REFERENCE   = "automatic"
FDR         = 0.1
N_RESULTS   = 20000
N_BURNIN    = 5000
SEED        = 1
FONT        = 8
TISSUE_COL  = {"Normal": "#4DBBD5", "Tumor": "#E64B35"}
EPI_MAJOR   = ["Hepatocyte", "Cholangiocyte"]
EPI_SUB     = ["Hepatocyte", "Proliferating hepatocyte", "Cholangiocyte", "Malignant", "Proliferating"]
COV         = ["sample", "patient", "tissue", "dataset"]
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": FONT})


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def save(fig, name):
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{FIG_DIR}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ count tables
def cell_labels():
    obs = ad.read_h5ad(H5AD, backed="r").obs.copy()
    for c in ["celltype", "subtype_fine", "compartment", "celltype_cnv", "tissue", "patient", "sample", "dataset"]:
        obs[c] = obs[c].astype(str)
    mal = obs["celltype_cnv"] == "Malignant"
    sub = obs["subtype_fine"].copy()
    sub[(sub == "Proliferating") & (obs["compartment"] == "TNK")] = "Proliferating T/NK"
    sub[(sub == "Proliferating") & (obs["compartment"] == "Myeloid")] = "Proliferating myeloid"
    sub[mal] = "Malignant"
    labels = {
        "major":         obs["celltype"],
        "major_noEpi":   obs["celltype"].where(~(obs["celltype"].isin(EPI_MAJOR) | mal)),
        "subtype":       sub,
        "subtype_noEpi": sub.where(~sub.isin(EPI_SUB)),
    }
    for run, keep in COMPARTMENT_SUB.items():
        labels[run] = sub.where(sub.isin(keep))
    # sanity check: every subtype_fine label should land in exactly one compartment run
    used = set().union(*COMPARTMENT_SUB.values()) | set(EPI_SUB)
    left = sorted(set(sub.unique()) - used)
    log("compartment values:", dict(obs["compartment"].value_counts()))
    log("subtypes in no compartment run:", left if left else "none")
    return obs, labels


def count_table(obs, lab, run):
    o = obs.assign(label=lab).dropna(subset=["label"])
    o = o[o["label"] != "nan"]
    tab = o.groupby(["sample", "label"]).size().unstack(fill_value=0)
    meta = obs.groupby("sample")[["patient", "tissue", "dataset"]].first()
    df = meta.join(tab).reset_index()
    if run in COMPARTMENT_SUB:
        types = list(tab.columns)
        df[types] = df[types].fillna(0)
        n = df[types].sum(axis=1)
        dropped = df.loc[n < MIN_COMP_CELLS, "sample"].tolist()
        df = df[n >= MIN_COMP_CELLS]
        paired = df.groupby("patient")["tissue"].nunique()
        unpaired = df.loc[df["patient"].isin(paired[paired < 2].index), "sample"].tolist()
        df = df[df["patient"].isin(paired[paired == 2].index)].reset_index(drop=True)
        log(f"{run}: dropped (< {MIN_COMP_CELLS} cells) {dropped}; dropped (unpaired) {unpaired}; "
            f"{df['patient'].nunique()} patient pairs kept")
    df.to_csv(f"{OUT_DIR}/counts_{run}.csv", index=False)
    return df


# ------------------------------------------------------------------ scCODA
def run_one(df, run):
    types = [c for c in df.columns if c not in COV]
    df = df[COV + types]
    log(f"{run}: {len(df)} samples, {len(types)} cell types, {int(df[types].values.sum())} cells")
    data = dat.from_pandas(df, covariate_columns=COV)
    formula = "C(tissue, Treatment('Normal'))" + (" + C(patient)" if PAIRED else "")
    np.random.seed(SEED); tf.random.set_seed(SEED)
    model = mod.CompositionalAnalysis(data, formula=formula, reference_cell_type=REFERENCE)
    log(f"{run}: sampling (HMC, {N_RESULTS} samples) ...")
    res = model.sample_hmc(num_results=N_RESULTS, num_burnin=N_BURNIN)
    res.set_fdr(est_fdr=FDR)

    eff = res.effect_df.reset_index()
    eff = eff[eff["Covariate"].str.contains("tissue")].copy()
    eff["credible"] = eff["Final Parameter"] != 0
    ref = getattr(model, "reference_cell_type", REFERENCE)
    ref_name = types[ref] if isinstance(ref, (int, np.integer)) else ref
    eff["reference_cell_type"] = ref_name
    eff["run"] = run
    eff = eff.sort_values("log2-fold change", ascending=False)
    eff.to_csv(f"{OUT_DIR}/scCODA_{run}.csv", index=False)
    log(f"{run}: reference cell type = {ref_name}")
    log("\n" + eff[["Cell Type", "Inclusion probability", "log2-fold change", "credible"]].to_string(index=False))
    plot_bars(eff, run, ref_name)
    plot_props(df, types, eff, run)
    return eff


def plot_bars(eff, run, ref_name):
    fig, ax = plt.subplots(figsize=(4.2, 0.2 * len(eff) + 1))
    y = np.arange(len(eff))[::-1]
    cols = np.where(eff["credible"], np.where(eff["log2-fold change"] > 0, "#E64B35", "#3C5488"), "#CCCCCC")
    ax.barh(y, eff["log2-fold change"], color=cols)
    ax.set_yticks(y); ax.set_yticklabels(eff["Cell Type"])
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("log2 fold change, Tumor vs Normal (scCODA)")
    ax.set_title(f"{run}: credible effects at FDR {FDR} (reference: {ref_name})", fontsize=FONT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, f"scCODA_log2FC_{run}")


def plot_props(df, types, eff, run):
    prop = df[types].div(df[types].sum(axis=1), axis=0) * 100
    prop[["patient", "tissue"]] = df[["patient", "tissue"]]
    order = list(eff["Cell Type"])
    ncol = 6
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(9, 1.6 * nrow + 0.3), squeeze=False)
    for i, ctp in enumerate(order):
        ax = axes[i // ncol][i % ncol]
        wide = prop.pivot_table(index="patient", columns="tissue", values=ctp)
        for _, r in wide.iterrows():
            ax.plot([0, 1], [r.get("Normal"), r.get("Tumor")], color="#BBBBBB", lw=0.5)
        for j, t in enumerate(["Normal", "Tumor"]):
            v = wide[t].dropna().values if t in wide else np.array([])
            ax.boxplot(v, positions=[j], widths=0.5, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor=TISSUE_COL[t], alpha=0.8, linewidth=0.5),
                       medianprops=dict(color="black", linewidth=0.8),
                       whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5))
            ax.scatter(np.full(len(v), j), v, s=4, color=TISSUE_COL[t], zorder=3)
        cred = bool(eff.loc[eff["Cell Type"] == ctp, "credible"].iloc[0])
        ax.set_title(ctp + (" *" if cred else ""), fontsize=FONT - 1, fontweight="bold" if cred else "normal")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["N", "T"], fontsize=FONT - 1)
        ax.tick_params(axis="y", labelsize=FONT - 2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for k in range(len(order), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    ylab = (f"% of {run} cells in sample" if run in COMPARTMENT_SUB
            else "% of cells in sample" + (" (epithelium removed)" if "noEpi" in run else ""))
    fig.supylabel(ylab, fontsize=FONT)
    fig.suptitle(f"{run}: * credible change (scCODA, FDR {FDR})", fontsize=FONT, y=1.0)
    fig.tight_layout()
    save(fig, f"scCODA_paired_props_{run}")


def plot_compare(allres, level):
    """With vs without epithelium, side by side: colour = log2FC, ring = credible."""
    runs = [r for r in (level, f"{level}_noEpi") if r in allres]
    if len(runs) < 2:
        return
    w = pd.concat([allres[r].set_index("Cell Type")[["log2-fold change", "credible"]].add_prefix(f"{r}|")
                   for r in runs], axis=1)
    w = w.sort_values(f"{runs[1]}|log2-fold change", ascending=False, na_position="last")
    fc = w[[f"{r}|log2-fold change" for r in runs]].to_numpy(float)
    cr = w[[f"{r}|credible" for r in runs]].fillna(False).to_numpy(bool)
    v = np.nanmax(np.abs(fc))
    fig, ax = plt.subplots(figsize=(3.6, 0.22 * len(w) + 1.2))
    for j in range(len(runs)):
        y = np.arange(len(w))
        ok = ~np.isnan(fc[:, j])
        s_ = ax.scatter(np.full(ok.sum(), j), y[ok], c=fc[ok, j], cmap="RdBu_r", vmin=-v, vmax=v, s=60,
                        edgecolors=np.where(cr[ok, j], "black", "none"), linewidths=1.2)
        ax.scatter(np.full((~ok).sum(), j), y[~ok], marker="x", s=12, color="#bbbbbb")
    ax.set_xticks(range(len(runs)), ["with\nepithelium", "without\nepithelium"])
    ax.set_yticks(range(len(w)), w.index, fontsize=FONT - 1)
    ax.set_xlim(-0.6, len(runs) - 0.4); ax.set_ylim(len(w) - 0.5, -0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(s_, ax=ax, shrink=0.4).set_label("log2FC Tumor vs Normal", fontsize=FONT - 1)
    ax.set_title(f"scCODA {level}: black ring = credible (FDR {FDR})\nx = not in this run",
                 fontsize=FONT, loc="left")
    save(fig, f"scCODA_compare_{level}")


def main():
    runs = os.environ.get("SCCODA_RUNS", ",".join(RUNS)).split(",")
    obs, labels = cell_labels()
    allres = {}
    for run in runs:
        df = count_table(obs, labels[run], run)
        allres[run] = run_one(df, run)
    summ = pd.concat(allres.values(), ignore_index=True)[
        ["run", "Cell Type", "log2-fold change", "Inclusion probability", "credible", "reference_cell_type"]]
    summ.to_csv(f"{OUT_DIR}/scCODA_summary.csv", index=False)
    for level in ("major", "subtype"):
        plot_compare(allres, level)
    log("Finished. Results in", OUT_DIR, "and", FIG_DIR)


if __name__ == "__main__":
    main()
