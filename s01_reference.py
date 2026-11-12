import sys
import figs                     
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

import config as C


def get_counts(a):
    src = C.SC_COUNTS
    if src.startswith("layers:"):
        X = a.layers[src.split(":", 1)[1]]
        out = sc.AnnData(X=X, obs=a.obs.copy(), var=a.var.copy())
    elif src == "raw":
        out = sc.AnnData(X=a.raw.X, obs=a.obs.copy(), var=a.raw.var.copy())
    else:
        out = sc.AnnData(X=a.X, obs=a.obs.copy(), var=a.var.copy())
    X = out.X if sp.issparse(out.X) else sp.csr_matrix(out.X)
    sample = X[:2000].data
    if sample.size and (np.any(sample < 0) or not np.allclose(sample, np.round(sample))):
        sys.exit(f"SC_COUNTS='{src}' does not look like raw integer counts - set it in config.py")
    out.X = X.astype(np.float32)
    return out


def train_kwargs(model_cls, epochs, accelerator=None, **extra):
    """
    Device selection. cell2location/scvi-tools >=1.0 use accelerator=, older versions use use_gpu=.
    accelerator: "auto" | "cpu" | "mps" | "gpu"   (default: config.ACCELERATOR)
    """
    import inspect, torch
    acc = (accelerator or getattr(C, "ACCELERATOR", "auto")).lower()
    if acc == "auto":
        acc = "gpu" if torch.cuda.is_available() else "cpu"
    if acc == "mps" and not (getattr(torch.backends, "mps", None) and torch.backends.mps.is_available()):
        print("MPS requested but not available - using CPU")
        acc = "cpu"
    if acc == "gpu" and not torch.cuda.is_available():
        print("CUDA GPU requested but not available - using CPU")
        acc = "cpu"
    params = inspect.signature(model_cls.train).parameters
    takes_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params.values())
    kw = dict(max_epochs=epochs, **extra)
    # Cell2location.train passes the device through **kwargs, so it is not in the signature
    if "accelerator" in params or (takes_kwargs and "use_gpu" not in params):
        kw["accelerator"] = acc
    elif "use_gpu" in params:
        kw["use_gpu"] = acc == "gpu"
    print(f"Training {model_cls.__name__} for {epochs} epochs on {acc.upper()}")
    return kw


def prepare_reference(plot=True):
    C.OUT_DIR.mkdir(exist_ok=True, parents=True)
    a = sc.read_h5ad(C.SC_REF_H5AD)
    for k in [C.CELLTYPE_KEY, C.REF_BATCH_KEY, C.MALIGNANT_KEY, C.TISSUE_KEY]:
        if k not in a.obs:
            sys.exit(f"'{k}' not in .obs. Available: {list(a.obs.columns)}")
    ref = get_counts(a)
    del a

    # ---- compartment-specific names for cycling cells ("Proliferating" in T/NK and Myeloid)
    lab = ref.obs[C.CELLTYPE_KEY].astype(str)
    ck = getattr(C, "COMPARTMENT_KEY", None)
    if ck and ck in ref.obs:
        comp = ref.obs[ck].astype(str)
        for cp, new in getattr(C, "PROLIF_RENAME", {}).items():
            lab[(lab == "Proliferating") & (comp == cp)] = new

    # ---- malignant split of epithelial labels
    call = ref.obs[C.MALIGNANT_KEY].astype(str).str.lower()
    epi = lab.isin(C.EPITHELIAL_TYPES)
    lab = lab.where(~(epi & (call == "malignant")), "Malignant")
    drop = epi & call.isin(["unresolved", "undetermined", "nan"])
    tum = ref.obs[C.TISSUE_KEY].astype(str).str.lower().str.startswith("t")
    drop |= epi & tum & (call != "malignant")        # non-malignant epithelium inside tumour tissue is ambiguous
    for src, dst in getattr(C, "EPI_NORMAL_RENAME", {}).items():
        lab[(lab == src) & epi & (call != "malignant")] = dst
    lab[drop] = "DROP"
    lab[lab.isin(getattr(C, "DROP_TYPES", []))] = "DROP"
    ref.obs["c2l_label"] = lab
    ref = ref[ref.obs["c2l_label"] != "DROP"].copy()
    counts = ref.obs["c2l_label"].value_counts()
    keep = counts[counts >= C.MIN_CELLS_PER_TYPE].index
    print("Reference labels (cells):\n", counts.to_string())
    if len(keep) < len(counts):
        print("Dropped (<%d cells): %s" % (C.MIN_CELLS_PER_TYPE, list(counts.index.difference(keep))))
    ref = ref[ref.obs["c2l_label"].isin(keep)].copy()

    # ---- gene filtering (cell2location recommended defaults)
    ref.var["SYMBOL"] = ref.var_names
    ref = ref[:, ~ref.var_names.str.upper().str.startswith("MT-")].copy()
    with figs.capture_shows(C.OUT_DIR / "ref_gene_filter"):
        selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    ref = ref[:, selected].copy()
    print(f"Reference: {ref.n_obs} cells x {ref.n_vars} genes, {ref.obs['c2l_label'].nunique()} types")

    # patients are nested within datasets, so the per-sample batch effect already absorbs the dataset effect
    covs = None
    if C.DATASET_KEY and C.DATASET_KEY in ref.obs and ref.obs[C.DATASET_KEY].nunique() > 1:
        nested = ref.obs.groupby(C.REF_BATCH_KEY, observed=True)[C.DATASET_KEY].nunique().max() == 1
        covs = None if nested else [C.DATASET_KEY]
    print(f"batch_key={C.REF_BATCH_KEY} ({ref.obs[C.REF_BATCH_KEY].nunique()} batches), covariates={covs}")
    RegressionModel.setup_anndata(ref, batch_key=C.REF_BATCH_KEY, labels_key="c2l_label",
                                  categorical_covariate_keys=covs)
    return ref


def signatures(ref):
    names = list(ref.uns["mod"]["factor_names"])
    key = "means_per_cluster_mu_fg"
    src = ref.varm[key] if key in ref.varm else ref.var
    inf = src[[f"{key}_{n}" for n in names]].copy()
    inf.columns = names
    return inf


def plot_signature_correlation(inf):
    """Pearson correlation of log signatures on the most variable genes: pairs close to 1 are hard to separate."""
    lg = np.log1p(inf)
    top = lg.var(axis=1).sort_values(ascending=False).index[:2000]
    cor = np.corrcoef(lg.loc[top].T.to_numpy())
    order = pd.Series(cor.mean(0), index=inf.columns).sort_values().index  # simple ordering
    idx = [list(inf.columns).index(o) for o in order]
    cor = cor[np.ix_(idx, idx)]
    n = len(order)
    fig, ax = plt.subplots(figsize=(0.38 * n + 3, 0.38 * n + 2))
    im = ax.imshow(cor, cmap="Purples", vmin=max(0, cor.min()), vmax=1)
    ax.set_xticks(range(n), order, rotation=90, fontsize=7)
    ax.set_yticks(range(n), order, fontsize=7)
    for i in range(n):
        for j in range(n):
            if i != j and cor[i, j] >= 0.9:
                ax.text(j, i, f"{cor[i, j]:.2f}", ha="center", va="center", fontsize=5, color="white")
    fig.colorbar(im, ax=ax, shrink=0.6).set_label("Pearson r (log signature, 2,000 most variable genes)", fontsize=7)
    ax.set_title("Reference signature similarity (values shown where r >= 0.90)", fontsize=9, loc="left")
    figs.save(fig, C.OUT_DIR / "ref_signature_correlation")
    pd.DataFrame(cor, index=order, columns=order).to_csv(C.OUT_DIR / "ref_signature_correlation.csv")


def main(replot=False):
    ref = prepare_reference(plot=True)
    if replot:
        mod = RegressionModel.load(str(C.OUT_DIR / "ref_model"), ref)
    else:
        try:
            import scvi
            scvi.settings.seed = 0
        except Exception:
            pass
        mod = RegressionModel(ref)
        mod.train(**train_kwargs(RegressionModel, C.REF_EPOCHS))
        mod.save(str(C.OUT_DIR / "ref_model"), overwrite=True)

    with figs.capture_shows(C.OUT_DIR / "ref_training_history"):
        mod.plot_history(20)

    ref = mod.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
    # plot_QC makes two separate plots: 1 = reconstruction accuracy, 2 = signature vs mean expression
    with figs.capture_shows(C.OUT_DIR / "ref_QC"):
        mod.plot_QC()

    inf = signatures(ref)
    plot_signature_correlation(inf)
    if not replot:
        inf.to_csv(C.OUT_DIR / "reference_signatures.csv")
        print(f"Wrote {C.OUT_DIR/'reference_signatures.csv'}  ({inf.shape[0]} genes x {inf.shape[1]} types)")
    print(f"Figures (png/svg/pdf) in {C.OUT_DIR}")


if __name__ == "__main__":
    main(replot="--replot" in sys.argv)
