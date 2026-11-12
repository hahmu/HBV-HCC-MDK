import os, sys, time

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np
import pandas as pd
import scanpy as sc
import figs                     
import matplotlib.pyplot as plt

import cell2location
from cell2location.models import Cell2location

import config as C
from s01_reference import train_kwargs
from regions import assign_regions, region_outputs


def plot_sections(adata, cols, fname, max_cols=6):
    samples = adata.obs["sample"].unique()
    cols = list(cols)[:max_cols]
    fig, axes = plt.subplots(len(samples), len(cols), figsize=(3.2 * len(cols), 3.2 * len(samples)),
                             squeeze=False)
    for i, s in enumerate(samples):
        sub = adata[adata.obs["sample"] == s]
        sf = adata.uns["spatial"][s]["scalefactors"]
        img = adata.uns["spatial"][s]["images"].get("lowres")
        xy = sub.obsm["spatial"] * sf["tissue_lowres_scalef"]
        for j, c in enumerate(cols):
            ax = axes[i, j]
            if img is not None:
                ax.imshow(img, alpha=0.5)
            v = sub.obs[c].to_numpy()
            ax.scatter(xy[:, 0], xy[:, 1], c=v, s=4, cmap="magma", rasterized=True,
                       vmin=0, vmax=np.quantile(adata.obs[c], 0.992) + 1e-9)
            ax.set_axis_off()
            if i == 0:
                ax.set_title(c, fontsize=9)
            if j == 0:
                ax.text(-0.05, 0.5, s, transform=ax.transAxes, rotation=90, va="center", ha="right")
    fig.tight_layout()
    figs.save(fig, fname, dpi=200)


def load_inputs():
    vis = sc.read_h5ad(C.OUT_DIR / "visium_raw.h5ad")
    inf = pd.read_csv(C.OUT_DIR / "reference_signatures.csv", index_col=0)
    vis = vis[:, ~vis.var["MT"].astype(bool)].copy()
    shared = vis.var_names.intersection(inf.index)
    print(f"Genes: reference {inf.shape[0]}, Visium {vis.n_vars}, shared {len(shared)}")
    if len(shared) < 3000:
        print("WARNING: <3000 shared genes - check that both use gene SYMBOLS (not Ensembl IDs).")
    return vis[:, shared].copy(), inf.loc[shared]


def benchmark(epochs=30):
    """Time the same short training on CPU and on MPS; extrapolate to SPATIAL_EPOCHS."""
    import torch
    devices = ["cpu"] + (["mps"] if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available() else [])
    only = [a.split("=", 1)[1] for a in sys.argv if a.startswith("--device=")]
    if only:
        devices = only
    if torch.cuda.is_available():
        devices.append("gpu")
    vis, inf = load_inputs()
    res = {}
    for dev in devices:
        v = vis.copy()
        Cell2location.setup_anndata(v, batch_key="sample")
        mod = Cell2location(v, cell_state_df=inf, N_cells_per_location=C.N_CELLS_PER_LOCATION,
                            detection_alpha=C.DETECTION_ALPHA)
        t0 = time.time()
        try:
            mod.train(**train_kwargs(Cell2location, epochs, accelerator=dev, batch_size=None, train_size=1))
            dt = time.time() - t0
            h = getattr(mod, "history", {}) or {}
            hist = h["elbo_train"].to_numpy().ravel() if "elbo_train" in h else np.array([np.nan])
            res[dev] = dict(sec_per_epoch=dt / epochs, est_hours=dt / epochs * C.SPATIAL_EPOCHS / 3600,
                            elbo_last=float(hist[-1]), ok=bool(np.isfinite(hist).all()) or np.isnan(hist).all())
        except Exception as e:
            res[dev] = dict(error=f"{type(e).__name__}: {str(e)[:200]}")
    print("\n================ BENCHMARK ================")
    print(f"{epochs} epochs, {vis.n_obs} spots x {vis.n_vars} genes x {inf.shape[1]} cell types")
    for dev, r in res.items():
        if "error" in r:
            print(f"  {dev.upper():4s}  FAILED  {r['error']}")
        else:
            print(f"  {dev.upper():4s}  {r['sec_per_epoch']:.2f} s/epoch -> ~{r['est_hours']:.1f} h for "
                  f"{C.SPATIAL_EPOCHS} epochs   (loss {r['elbo_last']:.3g}, finite={r['ok']})")
    ok = {d: r for d, r in res.items() if "error" not in r and r["ok"]}
    if ok:
        best = min(ok, key=lambda d: ok[d]["sec_per_epoch"])
        print(f"\nFastest working device: {best.upper()}  ->  set ACCELERATOR = \"{best}\" in config.py")
        if "cpu" in ok and best != "cpu":
            gap = abs(ok[best]["elbo_last"] - ok["cpu"]["elbo_last"]) / abs(ok["cpu"]["elbo_last"])
            print(f"Loss after {epochs} epochs differs from CPU by {gap:.1%} "
                  f"({'fine' if gap < 0.05 else 'check: larger than expected'})")


def train_chunked(vis, inf):
    """
    Train in chunks of SPATIAL_CHUNK epochs. After every chunk the model is saved to
    spatial_model_ckpt/ and the loss history appended to spatial_loss_history.csv, so an
    interrupted run (sleep, closed Terminal, crash) resumes from the last checkpoint.
    Stops at SPATIAL_EPOCHS, or earlier once the loss has plateaued (after SPATIAL_MIN_EPOCHS).
    """
    ckpt = C.OUT_DIR / "spatial_model_ckpt"
    hist_f = C.OUT_DIR / "spatial_loss_history.csv"
    Cell2location.setup_anndata(vis, batch_key="sample")
    hist = np.array([])
    if (ckpt / "model.pt").exists() and hist_f.exists():
        hist = pd.read_csv(hist_f)["elbo_train"].to_numpy()
        # Repair histories written by the first version of this script, which appended the
        # model's cumulative history after every chunk (1k + 2k + ... + k*1k entries).
        ch, n_ = C.SPATIAL_CHUNK, len(hist)
        k = int(round((np.sqrt(8 * n_ / ch + 1) - 1) / 2))
        if k > 1 and ch * k * (k + 1) // 2 == n_ and np.isclose(hist[0], hist[n_ - k * ch], rtol=0.05):
            hist = hist[-k * ch:]
            pd.DataFrame({"elbo_train": hist}).to_csv(hist_f, index=False)
            print(f"Loss history repaired: {n_} duplicated entries -> {len(hist)} real epochs")
        try:
            mod = Cell2location.load(str(ckpt), vis)
            print(f"Resuming from checkpoint after {len(hist)} epochs")
        except Exception as e:
            print(f"Checkpoint could not be loaded ({type(e).__name__}: {e}) - starting from scratch")
            hist = np.array([])
            mod = None
    else:
        mod = None
    if mod is None:
        mod = Cell2location(vis, cell_state_df=inf,
                            N_cells_per_location=C.N_CELLS_PER_LOCATION,
                            detection_alpha=C.DETECTION_ALPHA)

    W = C.PLATEAU_WINDOW
    while len(hist) < C.SPATIAL_EPOCHS:
        n = min(C.SPATIAL_CHUNK, C.SPATIAL_EPOCHS - len(hist))
        t0 = time.time()
        kw = train_kwargs(Cell2location, n, batch_size=None, train_size=1)
        try:
            mod.train(**kw, enable_progress_bar=False)
        except TypeError:                     # older scvi: no enable_progress_bar
            mod.train(**kw)
        # mod.history can hold ALL epochs trained in this session, not only this chunk:
        # keep only the last n entries (this chunk)
        h = mod.history["elbo_train"].to_numpy().ravel()[-n:]
        hist = np.r_[hist, h]
        pd.DataFrame({"elbo_train": hist}).to_csv(hist_f, index=False)
        mod.save(str(ckpt), overwrite=True)
        dt = time.time() - t0
        left = (C.SPATIAL_EPOCHS - len(hist)) * dt / n / 3600
        print(time.strftime("[%H:%M]"), f"epoch {len(hist)}/{C.SPATIAL_EPOCHS}  loss {hist[-1]:.4g}  "
              f"({dt / n:.2f} s/epoch, <= {left:.1f} h left)  checkpoint saved", flush=True)
        if not np.isfinite(hist[-1]):
            sys.exit("Loss is not finite - stop and check the inputs / device")
        if len(hist) >= max(C.SPATIAL_MIN_EPOCHS, 2 * W):
            prev, last = hist[-2 * W:-W].mean(), hist[-W:].mean()
            gain = (prev - last) / abs(last)
            print(f"   loss improvement over the last {W} epochs: {gain:.3%}")
            if gain < C.PLATEAU_TOL:
                print(f"Loss has plateaued (< {C.PLATEAU_TOL:.1%}) - stopping at epoch {len(hist)}")
                break
    return mod, hist


def plot_loss(hist):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    x = np.arange(1, len(hist) + 1)
    axes[0].plot(x, hist, lw=0.6)
    axes[0].set_title("ELBO loss, all epochs")
    k = min(len(hist) // 2, 5000)
    axes[1].plot(x[-k:], hist[-k:], lw=0.6)
    axes[1].set_title(f"last {k} epochs")
    for a in axes:
        a.set_xlabel("epoch"); a.set_ylabel("-ELBO")
    fig.tight_layout()
    figs.save(fig, C.OUT_DIR / "spatial_training_history")


def main():
    vis, inf = load_inputs()
    print(f"{vis.n_obs} spots, {vis.obs['sample'].nunique()} sections, {vis.n_vars} genes, "
          f"{inf.shape[1]} cell types")

    mod, hist = train_chunked(vis, inf)
    mod.save(str(C.OUT_DIR / "spatial_model"), overwrite=True)      # final model, before any plotting
    plot_loss(hist)

    vis = mod.export_posterior(vis, sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs})
    with figs.capture_shows(C.OUT_DIR / "spatial_QC"):
        mod.plot_QC()

    names = list(vis.uns["mod"]["factor_names"])
    q05 = pd.DataFrame(vis.obsm["q05_cell_abundance_w_sf"].to_numpy(), index=vis.obs_names, columns=names)
    vis.obs[names] = q05
    props = q05.div(q05.sum(1) + 1e-9, axis=0)
    for n in names:
        vis.obs[f"{n}_prop"] = props[n]

    
    assign_regions(vis, names)

    vis.write_h5ad(C.OUT_DIR / "visium_c2l.h5ad", compression="gzip")
    q05.groupby(vis.obs["sample"], observed=True).mean().to_csv(C.OUT_DIR / "mean_abundance_per_section.csv")
    region_outputs(vis)

    top = q05.mean().sort_values(ascending=False).index
    main_types = ["Malignant"] + [t for t in top if t != "Malignant"] if "Malignant" in names else list(top)
    for k in range(0, len(main_types), 6):
        try:
            plot_sections(vis, main_types[k:k + 6], C.OUT_DIR / f"map_celltypes_{k // 6 + 1}")
        except Exception as e:
            print(f"map_celltypes_{k // 6 + 1} skipped: {e}")
    print("Done. Check spatial_training_history, spatial_QC_* and map_celltypes_* before s03.")


if __name__ == "__main__":
    if "--benchmark" in sys.argv:
        benchmark()
    else:
        main()
