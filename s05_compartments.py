import os
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.decomposition import NMF
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
import scipy.sparse as sp
import figs
import matplotlib.pyplot as plt

import config as C

SMOOTH = os.environ.get("SMOOTH", "0") == "1"
OUT = C.OUT_DIR / ("compartments_smoothed" if SMOOTH else "compartments")
N_FACT = [int(k) for k in os.environ.get("N_FACT", "5,6,7,8,10").split(",")]
N_RESAMPLE, SUBSET = 5, 0.8
TOP_FIG_K = None            # None = make maps for every K


def load():
    a = ad.read_h5ad(C.OUT_DIR / "visium_c2l.h5ad")
    names = list(a.uns["mod"]["factor_names"])
    A = a.obsm["q05_cell_abundance_w_sf"]
    A = pd.DataFrame(np.asarray(A), index=a.obs_names, columns=names)
    A = A.clip(lower=0)
    print(f"{A.shape[0]} spots x {A.shape[1]} cell types, {a.obs['sample'].nunique()} sections")
    return a, A


def neighbour_matrix(obs):
    """Row-normalised spot + 6-hex-neighbour averaging matrix (same as s03)."""
    xy = np.c_[obs["array_col"].to_numpy(float), obs["array_row"].to_numpy(float) * np.sqrt(3)]
    pairs = cKDTree(xy).query_pairs(r=2.1, output_type="ndarray")
    n = len(obs)
    rows = np.r_[pairs[:, 0], pairs[:, 1], np.arange(n)]
    cols = np.r_[pairs[:, 1], pairs[:, 0], np.arange(n)]
    W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return sp.diags(1 / np.asarray(W.sum(1)).ravel()) @ W


def smooth(a, A):
    out = A.copy()
    for s in a.obs["sample"].unique():
        m = (a.obs["sample"] == s).to_numpy()
        out.loc[m] = neighbour_matrix(a.obs.loc[m]) @ A.loc[m].to_numpy()
    return out


def fit(X, k, seed=0):
    m = NMF(n_components=k, init="nndsvda", random_state=seed, max_iter=2000, tol=1e-6,
            alpha_W=0.01, alpha_H="same", l1_ratio=0.0)
    W = m.fit_transform(X)
    return W, m.components_, m.reconstruction_err_


def stability(Hs):
    """Mean correlation of matched compartments between fits on spot subsets (1 = identical solutions)."""
    ref = Hs[0] / (Hs[0].sum(0, keepdims=True) + 1e-12)
    out = []
    for H in Hs[1:]:
        h = H / (H.sum(0, keepdims=True) + 1e-12)
        cor = np.corrcoef(ref, h)[:len(ref), len(ref):]
        r, c = linear_sum_assignment(-np.nan_to_num(cor))
        out.append(np.nanmean(cor[r, c]))
    return float(np.mean(out)) if out else np.nan


def name_compartments(frac):
    names = []
    for k in range(frac.shape[0]):
        top = frac.iloc[k].sort_values(ascending=False)
        top = top[top >= 0.25].index[:3].tolist() or [frac.iloc[k].idxmax()]
        names.append(f"C{k + 1}: " + " / ".join(top))
    return names


def dotplot(frac, path):
    order_c = frac.idxmax(axis=0).map({n: i for i, n in enumerate(frac.index)}).sort_values().index
    f = frac[order_c]
    fig, ax = plt.subplots(figsize=(0.32 * f.shape[1] + 4.5, 0.38 * f.shape[0] + 2))
    yy, xx = np.meshgrid(range(f.shape[0]), range(f.shape[1]), indexing="ij")
    v = f.to_numpy()
    s = ax.scatter(xx.ravel(), yy.ravel(), s=8 + 160 * v.ravel(), c=v.ravel(), cmap="Purples", vmin=0, vmax=1,
                   edgecolors="#555555", linewidths=0.3)
    ax.set_xticks(range(f.shape[1]), f.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(f.shape[0]), f.index, fontsize=7)
    ax.set_xlim(-0.6, f.shape[1] - 0.4); ax.set_ylim(f.shape[0] - 0.5, -0.5)
    ax.grid(color="#eeeeee", lw=0.5); ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.colorbar(s, ax=ax, shrink=0.4, pad=0.01).set_label("fraction of the cell type's\nabundance in compartment", fontsize=7)
    ax.set_title("Tissue compartments (NMF of cell2location abundance)", fontsize=9, loc="left")
    figs.save(fig, path)


def maps(a, spot, comp_names, path_stem, per_fig=4):
    samples = list(a.obs["sample"].unique())
    for start in range(0, len(comp_names), per_fig):
        cols = comp_names[start:start + per_fig]
        fig, axes = plt.subplots(len(samples), len(cols), figsize=(2.6 * len(cols), 2.5 * len(samples)), squeeze=False)
        for i, s in enumerate(samples):
            m = (a.obs["sample"] == s).to_numpy()
            sf = a.uns["spatial"][s]["scalefactors"]["tissue_lowres_scalef"]
            img = a.uns["spatial"][s]["images"].get("lowres")
            xy = a.obsm["spatial"][m] * sf
            for j, c in enumerate(cols):
                ax = axes[i, j]
                if img is not None:
                    ax.imshow(img, alpha=0.45)
                v = spot.loc[m, c].to_numpy()
                ax.scatter(xy[:, 0], xy[:, 1], c=v, s=2.5, cmap="magma", rasterized=True,
                           vmin=0, vmax=np.quantile(spot[c], 0.99) + 1e-9)
                ax.set_axis_off()
                if i == 0:
                    ax.set_title(c.replace(": ", ":\n", 1), fontsize=7)
                if j == 0:
                    ax.text(-0.05, 0.5, s, transform=ax.transAxes, rotation=90, va="center", ha="right", fontsize=8)
        fig.tight_layout()
        figs.save(fig, f"{path_stem}_{start // per_fig + 1}", dpi=150)


def axes_in_compartments(H, names, comp_names):
    """For each LIANA axis: share of sender and of receiver abundance in each compartment."""
    Hd = pd.DataFrame(H, index=comp_names, columns=names)
    rows, mat = [], {}
    for it in C.INTERACTIONS:
        name, snd, rcv = it[0], it[1], it[2]
        snd = [s for s in ([snd] if isinstance(snd, str) else snd) if s in names]
        rcv = [r for r in ([rcv] if isinstance(rcv, str) else rcv) if r in names]
        if not snd or not rcv:
            continue
        s = Hd[snd].sum(axis=1); s = s / s.sum()
        r = Hd[rcv].sum(axis=1); r = r / r.sum()
        shared = np.minimum(s, r)                          # both concentrated in the same compartment
        mat[name] = shared
        best = shared.idxmax()
        rows.append(dict(interaction=name, best_compartment=best, sender_share=round(s[best], 3),
                         receiver_share=round(r[best], 3), shared_score=round(shared[best], 3),
                         sender_top=s.idxmax(), receiver_top=r.idxmax(),
                         same_top_compartment=s.idxmax() == r.idxmax()))
    return pd.DataFrame(rows), pd.DataFrame(mat).T


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a, A = load()
    if SMOOTH:
        print("SMOOTH=1: abundance averaged over each spot and its 6 neighbours")
        A = smooth(a, A)
    X = A.to_numpy()
    rng = np.random.default_rng(0)
    stab = []
    for k in N_FACT:
        W, H, err = fit(X, k)
        subs = [fit(X[rng.random(len(X)) < SUBSET], k)[1] for _ in range(N_RESAMPLE)]
        st = stability([H] + subs)
        var_expl = 1 - err ** 2 / (X ** 2).sum()
        stab.append(dict(n_fact=k, smoothed=SMOOTH, reconstruction_error=round(err, 3),
                         variance_explained=round(var_expl, 4), stability_mean_r=round(st, 4)))
        print(f"K = {k}: variance explained {var_expl:.3f}, stability on 80 % spot subsets r = {st:.3f}")

        d = OUT / f"K{k}"
        d.mkdir(exist_ok=True)
        frac = pd.DataFrame(H / (H.sum(0, keepdims=True) + 1e-12), columns=A.columns)
        comp_names = name_compartments(frac)
        frac.index = comp_names
        frac.round(4).to_csv(d / "celltype_fractions.csv")
        spot = pd.DataFrame(W * H.sum(1), index=A.index, columns=comp_names)      # abundance per compartment
        spot.assign(sample=a.obs["sample"].values).to_csv(d / "spot_compartments.csv.gz")
        dotplot(frac, d / "celltype_fractions_dotplot")

        tab, mat = axes_in_compartments(H, list(A.columns), comp_names)
        tab.to_csv(d / "axes_in_compartments.csv", index=False)
        if len(mat):
            fig, ax = plt.subplots(figsize=(0.45 * mat.shape[1] + 3.5, 0.35 * mat.shape[0] + 1.8))
            im = ax.imshow(mat.to_numpy(), cmap="Purples", vmin=0, vmax=max(0.3, mat.to_numpy().max()), aspect="auto")
            ax.set_xticks(range(mat.shape[1]), [c.split(":")[0] for c in mat.columns], fontsize=8)
            ax.set_yticks(range(mat.shape[0]), mat.index, fontsize=8)
            fig.colorbar(im, ax=ax, shrink=0.6).set_label("min(sender share, receiver share)", fontsize=7)
            ax.set_title("LIANA axes: sender and receiver in the same compartment\n"
                         "(compartment names in celltype_fractions.csv)", fontsize=8, loc="left")
            figs.save(fig, d / "axes_in_compartments")
        if TOP_FIG_K is None or k == TOP_FIG_K:
            maps(a, spot, comp_names, d / "compartment_maps")
        print("   compartments:", *comp_names, sep="\n     ")
        if len(tab):
            print(tab[["interaction", "best_compartment", "sender_share", "receiver_share", "same_top_compartment"]]
                  .to_string(index=False))
    f = OUT / "compartment_stability.csv"
    new = pd.DataFrame(stab)
    if f.exists():
        old = pd.read_csv(f)
        if "variance_explained" in old.columns:          # drop rows from the old (seed-only) version
            new = pd.concat([old[~old.n_fact.isin(new.n_fact)], new])
    new.sort_values("n_fact").to_csv(f, index=False)
    print(pd.read_csv(f).to_string(index=False))
    print(f"\nOutputs in {OUT}. Choose K where compartments are stable (r close to 1) and biologically "
          f"distinct.")


if __name__ == "__main__":
    main()
