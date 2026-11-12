import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, ListedColormap
from matplotlib.patches import Patch

import config as C

K = int(os.environ.get("K", "10"))
SMOOTHED = os.environ.get("SMOOTHED", "1") == "1"
FIX_SECTIONS = [s for s in os.environ.get("SECTIONS", "").split(",") if s]
AXES = [a for a in os.environ.get("AXES", "").split(",") if a]
N_HBV, N_NONVIRAL = int(os.environ.get("N_HBV", "1")), 1
NICHE_MODE = os.environ.get("NICHE_MODE", "none").lower()
PER_PATIENT_AXES = [a for a in os.environ.get("PER_PATIENT_AXES", "MDK-SDC2").split(",") if a]
FINAL_AXES = ["SPP1-CD44", "CXCL12-CXCR4", "HLA-I-CD8", "CCL21-CCR7", "MDK-SDC2", "MDK-LRP1", "MDK-SDC3"]
SUMMARY_AXES = [a for a in os.environ.get("SUMMARY_AXES", ",".join(FINAL_AXES)).split(",") if a]
TUMOUR_TYPES = [t for t in os.environ.get("TUMOUR_TYPES", "Malignant,Hepatocyte").split(",") if t]
MIN_FRAC = float(os.environ.get("MIN_FRAC", "0.5"))
NICHE_COMP = os.environ.get("NICHE_COMP", "").strip()
if NICHE_COMP:
    NICHE_MODE = "fixed"
EDGE_LIM = 6                                    # colour range of the edge map, in spots
NONVIRAL_PREFIX = "A009"
OUT = C.OUT_DIR / "figures_final"
INT = C.OUT_DIR / "interactions"
COMP = C.OUT_DIR / ("compartments_smoothed" if SMOOTHED else "compartments") / f"K{K}"
SPOT_SCALE = float(os.environ.get("SPOT_SCALE", "1.6"))   # 1.0 = true Visium spot size; 1.6 ~ cell2location tutorial look
QMAX = 0.992                                   # same colour cap as cell2location plot_spatial
SENDER_COL, RECEIVER_COL = "#00E5FF", "#FF2DA0" # cyan + magenta -> white where both are high
TUTORIAL_COLS = ["#FFD600", "#FF7A00", "#2E7DFF", "#00C853", "#C51162", "#B0BEC5", "#FFFFFF",
                 "#AA00FF", "#00E5FF", "#8D6E63"]

try:
    import figs                                 # your helper (png + pdf); fall back if absent
    def save(fig, path, dpi=200): figs.save(fig, path, dpi=dpi)
except ImportError:
    def save(fig, path, dpi=200):
        fig.savefig(f"{path}.png", dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        fig.savefig(f"{path}.pdf", bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)


# ------------------------------------------------------------------ helpers (same maths as s03)
def neighbour_matrix(obs):
    xy = np.c_[obs["array_col"].to_numpy(float), obs["array_row"].to_numpy(float) * np.sqrt(3)]
    pairs = cKDTree(xy).query_pairs(r=2.1, output_type="ndarray")
    n = len(obs)
    rows = np.r_[pairs[:, 0], pairs[:, 1], np.arange(n)]
    cols = np.r_[pairs[:, 1], pairs[:, 0], np.arange(n)]
    W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return sp.diags(1 / np.asarray(W.sum(1)).ravel()) @ W


def patient_of(s):
    """B0013T1 / B013T5 -> B013, C004T1P -> C004, A009N -> A009"""
    import re
    m = re.match(r"([A-Za-z])0*(\d+)", s)
    return f"{m.group(1)}{int(m.group(2)):03d}" if m else s


def edge_distance(obs, ab):
    """signed distance (in spot spacings) to the tumour edge: < 0 inside tumour, > 0 outside."""
    xy = np.c_[obs["array_col"].to_numpy(float), obs["array_row"].to_numpy(float) * np.sqrt(3)]
    tt = [t for t in TUMOUR_TYPES if t in ab.columns]
    tum = ab[tt].sum(axis=1).to_numpy()
    T = tum / (ab.sum(axis=1).to_numpy() + 1e-12) >= MIN_FRAC
    d = np.full(len(xy), np.nan)
    if T.any() and (~T).any():
        d[T] = -cKDTree(xy[~T]).query(xy[T])[0] / 2.0      # neighbour spacing = 2 in these units
        d[~T] = cKDTree(xy[T]).query(xy[~T])[0] / 2.0
    return d, T


def gene_vec(X, var_index, spec):
    parts = []
    for part in spec.split("+"):
        alts = [g for g in part.split("|") if g in var_index]
        if not alts:
            return None
        M = X[:, [var_index.get_loc(g) for g in alts]]
        parts.append((M.toarray() if sp.issparse(M) else np.asarray(M)).mean(1))
    M = np.column_stack(parts)
    v = np.exp(np.log(M + 1e-3).mean(1)) - 1e-3 if M.shape[1] > 1 else M[:, 0]
    return np.clip(v, 0, None)


def as_list(x):
    return [x] if isinstance(x, str) else list(x)


def scale(v, ref=None):
    ref = v if ref is None else ref
    hi = np.quantile(ref, QMAX) if np.any(ref > 0) else 1
    return np.clip(v / (hi + 1e-12), 0, 1)


def screen_blend(layers):
    """layers: list of (values 0-1, colour). Screen blending on black, like overlapping light."""
    out = np.zeros((len(layers[0][0]), 3))
    for v, col in layers:
        out = 1 - (1 - out) * (1 - v[:, None] * np.array(to_rgb(col))[None, :])
    return np.clip(out, 0, 1)


class Section:
    def __init__(self, c2l, s):
        self.s = s
        self.m = (c2l.obs["sample"] == s).to_numpy()
        sp_ = c2l.uns["spatial"][s]
        self.sf = sp_["scalefactors"]["tissue_lowres_scalef"]
        self.img = sp_["images"].get("lowres")
        self.xy = c2l.obsm["spatial"][self.m] * self.sf
        d = sp_["scalefactors"].get("spot_diameter_fullres", 89.0) * self.sf
        self.d = d

    def axis(self, ax, title=None, dark=True):
        ax.set_facecolor("black" if dark else "white")
        if self.img is not None:
            img = self.img.mean(-1) if (dark and self.img.ndim == 3) else self.img
            ax.imshow(img, cmap="gray" if dark else None, alpha=0.35 if dark else 1)
        pad = 3 * self.d
        ax.set_xlim(self.xy[:, 0].min() - pad, self.xy[:, 0].max() + pad)
        ax.set_ylim(self.xy[:, 1].max() + pad, self.xy[:, 1].min() - pad)
        ax.set_xticks([]); ax.set_yticks([])
        for spn in ax.spines.values():
            spn.set_visible(False)
        if title:
            ax.set_title(title, fontsize=8, color="black")

    def spots(self, ax, c, **kw):
        # spot size in points^2 so circles roughly match the Visium spot diameter on screen
        fig = ax.figure
        ax.apply_aspect()                                   # measure the box actually drawn
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        xr = np.ptp(ax.get_xlim())
        pts = self.d / xr * bbox.width * 72
        return ax.scatter(self.xy[:, 0], self.xy[:, 1], c=c, s=max(pts * SPOT_SCALE, 0.6) ** 2,
                          linewidths=0, rasterized=True, **kw)


# ------------------------------------------------------------------ figures
def fig_compartments(c2l, comp, secs, sections_all):
    cnames = list(comp.columns)
    short = [c.split(":")[0] for c in cnames]
    cmap = ListedColormap(TUTORIAL_COLS[:len(cnames)])
    dom = comp.to_numpy().argmax(1)
    strength = comp.to_numpy().max(1) / (comp.to_numpy().sum(1) + 1e-12)

    ncol = 4
    nrow = int(np.ceil(len(sections_all) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * nrow + 0.8), squeeze=False)
    fig.patch.set_facecolor("white")
    for ax, s in zip(axes.ravel(), sections_all):
        S = secs[s]; S.axis(ax, s + (" (non-viral)" if s.startswith(NONVIRAL_PREFIX) else ""))
        cols = cmap(dom[S.m])
        cols[:, 3] = 0.35 + 0.65 * strength[S.m]
        S.spots(ax, cols)
    for ax in axes.ravel()[len(sections_all):]:
        ax.axis("off")
    fig.legend(handles=[Patch(color=TUTORIAL_COLS[i], label=cnames[i]) for i in range(len(cnames))],
               loc="lower center", ncol=2, fontsize=7, frameon=False)
    fig.suptitle(f"Niche compartments (NMF K={K}{', smoothed' if SMOOTHED else ''}): dominant compartment "
                 "per spot (opacity = how dominant)", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.12, 1, 0.97))
    save(fig, OUT / "A_compartments_overview")

    # tutorial-style overlay: each compartment its own colour, screen-blended
    for s in sections_all[:4]:
        S = secs[s]
        fig, ax = plt.subplots(figsize=(6.5, 6.5)); fig.patch.set_facecolor("black")
        S.axis(ax)
        layers = [(scale(comp.loc[S.m, c].to_numpy(), comp[c].to_numpy()), TUTORIAL_COLS[i])
                  for i, c in enumerate(cnames)]
        S.spots(ax, screen_blend(layers))
        ax.legend(handles=[Patch(color=TUTORIAL_COLS[i], label=short[i] + " " + cnames[i].split(": ", 1)[-1])
                           for i in range(len(cnames))], fontsize=6, loc="upper left", bbox_to_anchor=(1, 1),
                  frameon=False, labelcolor="white")
        ax.set_title(s, color="white", fontsize=9)
        save(fig, OUT / f"A_compartments_overlay_{s}")


def pick_sections(lr, name, sections_all):
    if FIX_SECTIONS:
        return [s for s in FIX_SECTIONS if s in sections_all]
    d = lr[(lr.interaction == name) & (lr.region == "All")].sort_values("z", ascending=False)
    d = d[d["sample"].isin(sections_all)]
    hbv_all = [s for s in d["sample"] if not s.startswith(NONVIRAL_PREFIX)]   # sorted by LR z
    if name in PER_PATIENT_AXES:                        # best section of every HBV patient
        seen, hbv = set(), []
        for s in hbv_all:
            if patient_of(s) not in seen:
                seen.add(patient_of(s)); hbv.append(s)
    else:
        hbv = hbv_all[:N_HBV]
    nv = [s for s in d["sample"] if s.startswith(NONVIRAL_PREFIX)][:N_NONVIRAL]
    return nv + hbv or sections_all[:2]                 # non-viral control first, then HBV


def fig_axis(it, c2l, abund, X, vidx, comp, secs, lr, nb, spot_scores, sections_all):
    name, snd, rcv, lig, rec = it
    snd = [x for x in as_list(snd) if x in abund]; rcv = [x for x in as_list(rcv) if x in abund]
    if not snd or not rcv:
        print(f"[{name}] skipped: sender/receiver labels missing"); return
    Lall, Rall = gene_vec(X, vidx, lig), gene_vec(X, vidx, rec)
    S_all = abund[snd].sum(axis=1).to_numpy(); R_all = abund[rcv].sum(axis=1).to_numpy()
    sc_col = f"{name}|celltype"
    score_all = spot_scores[sc_col].reindex(c2l.obs_names).to_numpy() if sc_col in spot_scores else None

    # niche compartment where sender and receiver are both concentrated (min of the two shares)
    frac = comp.div(comp.sum(axis=1) + 1e-12, axis=0)          # compartment proportions per spot
    sh_s = frac.mul(S_all, axis=0).sum(); sh_s /= sh_s.sum()    # where the sender's abundance sits
    sh_r = frac.mul(R_all, axis=0).sum(); sh_r /= sh_r.sum()
    base = frac.mean()                                           # share of tissue the compartment covers
    best_c = np.minimum(sh_s / base, sh_r / base).idxmax()      # both sender and receiver enriched there
    best_r = (sh_r / base).idxmax()                              # where the receiver is most enriched
    niche_c = best_r if NICHE_MODE == "receiver" else best_c
    if NICHE_MODE == "fixed":
        hits = [c for c in comp.columns if c.split(":")[0].strip() == NICHE_COMP]
        if not hits:
            raise SystemExit(f"NICHE_COMP={NICHE_COMP} not in compartments {list(comp.columns)}")
        niche_c = hits[0]
    last = {"edge": "Distance to tumour edge\norange = tumour, white = edge, purple = outside",
            "receiver": f"Receiver niche {niche_c.split(':')[0]}\n{niche_c.split(': ', 1)[-1]}",
            "both": f"Niche {niche_c.split(':')[0]}\n{niche_c.split(': ', 1)[-1]}",
            "fixed": f"Compartment {niche_c.split(':')[0]}\n{niche_c.split(': ', 1)[-1]}"}.get(NICHE_MODE)

    sel = pick_sections(lr, name, sections_all)
    cols = ["H&E", f"Sender (cyan) + receiver (magenta)\nwhite = both", f"Ligand {lig}\n(spot + 6 neighbours)",
            f"Receptor {rec}", "LR score (cell-type weighted)"] + ([last] if last else [])
    fig, axes = plt.subplots(len(sel), len(cols), figsize=(2.55 * len(cols), 2.6 * len(sel) + 0.9), squeeze=False)
    fig.patch.set_facecolor("white")
    for i, s in enumerate(sel):
        S = secs[s]; m = S.m
        W = neighbour_matrix(c2l.obs.loc[m])
        q = lr[(lr.interaction == name) & (lr["sample"] == s) & (lr.region == "All")]
        z = nb[(nb.interaction == name) & (nb["sample"] == s)]
        lab = s + (" (non-viral)" if s.startswith(NONVIRAL_PREFIX) else "")
        stats = []
        if len(z): stats.append(f"coloc z={z.coloc_z.iloc[0]:.1f}, q={z.q.iloc[0]:.2g}")
        if len(q): stats.append(f"LR z={q.z.iloc[0]:.1f}, q={q.q.iloc[0]:.2g}")

        ax = axes[i, 0]; S.axis(ax, dark=False)
        ax.set_ylabel(lab, fontsize=8)
        ax.text(0.02, 0.02, "\n".join(stats), transform=ax.transAxes, fontsize=6, va="bottom",
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))

        ax = axes[i, 1]; S.axis(ax)
        S.spots(ax, screen_blend([(scale(W @ S_all[m], S_all), SENDER_COL), (scale(R_all[m], R_all), RECEIVER_COL)]))

        for j, (vec, cmap, niche) in enumerate([(Lall, "magma", True), (Rall, "magma", False),
                                                 (score_all, "viridis", False)], start=2):
            ax = axes[i, j]; S.axis(ax)
            if vec is None:
                ax.text(0.5, 0.5, "not available", color="white", ha="center", transform=ax.transAxes, fontsize=7)
                continue
            v = vec[m]; v = W @ v if niche else v
            sc_ = S.spots(ax, v, cmap=cmap, vmin=0, vmax=np.nanquantile(vec, QMAX) + 1e-9)

        if NICHE_MODE in ("receiver", "both", "fixed"):
            ax = axes[i, 5]; S.axis(ax)
            cv = comp.loc[m, niche_c].to_numpy()
            S.spots(ax, cv, cmap="cividis", vmin=0, vmax=np.quantile(comp[niche_c], QMAX) + 1e-9)
        elif NICHE_MODE == "edge":
            ax = axes[i, 5]; S.axis(ax)
            dist, T = edge_distance(c2l.obs.loc[m], abund.loc[m])
            if np.isfinite(dist).any():
                S.spots(ax, np.clip(dist, -EDGE_LIM, EDGE_LIM), cmap="PuOr", vmin=-EDGE_LIM, vmax=EDGE_LIM)
                ax.text(0.02, 0.02, f"tumour spots {T.mean():.0%}", transform=ax.transAxes, fontsize=6,
                        color="white", va="bottom")
            else:
                ax.text(0.5, 0.5, "no tumour edge\n(" + ("all" if T.all() else "no") + " spots tumour)",
                        color="white", ha="center", va="center", transform=ax.transAxes, fontsize=7)
    for j, c in enumerate(cols):
        axes[0, j].set_title(c, fontsize=7.5)
    lbl = lambda xs: ", ".join(xs) if len(xs) <= 3 else ", ".join(xs[:3]) + f" +{len(xs) - 3}"
    fig.suptitle(f"{name}   sender: {lbl(snd)}   ->   receiver: {lbl(rcv)}", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = OUT / f"B_{name.replace('/', '_').replace(' ', '_')}"
    fig.savefig(f"{path}.svg", bbox_inches="tight", facecolor="white", dpi=200)   # for build_spatial_composite.py
    save(fig, path)
    print(f"[{name}] sections {sel}; last column = {NICHE_MODE}"
          + (f" ({niche_c})" if NICHE_MODE in ("receiver", "both", "fixed") else "")
          + f"; sender+receiver niche = {best_c}; receiver niche = {best_r}")


def fig_summary(lr, nb, names, sections_all):
    lz = lr[lr.region == "All"].pivot_table(index="interaction", columns="sample", values="z")
    lq = lr[lr.region == "All"].pivot_table(index="interaction", columns="sample", values="q")
    cz = nb.pivot_table(index="interaction", columns="sample", values="coloc_z")
    cq = nb.pivot_table(index="interaction", columns="sample", values="q")
    order_s = [s for s in sections_all if s.startswith(NONVIRAL_PREFIX)] + \
              [s for s in sections_all if not s.startswith(NONVIRAL_PREFIX)]
    order_i = [n for n in names if n in lz.index or n in cz.index]
    fig, axes = plt.subplots(1, 2, figsize=(0.42 * len(order_s) * 2 + 4, 0.42 * len(order_i) + 2.2))
    for ax, (Z, Q, t) in zip(axes, [(cz, cq, "Sender-receiver co-localisation (z)"), (lz, lq, "LR activity (z)")]):
        Z = Z.reindex(index=order_i, columns=order_s); Q = Q.reindex(index=order_i, columns=order_s)
        lim = np.nanquantile(np.abs(Z.to_numpy()), 0.95) if Z.notna().any().any() else 1
        im = ax.imshow(Z.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        for (r, c), qv in np.ndenumerate(Q.to_numpy()):
            if qv < 0.05:
                ax.text(c, r, "*", ha="center", va="center", fontsize=9)
        ax.set_xticks(range(len(order_s)), order_s, rotation=90, fontsize=7)
        ax.set_yticks(range(len(order_i)), order_i, fontsize=7.5)
        nv = sum(s.startswith(NONVIRAL_PREFIX) for s in order_s)
        if nv: ax.axvline(nv - 0.5, color="k", lw=1)
        n_other = sum(not n.startswith("MDK") for n in order_i)
        if 0 < n_other < len(order_i): ax.axhline(n_other - 0.5, color="k", lw=0.8, ls="--")
        ax.set_title(t + "\n* q < 0.05;  left of line = A009 non-viral", fontsize=8, loc="left")
        fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    save(fig, OUT / "C_axes_summary")


# ------------------------------------------------------------------ main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c2l = sc.read_h5ad(C.OUT_DIR / "visium_c2l.h5ad")
    names = list(c2l.uns["mod"]["factor_names"])
    abund = c2l.obs[names].copy()
    expr = sc.read_h5ad(C.OUT_DIR / "visium_raw.h5ad")[c2l.obs_names].copy()
    sc.pp.normalize_total(expr, target_sum=1e4); sc.pp.log1p(expr)
    X = expr.X.tocsr() if sp.issparse(expr.X) else expr.X
    vidx = expr.var_names

    comp = pd.read_csv(COMP / "spot_compartments.csv.gz", index_col=0).drop(columns="sample", errors="ignore")
    comp = comp.reindex(c2l.obs_names).fillna(0)
    lr = pd.read_csv(INT / "lr_activity_by_section_region.csv")
    nb = pd.read_csv(INT / "neighbourhood_enrichment.csv")
    spot_scores = pd.read_csv(INT / "spot_scores.csv.gz", index_col=0)

    sections_all = list(c2l.obs["sample"].unique())
    secs = {s: Section(c2l, s) for s in sections_all}
    inter = [it for it in C.INTERACTIONS if not AXES or it[0] in AXES]
    missing = [it[0] for it in inter if it[0] not in set(lr.interaction)]
    if missing:
        print("WARNING: not in s03 output (rerun s03 after editing config):", missing)

    fig_compartments(c2l, comp, secs, sections_all)
    for it in inter:
        fig_axis(it, c2l, abund, X, vidx, comp, secs, lr, nb, spot_scores, sections_all)
    summ = [a for a in SUMMARY_AXES if a in set(lr.interaction) | set(nb.interaction)]
    miss = [a for a in SUMMARY_AXES if a not in summ]
    if miss:
        print("C_axes_summary: not in s03 output, skipped (add to config.INTERACTIONS + rerun s03):", miss)
    fig_summary(lr, nb, summ, sections_all)
    print(f"Figures in {OUT}")


if __name__ == "__main__":
    main()
