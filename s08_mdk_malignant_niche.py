#!/usr/bin/env python
"""
s08_mdk_malignant_niche.py
Which cell types / NMF compartments share the niche with MDK-high malignant cells?
(GSE281759 Visium + cell2location; conda env with scanpy)

Why not just use cell2location: its malignant-specific MDK is abundance x a fixed reference
signature, so it cannot separate MDK-high from MDK-low malignant cells. Instead:
 1  per spot: cell-type proportions from cell2location (q05 abundance / spot total)
 2  malignant-rich spots = malignant proportion >= --min-mal (default 0.5)
 3  MDK = observed log-normalised expression in the spot
 4  in malignant-rich spots of each section, for every cell type (and NMF compartment):
      partial Spearman of MDK with the proportion IN THE SPOT and in its RING of neighbours
      (adjusted for malignant proportion, cholangiocyte proportion and log total counts:
       MDK is also made by cholangiocytes, and detection scales with depth)
      + MDK-high vs MDK-low tertile contrast (mean proportion in spot + ring)
 5  combine: mean per patient, then across HBV patients (B013, C008, C004); A009 (non-viral)
    reported separately. Consistency (n sections / patients with the same sign) is the main
    evidence: neighbouring spots are not independent, so per-section p values are optimistic.
Outputs: tables, heatmap, niche-composition bar plot, maps of the best section.

Usage
  python s08_mdk_malignant_niche.py --h5ad c2l_out/sp_with_c2l.h5ad --section-col sample \
      [--nmf c2l_out/nmf_K10_smoothed_spot_weights.csv] --out c2l_out/mdk_niche
  --nmf: optional CSV, rows = spot IDs (same as h5ad obs_names), columns = compartments
"""
import argparse, os, re, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PATIENT = lambda s: re.match(r"([A-Z]\d+)", s).group(1) if re.match(r"([A-Z]\d+)", s) else s
NONVIRAL = {"A009"}
plt.rcParams.update({"font.size": 8})


def abundance(ad):
    for k in ["q05_cell_abundance_w_sf", "means_cell_abundance_w_sf", "q50_cell_abundance_w_sf"]:
        if k in ad.obsm:
            A = pd.DataFrame(np.asarray(ad.obsm[k]), index=ad.obs_names)
            cols = list(ad.obsm[k].columns) if hasattr(ad.obsm[k], "columns") else None
            if cols is None and "mod" in ad.uns and "factor_names" in ad.uns["mod"]:
                cols = list(ad.uns["mod"]["factor_names"])              # s02 stores names here
            if cols is None: sys.exit(f"obsm['{k}'] has no cell-type names")
            A.columns = [re.sub(r"^(q05|means|q50)cell_abundance_w_sf_", "", str(c)) for c in cols]
            A = A.clip(lower=0)
            print(f"cell2location abundance: obsm['{k}'] ({A.shape[1]} cell types)")
            return A
    cols = [c for c in ad.obs.columns if c.startswith("q05cell_abundance_w_sf_")]
    if cols:
        A = ad.obs[cols].copy(); A.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in cols]; return A
    sys.exit("No cell2location abundance in obsm/obs")


def lognorm_gene(ad, gene):
    X = ad.layers["counts"] if "counts" in ad.layers else (ad.raw.X if ad.raw is not None else ad.X)
    var = ad.raw.var_names if ("counts" not in ad.layers and ad.raw is not None) else ad.var_names
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1
    i = list(var).index(gene)
    x = X[:, i]; x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    return np.log1p(x / lib * 1e4), np.log1p(lib)


def ring(xy):
    d, idx = cKDTree(xy).query(xy, k=min(7, len(xy)))
    step = np.median(d[:, 1]); idx = idx[:, 1:].astype(float); idx[d[:, 1:] > 1.5 * step] = np.nan
    return idx


def ring_mean(idx, v):
    v = np.asarray(v, float); out = np.full(len(idx), np.nan)
    ok = ~np.isnan(idx)
    vals = np.where(ok, v[np.nan_to_num(idx).astype(int)], np.nan)
    with np.errstate(invalid="ignore"): out = np.nanmean(vals, 1)
    return out


def partial_rho(x, y, Z):
    rk = stats.rankdata
    Zr = np.column_stack([np.ones(len(x))] + [rk(z) for z in Z])
    res = lambda a: rk(a) - Zr @ np.linalg.lstsq(Zr, rk(a), rcond=None)[0]
    r, p = stats.pearsonr(res(x), res(y))
    return r, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--section-col", default="sample")
    ap.add_argument("--nmf")
    ap.add_argument("--gene", default="MDK")
    ap.add_argument("--malignant", default="Malignant")
    ap.add_argument("--chol", default="Cholangiocyte")
    ap.add_argument("--min-mal", type=float, default=0.5)
    ap.add_argument("--min-spots", type=int, default=50)
    ap.add_argument("--exclude", nargs="*", default=["C004T1R"])
    ap.add_argument("--out", default="mdk_niche")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ad = sc.read_h5ad(a.h5ad)
    ad = ad[~ad.obs[a.section_col].astype(str).isin(a.exclude)].copy()
    A = abundance(ad)
    P = A.div(A.sum(1).replace(0, np.nan), axis=0).fillna(0)
    if a.malignant not in P: sys.exit(f"'{a.malignant}' not in {list(P.columns)}")
    feats = {f"celltype:{c}": P[c].values for c in P.columns if c != a.malignant}
    if a.nmf:
        N = pd.read_csv(a.nmf, index_col=0)
        N = N.select_dtypes("number").reindex(ad.obs_names)          # drops the 'sample' column of s05
        print(f"NMF: {N.shape[1]} compartments, {N.notna().all(1).mean():.0%} of spots matched")
        N = N.div(N.sum(1).replace(0, np.nan), axis=0)
        feats.update({f"NMF:{c}": N[c].values for c in N.columns})
    mdk, loglib = lognorm_gene(ad, a.gene)
    mal = P[a.malignant].values
    chol = P[a.chol].values if a.chol in P else np.zeros(ad.n_obs)
    sec = ad.obs[a.section_col].astype(str).values
    xy = np.asarray(ad.obsm["spatial"], float)

    rows, comp = [], []
    for s in pd.unique(sec):
        m = sec == s
        idx = ring(xy[m])
        rich = mal[m] >= a.min_mal
        n = rich.sum()
        print(f"[{s}] {m.sum()} spots, malignant-rich {n}, {a.gene}+ among them {(mdk[m][rich] > 0).mean():.0%}")
        if n < a.min_spots: continue
        g = mdk[m][rich]
        Z = [mal[m][rich], chol[m][rich], loglib[m][rich]]
        t = pd.qcut(stats.rankdata(g, "ordinal"), 3, labels=False)        # tertiles (ties broken by order)
        hi, lo = t == 2, t == 0
        for f, v in feats.items():
            v_spot = v[m]; v_ring = ring_mean(idx, v_spot)
            vs, vr = v_spot[rich], v_ring[rich]
            ok = ~np.isnan(vr)
            r1, p1 = partial_rho(g, vs, Z)
            r2, p2 = partial_rho(g[ok], vr[ok], [z[ok] for z in Z])
            niche = (vs + np.nan_to_num(vr)) / 2
            rows.append(dict(section=s, patient=PATIENT(s), hbv=PATIENT(s) not in NONVIRAL, feature=f, n_spots=n,
                             rho_spot=r1, p_spot=p1, rho_ring=r2, p_ring=p2,
                             niche_MDKhigh=niche[hi].mean(), niche_MDKlow=niche[lo].mean(),
                             log2_hi_vs_lo=np.log2((niche[hi].mean() + 1e-3) / (niche[lo].mean() + 1e-3))))
        comp.append(pd.DataFrame({"section": s, "tertile": np.where(hi, "MDK-high", np.where(lo, "MDK-low", "mid")),
                                  **{f: (feats[f][m][rich] + np.nan_to_num(ring_mean(idx, feats[f][m])[rich])) / 2
                                     for f in feats}}))
    R = pd.DataFrame(rows)
    R.to_csv(f"{a.out}/mdk_niche_by_section.csv", index=False)

    # combine: per patient mean, then HBV patients
    pp = R.groupby(["feature", "patient", "hbv"])[["rho_spot", "rho_ring", "log2_hi_vs_lo"]].mean().reset_index()
    out = []
    for f, d in pp.groupby("feature"):
        h = d[d.hbv]; nv = d[~d.hbv]; secs = R[(R.feature == f) & R.hbv]
        out.append(dict(feature=f,
                        HBV_mean_rho_ring=h.rho_ring.mean(), HBV_mean_rho_spot=h.rho_spot.mean(),
                        HBV_mean_log2_hi_vs_lo=h.log2_hi_vs_lo.mean(),
                        HBV_patients_pos=int((h.rho_ring > 0).sum()), HBV_patients=len(h),
                        HBV_sections_pos_sig=int(((secs.rho_ring > 0) & (secs.p_ring < 0.05)).sum()),
                        HBV_sections_neg_sig=int(((secs.rho_ring < 0) & (secs.p_ring < 0.05)).sum()),
                        HBV_sections=len(secs),
                        A009_rho_ring=nv.rho_ring.mean() if len(nv) else np.nan))
    S = pd.DataFrame(out).sort_values("HBV_mean_rho_ring", ascending=False)
    S.to_csv(f"{a.out}/mdk_niche_summary.csv", index=False)
    pd.set_option("display.width", 220)
    print("\nCell types / compartments around MDK-high malignant spots (ring = neighbouring spots; "
          "adjusted for malignant %, cholangiocyte %, depth)")
    print(S.round(3).to_string(index=False))

    # figures
    order = S.feature.tolist()
    H = R.pivot_table(index="feature", columns="section", values="rho_ring").reindex(order)
    H = H[sorted(H.columns, key=lambda c: (PATIENT(c) in NONVIRAL, c))]
    fig, ax = plt.subplots(figsize=(0.45 * H.shape[1] + 3.5, 0.2 * H.shape[0] + 1.5))
    v = np.nanmax(np.abs(H.values))
    im = ax.imshow(H.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(H.shape[1]), H.columns, rotation=90); ax.set_yticks(range(H.shape[0]), H.index, fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.5).set_label(f"partial rho, {a.gene} vs neighbourhood proportion")
    ax.set_title(f"Who surrounds {a.gene}-high malignant spots (A009 = non-viral, right)", loc="left")
    fig.tight_layout(); fig.savefig(f"{a.out}/mdk_niche_heatmap.png", dpi=200); plt.close(fig)

    C = pd.concat(comp); C = C[C.tertile != "mid"]
    C = C[C.section.map(PATIENT).map(lambda p: p not in NONVIRAL)]
    top = [f for f in order if f.startswith("celltype:")][:8] + [f for f in order if f.startswith("celltype:")][-5:]
    M = C.groupby(["section", "tertile"])[top].mean().groupby("tertile").mean().T
    fig, ax = plt.subplots(figsize=(5, 0.28 * len(top) + 1))
    y = np.arange(len(top))
    ax.barh(y + 0.2, M["MDK-high"], 0.4, label=f"{a.gene}-high malignant spots", color="#B2182B")
    ax.barh(y - 0.2, M["MDK-low"], 0.4, label=f"{a.gene}-low malignant spots", color="#4393C3")
    ax.set_yticks(y, [t.split(":", 1)[1] for t in top]); ax.invert_yaxis(); ax.legend(fontsize=7)
    ax.set_xlabel("mean proportion in spot + neighbours (HBV sections)")
    fig.tight_layout(); fig.savefig(f"{a.out}/mdk_niche_composition.png", dpi=200); plt.close(fig)

    best = R[R.hbv].groupby("section").n_spots.first().idxmax()
    f1 = [f for f in order if f.startswith("celltype:")][0]
    m = sec == best
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    for k, (v, t, cm) in enumerate([(mal[m], f"{a.malignant} proportion", "Reds"),
                                    (np.where(mal[m] >= a.min_mal, mdk[m], np.nan), f"{a.gene} in malignant-rich spots", "magma"),
                                    (feats[f1][m], f1.split(":", 1)[1], "Blues")]):
        s_ = ax[k].scatter(xy[m, 0], -xy[m, 1], c=v, s=4, cmap=cm); ax[k].set_title(t)
        ax[k].set_aspect("equal"); ax[k].axis("off"); fig.colorbar(s_, ax=ax[k], shrink=0.6)
    fig.suptitle(best); fig.tight_layout(); fig.savefig(f"{a.out}/mdk_niche_map_{best}.png", dpi=170); plt.close(fig)
    print(f"\nDone. Results in {a.out}")


if __name__ == "__main__":
    main()
