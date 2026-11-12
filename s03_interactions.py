import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as C

OUT = C.OUT_DIR / "interactions"
RNG = np.random.default_rng(0)
REGIONS = ["Tumour", "Border", "Non-tumour"]


# ------------------------------------------------------------------ helpers
def neighbour_matrix(obs):
    """Row-normalised spot + 6-hex-neighbour averaging matrix (Visium array coordinates)."""
    xy = np.c_[obs["array_col"].to_numpy(float), obs["array_row"].to_numpy(float) * np.sqrt(3)]
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r=2.1, output_type="ndarray")
    n = len(obs)
    rows = np.r_[pairs[:, 0], pairs[:, 1], np.arange(n)]
    cols = np.r_[pairs[:, 1], pairs[:, 0], np.arange(n)]
    W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    deg = np.asarray(W.sum(1)).ravel()
    return sp.diags(1 / deg) @ W


def _genes(spec):
    return [g for part in spec.split("+") for g in part.split("|")]


def gene_vec(X, var_index, spec):
    """'A+B' = complex (geometric mean, all subunits needed); 'A|B' = alternatives (mean of present genes)."""
    parts = []
    for part in spec.split("+"):
        alts = [g for g in part.split("|") if g in var_index]
        if not alts:
            return None, [part]
        M = X[:, [var_index.get_loc(g) for g in alts]]
        M = M.toarray() if sp.issparse(M) else np.asarray(M)
        parts.append(M.mean(1))
    M = np.column_stack(parts)
    v = np.exp(np.log(M + 1e-3).mean(1)) - 1e-3 if M.shape[1] > 1 else M[:, 0]
    return np.clip(v, 0, None), []


def resolve_interactions(names, var_index, props, abund):
    """Sum fractions of listed sender/receiver labels into new columns; report missing labels/genes."""
    inter = []
    for name, snd, rcv, lig, rec in C.INTERACTIONS:
        it = dict(name=name, ligand=lig, receptor=rec)
        issues = []
        for role, labs in (("sender", snd), ("receiver", rcv)):
            labs = [labs] if isinstance(labs, str) else list(labs)
            present = [l for l in labs if l in names]
            absent = [l for l in labs if l not in names]
            if absent:
                issues.append(f"{role} labels not in reference: {absent}")
            if present:
                key = "+".join(present)
                props[key] = props[present].sum(1)
                abund[key] = abund[present].sum(1)
                it[role] = key
            else:
                it[role] = None
        for spec in (lig, rec):
            for part in spec.split("+"):
                alts = part.split("|")
                miss = [g for g in alts if g not in var_index]
                if len(miss) == len(alts):
                    issues.append(f"gene(s) {part} absent from Visium")
        it["ct_ok"] = it["sender"] is not None and it["receiver"] is not None
        if issues:
            print(f"[{name}] " + "; ".join(issues))
        inter.append(it)
    return inter


def perm_test(stat_fn, x, n_perm):
    obs = stat_fn(x)
    perms = np.column_stack([RNG.permutation(x) for _ in range(n_perm)])
    null = stat_fn(perms)
    z = (obs - null.mean(-1)) / (null.std(-1) + 1e-12)
    p = ((null >= obs[..., None]).sum(-1) + 1) / (n_perm + 1)
    return obs, z, p


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    o = np.argsort(q)
    r = q[o] * len(q) / (np.arange(len(q)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    tmp = np.empty_like(q); tmp[o] = np.minimum(r, 1)
    out[ok] = tmp
    return out


# ------------------------------------------------------------------ main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c2l = sc.read_h5ad(C.OUT_DIR / "visium_c2l.h5ad")
    expr = sc.read_h5ad(C.OUT_DIR / "visium_raw.h5ad")[c2l.obs_names].copy()  # all genes incl. LR
    sc.pp.normalize_total(expr, target_sum=1e4)
    sc.pp.log1p(expr)
    X = expr.X.tocsr() if sp.issparse(expr.X) else expr.X
    vidx = expr.var_names

    names = list(c2l.uns["mod"]["factor_names"])
    props = c2l.obs[[f"{n}_prop" for n in names]].to_numpy()
    props = pd.DataFrame(props, index=c2l.obs_names, columns=names)
    abund = c2l.obs[names]

    props, abund = props.copy(), abund.copy()
    inter = resolve_interactions(names, vidx, props, abund)
    print("Reference labels available:", names)

    lr_rows, nb_rows, coloc = [], [], []
    spot_scores = pd.DataFrame(index=c2l.obs_names)
    for s in c2l.obs["sample"].unique():
        idx = np.where(c2l.obs["sample"] == s)[0]
        obs = c2l.obs.iloc[idx]
        W = neighbour_matrix(obs)
        Xs = X[idx]
        P = props.iloc[idx]
        A = abund.iloc[idx]
        region = obs["region_auto"].to_numpy()
        label = obs["region_label"].iloc[0]

        # A1. within-spot co-localisation
        rho = spearmanr(A[names].to_numpy()).statistic if len(names) > 2 else np.array([[1]])
        coloc.append(pd.DataFrame(rho, index=names, columns=names).assign(sample=s))

        for it in inter:
            base = dict(sample=s, section_label=label, interaction=it["name"],
                        sender=it["sender"], receiver=it["receiver"])
            # A2. neighbourhood enrichment sender(niche) x receiver(spot)
            if it["ct_ok"]:
                snd = P[it["sender"]].to_numpy(); rcv = P[it["receiver"]].to_numpy()
                stat = lambda v: (W @ v * rcv[:, None]).mean(0) if v.ndim == 2 else np.atleast_1d((W @ v * rcv).mean())
                o, z, p = perm_test(stat, snd, C.N_PERMUTATIONS)
                nb_rows.append({**base, "coloc_obs": float(o[0]), "coloc_z": float(z[0]), "coloc_p": float(p[0])})

            # B. LR activity
            L, missL = gene_vec(Xs, vidx, it["ligand"])
            R, missR = gene_vec(Xs, vidx, it["receptor"])
            if L is None or R is None:
                continue
            if it["ct_ok"]:
                wct = np.sqrt((W @ P[it["sender"]].to_numpy()) * P[it["receiver"]].to_numpy())
            else:
                wct = np.ones(len(idx))

            def lr_stat(lv, mask=None):
                ln = W @ lv
                sc_ = np.sqrt(ln * (R if lv.ndim == 1 else R[:, None])) * (wct if lv.ndim == 1 else wct[:, None])
                if mask is not None:
                    sc_ = sc_[mask]
                return np.atleast_1d(sc_.mean(0))

            spot = np.sqrt((W @ L) * R)
            spot_scores.loc[obs.index, f"{it['name']}|LR"] = spot
            spot_scores.loc[obs.index, f"{it['name']}|celltype"] = spot * wct
            sanity_s = spearmanr(L, A[it["sender"]]).statistic if it["ct_ok"] else np.nan
            sanity_r = spearmanr(R, A[it["receiver"]]).statistic if it["ct_ok"] else np.nan
            for reg in ["All"] + REGIONS:
                mask = None if reg == "All" else region == reg
                if mask is not None and mask.sum() < 20:
                    continue
                o, z, p = perm_test(lambda v: lr_stat(v, mask), L, C.N_PERMUTATIONS)
                lr_rows.append({**base, "region": reg,
                                "n_spots": int(len(idx) if mask is None else mask.sum()),
                                "score": float(o[0]), "z": float(z[0]), "p": float(p[0]),
                                "ligand_pct_spots": float((L > 0).mean() * 100),
                                "receptor_pct_spots": float((R > 0).mean() * 100),
                                "rho_ligand_vs_sender": sanity_s, "rho_receptor_vs_receiver": sanity_r})
        print(f"{s}: {len(idx)} spots, regions {pd.Series(region).value_counts().to_dict()}")

    # ---- write tables
    lr = pd.DataFrame(lr_rows)
    if len(lr):
        lr["q"] = bh(lr["p"])
        lr.to_csv(OUT / "lr_activity_by_section_region.csv", index=False)
        summ = (lr[lr.region != "All"].pivot_table(index="interaction", columns="region", values="score",
                                                   aggfunc="mean").reindex(columns=REGIONS))
        summ["Tumour_vs_NonTumour_log2FC"] = np.log2((summ["Tumour"] + 1e-6) / (summ["Non-tumour"] + 1e-6))
        sig = lr[lr.region == "All"].groupby("interaction").apply(
            lambda d: pd.Series({"sections_q<0.05": int((d.q < 0.05).sum()), "sections": len(d),
                                 "median_z": d.z.median()}))
        summ = summ.join(sig).sort_values("median_z", ascending=False)
        summ.to_csv(OUT / "lr_summary.csv")
        print("\nLR summary (mean cell-type-weighted score by region):\n", summ.round(3).to_string())
    nb = pd.DataFrame(nb_rows)
    if len(nb):
        nb["q"] = bh(nb["coloc_p"])
        nb.to_csv(OUT / "neighbourhood_enrichment.csv", index=False)
    co = pd.concat(coloc)
    mean_co = co.groupby(level=0).mean(numeric_only=True).loc[names, names]
    mean_co.to_csv(OUT / "celltype_colocalisation_spearman_mean.csv")
    spot_scores.to_csv(OUT / "spot_scores.csv.gz")

    # ---- figures
    fig, ax = plt.subplots(figsize=(0.45 * len(names) + 3, 0.45 * len(names) + 2))
    im = ax.imshow(mean_co.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6)
    ax.set_xticks(range(len(names)), names, rotation=90, fontsize=8)
    ax.set_yticks(range(len(names)), names, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.6, label="mean Spearman across sections")
    ax.set_title("Cell-type co-localisation (cell2location abundance)")
    fig.tight_layout(); fig.savefig(OUT / "celltype_colocalisation.png", dpi=150); plt.close(fig)

    if len(lr):
        reg = lr[lr.region.isin(REGIONS)]
        inames = summ.index.tolist()
        fig, axes = plt.subplots(1, len(inames), figsize=(2.3 * len(inames), 3.2), squeeze=False)
        colors = {"Tumour": "#c0392b", "Border": "#e1a33b", "Non-tumour": "#3b7dc0"}
        for ax, nm in zip(axes[0], inames):
            d = reg[reg.interaction == nm]
            for k, r in enumerate(REGIONS):
                v = d.loc[d.region == r, "score"].to_numpy()
                ax.scatter(np.full(len(v), k) + RNG.uniform(-0.12, 0.12, len(v)), v, s=14, color=colors[r])
                if len(v):
                    ax.hlines(np.median(v), k - 0.3, k + 0.3, color="k", lw=1.2)
            ax.set_xticks(range(3), ["T", "B", "NT"]); ax.set_title(nm, fontsize=8)
        axes[0, 0].set_ylabel("cell-type-weighted LR score\n(per section)")
        fig.tight_layout(); fig.savefig(OUT / "lr_by_region.png", dpi=150); plt.close(fig)

        # spatial maps of the cell-type-weighted score, all sections x interactions
        samples = c2l.obs["sample"].unique()
        fig, axes = plt.subplots(len(samples), len(inames), figsize=(2.6 * len(inames), 2.6 * len(samples)),
                                 squeeze=False)
        for i, s in enumerate(samples):
            m = (c2l.obs["sample"] == s).to_numpy()
            sf = c2l.uns["spatial"][s]["scalefactors"]["tissue_lowres_scalef"]
            img = c2l.uns["spatial"][s]["images"].get("lowres")
            xy = c2l.obsm["spatial"][m] * sf
            for j, nm in enumerate(inames):
                ax = axes[i, j]
                col = f"{nm}|celltype"
                if img is not None:
                    ax.imshow(img, alpha=0.45)
                if col in spot_scores:
                    v = spot_scores.loc[m, col].to_numpy()
                    ax.scatter(xy[:, 0], xy[:, 1], c=v, s=3, cmap="viridis",
                               vmax=np.nanquantile(spot_scores[col], 0.99) + 1e-9)
                ax.set_axis_off()
                if i == 0: ax.set_title(nm, fontsize=8)
                if j == 0: ax.text(-0.05, 0.5, s, transform=ax.transAxes, rotation=90, va="center", ha="right")
        fig.tight_layout(); fig.savefig(OUT / "lr_spatial_maps.png", dpi=120); plt.close(fig)
    print(f"\nOutputs in {OUT.resolve()}")


if __name__ == "__main__":
    main()
