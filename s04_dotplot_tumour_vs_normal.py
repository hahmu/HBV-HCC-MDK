import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import wilcoxon, mannwhitneyu
import figs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import config as C
from s03_interactions import neighbour_matrix, gene_vec, resolve_interactions, bh

OUT = C.OUT_DIR / "dotplot"
MIN_FRAC = 0.05           
RICH_Q = 0.75             
CMAP = "Purples"          
GROUPS = ["Tumour", "Normal"]


def compute(c2l, X, vidx):
    names = list(c2l.uns["mod"]["factor_names"])
    props = pd.DataFrame(c2l.obs[[f"{n}_prop" for n in names]].to_numpy(), index=c2l.obs_names, columns=names)
    abund = c2l.obs[names].copy()
    inter = [it for it in resolve_interactions(names, vidx, props, abund) if it["ct_ok"]]

    labels = c2l.obs.groupby("sample", observed=True)["region_label"].first()
    have_sections = {"Tumour", "Non-tumour"} <= set(labels)
    mode = C.COMPARE_BY if C.COMPARE_BY != "auto" else ("section" if have_sections else "spot")
    print(f"Comparison mode: {mode}  (section labels: {labels.value_counts().to_dict()})")

    rows = []
    for s in c2l.obs["sample"].unique():
        idx = np.where(c2l.obs["sample"] == s)[0]
        obs = c2l.obs.iloc[idx]
        W = neighbour_matrix(obs)
        Xs = X[idx]
        P = props.iloc[idx]
        patient = obs["patient"].iloc[0] if "patient" in obs else s
        if mode == "section":
            g = {"Tumour": "Tumour", "Non-tumour": "Normal"}.get(labels[s])
            masks = {g: np.ones(len(idx), bool)} if g else {}
        else:
            reg = obs["region_auto"].to_numpy()
            masks = {"Tumour": reg == "Tumour", "Normal": reg == "Non-tumour"}
        for it in inter:
            L, _ = gene_vec(Xs, vidx, it["ligand"])
            R, _ = gene_vec(Xs, vidx, it["receptor"])
            if L is None or R is None:
                continue
            Ln = W @ L
            snd = W @ P[it["sender"]].to_numpy()
            rcv = P[it["receiver"]].to_numpy()
            lr = np.sqrt(Ln * R)
            w = lr * np.sqrt(snd * rcv)
            active = (Ln > 0) & (R > 0) & (snd >= MIN_FRAC) & (rcv >= MIN_FRAC)
            snd_rich = snd >= np.quantile(snd, RICH_Q)
            rcv_rich = rcv >= np.quantile(rcv, RICH_Q)
            for grp, mk in masks.items():
                if mk.sum() < 20:
                    continue
                rows.append(dict(sample=s, patient=patient, group=grp, interaction=it["name"],
                                 n_spots=int(mk.sum()),
                                 weighted_score=float(w[mk].mean()), lr_score=float(lr[mk].mean()),
                                 pct_active=float(active[mk].mean() * 100),
                                 ligand_in_sender_rich=float(L[mk & snd_rich].mean()) if (mk & snd_rich).any() else np.nan,
                                 receptor_in_receiver_rich=float(R[mk & rcv_rich].mean()) if (mk & rcv_rich).any() else np.nan,
                                 sender_frac=float(snd[mk].mean()), receiver_frac=float(rcv[mk].mean()),
                                 sender=it["sender"], receiver=it["receiver"],
                                 ligand=it["ligand"], receptor=it["receptor"]))
    return pd.DataFrame(rows), [it["name"] for it in inter], mode


def stats(df, order):
    out = []
    for nm in order:
        d = df[df.interaction == nm]
        t = d[d.group == "Tumour"].groupby("patient")["weighted_score"].mean()
        n = d[d.group == "Normal"].groupby("patient")["weighted_score"].mean()
        pairs = t.index.intersection(n.index)
        ts = d[d.group == "Tumour"].set_index("sample")["weighted_score"]
        ns = d[d.group == "Normal"].set_index("sample")["weighted_score"]
        sp_ = ts.index.intersection(ns.index)
        sec_up, sec_n = int((ts[sp_] > ns[sp_]).sum()), len(sp_)
        p, test = np.nan, "n/a"
        if len(pairs) >= 3:
            diff = t[pairs] - n[pairs]
            p = wilcoxon(diff).pvalue if np.any(diff != 0) else 1.0
            test = "paired Wilcoxon"
        elif len(t) >= 2 and len(n) >= 2:
            p = mannwhitneyu(t, n).pvalue
            test = "Mann-Whitney"
        tm, nmn = t.mean(), n.mean()
        out.append(dict(interaction=nm, tumour_mean=tm, normal_mean=nmn,
                        log2FC=np.log2((tm + 1e-6) / (nmn + 1e-6)),
                        n_pairs=len(pairs), n_tumour=len(t), n_normal=len(n),
                        tumour_higher_in_pairs=int((t[pairs] > n[pairs]).sum()) if len(pairs) else np.nan,
                        sections_T_gt_N=sec_up, sections_paired=sec_n,
                        test=test, p=p))
    st = pd.DataFrame(out)
    st["q"] = bh(st["p"])
    return st


def p_label(p):
    if np.isnan(p):
        return "n.s."
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def size_of(pct):
    return 15 + 3.0 * np.asarray(pct)       # 0% -> 15, 100% -> 315 pt^2


def dot_panel(ax, mat_col, mat_size, rows, cols, vmax_by_row=True):
    col = mat_col.to_numpy(float)
    if vmax_by_row:
        col = col / (np.nanmax(col, axis=1, keepdims=True) + 1e-12)
    yy, xx = np.meshgrid(range(len(rows)), range(len(cols)), indexing="ij")
    sc_ = ax.scatter(xx.ravel(), yy.ravel(), s=size_of(mat_size.to_numpy(float).ravel()),
                     c=col.ravel(), cmap=CMAP, vmin=0, vmax=1, edgecolors="#555555", linewidths=0.5)
    ax.set_xticks(range(len(cols)), cols)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.grid(color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.tick_params(length=0)
    return sc_


def size_legend(fig, anchor):
    hs = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(size_of(v)), markerfacecolor="#bdbdbd",
                 markeredgecolor="#555555", markeredgewidth=0.5, label=f"{v}%") for v in (10, 30, 60, 90)]
    fig.legend(handles=hs, title="% active spots", loc="upper left", bbox_to_anchor=anchor,
               frameon=False, labelspacing=1.2, fontsize=8, title_fontsize=8)


def ylabel_rows(df, order):
    first = df.drop_duplicates("interaction").set_index("interaction")
    return [f"{nm}\n{first.loc[nm,'sender'].replace('+','/')} → {first.loc[nm,'receiver'].replace('+','/')}"
            for nm in order]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c2l = sc.read_h5ad(C.OUT_DIR / "visium_c2l.h5ad")
    expr = sc.read_h5ad(C.OUT_DIR / "visium_raw.h5ad")[c2l.obs_names].copy()
    sc.pp.normalize_total(expr, target_sum=1e4)
    sc.pp.log1p(expr)
    X = expr.X.tocsr() if sp.issparse(expr.X) else expr.X

    df, order, mode = compute(c2l, X, expr.var_names)
    if df.empty:
        raise SystemExit("No interaction could be scored - check the label/gene warnings above.")
    df.to_csv(OUT / "section_level_scores.csv", index=False)
    st = stats(df, order)
    st.to_csv(OUT / "tumour_vs_normal_stats.csv", index=False)
    print(st.round(4).to_string(index=False))
    rows = ylabel_rows(df, order)
    normal_name = "Non-tumour sections" if mode == "section" else "Non-tumour spots"

    # ---------------- 1. pooled tumour vs normal
    g = df.groupby(["interaction", "group"])
    col = g["weighted_score"].mean().unstack().reindex(index=order, columns=GROUPS)
    siz = g["pct_active"].mean().unstack().reindex(index=order, columns=GROUPS)
    st_i = st.set_index("interaction")
    h = 0.8 * len(order) + 1.9
    fig = plt.figure(figsize=(10.5, h))
    gs = fig.add_gridspec(1, 4, width_ratios=[2.4, 1.0, 2.2, 1.3], wspace=0.05)
    ax_lab, ax, ax_txt, ax_leg = (fig.add_subplot(gs[0, i]) for i in range(4))
    s_ = dot_panel(ax, col, siz, rows, GROUPS)
    ax.set_yticklabels([])
    ax.set_xlim(-0.7, 1.7)
    # row labels in their own panel so long sender -> receiver names never collide
    for a_ in (ax_lab, ax_txt, ax_leg):
        a_.set_ylim(len(order) - 0.5, -0.5); a_.set_axis_off()
    for i, r in enumerate(rows):
        nm, cells = r.split("\n")
        ax_lab.text(1, i - 0.13, nm, ha="right", va="center", fontsize=10, fontweight="bold")
        ax_lab.text(1, i + 0.2, cells, ha="right", va="center", fontsize=8, color="#555555")
    # stats columns
    ax_txt.set_xlim(0, 1)
    for x, t in ((0.04, "log2FC T/N"), (0.40, "p (q)"), (0.78, "T>N patients\n(sections)")):
        ax_txt.text(x, 1.01, t, fontsize=8, color="#666666", va="bottom", transform=ax_txt.transAxes)
    for i, nm in enumerate(order):
        r = st_i.loc[nm]
        ptxt = "n/a" if np.isnan(r.p) else f"{r.p:.2g} ({r.q:.2g}) {p_label(r.p)}"
        pairs = f"{int(r.tumour_higher_in_pairs)}/{int(r.n_pairs)}" if r.n_pairs else f"{r.n_tumour}T/{r.n_normal}N"
        if r.sections_paired:
            pairs += f" ({int(r.sections_T_gt_N)}/{int(r.sections_paired)})"
        ax_txt.text(0.04, i, f"{r.log2FC:+.2f}", fontsize=9, va="center", fontweight="bold",
                    color="#7a1f1f" if r.log2FC > 0 else "#1f3f7a")
        ax_txt.text(0.40, i, ptxt, fontsize=9, va="center")
        ax_txt.text(0.78, i, pairs, fontsize=9, va="center")
    # legends: colour bar + size key
    cax = ax_leg.inset_axes([0.08, 0.55, 0.12, 0.40])
    cb = fig.colorbar(s_, cax=cax)
    cb.set_label("weighted score\n(row-scaled)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for k, v in enumerate((10, 30, 60, 90)):
        y = 0.08 + k * 0.1
        ax_leg.scatter(0.15, y, s=size_of(v), color="#bdbdbd", edgecolors="#555555", lw=0.5,
                       transform=ax_leg.transAxes, clip_on=False)
        ax_leg.text(0.35, y, f"{v}%", fontsize=8, va="center", transform=ax_leg.transAxes)
    ax_leg.text(0.02, 0.48, "% active spots", fontsize=8, transform=ax_leg.transAxes)
    n_min = int(st.n_pairs.max()) if st.n_pairs.max() > 0 else 0
    floor = f" With {n_min} pairs the smallest possible paired-Wilcoxon p is {2 / 2 ** n_min:.3g}." if n_min else ""
    fig.suptitle(f"Spatial interaction activity: tumour vs normal ({'paired sections' if mode == 'section' else 'spots within sections'})",
                 x=0.02, ha="left", fontsize=11)
    fig.text(0.02, 0.01, f"Normal = {normal_name}. Colour: mean cell-type-weighted LR score, scaled to each row's max. "
             f"Stats on per-section means ({st.test.iloc[0]}, BH-adjusted q).\n{floor.strip()}",
             fontsize=7, color="#666666")
    fig.subplots_adjust(left=0.01, right=0.99, top=1 - 0.75 / h, bottom=0.75 / h)
    figs.save(fig, OUT / "dotplot_tumour_vs_normal")

    # ---------------- 2. per patient
    pats = sorted(df.patient.unique())
    cols = [f"{p} {gr[0]}" for p in pats for gr in GROUPS]
    piv_c = df.pivot_table(index="interaction", columns=["patient", "group"], values="weighted_score")
    piv_s = df.pivot_table(index="interaction", columns=["patient", "group"], values="pct_active")
    idx = pd.MultiIndex.from_tuples([(p, gr) for p in pats for gr in GROUPS])
    piv_c = piv_c.reindex(index=order, columns=idx)
    piv_s = piv_s.reindex(index=order, columns=idx).fillna(0)
    fig, ax = plt.subplots(figsize=(0.55 * len(cols) + 4.5, 0.75 * len(order) + 1.8))
    s_ = dot_panel(ax, piv_c, piv_s, rows, cols)
    ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8)
    for k in range(1, len(pats)):
        ax.axvline(2 * k - 0.5, color="#bbbbbb", lw=0.8)
    cb = fig.colorbar(s_, ax=ax, shrink=0.45, pad=0.02)
    cb.set_label("mean weighted score\n(scaled to row max)", fontsize=8)
    ax.set_title("Per patient (T = tumour, N = normal); empty = section/region absent", fontsize=10, loc="left")
    size_legend(fig, (0.88, 0.98))
    fig.tight_layout(rect=(0, 0, 0.87, 1))
    figs.save(fig, OUT / "dotplot_per_patient")

    # ---------------- 3. components: is the change driven by ligand, receptor or cell abundance?
    comp = []
    for nm in order:
        d = df[df.interaction == nm]
        lig, rec = d.ligand.iloc[0], d.receptor.iloc[0]
        for kind, colname, lab in (("ligand", "ligand_in_sender_rich", f"{lig} (sender-rich)"),
                                   ("receptor", "receptor_in_receiver_rich", f"{rec} (receiver-rich)"),
                                   ("sender", "sender_frac", f"{d.sender.iloc[0]} fraction"),
                                   ("receiver", "receiver_frac", f"{d.receiver.iloc[0]} fraction")):
            v = d.groupby("group")[colname].mean().reindex(GROUPS)
            comp.append(dict(interaction=nm, row=f"{nm} · {lab}".replace("+", "/"), Tumour=v["Tumour"], Normal=v["Normal"]))
    cp = pd.DataFrame(comp).set_index("row")
    m = cp[GROUPS].to_numpy(float)
    l2 = np.log2((m[:, 0] + 1e-3) / (m[:, 1] + 1e-3))
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(cp) + 1.4))
    lim = np.nanmax(np.abs(l2)) if np.isfinite(l2).any() else 1
    im = ax.imshow(l2[:, None], cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_yticks(range(len(cp)), cp.index, fontsize=7)
    ax.set_xticks([0], ["log2 Tumour/Normal"])
    for i, v in enumerate(l2):
        ax.text(0, i, f"{v:+.2f}", ha="center", va="center", fontsize=7,
                color="white" if abs(v) > 0.6 * lim else "#222222")
    for k in range(4, len(cp), 4):
        ax.axhline(k - 0.5, color="white", lw=2)
    fig.colorbar(im, ax=ax, shrink=0.5).set_label("log2FC", fontsize=8)
    ax.set_title("What drives each axis: ligand / receptor expression vs cell abundance", fontsize=9, loc="left")
    fig.tight_layout()
    figs.save(fig, OUT / "dotplot_components")
    print(f"\nOutputs in {OUT.resolve()}")


if __name__ == "__main__":
    main()
