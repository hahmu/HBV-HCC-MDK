import os, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge, FancyArrowPatch
from matplotlib.lines import Line2D

BASE    = os.environ.get("HBV_BASE", "./HBV_HCC")
LDIR    = f"{BASE}/results/liana"
OUT_DIR = f"{LDIR}/figures_story"
os.makedirs(OUT_DIR, exist_ok=True)

MAG, SPEC = 0.01, 0.05
EXCLUDE_GENES = {"B2M", "ALB", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3",
                 "APOE", "APOH", "HP", "HPX", "TF", "TTR", "FGA", "FGB", "FGG", "SERPINA1",
                 "AHSG", "ORM1", "ORM2", "RBP4", "AMBP", "CP", "VTN", "FTL", "FTH1"}
KT_CMAP = LinearSegmentedColormap.from_list("ktplot", ["#104E8B", "#FFDAB9", "#8B0A50"])  # dodgerblue4-peachpuff-deeppink4

senders, major receivers)
CD8   = ["CD8 Tem", "CD8 Trm", "CD8 Temra"]
TAM   = ["SPP1+ TAM", "FOLR2+ TAM", "Kupffer"]
MONO  = ["Classical monocyte", "Non-classical monocyte"]
HSC   = ["cyHSC", "Mesenchymal-like HSC", "myHSC"]
EC    = ["Tumour EC", "Arterial EC", "Venous EC", "LSEC", "Proliferating EC"]
TNK   = CD8 + ["CD4 T", "Treg", "MAIT", "gdT", "NK", "Liver-resident NK"]
NKR   = ["NK", "Liver-resident NK"] + CD8
CCR7R = ["mregDC", "Treg", "B cell", "CD4 T", "Tfh-like CD4 T"]

AXES = [
    ("SPP1-CD44",       ["SPP1"],   ["CD44"],   ["SPP1+ TAM"],        TNK + MONO + ["mregDC", "cDC2"], ["Myeloid"],     ["T_cell", "Myeloid", "NK_cell"]),
    ("CXCL12-CXCR4",    ["CXCL12"], ["CXCR4"],  EC + HSC,             TNK + TAM + MONO, ["Endothelial", "Fibroblast_HSC"], ["T_cell", "NK_cell", "Myeloid"]),
    ("HLA-I-CD8",       ["HLA-A", "HLA-B", "HLA-C"], ["CD8A", "CD8B"], CD8 + ["Malignant", "Tumour EC"], CD8,
                                                                       ["T_cell", "Endothelial", "Malignant"], ["T_cell"]),
    ("CCL21-CCR7",      ["CCL21"],  ["CCR7"],   ["Lymphatic EC"],     CCR7R,          ["Endothelial"],                   ["T_cell", "B_cell", "Myeloid"]),
    ("MDK-SDC2",        ["MDK"],    ["SDC2"],   ["Malignant"],        HSC + ["vSMC"], ["Malignant"],                     ["Fibroblast_HSC"]),
    ("MDK-LRP1",        ["MDK"],    ["LRP1"],   ["Malignant"],        ["myHSC", "Mesenchymal-like HSC"] + TAM, ["Malignant"], ["Fibroblast_HSC", "Myeloid"]),
    ("MDK-SDC3",        ["MDK"],    ["SDC3"],   ["Malignant"],        TAM + ["Proliferating myeloid"], ["Malignant"],   ["Myeloid"]),
]
AXIS_COLOURS = dict(zip([a[0] for a in AXES],
                        ["#00A087", "#7E6148", "#3C5488", "#4DBBD5", "#E64B35", "#F39B7F", "#DC0000"]))
FONT = 8
plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": FONT})


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{OUT_DIR}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def load(cond):
    f = f"{LDIR}/liana_{cond}.csv"
    if not os.path.exists(f):
        log(f"{f} not found - {cond} skipped"); return None
    d = pd.read_csv(f)
    bad = lambda c: any(g in EXCLUDE_GENES for g in str(c).split("_"))
    return d[~d["ligand_complex"].map(bad) & ~d["receptor_complex"].map(bad)].copy()


def significant(d):
    return d[(d["magnitude_rank"] <= MAG) & (d["specificity_rank"] <= SPEC)]


# ------------------------------------------------------------------ 1. count heatmap (ktplot style)
def tested_types(d):
    return set(d["source"]) | set(d["target"])


def heatmap_counts(d, cond, order=None, vmax=None, all_types=None):
    """all_types: union of labels over both tissues; labels not tested in this tissue are drawn grey."""
    s = significant(d)
    s.to_csv(f"{OUT_DIR}/lr_{cond}_significant_filtered.csv", index=False)
    types = sorted(tested_types(d))
    m = s.groupby(["source", "target"]).size().unstack(fill_value=0).reindex(index=types, columns=types, fill_value=0)
    m.to_csv(f"{OUT_DIR}/counts_{cond}.csv")
    if order is None:
        sym = np.log1p(m.values + m.values.T)
        order = [types[i] for i in leaves_list(linkage(sym, "average", metric="euclidean"))] if len(types) > 2 else types
    full = list(all_types) if all_types is not None else types
    order = [t for t in order if t in full] + sorted(t for t in full if t not in order)
    m = m.reindex(index=order, columns=order).astype(float)          # NaN = not tested in this tissue
    n = len(order)
    fig = plt.figure(figsize=(0.26 * n + 3.2, 0.26 * n + 2.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[n, 2.2], height_ratios=[2.2, n], wspace=0.03, hspace=0.03)
    ax = fig.add_subplot(gs[1, 0]); axt = fig.add_subplot(gs[0, 0], sharex=ax); axr = fig.add_subplot(gs[1, 1], sharey=ax)
    cmap = KT_CMAP.copy(); cmap.set_bad("#D9D9D9")
    im = ax.imshow(np.ma.masked_invalid(m.values), cmap=cmap, aspect="auto", vmin=0, vmax=vmax)
    untested = [t for t in order if t not in types]
    ax.set_xticks(range(n), order, rotation=90); ax.set_yticks(range(n), order)
    ax.set_xlabel("Receiver (target)"); ax.set_ylabel("Sender (source)")
    axt.bar(range(n), np.nansum(m.values, 0), color="#8B0A50", width=0.8); axt.set_ylabel("received", fontsize=FONT - 1)
    axr.barh(range(n), np.nansum(m.values, 1), color="#104E8B", height=0.8); axr.set_xlabel("sent", fontsize=FONT - 1)
    for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        if lab.get_text() in untested:
            lab.set_color("#999999")
    for a in (axt, axr):
        for s_ in ("top", "right"):
            a.spines[s_].set_visible(False)
        a.tick_params(labelsize=FONT - 2)
    plt.setp(axt.get_xticklabels(), visible=False); plt.setp(axr.get_yticklabels(), visible=False)
    cax = fig.add_axes([0.92, 0.82, 0.015, 0.12])
    fig.colorbar(im, cax=cax).set_label("interactions", fontsize=FONT - 1)
    fig.suptitle(f"{cond}: significant interactions (magnitude rank ≤ {MAG}, specificity rank ≤ {SPEC}, "
                 f"ambient genes removed; total {len(s)})"
                 + (f"\ngrey = not tested in this tissue (< min_cells): {', '.join(untested)}" if untested else ""),
                 x=0.02, ha="left", fontsize=FONT + 1)
    save(fig, f"heatmap_counts_{cond}")
    log(f"{cond}: {len(s)} significant interactions -> heatmap_counts_{cond}")
    return order


# ------------------------------------------------------------------ dot plot helper (ktplot / CellPhoneDB style)
def lr_dotplot(df, rows, cols, title, name, col_groups=None):
    """df needs columns row, col, lr_means, specificity_rank, magnitude_rank."""
    if df.empty:
        log(f"{name}: nothing to plot"); return
    df = df.assign(**{v: pd.to_numeric(df[v], errors="coerce") for v in ("lr_means", "specificity_rank", "magnitude_rank")})
    piv = lambda v: (df.pivot_table(index="row", columns="col", values=v, aggfunc="first")
                     .reindex(index=rows, columns=cols).astype(float))
    mean, spec, mag = piv("lr_means"), piv("specificity_rank"), piv("magnitude_rank")
    yy, xx = np.meshgrid(range(len(rows)), range(len(cols)), indexing="ij")
    ok = ~np.isnan(mean.values)
    size = 8 + 22 * np.clip(-np.log10(np.clip(spec.values, 1e-6, 1)), 0, 6)
    fig, ax = plt.subplots(figsize=(0.24 * len(cols) + 3.5, 0.24 * len(rows) + 2.2))
    sc_ = ax.scatter(xx[ok], yy[ok], s=size[ok], c=mean.values[ok], cmap=KT_CMAP,
                     edgecolors=np.where(mag.values[ok] <= MAG, "black", "none"), linewidths=0.7)
    ax.set_xticks(range(len(cols)), cols, rotation=90, fontsize=FONT - 1)
    ax.set_yticks(range(len(rows)), rows, fontsize=FONT - 1)
    ax.set_xlim(-0.7, len(cols) - 0.3); ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.grid(color="#eeeeee", lw=0.5); ax.set_axisbelow(True)
    if col_groups:
        for x in col_groups:
            ax.axvline(x - 0.5, color="#999999", lw=0.8)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cax = ax.inset_axes([1.02, 0.62, 0.015, 0.36])
    fig.colorbar(sc_, cax=cax).set_label("mean L-R expression", fontsize=FONT - 1)
    hs = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(8 + 22 * k), color="#999999", label=f"1e-{k}")
          for k in (1, 3, 5)] + [Line2D([], [], ls="", marker="o", markersize=7, markerfacecolor="white",
                                        markeredgecolor="black", label=f"magnitude ≤ {MAG}")]
    ax.legend(handles=hs, title="specificity rank", loc="upper left", bbox_to_anchor=(1.0, 0.52),
              frameon=False, fontsize=FONT - 1, title_fontsize=FONT - 1, labelspacing=1.1)
    ax.set_title(title, fontsize=FONT + 1, loc="left")
    save(fig, name)


# ------------------------------------------------------------------ 2. overview dot plot
def dotplot_top(d, cond, senders=("myHSC", "SPP1+ TAM", "FOLR2+ TAM", "Kupffer", "Tumour EC",
                                   "Mesenchymal-like HSC", "Malignant", "CD8 Trm", "Treg"), per_sender=5):
    s = significant(d)
    s = s[s["source"].isin(senders)]
    if s.empty:
        return
    s = s.assign(lr=s["ligand_complex"] + " → " + s["receptor_complex"])
    top = (s.sort_values("magnitude_rank").drop_duplicates(["source", "lr"])
           .groupby("source").head(per_sender))
    keep = s[s["lr"].isin(top["lr"]) & s["source"].isin(top["source"])]
    keep = keep.sort_values("magnitude_rank").groupby(["source", "lr"]).head(4)   # 4 best receivers per LR
    keep = keep.assign(row=keep["lr"], col=keep["source"] + " → " + keep["target"])
    rows = list(dict.fromkeys(top.sort_values(["source"])["lr"]))
    cols = list(dict.fromkeys(keep.sort_values(["source", "target"])["col"]))
    lr_dotplot(keep, rows, cols, f"{cond}: top {per_sender} specific interactions of key senders",
               f"dotplot_top_{cond.lower()}")


# ------------------------------------------------------------------ 
def axis_rows(d, cond):
    out = []
    for name, ligs, recs, snd, rcv, *_ in AXES:
        x = d[d["ligand_complex"].isin(ligs) & d["receptor_complex"].isin(recs) &
              d["source"].isin(snd) & d["target"].isin(rcv)].copy()
        x["axis"], x["condition"] = name, cond
        out.append(x)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def dotplot_axes(dt, dn):
    ax_t = axis_rows(dt, "Tumor")
    ax_n = axis_rows(dn, "Normal") if dn is not None else pd.DataFrame(columns=ax_t.columns)
    both = pd.concat([ax_t, ax_n], ignore_index=True)
    both.to_csv(f"{OUT_DIR}/axes_subtype_Tumor_Normal.csv", index=False)
    # columns: sender -> receiver combinations that are significant in tumour for that axis (top 6 per axis)
    sig_t = ax_t.sort_values("magnitude_rank")          # top pairs per axis, weak ones included
    pairs = (sig_t.assign(pair=sig_t["source"] + " → " + sig_t["target"])
             .drop_duplicates(["axis", "pair"]).groupby("axis").head(6))
    cols, groups = [], []
    for cond in ("Tumor", "Normal"):
        groups.append(len(cols))
        cols += [f"{p} | {cond[0]}" for p in dict.fromkeys(pairs["pair"])]
    both = both.assign(pair=both["source"] + " → " + both["target"],
                       lr=both["ligand_complex"] + " → " + both["receptor_complex"])
    both = both[both["pair"].isin(pairs["pair"])]
    both["col"] = both["pair"] + " | " + both["condition"].str[0]
    both["row"] = both["axis"] + ":  " + both["lr"]
    rows = [r for a in [x[0] for x in AXES] for r in dict.fromkeys(both.loc[both["axis"] == a, "row"])]
    lr_dotplot(both, rows, cols, "Selected axes at subtype level: T = tumour, N = adjacent normal "
               "(empty = cell types or genes below LIANA thresholds)", "dotplot_axes_Tumor_vs_Normal",
               col_groups=groups[1:])
    return ax_t, ax_n


# ------------------------------------------------------------------
def chord(ax_df, cond, order):
    if ax_df.empty:
        log(f"chord {cond}: no axis edges"); return
    e = (ax_df.assign(w=-np.log10(ax_df["magnitude_rank"].clip(lower=1e-10)))
         .sort_values("w", ascending=False).drop_duplicates(["axis", "source", "target"])
         .groupby("axis").head(6))
    used = set(e["source"]) | set(e["target"])
    nodes = [t for t in order if t in used] + sorted(used - set(order))   # types absent from the tumour order
    n = len(nodes)
    ang = {t: 2 * np.pi * i / n for i, t in enumerate(nodes)}
    pal = plt.get_cmap("tab20")
    ncol = {t: pal(i % 20) for i, t in enumerate(nodes)}
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    half = np.degrees(np.pi / n) * 0.85
    for t, a in ang.items():
        ax.add_patch(Wedge((0, 0), 1.0, np.degrees(a) - half, np.degrees(a) + half, width=0.06, color=ncol[t]))
        r, deg = 1.1, np.degrees(a)
        ax.text(r * np.cos(a), r * np.sin(a), t, rotation=deg if np.cos(a) >= 0 else deg + 180,
                rotation_mode="anchor", ha="left" if np.cos(a) >= 0 else "right", va="center", fontsize=FONT)
    wmax = e["w"].max()
    for _, r in e.iterrows():
        a1, a2 = ang[r["source"]], ang[r["target"]]
        p1 = 0.93 * np.array([np.cos(a1), np.sin(a1)]); p2 = 0.93 * np.array([np.cos(a2), np.sin(a2)])
        lw = 0.6 + 3.4 * r["w"] / wmax
        if r["source"] == r["target"]:
            c = 0.6 * p1; perp = np.array([-p1[1], p1[0]]) * 0.18
            path = Path([p1, c + perp, c - perp, p1], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            ax.add_patch(PathPatch(path, fc="none", ec=AXIS_COLOURS[r["axis"]], lw=lw, alpha=0.75))
        else:
            ax.add_patch(FancyArrowPatch(p1, p2, connectionstyle="arc3,rad=0.25", arrowstyle="-|>",
                                         mutation_scale=6 + 4 * r["w"] / wmax, lw=lw,
                                         ls="-" if r["magnitude_rank"] <= 0.05 else (0, (3, 2)),
                                         color=AXIS_COLOURS[r["axis"]], alpha=0.75))
    hs = [Line2D([], [], color=c, lw=3, label=k) for k, c in AXIS_COLOURS.items() if k in set(e["axis"])]
    ax.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=False, fontsize=FONT)
    ax.set_title(f"{cond}: selected axes, top 6 sender-receiver pairs each\n"
                 f"edge width = -log10 magnitude rank; dashed = weak (rank > 0.05)",
                 fontsize=FONT + 1)
    save(fig, f"chord_axes_{cond}")


# ------------------------------------------------------------------ 
def dotplot_de():
    f = f"{LDIR}/de_informed/lr_de_informed_all.csv"
    if not os.path.exists(f):
        log("DE-informed table not found - skipped"); return
    de = pd.read_csv(f)
    out = []
    for name, ligs, recs, _, _, msnd, mrcv in AXES:
        recs = [r for r in recs if r != "CD8B"]
        x = de[de["ligand_complex"].isin(ligs) & de["receptor_complex"].isin(recs) &
               de["source"].isin(msnd) & de["target"].isin(mrcv)].copy()
        out.append(x.assign(axis=name))
    d = pd.concat(out, ignore_index=True)
    d.to_csv(f"{OUT_DIR}/axes_DE_informed.csv", index=False)
    if d.empty:
        return
    d["label"] = (d["axis"] + ":  " + d["ligand_complex"] + " → " + d["receptor_complex"] + "   "
                  + d["source"] + " → " + d["target"])
    st = d[["ligand_stat", "receptor_stat"]].to_numpy(float)
    pa = d[["ligand_padj", "receptor_padj"]].to_numpy(float)
    size = 10 + 16 * np.clip(-np.log10(np.clip(np.nan_to_num(pa, nan=1), 1e-300, 1)), 0, 10)
    v = np.nanmax(np.abs(st))
    fig, ax = plt.subplots(figsize=(6.8, 0.22 * len(d) + 1.4))
    yy, xx = np.meshgrid(range(len(d)), [0, 1], indexing="ij")
    s_ = ax.scatter(xx.ravel(), yy.ravel(), c=st.ravel(), s=size.ravel(), cmap="RdBu_r", vmin=-v, vmax=v,
                    edgecolors=np.where(pa.ravel() < 0.05, "black", "none"), linewidths=0.7)
    ax.set_xticks([0, 1], ["ligand\n(sender)", "receptor\n(receiver)"])
    ax.set_yticks(range(len(d)), d["label"], fontsize=FONT - 1)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(len(d) - 0.5, -0.5)
    for s__ in ax.spines.values():
        s__.set_visible(False)
    fig.colorbar(s_, ax=ax, shrink=0.35).set_label("DESeq2 stat, Tumor vs Normal", fontsize=FONT - 1)
    ax.set_title("Selected axes in the paired pseudobulk test (major types)\n"
                 "dot size = -log10 padj; black outline = padj < 0.05", fontsize=FONT + 1, loc="left")
    save(fig, "dotplot_axes_DE")


# ------------------------------------------------------------------ 
LINEAGE_COLOURS = {"Malignant": "#B2182B", "Epithelial": "#E08214", "Endothelial": "#4DBBD5",
                   "Stromal": "#00A087", "Myeloid": "#7E6148", "T/NK": "#3C5488", "B/Plasma": "#8491B4",
                   "Other": "#999999"}


def lineage(t):
    t_ = t.lower()
    if t == "Malignant":
        return "Malignant"
    if t_ in ("hepatocyte", "cholangiocyte"):
        return "Epithelial"
    if " ec" in f" {t_}" or t_ == "lsec" or "endothel" in t_:
        return "Endothelial"
    if "hsc" in t_ or t_ in ("vsmc",) or "fibro" in t_:
        return "Stromal"
    if any(k in t_ for k in ("tam", "kupffer", "monocyte", "dc", "myeloid", "mast", "macrophage")):
        return "Myeloid"
    if any(k in t_ for k in ("b cell", "plasma")):
        return "B/Plasma"
    if any(k in t_ for k in ("cd4", "cd8", "treg", "tfh", "mait", "gdt", "nk", "ilc", "t/nk", "proliferating t")):
        return "T/NK"
    return "Other"


def chord_top50(d, cond, n_top=50, n_label=10):
    s = significant(d).sort_values(["magnitude_rank", "specificity_rank"]).head(n_top).reset_index(drop=True)
    if s.empty:
        return
    floor = s.loc[s["magnitude_rank"] > 0, "magnitude_rank"].min()
    floor = floor / 10 if pd.notna(floor) else 1e-10
    s["w"] = -np.log10(s["magnitude_rank"].clip(lower=floor))
    s["lineage"] = s["source"].map(lineage)
    s["n"] = np.arange(1, len(s) + 1)
    s.to_csv(f"{OUT_DIR}/chord_top{n_top}_{cond}.csv", index=False)
    # nodes grouped by lineage, arc length = number of chord ends
    ends = pd.concat([s["source"], s["target"]]).value_counts()
    lin_order = list(LINEAGE_COLOURS)
    nodes = sorted(ends.index, key=lambda t: (lin_order.index(lineage(t)), t))
    gap = 0.02 * 2 * np.pi
    unit = (2 * np.pi - gap * len(nodes)) / ends.sum()
    start, cursor, a0 = {}, {}, np.pi / 2
    for t in nodes:
        start[t] = a0; cursor[t] = a0; a0 -= ends[t] * unit + gap        # clockwise
    fig, ax = plt.subplots(figsize=(8, 9.5))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.95, 1.5)
    for t in nodes:
        a_end = start[t] - ends[t] * unit
        ax.add_patch(Wedge((0, 0), 1.0, np.degrees(a_end), np.degrees(start[t]), width=0.05,
                           color=LINEAGE_COLOURS[lineage(t)]))
        am = (start[t] + a_end) / 2
        ax.text(1.07 * np.cos(am), 1.07 * np.sin(am), t, fontsize=FONT - 1, va="center",
                rotation=np.degrees(am) if np.cos(am) >= 0 else np.degrees(am) + 180, rotation_mode="anchor",
                ha="left" if np.cos(am) >= 0 else "right")
    pt = lambda a, r=0.94: r * np.array([np.cos(a), np.sin(a)])
    wmax = s["w"].max()
    for _, r in s.iloc[::-1].iterrows():                                # strongest drawn last (on top)
        a1 = cursor[r["source"]] - unit / 2; cursor[r["source"]] -= unit
        a2 = cursor[r["target"]] - unit / 2; cursor[r["target"]] -= unit
        p1, p2 = pt(a1), pt(a2)
        col = LINEAGE_COLOURS[r["lineage"]]
        lw = 0.5 + 3.5 * r["w"] / wmax
        if r["source"] == r["target"]:
            ax.add_patch(PathPatch(Path([p1, 0.55 * p1 + 0.1 * np.array([-p1[1], p1[0]]),
                                         0.55 * p1 - 0.1 * np.array([-p1[1], p1[0]]), p1],
                                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                   fc="none", ec=col, lw=lw, alpha=0.7))
            mid = 0.62 * p1
        else:
            ax.add_patch(FancyArrowPatch(p1, p2, connectionstyle="arc3,rad=0.3", arrowstyle="-|>",
                                         mutation_scale=5 + 4 * r["w"] / wmax, lw=lw, color=col, alpha=0.7))
            mid = (p1 + p2) / 2 * 0.75
        if r["n"] <= n_label:
            ax.text(*mid, str(r["n"]), fontsize=FONT - 1, ha="center", va="center", fontweight="bold",
                    bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec=col, lw=0.6))
    hs = [Line2D([], [], color=c, lw=3, label=k) for k, c in LINEAGE_COLOURS.items() if k in set(s["lineage"])]
    ax.legend(handles=hs, title="sender lineage", loc="upper left", bbox_to_anchor=(0.98, 1.0), frameon=False,
              fontsize=FONT - 1, title_fontsize=FONT - 1)
    lines = [f"{r['n']}. {r['ligand_complex']}–{r['receptor_complex']}, {r['source']} → {r['target']}"
             for _, r in s.head(n_label).iterrows()]
    half = (len(lines) + 1) // 2
    for j, block in enumerate((lines[:half], lines[half:])):
        ax.text(-1.5 + 1.55 * j, -1.42, "\n".join(block), fontsize=FONT - 1.5, va="top", ha="left")
    ax.set_title(f"{'Tumour' if cond == 'Tumor' else 'Non-tumour'}: top {n_top} interactions "
                 f"(width = −log10 magnitude rank)", fontsize=FONT + 1)
    save(fig, f"chord_top{n_top}_{cond}")


def main():
    dt, dn = load("Tumor"), load("Normal")
    if dt is None:
        raise SystemExit("liana_Tumor.csv is needed")
    all_types = sorted(tested_types(dt) | (tested_types(dn) if dn is not None else set()))
    # shared colour scale: same vmax for tumour and normal (max cell-pair count in either)
    vmax = max(significant(x).groupby(["source", "target"]).size().max() for x in (dt, dn) if x is not None)
    log(f"shared colour scale: 0-{vmax} interactions per sender-receiver pair")
    order = heatmap_counts(dt, "Tumor", vmax=vmax, all_types=all_types)
    if dn is not None:
        heatmap_counts(dn, "Normal", order=order, vmax=vmax, all_types=all_types)
    dotplot_top(dt, "Tumor")
    ax_t, ax_n = dotplot_axes(dt, dn)
    chord(ax_t, "Tumor", order)
    if dn is not None and len(ax_n):
        chord(ax_n, "Normal", order)
    chord_top50(dt, "Tumor")
    if dn is not None:
        chord_top50(dn, "Normal")
    table_s3(dt, dn, all_types)
    dotplot_de()
    log("Finished. Figures and tables in", OUT_DIR)


if __name__ == "__main__":
    main()
