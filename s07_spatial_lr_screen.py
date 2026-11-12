import os, re, time, json, urllib.request, urllib.parse
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests

ADATA      = os.environ.get("ADATA", "c2l_out/sp_with_c2l.h5ad")
OUT        = os.environ.get("OUT", "c2l_out/spatial_lr_screen")
SAMPLE_KEY = os.environ.get("SAMPLE_KEY", "sample")
ABUND_KEY  = os.environ.get("ABUND_KEY", "q05_cell_abundance_w_sf")
N_PERMS    = int(os.environ.get("N_PERMS", 500))
NZ_PROP    = float(os.environ.get("NZ_PROP", 0.05))   
DO_PUBMED  = os.environ.get("PUBMED", "1") == "1"
PUBMED_TOP = int(os.environ.get("PUBMED_TOP", 300))
EXCLUDE    = {"C004T1R"}                              
HBV_PTS    = {"B013", "C004", "C008"}
NONVIRAL   = {"A009"}
os.makedirs(OUT, exist_ok=True)


def patient_of(section):
    """A009N -> A009, B0013T1 -> B013, B013T5 -> B013, C004T1P -> C004."""
    m = re.match(r"^([A-Z])(\d+)", section)
    return f"{m.group(1)}{int(m.group(2)):03d}"


def run_section(ad):
    ad = ad.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    # bandwidth = 1.5 x median nearest-spot distance (first ring of Visium neighbours),
    # so the result does not depend on the pixel scale of obsm['spatial']
    xy = ad.obsm["spatial"]
    d, _ = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)
    bw = 1.5 * np.median(d[:, 1])
    li.ut.spatial_neighbors(ad, bandwidth=bw, cutoff=0.1, kernel="gaussian", set_diag=True)
    lr = li.mt.bivariate(
        ad, resource_name="consensus",
        local_name="cosine", global_name="morans",
        n_perms=N_PERMS, nz_prop=NZ_PROP,
        mask_negatives=False, add_categories=False,
        use_raw=False, verbose=False,
    )
    if lr is None:                       # older liana versions store it in obsm
        lr = ad.obsm["local_scores"]
    return lr


def attribution(lr, ad, top=3):
    """Spearman of per-spot local score with each cell type's abundance."""
    if ABUND_KEY not in ad.obsm:
        return {}
    ab = pd.DataFrame(ad.obsm[ABUND_KEY], index=ad.obs_names)
    ab.columns = [re.sub(r"^q05cell_abundance_w_sf_", "", str(c)) for c in ab.columns]
    ab = ab.loc[lr.obs_names]
    X = lr.X.toarray() if hasattr(lr.X, "toarray") else np.asarray(lr.X)
    ab_rank = ab.rank().values
    ab_rank = (ab_rank - ab_rank.mean(0)) / (ab_rank.std(0) + 1e-12)
    out = {}
    for j, pair in enumerate(lr.var_names):
        r = pd.Series(X[:, j]).rank().values
        r = (r - r.mean()) / (r.std() + 1e-12)
        rho = pd.Series((ab_rank * r[:, None]).mean(0), index=ab.columns)
        out[pair] = rho
    return out


# ---------------- per-section screen ----------------
adata = sc.read_h5ad(ADATA)
sections = [s for s in adata.obs[SAMPLE_KEY].astype(str).unique() if s not in EXCLUDE]
rows, rho_all = [], {}
for s in sorted(sections):
    pt = patient_of(s)
    print(f"[{s}] patient {pt} ...", flush=True)
    ad = adata[adata.obs[SAMPLE_KEY].astype(str) == s].copy()
    lr = run_section(ad)
    v = lr.var.copy()
    v["pair"] = lr.var_names
    v["section"], v["patient"] = s, pt
    v["group"] = "HBV" if pt in HBV_PTS else ("NonViral" if pt in NONVIRAL else "Other")
    rows.append(v)
    rho_all[s] = attribution(lr, ad)

per = pd.concat(rows, ignore_index=True)
pcol = "morans_pvals" if "morans_pvals" in per.columns else [c for c in per.columns if "pval" in c][0]
per.to_csv(f"{OUT}/per_section_morans.csv", index=False)

# ---------------- combine across HBV sections ----------------
p_floor = 1.0 / (N_PERMS + 1)
per["p"] = per[pcol].clip(lower=p_floor, upper=1.0)
# signed z: positive = spatial co-localisation of ligand and receptor
per["z"] = np.sign(per["morans"]) * stats.norm.isf(per["p"] / 2)
per["sig_pos"] = (per["p"] < 0.05) & (per["morans"] > 0)

hbv_secs = sorted(per.loc[per.group == "HBV", "section"].unique())
k = len(hbv_secs)
hbv = per[per.group == "HBV"]
# pairs filtered out of a section (nz_prop) contribute z = 0 (conservative)
agg = hbv.groupby("pair").agg(
    ligand=("ligand", "first"), receptor=("receptor", "first"),
    sum_z=("z", "sum"), mean_I=("morans", "mean"),
    n_sig=("sig_pos", "sum"), n_tested=("section", "nunique"),
)
agg["stouffer_z"] = agg["sum_z"] / np.sqrt(k)
agg["p_comb"] = stats.norm.sf(agg["stouffer_z"])          # one-sided: co-localised
agg["fdr"] = multipletests(agg["p_comb"], method="fdr_bh")[1]
agg["n_HBV_sections"] = k
agg["n_patients_sig"] = hbv[hbv.sig_pos].groupby("pair")["patient"].nunique()
agg["n_patients_sig"] = agg["n_patients_sig"].fillna(0).astype(int)

nv = per[per.group == "NonViral"].groupby("pair").agg(
    A009_mean_I=("morans", "mean"), A009_n_sig=("sig_pos", "sum"))
agg = agg.join(nv, how="left")

# ---------------- cell-type attribution (HBV sections, mean rho) ----------------
def top_ct(pair, n=3):
    rhos = [rho_all[s][pair] for s in hbv_secs if pair in rho_all.get(s, {})]
    if not rhos:
        return ""
    m = pd.concat(rhos, axis=1).mean(1).sort_values(ascending=False)
    return "; ".join(f"{ct} ({r:.2f})" for ct, r in m.head(n).items())

agg = agg.sort_values(["fdr", "stouffer_z"], ascending=[True, False])
sig = agg[(agg.fdr < 0.05) & (agg.n_patients_sig >= 2)].copy()
sig["top_celltypes"] = [top_ct(p) for p in sig.index]

# ---------------- PubMed novelty count ----------------
def pubmed_count(lig, rec):
    def block(sym):
        subs = [x for x in re.split(r"[_^]", sym) if x]
        return "(" + " OR ".join(f'{x}[tiab]' for x in subs) + ")"
    term = f'("hepatocellular carcinoma"[tiab] OR HCC[tiab]) AND {block(lig)} AND {block(rec)}'
    q = {"db": "pubmed", "term": term, "retmode": "json", "rettype": "count"}
    if os.environ.get("NCBI_EMAIL"):   q["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"): q["api_key"] = os.environ["NCBI_API_KEY"]
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(q)
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            return int(json.load(r)["esearchresult"]["count"])
    except Exception:
        return np.nan

if DO_PUBMED and len(sig):
    delay = 0.12 if os.environ.get("NCBI_API_KEY") else 0.4
    counts = []
    for p, row in sig.head(PUBMED_TOP).iterrows():
        counts.append(pubmed_count(row.ligand, row.receptor)); time.sleep(delay)
    sig.loc[sig.index[:len(counts)], "pubmed_HCC_hits"] = counts

agg.to_csv(f"{OUT}/HBV_combined_all_pairs.csv")
sig.to_csv(f"{OUT}/HBV_significant_pairs.csv")
print(f"\n{len(sig)} pairs FDR<0.05 and significant in >=2 HBV patients "
      f"({k} HBV sections). Output: {OUT}/")
