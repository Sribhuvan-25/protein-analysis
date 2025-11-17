# cBioPortal.py
# -*- coding: utf-8 -*-
"""
NEPC vs CRPC expression pull from cBioPortal (prad_msk_mdanderson_2023)
- Robust clinical fetch, subtype mapping, strict study-matched RNA profile
- Boxplots (Matplotlib)
- Heatmap (Matplotlib ONLY, fixed-axes layout: no seaborn, no auto-layout)
"""

import os, math, argparse, warnings
from typing import List, Tuple
import requests
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --------- Globals ---------
CBIO = "https://www.cbioportal.org/api"
STUDY_ID = "prad_msk_mdanderson_2023"
REQ_TIMEOUT = 45
CHUNK = 250
NEPC_TOKENS = ["neuroendocrine", "nepc", "small cell", "ne-like", "ne like", "scpc"]
CRPC_LABEL, NEPC_LABEL = "Adeno", "Neuroendocrine"

# Disable all figure auto-layout globally (we’ll place axes explicitly for the heatmap)
plt.rcParams['figure.constrained_layout.use'] = False
plt.rcParams['figure.autolayout'] = False

# ---------------- HTTP helpers ----------------
def _headers(accept="application/json"):
    h = {"Accept": accept}
    tok = os.getenv("CBIOPORTAL_TOKEN")
    if tok: h["X-API-TOKEN"] = tok
    return h

def _get_json(url: str, params=None):
    r = requests.get(url, headers=_headers(), params=params, timeout=REQ_TIMEOUT)
    if r.status_code != 200:
        raise requests.HTTPError(f"GET {url} ({r.status_code}): {r.text[:300]}")
    return r.json()

def _post_json(url: str, payload: dict):
    r = requests.post(url, headers=_headers(), json=payload, timeout=REQ_TIMEOUT)
    if r.status_code != 200:
        raise requests.HTTPError(f"POST {url} ({r.status_code}): {r.text[:500]}")
    return r.json()

# ---------------- discovery ----------------
def list_samples_scoped(study_id: str) -> pd.DataFrame:
    url = f"{CBIO}/studies/{study_id}/samples"
    return pd.DataFrame(_get_json(url))

def list_molecular_profiles(study_id: str) -> pd.DataFrame:
    url = f"{CBIO}/molecular-profiles"
    return pd.DataFrame(_get_json(url, params={"studyId": study_id, "pageSize": 10000}))

def choose_mrna_profile(df_prof: pd.DataFrame, study_id: str) -> Tuple[str, str]:
    """
    Strictly select an RNA expression profile BELONGING to this study.
    Prefer raw/quantified RNA-seq (FPKM/RSEM/UQ) when available, else z-score.
    Returns (molecular_profile_id, profile_type) where profile_type in {"raw_fpkm","zscore"}.
    """
    if df_prof.empty:
        raise RuntimeError("No molecular profiles found.")

    # hard filter to this study when possible
    if "studyId" in df_prof.columns:
        df_prof = df_prof[df_prof["studyId"].astype(str) == study_id].copy()

    if df_prof.empty:
        # defensive fallback if the API didn't include studyId
        urlish = study_id.lower()
        keep = []
        for _, r in df_prof.iterrows():
            pid = str(r.get("molecularProfileId", "")).lower()
            name = str(r.get("name", "")).lower()
            if urlish in pid or urlish in name:
                keep.append(r)
        df_prof = pd.DataFrame(keep)

    if df_prof.empty:
        raise RuntimeError(f"No molecular profiles matched this study: {study_id}")

    def norm(s): return str(s or "").lower()

    def is_rna(row) -> bool:
        pid   = norm(row.get("molecularProfileId"))
        mtype = norm(row.get("molecularAlterationType"))
        dtype = norm(row.get("datatype"))
        name  = norm(row.get("name"))
        blob  = " ".join([pid, mtype, dtype, name])

        # Accept typical RNA/mRNA expression profiles; exclude protein/RPPA explicitly
        if any(tok in blob for tok in ["rppa", "protein"]):
            return False
        return any(tok in blob for tok in [
            "rna", "mrna", "expression", "gex", "rna_seq", "rsem", "fpkm", "counts"
        ])

    rna = df_prof[df_prof.apply(is_rna, axis=1)].copy()
    if rna.empty:
        raise RuntimeError("No RNA-like profiles in this study.")

    def has_fpkm_like(row) -> bool:
        blob = " ".join([
            norm(row.get("molecularProfileId")),
            norm(row.get("molecularAlterationType")),
            norm(row.get("datatype")),
            norm(row.get("name")),
        ])
        return any(tok in blob for tok in ["fpkm", "rsem", "rna_seq", "uq", "counts"])

    def has_zscore(row) -> bool:
        blob = " ".join([
            norm(row.get("molecularProfileId")),
            norm(row.get("molecularAlterationType")),
            norm(row.get("datatype")),
            norm(row.get("name")),
        ])
        return any(tok in blob for tok in ["zscore", "z_score", "z-score"])

    def score(row) -> int:
        s = 0
        if has_fpkm_like(row): s += 300     # prefer raw/quant RNA
        if has_zscore(row):    s += 200     # fallback to z-score
        # minor tiebreakers
        pid = norm(row.get("molecularProfileId"))
        name = norm(row.get("name"))
        if "seq" in pid or "seq" in name: s += 20
        return s

    rna["__score__"] = rna.apply(score, axis=1)
    top = rna.sort_values("__score__", ascending=False).iloc[0]

    profile_type = "raw_fpkm" if has_fpkm_like(top) and not has_zscore(top) else "zscore"
    print(f"[INFO] RNA profile (study-matched): {top['molecularProfileId']} (type={profile_type}, score={top['__score__']})")
    return top["molecularProfileId"], profile_type

# ---------------- clinical helpers ----------------
def get_clinical_wide(study_id: str) -> pd.DataFrame:
    url = f"{CBIO}/studies/{study_id}/clinical-data"
    df = pd.DataFrame(_get_json(url, params={"clinicalDataType": "SAMPLE", "pageSize": 100000}))
    if df.empty: raise RuntimeError("Study clinical table is empty.")
    id_col = next((c for c in ["entityId","sampleId","uniqueSampleKey"] if c in df.columns), None)
    if id_col is None: raise RuntimeError("No id column in clinical table.")
    if "clinicalAttributeId" not in df.columns or "value" not in df.columns:
        raise RuntimeError("Unexpected clinical schema.")
    wide = df.pivot_table(index=id_col, columns="clinicalAttributeId", values="value", aggfunc="first")
    wide.index.name = "sampleId"
    return wide.reset_index()

def pick_best_morphology_column(df_clin: pd.DataFrame) -> str:
    cand_names = ["MORPHOLOGY","HISTOLOGY","SUBTYPE","PHENOTYPE","HISTOLOGICAL_TYPE",
                  "HISTOLOGICAL_SUBTYPE","HISTOLOGIC_DIAGNOSIS","DIAGNOSIS","TUMOR_TYPE"]
    cols_noid = [c for c in df_clin.columns if c != "sampleId"]
    best, best_hits = None, -1
    for c in cols_noid:
        s = df_clin[c].dropna().astype(str).str.lower()
        hits = s.apply(lambda x: any(tok in x for tok in NEPC_TOKENS)).sum()
        if hits > best_hits:
            best, best_hits = c, hits
    if best_hits > 0:
        print(f"[INFO] Column with NEPC-like values: {best} (hits={best_hits})")
        return best
    scored = []
    for c in cols_noid:
        up = c.upper()
        hint = any(h in up for h in cand_names)
        nn = df_clin[c].notna().sum()
        scored.append(((100 if hint else 0)+nn, c))
    scored.sort(reverse=True)
    return scored[0][1]

def map_nepc_vs_crpc(series: pd.Series) -> pd.Series:
    def f(x):
        if pd.isna(x): return np.nan
        s = str(x).strip().lower()
        if any(tok in s for tok in NEPC_TOKENS): return NEPC_LABEL
        return CRPC_LABEL
    return series.map(f)

# ---------------- expression ----------------
def get_entrez_ids(gene_symbols: List[str]) -> dict:
    url = f"{CBIO}/genes"
    out = {}
    for symbol in gene_symbols:
        try:
            data = _get_json(url, params={"keyword": symbol, "pageSize": 10})
            for item in data or []:
                if item.get("hugoGeneSymbol","").upper() == symbol.upper():
                    eid = item.get("entrezGeneId")
                    if eid and eid > 0: out[symbol] = str(eid); break
            if symbol not in out and data:
                eid = data[0].get("entrezGeneId")
                if eid and eid > 0: out[symbol] = str(eid)
        except Exception as e:
            warnings.warn(f"Entrez lookup failed for {symbol}: {e}")
    return out

def fetch_expression(mol_profile_id: str, sample_ids: List[str], genes: List[str], study_id: str = None) -> pd.DataFrame:
    if study_id is None:
        study_id = mol_profile_id.split("_")[0] if "_" in mol_profile_id else STUDY_ID
    symbol_to_entrez = get_entrez_ids(genes)
    entrez_ids = [symbol_to_entrez.get(g,"") for g in genes if symbol_to_entrez.get(g)]
    frames = []
    url_profile = f"{CBIO}/molecular-profiles/{mol_profile_id}/molecular-data/fetch"

    for i in range(0, len(sample_ids), CHUNK):
        chunk_ids = sample_ids[i:i+CHUNK]
        formatted_sample_ids = [f"{study_id}:{sid}" for sid in chunk_ids]
        strategies = [
            {"geneSymbols": genes, "sampleIds": chunk_ids},
            {"geneSymbols": genes, "sampleIds": formatted_sample_ids},
        ]
        if entrez_ids:
            strategies += [
                {"entrezGeneIds": entrez_ids, "sampleIds": chunk_ids},
                {"entrezGeneIds": entrez_ids, "sampleIds": formatted_sample_ids},
            ]
        ok = False; last_error=None
        for payload in strategies:
            try:
                res = _post_json(url_profile, payload)
                if isinstance(res, list): df = pd.DataFrame(res)
                elif isinstance(res, dict): df = pd.DataFrame(res.get("data", [res]))
                else: df = pd.DataFrame()
                if not df.empty: frames.append(df); ok=True; break
            except Exception as e:
                last_error = e; continue
        if not ok:
            raise RuntimeError(f"Expression fetch failed for chunk {i//CHUNK+1}: {last_error}")

    if not frames: raise RuntimeError("No expression rows returned.")
    df = pd.concat(frames, ignore_index=True)
    if "sampleId" in df.columns:
        df["sampleId"] = df["sampleId"].astype(str).str.split(":", n=1).str[-1]
    gcol = "geneSymbol" if "geneSymbol" in df.columns else ("geneId" if "geneId" in df.columns else None)
    if gcol is None and "entrezGeneId" in df.columns:
        # map back
        inv = {str(v):k for k,v in symbol_to_entrez.items()}
        df["geneSymbol"] = df["entrezGeneId"].astype(str).map(inv).fillna(df["entrezGeneId"].astype(str))
        gcol = "geneSymbol"
    wide = df.pivot_table(index="sampleId", columns=gcol, values="value", aggfunc="first")
    return wide.reset_index()

# ---------------- stats/plots ----------------
def welch_t(x, y) -> Tuple[float,float]:
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2: return (np.nan, np.nan)
    return stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")

def cohens_d(x, y) -> float:
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2: return np.nan
    nx, ny = len(x), len(y)
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    pooled = ((nx-1)*vx + (ny-1)*vy) / (nx+ny-2)
    return (x.mean() - y.mean()) / math.sqrt(pooled) if pooled > 0 else np.nan

def boxplot_gene(df_long: pd.DataFrame, gene: str, out_png: str,
                 use_log_scale: bool = False, profile_type: str = "zscore",
                 disable_auto_log: bool = False):
    sub = df_long[df_long["gene"]==gene].dropna(subset=["value","subtype"]).copy()
    if sub.empty:
        warnings.warn(f"No data to plot for {gene}")
        return np.nan, np.nan, 0, 0
    ne = sub[sub["subtype"]==NEPC_LABEL]["value"].values
    cr = sub[sub["subtype"]==CRPC_LABEL]["value"].values

    all_vals = np.concatenate([ne, cr]) if len(ne)+len(cr) else np.array([0.0])
    has_negative = np.any(all_vals < 0)
    is_zscore = has_negative or (profile_type == "zscore" and abs(np.nanmean(all_vals)) < 5)

    should_log = use_log_scale or (profile_type=="raw_fpkm" and not is_zscore and not disable_auto_log)
    if should_log:
        if is_zscore:
            offset = abs(np.nanmin(all_vals)) + 1 if np.nanmin(all_vals) <= 0 else 0
            ne = np.log2(ne + offset); cr = np.log2(cr + offset)
            ylab = " (log2 transformed)"
        else:
            ne = np.log2(ne + 1); cr = np.log2(cr + 1)
            ylab = " (log2 FPKM)"
    else:
        ylab = " (Z-score)" if is_zscore else " (FPKM)"

    _, p = welch_t(ne, cr); d = cohens_d(ne, cr)

    plt.close('all')
    fig, ax = plt.subplots(figsize=(6,6), dpi=180)
    ax.boxplot([ne, cr],
               labels=[f"{NEPC_LABEL}\n(n={len(ne)})", f"{CRPC_LABEL}\n(n={len(cr)})"],
               showfliers=False, patch_artist=True, widths=0.6)
    rng = np.random.default_rng(7)
    ax.scatter(rng.normal(1,0.045,len(ne)), ne, s=26, alpha=0.75, zorder=3)
    ax.scatter(rng.normal(2,0.045,len(cr)), cr, s=26, alpha=0.75, zorder=3)
    ax.set_ylabel(f"{gene} Expression{ylab}", fontsize=12, fontweight="bold")

    sig = "ns"
    if pd.notna(p):
        if p < 1e-3: sig = "***"
        elif p < 1e-2: sig = "**"
        elif p < 5e-2: sig = "*"
    ax.set_title(f"{gene}: {NEPC_LABEL} vs {CRPC_LABEL}\nWelch p={p:.3g}, d={d:.2f}\n{sig}", fontsize=14, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.savefig(out_png, dpi=300)  # no bbox_inches='tight'
    plt.close(fig)
    return p, d, len(ne), len(cr)

# ---------- PURE MATPLOTLIB HEATMAP (fixed axes; no seaborn) ----------
def plot_heatmap_mycn_pp2a(df: pd.DataFrame, genes: List[str], out_png: str):
    """
    Build a heatmap with:
      - explicit axes rectangles (no tight/constrained layout)
      - a locked right-side vertical colorbar
      - a top subtype strip
      - gene labels on RIGHT
    """
    if "subtype" not in df.columns:
        raise ValueError("DataFrame must contain 'subtype' column for heatmap")

    desired_order = ['MYCN', 'PPP2R1B', 'PPP2R2B', 'PPP2R5C', 'PPP2CA', 'CIP2A', 'SET', 'ARPP19', 'PABIR1', 'PPME1']
    present = [g for g in desired_order if g in df.columns]
    missing = [g for g in desired_order if g not in df.columns]
    if len(present) < 2:
        raise ValueError(f"Insufficient genes for heatmap: {present}")

    expr = df[present].copy()
    # Z-score per gene
    Z = expr.apply(lambda c: (c - np.nanmean(c)) / (np.nanstd(c) if np.nanstd(c) > 0 else 1.0), axis=0)

    # order samples by subtype (NEPC first)
    subtype = df["subtype"].reindex(Z.index)
    order_map = {NEPC_LABEL: 0, CRPC_LABEL: 1}
    ordered_idx = subtype.map(order_map).sort_values(kind="mergesort").index
    Zs = Z.loc[ordered_idx]
    sub_sorted = subtype.loc[ordered_idx]

    # colors
    lut = {NEPC_LABEL: "#d62728", CRPC_LABEL: "#1f77b4"}
    strip_colors = sub_sorted.map(lut).to_list()
    n_ne = int((sub_sorted==NEPC_LABEL).sum())
    n_ad = int((sub_sorted==CRPC_LABEL).sum())

    # ---- build figure with absolute axes (in figure coords)
    plt.close('all')
    fig = plt.figure(figsize=(24, 10), dpi=300)

    # axes rects: [left, bottom, width, height] in 0..1
    ax_title = fig.add_axes([0.08, 0.94, 0.78, 0.04], frameon=False)
    ax_title.axis("off")
    ax_title.text(0.5, 0.5,
                  f"MYCN & PP2A Subunits: {NEPC_LABEL} (n={n_ne}) vs {CRPC_LABEL} (n={n_ad}) | "
                  f"Total: {len(sub_sorted)} samples | {len(present)} genes",
                  ha='center', va='center', fontsize=14, fontweight='bold')

    ax_strip = fig.add_axes([0.08, 0.89, 0.70, 0.035])
    # draw subtype strip
    rgb = np.array([[int(c.strip("#")[i:i+2], 16) for i in (0,2,4)] for c in strip_colors], float)/255.0
    ax_strip.imshow(rgb[np.newaxis, :, :], aspect="auto")
    ax_strip.set_axis_off()
    ax_strip.text(1.002, 0.5, "Subtype", transform=ax_strip.transAxes, va="center", ha="left", fontsize=10)

    ax_heat = fig.add_axes([0.08, 0.12, 0.70, 0.65])     # heatmap (reduced width to make room for legends)

    # draw heatmap (samples in columns, genes in rows)
    M = Zs[present].to_numpy().T  # genes x samples
    im = ax_heat.imshow(M, aspect="auto", interpolation="nearest",
                        cmap="RdBu_r", vmin=-2, vmax=2)
    
    # Explicitly prevent any automatic colorbar creation
    # Remove any existing colorbar that might have been auto-created
    if hasattr(im, 'colorbar'):
        im.colorbar = None

    # ticks/labels
    ax_heat.set_xlabel("Samples", fontsize=11, fontweight='bold')
    ax_heat.set_ylabel("Genes", fontsize=11, fontweight='bold')
    ax_heat.set_yticks(np.arange(len(present)))
    ax_heat.set_yticklabels(present)
    # move gene labels to RIGHT: mirror by turning off left spine/ticks and draw right ticks
    ax_heat.yaxis.tick_right()
    ax_heat.yaxis.set_label_position("right")
    ax_heat.tick_params(axis='y', which='major', pad=10, right=True, left=False, labelright=True, labelleft=False)
    ax_heat.set_xticks([])  # too dense; we hide them

    # legend (at top right) - create in separate axes to avoid conflicts
    legend_ax = fig.add_axes([0.83, 0.52, 0.12, 0.05], frameon=False)
    legend_ax.axis('off')
    legend_handles = [
        Patch(facecolor=lut[NEPC_LABEL], label=f'{NEPC_LABEL} (n={n_ne})'),
        Patch(facecolor=lut[CRPC_LABEL], label=f'{CRPC_LABEL} (n={n_ad})')
    ]
    legend_ax.legend(handles=legend_handles, loc='center', frameon=True,
                     fontsize=11, title='Subtype', title_fontsize=12, ncol=1)

    # Add colorbar for z-score legend UNDER the subtype legend
    cbar_ax = fig.add_axes([0.84, 0.35, 0.025, 0.12])  # positioned below subtype legend
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Z-score', fontsize=11, fontweight='bold', rotation=270, labelpad=15)
    cbar.ax.tick_params(labelsize=9)

    # Remove any automatically created colorbars
    axes_to_remove = []
    for ax in list(fig.axes):  # Use list() to avoid modification during iteration
        if ax == cbar_ax:  # Don't remove our intentional colorbar
            continue
        pos = ax.get_position()
        # Remove any narrow vertical axes that look like colorbars
        if pos.width < 0.05 and pos.height > 0.2:
            axes_to_remove.append(ax)
    for ax in axes_to_remove:
        ax.remove()

    # SAVE (no tight/constrained)
    fig.savefig(out_png, dpi=300, facecolor='white')
    plt.close(fig)

    return {
        "genes_requested": len(genes),
        "genes_included": len(present),
        "genes_missing": missing,
        "genes_present": present
    }

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="outputs_nepc", help="Output directory")
    ap.add_argument("--genes", nargs="+", default=["MYCN"],
                    help="Genes to plot")
    ap.add_argument("--clinical-column", default=None,
                    help="Override clinical column to map NEPC/CRPC (e.g., MORPHOLOGY)")
    ap.add_argument("--log-scale", action="store_true",
                    help="Force log transform for boxplots")
    ap.add_argument("--no-log-scale", action="store_true",
                    help="Disable automatic log transform for raw FPKM")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1) samples
    print(f"[INFO] Fetching samples for {STUDY_ID} …")
    df_samples = list_samples_scoped(STUDY_ID)
    if df_samples.empty: raise RuntimeError("No samples in study.")
    sample_ids = df_samples["sampleId"].astype(str).tolist()
    print(f"[INFO] N samples: {len(sample_ids)}")

    # 2) clinical (wide)
    print("[INFO] Fetching clinical table (study-scoped) …")
    clin_wide = get_clinical_wide(STUDY_ID)
    print(f"[INFO] Clinical columns (head): {list(clin_wide.columns)[:25]}")
    if args.clinical_column:
        chosen = args.clinical_column
        if chosen not in clin_wide.columns:
            raise RuntimeError(f"--clinical-column '{chosen}' not found.")
        print(f"[INFO] Using user-provided clinical column: {chosen}")
    else:
        chosen = pick_best_morphology_column(clin_wide)

    if "sampleId" not in clin_wide.columns:
        if clin_wide.index.name == "sampleId":
            clin_wide = clin_wide.reset_index()
        else:
            raise RuntimeError("No 'sampleId' column after clinical pivot.")

    clin = clin_wide[["sampleId", chosen]].copy()
    clin["subtype"] = map_nepc_vs_crpc(clin[chosen])

    # keep Neuroendocrine + Adenocarcinoma only for CRPC group
    before = len(clin)
    is_adeno = clin[chosen].astype(str).str.lower().str.contains("adenocarcinoma", na=False)
    is_neuro = clin["subtype"] == NEPC_LABEL
    clin = clin[is_neuro | is_adeno].copy()
    print(f"[INFO] Filtered out {before - len(clin)} non-Adeno samples from CRPC group")

    filtered_sample_ids = clin["sampleId"].astype(str).tolist()
    n_ne, n_cr = int((clin["subtype"]==NEPC_LABEL).sum()), int((clin["subtype"]==CRPC_LABEL).sum())
    print(f"[INFO] Subtype counts: {NEPC_LABEL}={n_ne}, {CRPC_LABEL}={n_cr}")

    # 3) profile + expression
    print("[INFO] Discovering molecular profiles …")
    prof_df = list_molecular_profiles(STUDY_ID)
    mol_profile, profile_type = choose_mrna_profile(prof_df, STUDY_ID)

    genes = list(dict.fromkeys([g.strip() for g in args.genes if g.strip()]))
    print(f"[INFO] Fetching expression for genes: {genes}")
    expr = fetch_expression(mol_profile, filtered_sample_ids, genes, STUDY_ID)

    fetched_genes = [g for g in genes if g in expr.columns]
    missing_fetched = [g for g in genes if g not in expr.columns]
    if missing_fetched:
        print(f"[WARN] Missing in expression: {', '.join(missing_fetched)}")

    df = expr.merge(clin[["sampleId","subtype"]], on="sampleId", how="left")

    # 4) long for boxplots
    long_rows = []
    for g in genes:
        if g not in df.columns:
            warnings.warn(f"{g} not present in expression matrix; skipping.")
            continue
        sub = df[["sampleId","subtype", g]].copy()
        sub["gene"] = g; sub.rename(columns={g:"value"}, inplace=True)
        long_rows.append(sub)
    if not long_rows: raise RuntimeError("No gene columns available to plot.")
    df_long = pd.concat(long_rows, ignore_index=True)

    # 5) boxplots + stats
    rows = []
    for g in genes:
        out_png = os.path.join(args.outdir, "boxplots", f"{g}_box_{NEPC_LABEL}_vs_{CRPC_LABEL}.png")
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        p, d, n1, n0 = boxplot_gene(df_long, g, out_png,
                                    use_log_scale=args.log_scale,
                                    profile_type=profile_type,
                                    disable_auto_log=args.no_log_scale)
        rows.append({"gene": g, "p_welch": p, "cohen_d": d,
                     f"n_{NEPC_LABEL}": n1, f"n_{CRPC_LABEL}": n0,
                     "profile": mol_profile, "clinical_column": chosen})
        if os.path.exists(out_png): print(f"[OK] {out_png}")
    stats_df = pd.DataFrame(rows)
    if stats_df["p_welch"].notna().any():
        m = stats_df["p_welch"].notna()
        _, q, _, _ = multipletests(stats_df.loc[m,"p_welch"].values, method="fdr_bh")
        stats_df.loc[m,"q_bh"] = q
    out_csv = os.path.join(args.outdir, "cbio_nepc_boxplot_stats.csv")
    stats_df.to_csv(out_csv, index=False)
    print(f"[DONE] Stats table: {out_csv}")

    # 6) heatmap
    if len(genes) > 1:
        available = [g for g in genes if g in df.columns]
        if len(available) >= 2:
            heatmap_png = os.path.join(args.outdir, "heatmaps", f"MYCN_PP2A_{NEPC_LABEL}_vs_{CRPC_LABEL}_heatmap.png")
            os.makedirs(os.path.dirname(heatmap_png), exist_ok=True)
            df_heatmap = df.set_index("sampleId") if "sampleId" in df.columns else df
            print(f"[INFO] Generating heatmap with {len(available)} genes for {len(df_heatmap)} samples…")
            summary = plot_heatmap_mycn_pp2a(df_heatmap, genes, heatmap_png)
            if os.path.exists(heatmap_png):
                print(f"[OK] {heatmap_png}")
                if summary:
                    print(f"[INFO] Heatmap: {summary['genes_included']}/{summary['genes_requested']} genes")
                    if summary['genes_missing']:
                        print(f"[WARN] Missing genes: {', '.join(summary['genes_missing'])}")
            else:
                print(f"[ERROR] Heatmap file was not created: {heatmap_png}")
        else:
            print(f"[WARN] Skipping heatmap: only {len(available)} gene(s) available")

if __name__ == "__main__":
    main()
