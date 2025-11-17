# nepc_analysis.py
# -*- coding: utf-8 -*-

"""
NEPC analysis (GSE35988 / GPL6480):
Anchor = MYCN ; Features = PP2A core subunits

Generates:
  1) Boxplot: MYCN expression (NEPC vs CRPC) with Welch t-test
  2) Correlations: MYCN vs each PP2A subunit (NEPC-only, CRPC-only, combined-colored)
     - Pearson (r, p) and Spearman (rho, p) + BH-FDR across tests
  3) Z-scored clustered heatmap for (MYCN + PP2A subunits) with subtype bar
  4) Tidy CSV with all stats
  5) One-page PDF summary assembling the key figures

Inputs (project root):
  - Data/Nepc/GSE35988-GPL6480_series_matrix.txt.gz
  - (optional) Data/Nepc/sample_labels_nepc.csv  => columns: sample_id,subtype ; subtype in {NEPC, CRPC}
  - (optional) Data/Nepc/GPL6480_family.soft.gz (speeds first run; otherwise script will still try)

Outputs:
  ./outputs_nepc/
    ├── boxplots/MYCN_box_NEPC_vs_CRPC.png
    ├── scatter/* (per-gene figures)
    ├── heatmaps/MYCN_PP2A_NEPC_vs_CRPC_heatmap.png
    ├── MYCN_PP2A_results_NEPC.csv
    └── NEPC_Mycn_PP2A_Summary.pdf
"""

import os, re, io, gzip, textwrap, warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Config
# -----------------------------
PROJECT_ROOT = "."
SERIES_MATRIX = os.path.join(PROJECT_ROOT, "Data/Nepc/GSE35988-GPL6480_series_matrix.txt.gz")

OUT_DIR  = os.path.join(PROJECT_ROOT, "outputs_nepc")
BOX_DIR  = os.path.join(OUT_DIR, "boxplots")
SCAT_DIR = os.path.join(OUT_DIR, "scatter")
HEAT_DIR = os.path.join(OUT_DIR, "heatmaps")
os.makedirs(BOX_DIR, exist_ok=True)
os.makedirs(SCAT_DIR, exist_ok=True)
os.makedirs(HEAT_DIR, exist_ok=True)

ANCHOR_GENE = "MYCN"
PP2A_GENES  = ["PPP2R5C", "PPP2R5D", "PPP2R1B", "PPP2CA", "PPP2R2D"]

GROUP_A = "NEPC"
GROUP_B = "CRPC"


# -----------------------------
# Utilities
# -----------------------------
def load_gpl_map(gpl_soft_path: str) -> dict:
    """Return dict: probe_id -> gene_symbol (first symbol if multiple; fallback empty)."""
    if not os.path.exists(gpl_soft_path):
        return {}
    probe_to_symbol = {}
    in_table = False
    header = []

    opener = gzip.open if gpl_soft_path.endswith(".gz") else open
    mode   = "rt" if gpl_soft_path.endswith(".gz") else "r"

    try:
        with opener(gpl_soft_path, mode, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("!platform_table_begin"):
                    in_table = True
                    header = f.readline().strip().split("\t")
                    continue
                if line.startswith("!platform_table_end"):
                    break
                if not in_table:
                    # Some GPLs are simple TSV with header "ID", "GENE_SYMBOL"/"Gene Symbol"
                    if line.startswith("ID\t"):
                        in_table = True
                        header = line.strip().split("\t")
                        continue
                    else:
                        continue
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.strip().split("\t")
                if len(parts) < len(header):
                    parts += [""] * (len(header) - len(parts))
                rec = dict(zip(header, parts))
                probe  = rec.get("ID") or rec.get("ID_REF") or ""
                symbol = rec.get("GENE_SYMBOL") or rec.get("Gene Symbol") or rec.get("Gene Symbol;") or ""
                if symbol:
                    symbol = re.split(r"[;/,\s]+", symbol.strip())[0]
                if probe and symbol and symbol not in {"---", "NA", ""}:
                    probe_to_symbol[probe] = symbol
    except Exception as e:
        warnings.warn(f"Error reading GPL file {gpl_soft_path}: {e}")
    return probe_to_symbol


def read_series_matrix(path: str) -> Tuple[pd.DataFrame, Dict[str, Dict[str,str]]]:
    """
    Returns:
      expr: DataFrame (rows = gene symbols if available / probe IDs fallback; cols = GSM IDs; values = log2 expr)
      meta: dict GSM -> {header:value} from the sample header lines
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Series matrix missing: {path}")

    # Parse GSM accessions
    gsm_ids: List[str] = []
    sample_meta: Dict[str, Dict[str, str]] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsm_ids = [x.strip().strip('"') for x in line.strip().split("\t")[1:]]
                break
    if not gsm_ids:
        raise ValueError("Could not find !Sample_geo_accession line.")

    for gsm in gsm_ids:
        sample_meta[gsm] = {}

    # Fill metadata + find table
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if raw.startswith("!Sample_"):
                parts = raw.strip().split("\t")
                tag = parts[0][len("!Sample_"):]
                vals = parts[1:]
                for gsm, val in zip(gsm_ids, vals):
                    sample_meta[gsm][tag] = val
            if raw.startswith("!series_matrix_table_begin"):
                break

        header = f.readline().strip().split("\t")
        cols = [x.strip().strip('"') for x in header[1:]]
        ids, data_rows = [], []
        for raw in f:
            if raw.startswith("!series_matrix_table_end"):
                break
            parts = raw.strip().split("\t")
            ids.append(parts[0].strip('"'))
            vals = [v.strip('"') for v in parts[1:]]
            def to_float(x):
                try: return float(x)
                except: return np.nan
            data_rows.append([to_float(x) for x in vals])

    expr = pd.DataFrame(data_rows, index=ids, columns=cols)

    # Map probes→symbols using GPL6480 (Agilent 4x44K)
    gpl_try = [
        os.path.join(PROJECT_ROOT, "Data/Nepc/GSE35988_family.soft.gz"),  # Try dataset-specific name first
        os.path.join(PROJECT_ROOT, "Data/Nepc/GPL6480_family.soft.gz"),
        os.path.join(PROJECT_ROOT, "Data/Nepc/GPL6480.soft"),
        os.path.join(PROJECT_ROOT, "GPL6480_family.soft.gz"),
        os.path.join(PROJECT_ROOT, "GPL6480.soft"),
    ]
    gpl_path = next((p for p in gpl_try if os.path.exists(p)), None)
    probe_map = load_gpl_map(gpl_path) if gpl_path else {}

    if probe_map:
        idx = expr.index.to_series().map(lambda x: probe_map.get(x, np.nan))
        has_sym = idx.notna()
        if has_sym.any():
            expr = expr.loc[has_sym]
            expr.index = idx.loc[has_sym].astype(str)
            expr = expr.groupby(expr.index).median()
            print(f"[INFO] Mapped {has_sym.sum()} probes to gene symbols")
        else:
            warnings.warn("GPL annotation found but no probes matched. Proceeding with probe-level IDs.")
    else:
        # Heuristic: if IDs already look like symbols, collapse by median
        def looks_like_symbol(s):
            return bool(re.match(r"^[A-Z0-9\-\.]{2,}$", s))
        frac = np.mean([looks_like_symbol(r) for r in expr.index[:min(2000, len(expr))]])
        if frac > 0.7:
            expr = expr.groupby(expr.index).median()
        else:
            warnings.warn("GPL6480 not found: rows likely probe IDs. Consider placing GPL6480_family.soft.gz in Data/Nepc/.")

    expr.index.name = "symbol_or_probe"
    return expr, sample_meta


def load_subtypes_nepc(sample_meta: Dict[str, Dict[str,str]]) -> pd.Series:
    """
    Return pd.Series: index = GSM IDs, values ∈ {'NEPC','CRPC'}
    Precedence:
      1) Data/Nepc/sample_labels_nepc.csv  (columns: sample_id, subtype)
      2) heuristic from title/characteristics/description
    """
    # 1) explicit labels
    labfile = os.path.join(PROJECT_ROOT, "Data/Nepc/sample_labels_nepc.csv")
    if os.path.exists(labfile):
        df = pd.read_csv(labfile)
        df.columns = [c.lower() for c in df.columns]
        mapping = {r["sample_id"]: r["subtype"] for _, r in df.iterrows()}
        ser = pd.Series({gsm: mapping.get(gsm, np.nan) for gsm in sample_meta})
        ser = ser.map(lambda x: GROUP_A if str(x).strip().upper() in {"NEPC","SCNC","SMALL CELL","NEUROENDOCRINE"} 
                                   else (GROUP_B if str(x).strip().upper() in {"CRPC","ADENO","ADENOCARCINOMA"} else np.nan))
        return ser

    # 2) heuristic
    def parse_one(meta: Dict[str,str]) -> str:
        text = " ".join(str(v) for v in meta.values()).lower()
        if re.search(r"\b(nepc|neuroendocrine|small[\s-]*cell)\b", text):
            return GROUP_A
        if re.search(r"\b(crpc|adenocarcinoma|adeno)\b", text):
            return GROUP_B
        return np.nan

    return pd.Series({gsm: parse_one(m) for gsm, m in sample_meta.items()})


def split_by_subtype(expr: pd.DataFrame, subtype: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame]:
    subtype = subtype.reindex(expr.columns)
    a_cols = list(subtype[subtype==GROUP_A].index)
    b_cols = list(subtype[subtype==GROUP_B].index)
    return expr[a_cols], expr[b_cols]


def welch_ttest(x: np.ndarray, y: np.ndarray):
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2: return (np.nan, np.nan)
    t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
    return t, p


def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2: return np.nan
    nx, ny = len(x), len(y)
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    pooled = ((nx-1)*vx + (ny-1)*vy) / (nx+ny-2)
    return (x.mean() - y.mean()) / np.sqrt(pooled) if pooled>0 else np.nan


def pearson_spearman(x: np.ndarray, y: np.ndarray) -> Dict[str,float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]
    if len(x) < 3:
        return dict(pearson_r=np.nan, pearson_p=np.nan, spearman_r=np.nan, spearman_p=np.nan, n=len(x))
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    return dict(pearson_r=pr, pearson_p=pp, spearman_r=sr, spearman_p=sp, n=len(x))


# -----------------------------
# Plotters
# -----------------------------
def plot_boxplot_groups(a: pd.Series, b: pd.Series, gene_name: str, outpng: str):
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    data = [a.values, b.values]  # order: NEPC, CRPC
    bp = ax.boxplot(
        data, labels=[GROUP_A, GROUP_B],
        showfliers=False, patch_artist=True, widths=0.5,
        boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.5),
        medianprops=dict(color='red', linewidth=2),
        whiskerprops=dict(color='black', linewidth=1.5),
        capprops=dict(color='black', linewidth=1.5)
    )
    # jittered points - all samples in black
    np.random.seed(42)
    # Filter out NaN/infinite values for plotting (only plot valid data points)
    a_finite = a.values[np.isfinite(a.values)]
    b_finite = b.values[np.isfinite(b.values)]
    # NEPC: black points
    if len(a_finite) > 0:
        ax.scatter(np.random.normal(1, 0.04, size=len(a_finite)), a_finite, 
                   s=30, color='black', alpha=0.6, edgecolors='black', linewidths=0.5, zorder=4)
    # CRPC: black points
    if len(b_finite) > 0:
        ax.scatter(np.random.normal(2, 0.04, size=len(b_finite)), b_finite, 
                   s=30, color='black', alpha=0.6, edgecolors='black', linewidths=0.5, zorder=3)
    ax.set_ylabel(f"{gene_name} (log2 expr)", fontsize=12, fontweight='bold')

    _, p = welch_ttest(a.values, b.values)
    d = cohens_d(a.values, b.values)
    p_text = f"p<0.001 (p={p:.2e})" if p < 1e-3 else (f"p={p:.4f}" if p < 0.01 else f"p={p:.3f}")
    ax.set_title(f"{gene_name}: {GROUP_A} vs {GROUP_B}\nWelch {p_text}, d={d:.2f}",
                 fontsize=13, fontweight='bold')
    # Show actual number of finite values plotted (excludes NaN)
    n_a_plot = len(a_finite) if len(a_finite) > 0 else 0
    n_b_plot = len(b_finite) if len(b_finite) > 0 else 0
    ax.set_xticklabels([f"{GROUP_A}\n(n={n_a_plot})", f"{GROUP_B}\n(n={n_b_plot})"], fontsize=11)

    # bracket + sig stars
    y_max, y_min = np.max([a.max(), b.max()]), np.min([a.min(), b.min()])
    y_range = y_max - y_min
    bracket_y = y_max + y_range * 0.05
    ax.plot([1,1,2,2], [bracket_y, bracket_y + y_range*0.02, bracket_y + y_range*0.02, bracket_y], color='black', lw=1.5)
    sig = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
    ax.text(1.5, bracket_y + y_range*0.03, sig, ha='center', va='bottom', fontsize=12, fontweight='bold')

    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.spines['left'].set_linewidth(1.5); ax.spines['bottom'].set_linewidth(1.5)
    fig.tight_layout(); fig.savefig(outpng, bbox_inches="tight"); plt.close(fig)
    return p


def plot_scatter(x: pd.Series, y: pd.Series, title: str, outpng: str, r_val=None, p_val=None):
    m = np.isfinite(x.values) & np.isfinite(y.values)
    xv, yv = x.values[m], y.values[m]
    fig, ax = plt.subplots(figsize=(4.0, 4.0), dpi=200)
    ax.scatter(xv, yv, s=18, alpha=0.85)
    if len(xv) >= 3:
        slope, intercept, *_ = stats.linregress(xv, yv)
        xs = np.linspace(xv.min(), xv.max(), 100); ys = slope*xs + intercept
        ax.plot(xs, ys, linewidth=2)
        if r_val is None or p_val is None:
            r_val, p_val = stats.pearsonr(xv, yv)
        ax.text(0.02, 0.96, f"R={r_val:.2f}\np={p_val:.2g}", transform=ax.transAxes,
                ha="left", va="top", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7"))
    ax.set_xlabel("MYCN"); ax.set_ylabel(title.split(" vs ")[0]); ax.set_title(title)
    ax.grid(True, alpha=0.2, linestyle=":"); fig.tight_layout()
    fig.savefig(outpng, bbox_inches="tight"); plt.close(fig)


def plot_heatmap_nepc(expr: pd.DataFrame, subtype: pd.Series, genes: List[str], outpng: str):
    import seaborn as sns
    present = [g for g in genes if g in expr.index]
    if not present:
        warnings.warn("No genes available for heatmap"); return
    X = expr.loc[present]
    Z = X.apply(lambda r: (r - np.nanmean(r)) / (np.nanstd(r) if np.nanstd(r)>0 else 1.0), axis=1).T
    subtype_aligned = subtype.reindex(Z.index)
    lut = {GROUP_A: "#d62728", GROUP_B: "#1f77b4"}
    row_colors = subtype_aligned.map(lut)

    g = sns.clustermap(
        Z.T, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        col_colors=row_colors, row_cluster=True, col_cluster=True,
        figsize=(12,6), dendrogram_ratio=0.15, colors_ratio=0.03,
        cbar_pos=(0.02,0.83,0.03,0.15),
        linewidths=0, method="average", metric="correlation",
        yticklabels=True, xticklabels=False,
        cbar_kws={"label":"Row Z-Score","orientation":"vertical","ticks":[-2,-1,0,1,2]}
    )
    g.fig.suptitle(f"MYCN & PP2A Subunits: {GROUP_A} vs {GROUP_B}", fontsize=13, fontweight='bold', y=0.98)

    from matplotlib.patches import Patch
    legend = [
        Patch(facecolor=lut[GROUP_A], label=f'{GROUP_A} (n={(subtype_aligned==GROUP_A).sum()})'),
        Patch(facecolor=lut[GROUP_B], label=f'{GROUP_B} (n={(subtype_aligned==GROUP_B).sum()})')
    ]
    g.ax_heatmap.legend(handles=legend, loc='upper left', bbox_to_anchor=(1.15,1), frameon=True, fontsize=9, title='Subtype', title_fontsize=10)

    g.ax_cbar.set_title('Color Key', fontsize=9, pad=10)
    plt.savefig(outpng, bbox_inches='tight', dpi=300); plt.close()


# -----------------------------
# PDF assembler
# -----------------------------
def build_summary_pdf(box_png: str, heat_png: str, scatter_pngs: List[str], caption: str, outpdf: str):
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(outpdf) as pdf:
        fig = plt.figure(figsize=(8.5,11), dpi=200)
        ax_heat = plt.axes([0.08, 0.74, 0.84, 0.22]); ax_heat.imshow(plt.imread(heat_png)); ax_heat.axis("off")
        ax_heat.set_title("PP2A Core Heatmap (Z-score)", fontsize=12, pad=6)

        ax_box  = plt.axes([0.08, 0.48, 0.38, 0.2]); ax_box.imshow(plt.imread(box_png)); ax_box.axis("off")

        pos = [
            [0.52, 0.52, 0.4, 0.16],
            [0.52, 0.34, 0.4, 0.16],
            [0.52, 0.16, 0.4, 0.16],
            [0.08, 0.28, 0.38, 0.16],
            [0.08, 0.10, 0.38, 0.16],
        ]
        for png, p in zip(scatter_pngs, pos):
            ax = plt.axes(p); ax.imshow(plt.imread(png)); ax.axis("off")

        ax_cap = plt.axes([0.08, 0.02, 0.84, 0.06]); ax_cap.axis("off")
        ax_cap.text(0, 0.9, f"{GROUP_A}: MYCN–PP2A summary (GSE35988 / GPL6480)", fontsize=11, weight="bold")
        ax_cap.text(0, 0.1, caption, fontsize=9, va="top", wrap=True)
        pdf.savefig(fig); plt.close(fig)


# -----------------------------
# Main
# -----------------------------
def main():
    print("[INFO] Reading series matrix…")
    expr, meta = read_series_matrix(SERIES_MATRIX)
    print(f"[INFO] matrix: {expr.shape[0]} rows x {expr.shape[1]} samples")

    subtype = load_subtypes_nepc(meta)  # 'NEPC' / 'CRPC' / NaN
    subtype = subtype.reindex(expr.columns)
    n_a = int((subtype==GROUP_A).sum()); n_b = int((subtype==GROUP_B).sum())
    if n_a==0 or n_b==0:
        warnings.warn(f"Subtype counts look imbalanced: {GROUP_A}={n_a}, {GROUP_B}={n_b}. "
                      f"If needed, add Data/Nepc/sample_labels_nepc.csv")
    else:
        print(f"[INFO] Subtype counts: {GROUP_A}={n_a}, {GROUP_B}={n_b}")

    # Ensure required genes
    need = [ANCHOR_GENE] + PP2A_GENES
    available = [g for g in need if g in expr.index]
    missing   = [g for g in need if g not in expr.index]
    if missing:
        warnings.warn(f"Missing genes in matrix: {missing}")
    genes = available
    if len(genes)==0:
        raise RuntimeError("None of the target genes found. Ensure GPL mapping worked (GPL6480).")

    # Split
    a_expr, b_expr = split_by_subtype(expr.loc[genes], subtype)  # NEPC, CRPC

    # 1) MYCN boxplot
    p_box = np.nan
    if ANCHOR_GENE in genes:
        p_box = plot_boxplot_groups(a_expr.loc[ANCHOR_GENE], b_expr.loc[ANCHOR_GENE], ANCHOR_GENE,
                                    os.path.join(BOX_DIR, "MYCN_box_NEPC_vs_CRPC.png"))
    else:
        warnings.warn("MYCN not found; boxplot skipped.")

    # 1b) PP2A subunit boxplots for key genes showing differential expression
    pp2a_boxplot_genes = ["PPP2R1B", "PPP2R5D", "PPP2R2D", "PPP2R5C"]
    for gene in pp2a_boxplot_genes:
        if gene in genes:
            plot_boxplot_groups(a_expr.loc[gene], b_expr.loc[gene], gene,
                               os.path.join(BOX_DIR, f"{gene}_box_NEPC_vs_CRPC.png"))
            print(f"[OK] Generated boxplot for {gene}")
        else:
            print(f"[WARN] {gene} not found; boxplot skipped.")

    # 2) Correlations: MYCN vs each PP2A subunit (grouped & colored)
    records, scatter_for_pdf = [], []
    for gene in PP2A_GENES:
        if gene not in genes or ANCHOR_GENE not in genes:
            continue
        # groupwise stats
        stats_a = pearson_spearman(a_expr.loc[ANCHOR_GENE].values, a_expr.loc[gene].values)
        stats_b = pearson_spearman(b_expr.loc[ANCHOR_GENE].values, b_expr.loc[gene].values)

        # combined colored with group-specific regression lines
        df = pd.DataFrame({"x": expr.loc[ANCHOR_GENE], "y": expr.loc[gene], "subtype": subtype}).dropna()
        fig, ax = plt.subplots(figsize=(4,4), dpi=200)
        
        # Calculate overall correlation for all samples
        if len(df) >= 3:
            r_all, p_all = stats.pearsonr(df["x"], df["y"])
        
        # Plot each group with its own regression line
        for lab, c in [(GROUP_A,"#d62728"), (GROUP_B,"#1f77b4")]:
            sub = df[df["subtype"]==lab]
            if len(sub) > 0:
                ax.scatter(sub["x"], sub["y"], s=18, label=lab, alpha=0.85, color=c)
                # Group-specific regression line if >= 3 points
                if len(sub) >= 3 and np.nanstd(sub["x"]) > 0:
                    slope, intercept, *_ = stats.linregress(sub["x"], sub["y"])
                    xs = np.linspace(sub["x"].min(), sub["x"].max(), 100)
                    ax.plot(xs, slope*xs + intercept, lw=2, color=c, alpha=0.7)
        
        # Display overall correlation stats
        if len(df) >= 3:
            ax.text(0.02, 0.96, f"All R={r_all:.2f}\np={p_all:.2g}", transform=ax.transAxes,
                    ha="left", va="top", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7"))
        
        ax.set_xlabel("MYCN"); ax.set_ylabel(gene); ax.set_title(f"{gene} vs MYCN — colored by subtype")
        ax.grid(True, alpha=0.2, linestyle=":"); ax.legend(frameon=True)
        out_col = os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_colored.png")
        fig.tight_layout(); fig.savefig(out_col, bbox_inches="tight"); plt.close(fig)
        scatter_for_pdf.append(out_col)

        # per-group panels (use precomputed r,p to match table)
        plot_scatter(a_expr.loc[ANCHOR_GENE], a_expr.loc[gene],
                     f"{gene} vs MYCN — {GROUP_A}",
                     os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_{GROUP_A}.png"),
                     r_val=stats_a['pearson_r'], p_val=stats_a['pearson_p'])
        plot_scatter(b_expr.loc[ANCHOR_GENE], b_expr.loc[gene],
                     f"{gene} vs MYCN — {GROUP_B}",
                     os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_{GROUP_B}.png"),
                     r_val=stats_b['pearson_r'], p_val=stats_b['pearson_p'])

        records.append(dict(gene=gene, group=GROUP_A, **stats_a))
        records.append(dict(gene=gene, group=GROUP_B, **stats_b))

    stats_df = pd.DataFrame.from_records(records)

    # 2b) BH-FDR across Pearson p-values within each group
    if stats_df.empty:
        print("[WARN] No correlation rows generated. Skipping FDR and PDF.")
        stats_df = pd.DataFrame(columns=["gene","group","pearson_r","pearson_p","spearman_r","spearman_p","n"])
    else:
        frames = []
        for grp in [GROUP_A, GROUP_B]:
            sub = stats_df[stats_df["group"]==grp].copy()
            m = np.isfinite(sub["pearson_p"])
            if m.sum()>0:
                rej, p_adj, *_ = multipletests(sub.loc[m,"pearson_p"].values, alpha=0.05, method="fdr_bh")
                sub.loc[m,"pearson_fdr"] = p_adj
                sub.loc[m,"pearson_sig"] = rej
            frames.append(sub)
        stats_df = pd.concat(frames, ignore_index=True)

    # 3) Heatmap
    heat_genes = [ANCHOR_GENE] + PP2A_GENES
    heat_genes = [g for g in heat_genes if g in genes]
    if heat_genes:
        plot_heatmap_nepc(expr, subtype, heat_genes,
                          os.path.join(HEAT_DIR, "MYCN_PP2A_NEPC_vs_CRPC_heatmap.png"))

    # 4) Save table + caption
    stats_df.to_csv(os.path.join(OUT_DIR, "MYCN_PP2A_results_NEPC.csv"), index=False)

    cap_lines = [f"MYCN {GROUP_A} vs {GROUP_B} Welch p = {p_box:.3g}."]
    for gene in PP2A_GENES:
        row_a = stats_df[(stats_df["gene"]==gene)&(stats_df["group"]==GROUP_A)]
        row_b = stats_df[(stats_df["gene"]==gene)&(stats_df["group"]==GROUP_B)]
        if len(row_a):
            pa, qA = row_a["pearson_p"].values[0], row_a.get("pearson_fdr", pd.Series([np.nan])).values[0]
            cap_lines.append(f"{gene} ({GROUP_A}): Pearson p={pa:.3g}, FDR={qA:.3g}")
        if len(row_b):
            pb, qB = row_b["pearson_p"].values[0], row_b.get("pearson_fdr", pd.Series([np.nan])).values[0]
            cap_lines.append(f"{gene} ({GROUP_B}): Pearson p={pb:.3g}, FDR={qB:.3g}")
    caption = " | ".join(cap_lines)

    # 5) Assemble PDF
    pdf_out = os.path.join(OUT_DIR, "NEPC_Mycn_PP2A_Summary.pdf")
    heat_png = os.path.join(HEAT_DIR, "MYCN_PP2A_NEPC_vs_CRPC_heatmap.png")
    box_png  = os.path.join(BOX_DIR, "MYCN_box_NEPC_vs_CRPC.png")
    scatter_pngs = [os.path.join(SCAT_DIR, f"{g}_vs_MYCN_colored.png")
                    for g in PP2A_GENES if os.path.exists(os.path.join(SCAT_DIR, f"{g}_vs_MYCN_colored.png"))][:5]
    if os.path.exists(box_png) and os.path.exists(heat_png) and scatter_pngs:
        build_summary_pdf(box_png, heat_png, scatter_pngs, caption, pdf_out)
        print(f"[DONE] Summary PDF: {pdf_out}")
    else:
        print("[WARN] Could not build summary PDF (missing figures).")

    print("[DONE] Outputs written to ./outputs_nepc")


if __name__ == "__main__":
    main()
