# main.py
# -*- coding: utf-8 -*-

"""
ETP-ALL analysis (GSE28703):
Anchor = MYCN ; Features = PP2A core subunits
Generates:
  1) Boxplot: MYCN expression (ETP vs nonETP) with Welch t-test and BH-FDR
  2) Correlations: MYCN vs each PP2A subunit (ETP-only, nonETP-only, combined-colored)
     - Pearson (r, p) and Spearman (rho, p) + BH-FDR across tests
  3) Z-scored heatmap for PP2A subunits (rows) across samples (columns), with subtype bar
  4) Tidy CSV with all stats
  5) One-page PDF summary assembling the key figures

Inputs (project root):
  - GSE28703_series_matrix.txt.gz
  - (optional) sample_labels.csv  => columns: sample_id,subtype ; subtype in {ETP, nonETP}

Outputs:
  ./outputs/
    ├── boxplots/MYCN_box_ETP_vs_nonETP.png
    ├── scatter/* (per-gene figures)
    ├── heatmaps/PP2A_core_heatmap.png
    ├── MYCN_PP2A_results.csv
    └── ETP_Mycn_PP2A_Summary.pdf
"""

import os, io, re, gzip, textwrap, warnings
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
SERIES_MATRIX = os.path.join(PROJECT_ROOT, "Data/GSE28703_series_matrix.txt.gz")

OUT_DIR     = os.path.join(PROJECT_ROOT, "outputs")
BOX_DIR     = os.path.join(OUT_DIR, "boxplots")
SCAT_DIR    = os.path.join(OUT_DIR, "scatter")
HEAT_DIR    = os.path.join(OUT_DIR, "heatmaps")
os.makedirs(BOX_DIR, exist_ok=True)
os.makedirs(SCAT_DIR, exist_ok=True)
os.makedirs(HEAT_DIR, exist_ok=True)

ANCHOR_GENE = "MYCN"
PP2A_GENES  = ["PPP2R5C", "PPP2R5D", "PPP2R1B", "PPP2CA", "PPP2R2D"]

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
    
    # Handle both .gz and plain text
    opener = gzip.open if gpl_soft_path.endswith('.gz') else open
    mode = "rt" if gpl_soft_path.endswith('.gz') else "r"
    
    try:
        with opener(gpl_soft_path, mode, encoding="utf-8", errors="ignore") as f:
            for line in f:
                # Check for SOFT format markers
                if line.startswith("!platform_table_begin") or line.startswith("^DATABASE"):
                    in_table = True
                    # Next line should be header
                    header_line = f.readline()
                    if header_line:
                        header = header_line.strip().split("\t")
                    continue
                if line.startswith("!platform_table_end") or (in_table and line.startswith("!")):
                    break
                # Check for tab-separated table format (like GPL13158-5065.txt)
                if line.startswith("ID\t"):
                    in_table = True
                    header = line.strip().split("\t")
                    continue
                if not in_table:
                    continue
                # Skip comment lines
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < len(header):
                    parts.extend([""] * (len(header) - len(parts)))
                rec = dict(zip(header, parts))
                probe = rec.get("ID") or rec.get("ID_REF") or ""
                symbol = rec.get("Gene Symbol") or rec.get("GENE_SYMBOL") or rec.get("Gene Symbol;") or ""
                if symbol:
                    # take first symbol if "A;B;C" or "A /// B"
                    symbol = re.split(r"[;,\s/]+", symbol.strip())[0]
                if probe and symbol and symbol not in ["---", "NA", ""]:
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

    # Parse header to extract sample IDs and metadata
    gsm_ids: List[str] = []
    sample_meta: Dict[str, Dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsm_ids = [x.strip().strip('"') for x in line.strip().split("\t")[1:]]
                break

    if not gsm_ids:
        # Try raw text (if someone already decompressed)
        if path.endswith(".gz"):
            alt = path[:-3]
            if os.path.exists(alt):
                return read_series_matrix(alt)
        raise ValueError("Could not find !Sample_geo_accession line.")

    # Initialize meta containers
    for gsm in gsm_ids:
        sample_meta[gsm] = {}

    # Re-scan to fill meta
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if raw.startswith("!Sample_"):
                parts = raw.strip().split("\t")
                tag = parts[0][len("!Sample_"):]
                vals = parts[1:]
                for gsm, val in zip(gsm_ids, vals):
                    sample_meta[gsm][tag] = val

            if raw.startswith("!series_matrix_table_begin"):
                break

        # Read expression table
        # First row is header: ID_REF <tab> GSM...
        header = f.readline().strip().split("\t")
        cols = [x.strip().strip('"') for x in header[1:]]  # Strip quotes from GSM IDs
        data_rows = []
        ids = []
        for raw in f:
            if raw.startswith("!series_matrix_table_end"):
                break
            parts = raw.strip().split("\t")
            ids.append(parts[0].strip('"'))  # Strip quotes from probe IDs
            vals = [v.strip('"') for v in parts[1:]]  # Strip quotes from values too
            # Convert to float if possible
            def to_float(x):
                try: return float(x)
                except: return np.nan
            data_rows.append([to_float(x) for x in vals])

    expr = pd.DataFrame(data_rows, index=ids, columns=cols)
    
    # Map probe IDs to gene symbols using GPL annotation
    GPL_SOFT = os.path.join(PROJECT_ROOT, "GPL13158_family.soft.gz")
    if not os.path.exists(GPL_SOFT):
        # Try alternative locations
        GPL_SOFT = os.path.join(PROJECT_ROOT, "Data", "GPL13158_family.soft.gz")
    if not os.path.exists(GPL_SOFT):
        # Try plain text version
        GPL_SOFT = os.path.join(PROJECT_ROOT, "GPL13158_family.soft")
    if not os.path.exists(GPL_SOFT):
        GPL_SOFT = os.path.join(PROJECT_ROOT, "Data", "GPL13158-5065.txt")
    
    probe_map = load_gpl_map(GPL_SOFT)
    
    if probe_map:
        # Map probe IDs to gene symbols
        idx = expr.index.to_series().map(lambda x: probe_map.get(x, np.nan))
        has_sym = idx.notna()
        if has_sym.any():
            expr = expr.loc[has_sym]
            expr.index = idx.loc[has_sym].astype(str)
            # Collapse duplicate symbols by median
            expr = expr.groupby(expr.index).median()
            print(f"[INFO] Mapped {has_sym.sum()} probes to gene symbols")
        else:
            warnings.warn("GPL annotation found but no probes matched. Rows may be probe IDs.")
    else:
        # Fallback: check if rows already look like gene symbols
        def looks_like_symbol(s):
            return bool(re.match(r"^[A-Z0-9\-\.]{2,}$", s))
        frac_symbols = np.mean([looks_like_symbol(r) for r in expr.index[: min(1000, len(expr))]])
        if frac_symbols > 0.7:
            # Already gene-level; if duplicates, median collapse
            expr = expr.groupby(expr.index).median()
        else:
            warnings.warn("GPL annotation not found; rows likely probe IDs. Download GPL13158_family.soft.gz.")
    
    expr.index.name = "symbol_or_probe"
    return expr, sample_meta


def load_subtypes(sample_meta: Dict[str, Dict[str,str]]) -> pd.Series:
    """
    Return pd.Series: index = GSM IDs, values ∈ {'ETP','nonETP'}
    Order of precedence:
      1) sample_labels.csv file in cwd
      2) infer from Sample_title / Sample_characteristics / Sample_description
    """
    # 1) override file
    labfile = os.path.join(PROJECT_ROOT, "sample_labels.csv")
    if os.path.exists(labfile):
        df = pd.read_csv(labfile)
        df.columns = [c.lower() for c in df.columns]
        m = dict(zip(df["sample_id"], df["subtype"]))
        ser = pd.Series({gsm: m.get(gsm, np.nan) for gsm in sample_meta})
        # normalize
        ser = ser.map(lambda x: "ETP" if str(x).strip().upper() in {"ETP","ETP-ALL","EARLY T-CELL PRECURSOR"} else ("nonETP" if pd.notna(x) else np.nan))
        return ser

    # 2) heuristic parse from metadata
    def parse_one(meta: Dict[str,str]) -> str:
        text = " ".join(str(v) for v in meta.values()).lower()
        if re.search(r"\b(etp|early t[- ]cell precursor)\b", text):
            return "ETP"
        # non-ETP can appear as "T-ALL", "NOT ETP" etc. Fall back to nonETP if we see ALL but not ETP
        if "t-all" in text or "t all" in text:
            return "nonETP"
        return np.nan

    out = {gsm: parse_one(meta) for gsm, meta in sample_meta.items()}
    ser = pd.Series(out)
    return ser


def split_by_subtype(expr: pd.DataFrame, subtype: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame]:
    subtype = subtype.reindex(expr.columns)
    etp_cols    = list(subtype[subtype=="ETP"].index)
    non_cols    = list(subtype[subtype=="nonETP"].index)
    return expr[etp_cols], expr[non_cols]


def welch_ttest(x: np.ndarray, y: np.ndarray) -> Tuple[float,float]:
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
def plot_boxplot_etp_vs_nonetp(etp: pd.Series, non: pd.Series, gene_name: str, outpng: str):
    """
    Generic boxplot for any gene comparing ETP vs nonETP with all data points shown.
    """
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    data = [non.values, etp.values]  # order: nonETP, ETP
    
    # Create boxplot with better styling to match reference image
    # showfliers=False because we'll plot ALL points manually
    bp = ax.boxplot(data, labels=["nonETP","ETP"], 
                     showfliers=False,
                     patch_artist=True,
                     widths=0.5,
                     boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.5),
                     medianprops=dict(color='red', linewidth=2),
                     whiskerprops=dict(color='black', linewidth=1.5),
                     capprops=dict(color='black', linewidth=1.5))
    
    # CRITICAL: Overlay ALL individual data points (not just outliers)
    # Add jitter to x-coordinates so points don't overlap
    np.random.seed(42)  # for reproducibility
    x_non = np.random.normal(1, 0.04, size=len(non))  # jitter around position 1
    x_etp = np.random.normal(2, 0.04, size=len(etp))  # jitter around position 2
    
    ax.scatter(x_non, non.values, alpha=0.6, s=30, color='black', zorder=3)
    ax.scatter(x_etp, etp.values, alpha=0.6, s=30, color='black', zorder=3)
    
    ax.set_ylabel(f"{gene_name} (log2 expr)", fontsize=12, fontweight='bold')
    
    # stats
    _, p = welch_ttest(etp.values, non.values)  # ETP vs nonETP
    d = cohens_d(etp.values, non.values)
    
    # Format p-value appropriately
    if p < 0.001:
        p_text = f"p<0.001 (p={p:.2e})"
    elif p < 0.01:
        p_text = f"p={p:.4f}"
    else:
        p_text = f"p={p:.3f}"
    
    ax.set_title(f"{gene_name}: ETP vs nonETP\nWelch {p_text}, d={d:.2f}", 
                 fontsize=13, fontweight='bold')
    
    # Add sample sizes to x-labels
    ax.set_xticklabels([f"nonETP\n(n={len(non)})", f"ETP\n(n={len(etp)})"], 
                       fontsize=11)
    
    # Statistical significance bracket
    y_max = max(np.max(non.values), np.max(etp.values))
    y_min = min(np.min(non.values), np.min(etp.values))
    y_range = y_max - y_min
    
    # Draw bracket for significance
    bracket_y = y_max + y_range * 0.05
    ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + y_range * 0.02, 
                            bracket_y + y_range * 0.02, bracket_y], 
            color='black', linewidth=1.5)
    
    # Add p-value text
    if p < 0.001:
        sig_text = '***'
    elif p < 0.01:
        sig_text = '**'
    elif p < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'
    
    ax.text(1.5, bracket_y + y_range * 0.03, sig_text, 
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Remove top and right spines for cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    fig.tight_layout()
    fig.savefig(outpng, bbox_inches="tight")
    plt.close(fig)
    return p


def plot_box_mycn(etp: pd.Series, non: pd.Series, outpng: str):
    """Legacy wrapper for MYCN boxplot - calls generic function"""
    return plot_boxplot_etp_vs_nonetp(etp, non, "MYCN", outpng)


def plot_scatter(x: pd.Series, y: pd.Series, title: str, outpng: str, r_val: float = None, p_val: float = None):
    """Simple scatter with OLS line and 95% CI band.
    If r_val and p_val are provided, use those instead of recalculating.
    """
    # mask
    m = np.isfinite(x.values) & np.isfinite(y.values)
    xv = x.values[m]; yv = y.values[m]
    fig, ax = plt.subplots(figsize=(4.0, 4.0), dpi=200)
    ax.scatter(xv, yv, s=18, alpha=0.85)
    # regression line
    if len(xv) >= 3:
        if r_val is None or p_val is None:
            # Calculate if not provided
            slope, intercept, r, p, _ = stats.linregress(xv, yv)
        else:
            # Use provided values (matches CSV/summary table)
            r = r_val
            p = p_val
            slope, intercept, _, _, _ = stats.linregress(xv, yv)
        xs = np.linspace(xv.min(), xv.max(), 100)
        ys = slope*xs + intercept
        ax.plot(xs, ys, linewidth=2)
        ax.text(0.02, 0.96, f"R={r:.2f}\np={p:.2g}", transform=ax.transAxes,
                ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7"))
    ax.set_xlabel("MYCN")
    ax.set_ylabel(title.split(" vs ")[0])  # gene name
    ax.set_title(title)
    ax.grid(True, alpha=0.2, linestyle=":")
    fig.tight_layout()
    fig.savefig(outpng, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap_z(expr_sub: pd.DataFrame, subtype: pd.Series, outpng: str):
    """Rows = genes, columns = samples; z-score per gene; subtype top bar."""
    genes = expr_sub.index.tolist()
    Z = expr_sub.copy()
    Z = Z.apply(lambda v: (v - np.nanmean(v))/ (np.nanstd(v) if np.nanstd(v)>0 else 1.0), axis=1)
    # order columns by subtype then by average expression to stabilize layout
    order = list(subtype.dropna().sort_values().index)
    Z = Z.loc[:, [c for c in order if c in Z.columns]]

    fig = plt.figure(figsize=(8.5, 3.3), dpi=220)
    ax = plt.axes([0.08, 0.25, 0.86, 0.7])
    im = ax.imshow(Z.values, aspect="auto", interpolation="nearest")
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes)
    ax.set_xticks(np.arange(len(Z.columns)))
    ax.set_xticklabels(Z.columns, rotation=90, fontsize=6)
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Z-score")

    # subtype bar
    ax2 = plt.axes([0.08, 0.18, 0.86, 0.03])
    # Create color bar: ETP=1 (darker), nonETP=0 (lighter) for Blues colormap
    bar = np.array([1 if subtype.get(c)=="ETP" else 0 for c in Z.columns])[None, :]
    im2 = ax2.imshow(bar, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("Subtype: dark blue=ETP, light blue=nonETP", fontsize=8)
    fig.savefig(outpng, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap_etp_vs_nonetp(expr: pd.DataFrame, subtype: pd.Series, genes: List[str], outpng: str):
    """
    Publication-quality clustered heatmap similar to Figure 1E style.
    Features: hierarchical clustering, dendrograms, annotation bars, z-score normalization.
    """
    import seaborn as sns
    from scipy.cluster import hierarchy
    
    # Filter to available genes
    present_genes = [g for g in genes if g in expr.index]
    if not present_genes:
        warnings.warn("No genes available for heatmap")
        return
    
    # Extract expression matrix (genes × samples)
    expr_sub = expr.loc[present_genes]
    
    # Z-score normalize per gene (row) across all samples
    Z = expr_sub.apply(lambda row: (row - np.nanmean(row)) / (np.nanstd(row) if np.nanstd(row) > 0 else 1.0), axis=1)
    Z = Z.T  # Transpose to samples × genes for seaborn clustermap
    
    # Create color mapping for subtypes
    subtype_aligned = subtype.reindex(Z.index)
    lut = {"ETP": "#d62728", "nonETP": "#1f77b4"}  # Red for ETP, Blue for nonETP
    row_colors = subtype_aligned.map(lut)
    
    # Create the clustered heatmap
    g = sns.clustermap(
        Z.T,  # Transpose back to genes × samples
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        col_colors=row_colors,
        row_cluster=True,  # Cluster genes
        col_cluster=True,  # Cluster samples
        figsize=(12, 6),
        dendrogram_ratio=0.15,
        colors_ratio=0.03,
        cbar_pos=(0.02, 0.83, 0.03, 0.15),
        linewidths=0,
        cbar_kws={
            "label": "Row Z-Score",
            "orientation": "vertical",
            "ticks": [-2, -1, 0, 1, 2]
        },
        yticklabels=True,
        xticklabels=False,
        method='average',
        metric='correlation'
    )
    
    # Customize the plot
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("", fontsize=11, fontweight='bold')
    
    # Add title
    g.fig.suptitle('MYCN and PP2A Subunit Expression: ETP-ALL vs nonETP-ALL', 
                   fontsize=13, fontweight='bold', y=0.98)
    
    # Improve gene labels
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), 
                                  fontsize=10, fontweight='bold')
    
    # Add subtype legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label=f'ETP (n={int((subtype_aligned == "ETP").sum())})'),
        Patch(facecolor='#1f77b4', label=f'nonETP (n={int((subtype_aligned == "nonETP").sum())})')
    ]
    g.ax_heatmap.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(1.15, 1),
        frameon=True,
        fontsize=9,
        title='Subtype',
        title_fontsize=10
    )
    
    # Add color key title
    g.ax_cbar.set_title('Color Key', fontsize=9, pad=10)
    
    plt.savefig(outpng, bbox_inches='tight', dpi=300)
    plt.close()

# -----------------------------
# PDF assembler
# -----------------------------
def build_summary_pdf(box_png: str, heat_png: str, scatter_pngs: List[str], caption: str, outpdf: str):
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(outpdf) as pdf:
        fig = plt.figure(figsize=(8.5, 11), dpi=200)

        # Heatmap on top
        ax_heat = plt.axes([0.08, 0.74, 0.84, 0.22])
        img = plt.imread(heat_png)
        ax_heat.imshow(img)
        ax_heat.axis("off")
        ax_heat.set_title("PP2A Core Heatmap (Z-score)", fontsize=12, pad=6)

        # Boxplot (left)
        ax_box = plt.axes([0.08, 0.48, 0.38, 0.2])
        ax_box.imshow(plt.imread(box_png))
        ax_box.axis("off")

        # Scatter grid (right/bottom)
        pos = [
            [0.52, 0.52, 0.4, 0.16],
            [0.52, 0.34, 0.4, 0.16],
            [0.52, 0.16, 0.4, 0.16],
            [0.08, 0.28, 0.38, 0.16],
            [0.08, 0.10, 0.38, 0.16],
        ]
        for png, p in zip(scatter_pngs, pos):
            ax = plt.axes(p)
            ax.imshow(plt.imread(png))
            ax.axis("off")

        # caption
        ax_cap = plt.axes([0.08, 0.02, 0.84, 0.06])
        ax_cap.axis("off")
        ax_cap.text(0, 0.9, "ETP-ALL: MYCN–PP2A summary (GSE28703)", fontsize=11, weight="bold")
        ax_cap.text(0, 0.1, caption, fontsize=9, va="top", wrap=True)
        pdf.savefig(fig)
        plt.close(fig)

# -----------------------------
# Main
# -----------------------------
def main():
    print("[INFO] Reading series matrix…")
    expr, meta = read_series_matrix(SERIES_MATRIX)
    print(f"[INFO] matrix: {expr.shape[0]} rows x {expr.shape[1]} samples")

    subtype = load_subtypes(meta)  # 'ETP' / 'nonETP' / NaN
    subtype = subtype.reindex(expr.columns)
    n_etp = int((subtype=="ETP").sum())
    n_non = int((subtype=="nonETP").sum())
    if n_etp==0 or n_non==0:
        warnings.warn(f"Subtypes appear imbalanced: ETP={n_etp}, nonETP={n_non}. Check sample_labels.csv.")
    else:
        print(f"[INFO] Subtype counts: ETP={n_etp}, nonETP={n_non}")

    # Ensure required genes exist (case-sensitive symbols)
    need = [ANCHOR_GENE] + PP2A_GENES
    available = [g for g in need if g in expr.index]
    missing   = [g for g in need if g not in expr.index]
    if missing:
        warnings.warn(f"Missing genes in matrix: {missing}")
    genes = available

    # Split by subtype
    etp_expr, non_expr = split_by_subtype(expr.loc[genes], subtype)

    # --- 1) MYCN boxplot
    p_box = np.nan
    if ANCHOR_GENE in genes:
        p_box = plot_box_mycn(etp_expr.loc[ANCHOR_GENE], non_expr.loc[ANCHOR_GENE],
                              os.path.join(BOX_DIR, "MYCN_box_ETP_vs_nonETP.png"))
    else:
        warnings.warn("MYCN not found; boxplot skipped.")
    
    # --- 1b) PP2A subunit boxplots for key genes showing differential expression
    pp2a_boxplot_genes = ["PPP2R1B", "PPP2R5D", "PPP2R2D", "PPP2R5C"]
    for gene in pp2a_boxplot_genes:
        if gene in genes:
            plot_boxplot_etp_vs_nonetp(etp_expr.loc[gene], non_expr.loc[gene], gene,
                                       os.path.join(BOX_DIR, f"{gene}_box_ETP_vs_nonETP.png"))
            print(f"[OK] Generated boxplot for {gene}")
        else:
            print(f"[WARN] {gene} not found; boxplot skipped.")

    # --- 2) Correlations: MYCN vs each PP2A subunit
    records = []
    scatter_for_pdf = []
    for gene in PP2A_GENES:
        if gene not in genes or ANCHOR_GENE not in genes:
            continue
        # ETP
        stats_etp = pearson_spearman(etp_expr.loc[ANCHOR_GENE].values, etp_expr.loc[gene].values)
        # nonETP
        stats_non = pearson_spearman(non_expr.loc[ANCHOR_GENE].values, non_expr.loc[gene].values)

        # colored combined (for intuition)
        series_anchor = expr.loc[ANCHOR_GENE]
        series_gene   = expr.loc[gene]
        comb_df = pd.DataFrame({
            "x": series_anchor,
            "y": series_gene,
            "subtype": subtype
        }).dropna()
        fig, ax = plt.subplots(figsize=(4,4), dpi=200)
        for lab, c in [("ETP","#d62728"), ("nonETP","#1f77b4")]:
            sub = comb_df[comb_df["subtype"]==lab]
            ax.scatter(sub["x"], sub["y"], s=18, label=lab, alpha=0.85)
        if len(comb_df)>=3:
            slope, intercept, r, p, _ = stats.linregress(comb_df["x"], comb_df["y"])
            xs = np.linspace(comb_df["x"].min(), comb_df["x"].max(), 100)
            ax.plot(xs, slope*xs + intercept, linewidth=2)
            ax.text(0.02, 0.96, f"All R={r:.2f}\np={p:.2g}", transform=ax.transAxes,
                    ha="left", va="top",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7"))
        ax.set_xlabel("MYCN")
        ax.set_ylabel(gene)
        ax.set_title(f"{gene} vs MYCN — colored by subtype")
        ax.grid(True, alpha=0.2, linestyle=":")
        ax.legend(frameon=True)
        out_col = os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_colored.png")
        fig.tight_layout(); fig.savefig(out_col, bbox_inches="tight"); plt.close(fig)
        scatter_for_pdf.append(out_col)

        # separate ETP
        plot_scatter(etp_expr.loc[ANCHOR_GENE], etp_expr.loc[gene],
                     f"{gene} vs MYCN — ETP", os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_ETP.png"),
                     r_val=stats_etp['pearson_r'], p_val=stats_etp['pearson_p'])
        # separate nonETP
        plot_scatter(non_expr.loc[ANCHOR_GENE], non_expr.loc[gene],
                     f"{gene} vs MYCN — nonETP", os.path.join(SCAT_DIR, f"{gene}_vs_MYCN_nonETP.png"),
                     r_val=stats_non['pearson_r'], p_val=stats_non['pearson_p'])

        records.append(dict(
            gene=gene, group="ETP", **stats_etp
        ))
        records.append(dict(
            gene=gene, group="nonETP", **stats_non
        ))

    stats_df = pd.DataFrame.from_records(records)

    # --- BH-FDR per group (ETP, nonETP) over Pearson p-values
    if stats_df.empty:
        print("[WARN] No correlation rows generated (likely due to missing symbols). Skipping FDR and PDF.")
        stats_df = pd.DataFrame(columns=["gene","group","pearson_r","pearson_p","spearman_r","spearman_p","n"])
    else:
        fdr_frames = []
        for grp in ["ETP","nonETP"]:
            sub = stats_df[(stats_df["group"]==grp)].copy()
            if len(sub):
                m = np.isfinite(sub["pearson_p"])
                if m.sum()>0:
                    rej, p_adj, *_ = multipletests(sub.loc[m, "pearson_p"].values, alpha=0.05, method="fdr_bh")
                    sub.loc[m, "pearson_fdr"] = p_adj
                    sub.loc[m, "pearson_sig"] = rej
            fdr_frames.append(sub)
        stats_df = pd.concat(fdr_frames, ignore_index=True)

    # --- 3) Enhanced heatmap: MYCN + PP2A subunits (ETP vs nonETP comparison)
    heatmap_genes = [ANCHOR_GENE] + PP2A_GENES  # Include MYCN
    heatmap_genes = [g for g in heatmap_genes if g in genes]
    if heatmap_genes:
        plot_heatmap_etp_vs_nonetp(expr, subtype, heatmap_genes, 
                                   os.path.join(HEAT_DIR, "MYCN_PP2A_ETP_vs_nonETP_heatmap.png"))
        # Also create the old version for backward compatibility
        plot_heatmap_z(expr.loc[[g for g in PP2A_GENES if g in genes]], subtype, 
                      os.path.join(HEAT_DIR, "PP2A_core_heatmap.png"))

    # --- Save table: add MYCN box p to caption info
    stats_df.to_csv(os.path.join(OUT_DIR, "MYCN_PP2A_results.csv"), index=False)

    # --- Compose caption
    cap_lines = []
    cap_lines.append(f"MYCN ETP vs nonETP Welch p = {p_box:.3g} (effect size d shown on figure).")
    for gene in PP2A_GENES:
        row_etp = stats_df[(stats_df["gene"]==gene) & (stats_df["group"]=="ETP")]
        row_non = stats_df[(stats_df["gene"]==gene) & (stats_df["group"]=="nonETP")]
        if len(row_etp):
            pe, pef = row_etp["pearson_p"].values[0], row_etp.get("pearson_fdr", pd.Series([np.nan])).values[0]
            cap_lines.append(f"{gene} (ETP): Pearson p={pe:.3g}, FDR={pef:.3g}")
        if len(row_non):
            pn, pnf = row_non["pearson_p"].values[0], row_non.get("pearson_fdr", pd.Series([np.nan])).values[0]
            cap_lines.append(f"{gene} (nonETP): Pearson p={pn:.3g}, FDR={pnf:.3g}")
    caption = " | ".join(cap_lines)

    # --- Assemble PDF
    pdf_out = os.path.join(OUT_DIR, "ETP_Mycn_PP2A_Summary.pdf")
    heat_png = os.path.join(HEAT_DIR, "PP2A_core_heatmap.png")
    box_png  = os.path.join(BOX_DIR, "MYCN_box_ETP_vs_nonETP.png")
    # pick up to 5 colored scatters (one per gene)
    scatter_pngs = [os.path.join(SCAT_DIR, f"{g}_vs_MYCN_colored.png") for g in PP2A_GENES if os.path.exists(os.path.join(SCAT_DIR, f"{g}_vs_MYCN_colored.png"))]
    scatter_pngs = scatter_pngs[:5]
    if os.path.exists(box_png) and os.path.exists(heat_png) and len(scatter_pngs)>0:
        build_summary_pdf(box_png, heat_png, scatter_pngs, caption, pdf_out)
        print(f"[DONE] Summary PDF: {pdf_out}")
    else:
        print("[WARN] Could not build summary PDF (missing figures).")
    print("[DONE] Outputs written to ./outputs")

if __name__ == "__main__":
    main()
