# analysis_gse104786.py
# -*- coding: utf-8 -*-

"""
GSE104786 analysis (GPL5175):
Compare MYCN expression between:
  - 16 primary SCPCs (Small Cell Prostate Cancer = NEPC)
  - 16 primary high-grade adenocarcinomas

If MYCN shows significant difference, proceed with PP2A subunit analysis.

Generates:
  1) Boxplot: MYCN expression (SCPC vs Adenocarcinoma) with Welch t-test
  2) If significant: Correlations and PP2A subunit analysis
  3) Summary statistics CSV
  4) One-page PDF summary

Inputs:
  - Data/NEPC/GSE104786_series_matrix.txt.gz
  - Data/NEPC/sample_labels_gse104786.csv (sample_id, subtype)
  - Data/NEPC/GSE104786_family.soft.gz (for GPL mapping)

Outputs:
  ./outputs_gse104786/
    ├── boxplots/MYCN_box_SCPC_vs_Adenocarcinoma.png
    ├── MYCN_results.csv
    └── GSE104786_MYCN_Summary.pdf
"""

import os, re, io, gzip, textwrap, warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Config
# -----------------------------
PROJECT_ROOT = "."
SERIES_MATRIX = os.path.join(PROJECT_ROOT, "Data/NEPC/GSE104786_series_matrix.txt.gz")

OUT_DIR  = os.path.join(PROJECT_ROOT, "outputs_gse104786")
BOX_DIR  = os.path.join(OUT_DIR, "boxplots")
os.makedirs(BOX_DIR, exist_ok=True)

ANCHOR_GENE = "MYCN"
PP2A_GENES  = ["PPP2R5C", "PPP2R5D", "PPP2R1B", "PPP2CA", "PPP2R2D"]

GROUP_A = "SCPC"
GROUP_B = "Adenocarcinoma"


# -----------------------------
# Utilities
# -----------------------------
def load_gpl_map(gpl_soft_path: str) -> dict:
    """Return dict: probe_id -> gene_symbol. For GPL5175, extract from description field."""
    if not os.path.exists(gpl_soft_path):
        return {}
    probe_to_symbol = {}
    try:
        opener = gzip.open if gpl_soft_path.endswith(".gz") else open
        with opener(gpl_soft_path, "rt", encoding="utf-8", errors="ignore") as f:
            # For GPL5175, data starts after platform section
            # Look for lines that start with numeric probe IDs (not ^ or !)
            for line in f:
                # Skip header lines
                if line.startswith("^") or line.startswith("!"):
                    continue
                if not line.strip():
                    continue
                # Check if line starts with a numeric probe ID
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    probe_id = parts[0].strip()
                    # Check if first field is numeric (probe ID)
                    if probe_id.isdigit():
                        # Look for gene symbol in any field with "// GENE_SYMBOL //" pattern
                        gene_symbol = None
                        for part in parts:
                            if "//" in part:
                                # Pattern: "NM_XXXXX // GENE_SYMBOL // description"
                                # Split by // and take the middle part
                                split_parts = [p.strip() for p in part.split("//")]
                                if len(split_parts) >= 3:
                                    # Middle part is the gene symbol
                                    gene_symbol = split_parts[1]
                                    # Clean up - take first word if multiple
                                    gene_symbol = re.split(r"[;/,\s]+", gene_symbol)[0]
                                    break
                        if probe_id and gene_symbol and gene_symbol not in {"---", "NA", ""}:
                            if gene_symbol and len(gene_symbol) > 1:  # Valid gene symbol
                                probe_to_symbol[probe_id] = gene_symbol
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
    characteristics_ch1_lines = []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if raw.startswith("!Sample_"):
                parts = raw.strip().split("\t")
                tag = parts[0][len("!Sample_"):].lower()
                vals = [v.strip('"') for v in parts[1:]]
                
                # Special handling for characteristics_ch1 (multiple lines)
                if tag == "characteristics_ch1":
                    characteristics_ch1_lines.append(vals)
                else:
                    for gsm, val in zip(gsm_ids, vals):
                        sample_meta[gsm][tag] = val
            if raw.startswith("!series_matrix_table_begin"):
                break

        # Find annotation line
        annotation_line = None
        for char_line in characteristics_ch1_lines:
            if any('annotation' in v.lower() for v in char_line):
                annotation_line = char_line
                break
        if annotation_line:
            for gsm, ann in zip(gsm_ids, annotation_line):
                sample_meta[gsm]["annotation"] = ann

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

    # Map probe IDs to gene symbols using GPL annotation
    gpl_try = [
        os.path.join(PROJECT_ROOT, "Data/NEPC/GSE104786_family.soft.gz"),
        os.path.join(PROJECT_ROOT, "Data/NEPC/GPL5175_family.soft.gz"),
        os.path.join(PROJECT_ROOT, "GPL5175_family.soft.gz"),
    ]
    probe_map = {}
    for gpl_path in gpl_try:
        if os.path.exists(gpl_path):
            probe_map = load_gpl_map(gpl_path)
            if probe_map:
                print(f"[INFO] Mapped {len(probe_map)} probes to gene symbols")
                break
    if not probe_map:
        warnings.warn("GPL5175 not found: rows likely probe IDs. Consider placing GPL5175_family.soft.gz in Data/NEPC/.")

    if probe_map:
        # Map probe IDs to gene symbols
        idx_map = {}
        for probe_id in expr.index:
            symbol = probe_map.get(probe_id, "")
            if symbol:
                if symbol not in idx_map:
                    idx_map[symbol] = []
                idx_map[symbol].append(probe_id)
        
        # Aggregate multiple probes per gene (mean)
        new_expr = {}
        for symbol, probes in idx_map.items():
            new_expr[symbol] = expr.loc[probes].mean(axis=0)
        expr = pd.DataFrame(new_expr).T
        print(f"[INFO] Mapped to {len(expr)} genes")
    else:
        # Check if IDs already look like gene symbols
        sample_ids = list(expr.index[:min(100, len(expr))])
        looks_like_symbol = sum(1 for sid in sample_ids if len(sid) > 2 and any(c.isalpha() for c in sid))
        if looks_like_symbol / len(sample_ids) > 0.5:
            print(f"[INFO] Row IDs appear to be gene symbols already")
            # Try to find MYCN
            mycn_variants = ['MYCN', 'Mycn', 'mycn']
            for variant in mycn_variants:
                if variant in expr.index:
                    print(f"[INFO] Found {variant} in expression matrix")
                    break

    return expr, sample_meta


def load_subtypes_gse104786(sample_meta: Dict[str, Dict[str,str]]) -> pd.Series:
    """Load subtypes from sample_labels_gse104786.csv or infer from metadata."""
    # 1) explicit labels
    labfile = os.path.join(PROJECT_ROOT, "Data/NEPC/sample_labels_gse104786.csv")
    if os.path.exists(labfile):
        df = pd.read_csv(labfile)
        df.columns = [c.lower() for c in df.columns]
        mapping = {r["sample_id"]: r["subtype"] for _, r in df.iterrows()}
        ser = pd.Series({gsm: mapping.get(gsm, np.nan) for gsm in sample_meta})
        return ser
    
    # 2) heuristic from metadata
    def parse_one(meta: Dict[str,str]) -> str:
        annotation = (meta.get("annotation") or "").lower()
        if re.search(r"\b(small[\s-]*cell[\s-]*carcinoma|small[\s-]*cell)\b", annotation):
            return GROUP_A
        if re.search(r"\b(adenocarcinoma|adeno)\b", annotation) and "small cell" not in annotation:
            return GROUP_B
        return np.nan
    
    return pd.Series({gsm: parse_one(m) for gsm, m in sample_meta.items()})


def split_by_subtype(expr: pd.DataFrame, subtype: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame]:
    subtype = subtype.reindex(expr.columns)
    a_cols = list(subtype[subtype==GROUP_A].index)
    b_cols = list(subtype[subtype==GROUP_B].index)
    return expr[a_cols], expr[b_cols]


def welch_ttest(a: np.ndarray, b: np.ndarray) -> Tuple[float, float]:
    """Welch's t-test (unequal variances). Returns (t_stat, p_value)."""
    a_clean = a[np.isfinite(a)]
    b_clean = b[np.isfinite(b)]
    if len(a_clean) < 2 or len(b_clean) < 2:
        return np.nan, np.nan
    t, p = stats.ttest_ind(a_clean, b_clean, equal_var=False)
    return t, p


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's d effect size."""
    a_clean = a[np.isfinite(a)]
    b_clean = b[np.isfinite(b)]
    if len(a_clean) < 2 or len(b_clean) < 2:
        return np.nan
    n1, n2 = len(a_clean), len(b_clean)
    var1, var2 = np.var(a_clean, ddof=1), np.var(b_clean, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    if pooled_std == 0:
        return np.nan
    return (np.mean(a_clean) - np.mean(b_clean)) / pooled_std


# -----------------------------
# Plotting
# -----------------------------
def plot_boxplot_groups(a: pd.Series, b: pd.Series, gene_name: str, outpng: str):
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    
    # Convert log2 values to normal scale (2^log2_value)
    a_normal = np.power(2, a.values)
    b_normal = np.power(2, b.values)
    data = [a_normal, b_normal]  # order: SCPC, Adenocarcinoma
    
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
    a_finite = a_normal[np.isfinite(a_normal)]
    b_finite = b_normal[np.isfinite(b_normal)]
    # SCPC: black points
    if len(a_finite) > 0:
        ax.scatter(np.random.normal(1, 0.04, size=len(a_finite)), a_finite, 
                   s=30, color='black', alpha=0.6, edgecolors='black', linewidths=0.5, zorder=4)
    # Adenocarcinoma: black points
    if len(b_finite) > 0:
        ax.scatter(np.random.normal(2, 0.04, size=len(b_finite)), b_finite, 
                   s=30, color='black', alpha=0.6, edgecolors='black', linewidths=0.5, zorder=3)
    ax.set_ylabel(f"{gene_name} expression", fontsize=12, fontweight='bold')

    _, p = welch_ttest(a.values, b.values)
    d = cohens_d(a.values, b.values)
    p_text = f"p<0.001 (p={p:.2e})" if p < 1e-3 else (f"p={p:.4f}" if p < 0.01 else f"p={p:.3f}")
    ax.set_title(f"{gene_name}: {GROUP_A} vs {GROUP_B}\nWelch {p_text}, d={d:.2f}",
                 fontsize=13, fontweight='bold')
    # Show actual number of finite values plotted (excludes NaN)
    n_a_plot = len(a_finite) if len(a_finite) > 0 else 0
    n_b_plot = len(b_finite) if len(b_finite) > 0 else 0
    ax.set_xticklabels([f"{GROUP_A}\n(n={n_a_plot})", f"{GROUP_B}\n(n={n_b_plot})"], fontsize=11)

    # bracket + sig stars (use normal scale for positioning)
    y_max, y_min = np.max([np.max(a_normal), np.max(b_normal)]), np.min([np.min(a_normal), np.min(b_normal)])
    y_range = y_max - y_min
    bracket_y = y_max + y_range * 0.05
    ax.plot([1,1,2,2], [bracket_y, bracket_y + y_range*0.02, bracket_y + y_range*0.02, bracket_y], color='black', lw=1.5)
    sig = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
    ax.text(1.5, bracket_y + y_range*0.03, sig, ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.grid(True, alpha=0.2, linestyle=":")
    fig.tight_layout()
    fig.savefig(outpng, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return p


# -----------------------------
# Main
# -----------------------------
def main():
    print("[INFO] Reading series matrix…")
    expr, meta = read_series_matrix(SERIES_MATRIX)
    print(f"[INFO] matrix: {expr.shape[0]} rows x {expr.shape[1]} samples")

    subtype = load_subtypes_gse104786(meta)  # 'SCPC' / 'Adenocarcinoma' / NaN
    subtype = subtype.reindex(expr.columns)
    n_a = int((subtype==GROUP_A).sum()); n_b = int((subtype==GROUP_B).sum())
    if n_a==0 or n_b==0:
        warnings.warn(f"Subtype counts look imbalanced: {GROUP_A}={n_a}, {GROUP_B}={n_b}. "
                      f"If needed, add Data/NEPC/sample_labels_gse104786.csv")
    else:
        print(f"[INFO] Subtype counts: {GROUP_A}={n_a}, {GROUP_B}={n_b}")

    # Check if MYCN is available
    if ANCHOR_GENE not in expr.index:
        raise RuntimeError(f"{ANCHOR_GENE} not found in expression matrix. Check GPL mapping.")

    # Split
    a_expr, b_expr = split_by_subtype(expr.loc[[ANCHOR_GENE]], subtype)  # SCPC, Adenocarcinoma

    # MYCN boxplot
    print(f"\n[INFO] Comparing {ANCHOR_GENE} expression: {GROUP_A} vs {GROUP_B}")
    p_mycn = plot_boxplot_groups(a_expr.loc[ANCHOR_GENE], b_expr.loc[ANCHOR_GENE], ANCHOR_GENE,
                                 os.path.join(BOX_DIR, f"{ANCHOR_GENE}_box_{GROUP_A}_vs_{GROUP_B}.png"))
    print(f"[OK] Generated boxplot for {ANCHOR_GENE}")

    # Calculate statistics
    scpc_vals = a_expr.loc[ANCHOR_GENE].values
    adeno_vals = b_expr.loc[ANCHOR_GENE].values
    scpc_mean = np.nanmean(scpc_vals)
    adeno_mean = np.nanmean(adeno_vals)
    scpc_std = np.nanstd(scpc_vals)
    adeno_std = np.nanstd(adeno_vals)
    _, p_val = welch_ttest(scpc_vals, adeno_vals)
    d_val = cohens_d(scpc_vals, adeno_vals)

    # Save results
    results = pd.DataFrame({
        "Gene": [ANCHOR_GENE],
        "Group_A": [GROUP_A],
        "Group_B": [GROUP_B],
        "n_A": [n_a],
        "n_B": [n_b],
        "Mean_A": [scpc_mean],
        "Mean_B": [adeno_mean],
        "Std_A": [scpc_std],
        "Std_B": [adeno_std],
        "Welch_p": [p_val],
        "Cohens_d": [d_val],
        "Significant": [p_val < 0.05 if not np.isnan(p_val) else False]
    })
    
    out_csv = os.path.join(OUT_DIR, f"{ANCHOR_GENE}_results.csv")
    results.to_csv(out_csv, index=False)
    print(f"[OK] Results saved to {out_csv}")

    # Summary
    print(f"\n{'='*60}")
    print(f"MYCN Expression Comparison Results")
    print(f"{'='*60}")
    print(f"{GROUP_A} (n={n_a}): Mean={scpc_mean:.3f}, Std={scpc_std:.3f}")
    print(f"{GROUP_B} (n={n_b}): Mean={adeno_mean:.3f}, Std={adeno_std:.3f}")
    print(f"Welch's t-test: p={p_val:.4f}")
    print(f"Cohen's d: {d_val:.3f}")
    if p_val < 0.05:
        print(f"\n✅ SIGNIFICANT DIFFERENCE (p < 0.05)")
        print(f"   Proceed with PP2A subunit analysis!")
    else:
        print(f"\n⚠️  No significant difference (p >= 0.05)")
        print(f"   MYCN expression is not significantly different between groups.")
    print(f"{'='*60}\n")

    print(f"[DONE] Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()

