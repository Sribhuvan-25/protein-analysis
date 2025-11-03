# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# Pipeline: MYCN–PP2A correlation & heatmaps in ETP-ALL vs non-ETP-ALL (GSE28703)
# Inputs required in project root:
#   - GSE28703_series_matrix.txt.gz
#   - GPL13158.annot.gz
# Outputs written to ./outputs/
# """

# import gzip
# import io
# import os
# import re
# from typing import Dict, List, Tuple

# import numpy as np
# import pandas as pd
# from scipy import stats
# import matplotlib.pyplot as plt
# import seaborn as sns

# # ----------------------------
# # 0) Config & genes of interest
# # ----------------------------
# PROJECT_ROOT = "."
# SERIES_MATRIX = os.path.join(PROJECT_ROOT, "GSE28703_series_matrix.txt.gz")
# GPL_ANNOT = os.path.join(PROJECT_ROOT, "GPL13158-5065.txt")

# OUT_DIR = os.path.join(PROJECT_ROOT, "outputs")
# SCATTER_DIR = os.path.join(OUT_DIR, "scatter")
# HEATMAP_DIR = os.path.join(OUT_DIR, "heatmaps")
# os.makedirs(SCATTER_DIR, exist_ok=True)
# os.makedirs(HEATMAP_DIR, exist_ok=True)

# # Core genes (MYCN + PP2A set)
# core_genes = [
#     "MYCN", "PPP2R2D", "PPP2R5C", "PPP2R5D", "PPP2R1B", "PPP2CA"
# ]

# # Extended/set for biological context (feel free to add)
# context_genes = [
#     "SET", "KIAA1524",  # CIP2A (synonym: KIAA1524)
#     "EZH2", "SOX2", "AR", "SYP"
# ]

# extended_genes = list(dict.fromkeys(core_genes + context_genes))  # unique preserve order

# RANDOM_SEED = 7
# np.random.seed(RANDOM_SEED)

# # ----------------------------
# # 1) Utilities
# # ----------------------------
# def read_gzip_text(path: str) -> str:
#     if not os.path.exists(path):
#         raise FileNotFoundError(f"Missing required file: {path}")
#     with gzip.open(path, "rb") as f:
#         return f.read().decode("utf-8", errors="replace")

# def read_plain_text(path: str) -> str:
#     """Read plain text file (not gzipped)."""
#     if not os.path.exists(path):
#         raise FileNotFoundError(f"Missing required file: {path}")
#     with open(path, "r", encoding="utf-8", errors="replace") as f:
#         return f.read()

# def parse_series_matrix(raw_text: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Returns:
#         expr_df: probe x sample expression (float)
#         meta_df: sample metadata table (wide)
#     """
#     # Extract table block
#     start = raw_text.find("!series_matrix_table_begin")
#     end = raw_text.find("!series_matrix_table_end")
#     if start == -1 or end == -1:
#         raise ValueError("Series matrix table markers not found.")
#     table_text = raw_text[start:end].splitlines()[1:]  # skip the begin line

#     # Read as TSV; first column is ID_REF (probe)
#     expr_df = pd.read_csv(io.StringIO("\n".join(table_text)), sep="\t", index_col=0)
#     expr_df.index.name = "ID_REF"
#     # Convert to numeric
#     expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

#     # Parse sample metadata: lines starting with !Sample_
#     meta_lines = [l for l in raw_text.splitlines() if l.startswith("!Sample_")]
#     # Build key -> list of values
#     meta_dict: Dict[str, List[str]] = {}
#     for line in meta_lines:
#         # Handle two formats:
#         # 1. Tab-separated: !Sample_title\t"value1"\t"value2"...
#         # 2. With equals: !Sample_title = value1 \t value2 ...
#         if "\t" in line:
#             # Tab-separated format
#             parts = line.split("\t")
#             if len(parts) < 2:
#                 continue
#             key = parts[0].replace("!Sample_", "")
#             values = parts[1:]  # All values after the field name
#             # Remove quotes if present
#             values = [v.strip('"') for v in values]
#             meta_dict[key] = values
#         elif " = " in line:
#             # Equals format (less common)
#             key, rest = line.split(" = ", 1)
#             values = rest.split("\t")
#             meta_dict.setdefault(key.replace("!Sample_", ""), []).extend(values)

#     # Build dataframe directly from meta_dict
#     meta_df = pd.DataFrame(meta_dict)
#     # Ensure columns align with expr_df columns (samples)
#     # Some keys may have different lengths; align by number of samples.
#     n_samples = expr_df.shape[1]
#     for c in meta_df.columns:
#         meta_df[c] = meta_df[c].astype(str)
#         if len(meta_df[c]) != n_samples:
#             # Pad or trim if necessary (rare)
#             vals = list(meta_df[c].values)
#             if len(vals) < n_samples:
#                 vals += [""] * (n_samples - len(vals))
#             else:
#                 vals = vals[:n_samples]
#             meta_df[c] = vals

#     meta_df.index = expr_df.columns  # index rows by sample id
#     meta_df.index.name = "sample_id"

#     return expr_df, meta_df

# def infer_etp_labels(meta_df: pd.DataFrame) -> pd.Series:
#     """
#     Infer 'ETP' vs 'non-ETP' from any metadata columns containing subtype.
#     Updated to search for "early T-cell precursor" pattern instead of just "ETP".
#     """
#     subtype = pd.Series(index=meta_df.index, dtype="object")

#     # Search all textual columns for ETP indicators
#     for col in meta_df.columns:
#         col_vals = meta_df[col].astype(str)
        
#         # Match "early T-cell precursor" (ETP samples)
#         # Pattern allows for variations: "early T-cell precursor", "early T cell precursor", etc.
#         is_etp = col_vals.str.contains(
#             r"early\s+T[-\s]?cell\s+precursor", 
#             case=False, 
#             na=False, 
#             regex=True
#         )
        
#         # Also check for explicit "non-ETP" or "non ETP" markers
#         is_non_etp = col_vals.str.contains(
#             r"non[-\s]?ETP", 
#             case=False, 
#             na=False, 
#             regex=True
#         )
        
#         # Mark ETP where found and not already set
#         subtype[is_etp & subtype.isna()] = "ETP"
#         # Mark nonETP explicitly
#         subtype[is_non_etp & subtype.isna()] = "nonETP"

#     # Everything else becomes nonETP by default (if unlabeled)
#     subtype = subtype.fillna("nonETP")

#     return subtype

# def load_probe_to_gene_map(gpl_text: str) -> pd.DataFrame:
#     """
#     Parse GPL annotation to extract 'ID' and 'Gene Symbol' columns.
#     """
#     # GPL files contain a table after a header starting with "^ANNOTATION_TABLE"
#     # Fallback: detect the section starting at a line with "ID\t"
#     lines = gpl_text.splitlines()
#     # Find first header line that starts with "ID\t"
#     header_idx = None
#     for i, line in enumerate(lines):
#         if line.startswith("ID\t"):
#             header_idx = i
#             break
#     if header_idx is None:
#         raise ValueError("Could not find GPL annotation table header (ID\\t...).")

#     table_text = "\n".join(lines[header_idx:])
#     annot = pd.read_csv(io.StringIO(table_text), sep="\t", dtype=str)
#     # Normalize column names
#     annot.columns = [c.strip() for c in annot.columns]
#     # Prefer common names
#     col_id = "ID"
#     # Try different variants for gene symbol column
#     symbol_col_candidates = ["Gene Symbol", "Gene Symbol [gene_assignment]", "Gene Symbol "]
#     symbol_col = None
#     for c in symbol_col_candidates:
#         if c in annot.columns:
#             symbol_col = c
#             break
#     if symbol_col is None:
#         # Fallback: look for any column containing 'Gene Symbol'
#         matches = [c for c in annot.columns if "Gene Symbol" in c]
#         if matches:
#             symbol_col = matches[0]
#         else:
#             raise ValueError("Could not find a 'Gene Symbol' column in GPL annotation.")

#     df = annot[[col_id, symbol_col]].rename(columns={col_id: "ID_REF", symbol_col: "symbol"})
#     # Clean symbol: take first if multiple separated by '///'
#     df["symbol"] = df["symbol"].astype(str).str.split("///").str[0].str.strip()
#     df = df.replace({"": np.nan})
#     df = df.dropna(subset=["symbol"]).drop_duplicates()
#     return df

# def collapse_probes_to_genes(expr_df: pd.DataFrame, p2g: pd.DataFrame) -> pd.DataFrame:
#     """
#     expr_df: probe x sample
#     p2g: columns [ID_REF, symbol]
#     Returns gene x sample, collapsed by median across probes per gene.
#     """
#     merged = p2g.merge(expr_df, left_on="ID_REF", right_index=True, how="inner")
#     # Group by symbol and take median across probes
#     gene_expr = merged.drop(columns=["ID_REF"]).groupby("symbol", as_index=True).median(numeric_only=True)
#     return gene_expr

# def fisher_z_test(r1: float, n1: int, r2: float, n2: int) -> Tuple[float, float]:
#     """
#     Fisher z-test for difference between two independent correlations (Pearson).
#     """
#     # Guard against r at exactly +/-1
#     r1 = np.clip(r1, -0.999999, 0.999999)
#     r2 = np.clip(r2, -0.999999, 0.999999)
#     z1 = np.arctanh(r1)
#     z2 = np.arctanh(r2)
#     z = (z1 - z2) / np.sqrt(1/(n1-3) + 1/(n2-3))
#     p = 2 * stats.norm.sf(abs(z))
#     return float(z), float(p)

# def compute_correlations(df: pd.DataFrame, labels: pd.Series, genes: List[str]) -> pd.DataFrame:
#     """
#     df: sample x gene (includes MYCN and target genes)
#     labels: sample -> 'ETP' or 'nonETP'
#     genes: list of target genes (including PP2A & context genes, NOT including MYCN)
#     """
#     assert "MYCN" in df.columns, "MYCN not present in the expression matrix."

#     # Split groups
#     etp_mask = labels == "ETP"
#     non_mask = labels == "nonETP"
#     df_etp = df.loc[etp_mask]
#     df_non = df.loc[non_mask]

#     rows = []
#     for g in genes:
#         if g not in df.columns:
#             rows.append({
#                 "gene": g,
#                 "r_pearson_etp": np.nan, "p_pearson_etp": np.nan,
#                 "r_spearman_etp": np.nan, "p_spearman_etp": np.nan,
#                 "r_pearson_non": np.nan, "p_pearson_non": np.nan,
#                 "r_spearman_non": np.nan, "p_spearman_non": np.nan,
#                 "z_diff": np.nan, "p_diff": np.nan,
#                 "notes": "missing gene"
#             })
#             continue

#         # ETP
#         rpe_e, ppe_e = stats.pearsonr(df_etp["MYCN"], df_etp[g]) if len(df_etp) > 3 else (np.nan, np.nan)
#         rsp_e, psp_e = stats.spearmanr(df_etp["MYCN"], df_etp[g]) if len(df_etp) > 3 else (np.nan, np.nan)
#         # nonETP
#         rpe_n, ppe_n = stats.pearsonr(df_non["MYCN"], df_non[g]) if len(df_non) > 3 else (np.nan, np.nan)
#         rsp_n, psp_n = stats.spearmanr(df_non["MYCN"], df_non[g]) if len(df_non) > 3 else (np.nan, np.nan)

#         # Fisher z
#         if np.isfinite(rpe_e) and np.isfinite(rpe_n):
#             z, pz = fisher_z_test(rpe_e, len(df_etp), rpe_n, len(df_non))
#         else:
#             z, pz = (np.nan, np.nan)

#         rows.append({
#             "gene": g,
#             "r_pearson_etp": rpe_e, "p_pearson_etp": ppe_e,
#             "r_spearman_etp": rsp_e, "p_spearman_etp": psp_e,
#             "r_pearson_non": rpe_n, "p_pearson_non": ppe_n,
#             "r_spearman_non": rsp_n, "p_spearman_non": psp_n,
#             "z_diff": z, "p_diff": pz,
#             "notes": ""
#         })

#     out = pd.DataFrame(rows)
#     return out

# def annotate_stats(ax: plt.Axes, text: str, loc: str = "upper left", pad: float = 0.02):
#     """
#     Place a small text box with statistics on a plot.
#     """
#     align = {"upper left": (0.02, 0.98, "left", "top"),
#              "upper right": (0.98, 0.98, "right", "top"),
#              "lower left": (0.02, 0.02, "left", "bottom"),
#              "lower right": (0.98, 0.02, "right", "bottom")}
#     x, y, ha, va = align.get(loc, align["upper left"])
#     ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
#             bbox=dict(boxstyle="round", alpha=0.2), fontsize=10)

# def scatter_reg_line(ax, x, y):
#     # Simple linear regression line & CI via seaborn regplot drawn onto given axis (no facet)
#     sns.regplot(x=x, y=y, ax=ax, scatter=False, ci=95, truncate=True)

# def plot_scatter_separate(df: pd.DataFrame, labels: pd.Series, gene: str, outdir: str):
#     """
#     Two separate figures: ETP-only and nonETP-only, with R and p on each.
#     """
#     for group in ["ETP", "nonETP"]:
#         mask = labels == group
#         sub = df.loc[mask, ["MYCN", gene]].dropna()
#         fig, ax = plt.subplots(figsize=(4, 4))
#         ax.scatter(sub["MYCN"], sub[gene], s=25)
#         if len(sub) > 3:
#             r, p = stats.pearsonr(sub["MYCN"], sub[gene])
#             annotate_stats(ax, f"R = {r:.2f}\np = {p:.2e}", loc="upper left")
#             scatter_reg_line(ax, sub["MYCN"], sub[gene])
#         ax.set_xlabel("MYCN")
#         ax.set_ylabel(gene)
#         ax.set_title(f"{gene} vs MYCN — {group}")
#         fig.tight_layout()
#         fig.savefig(os.path.join(outdir, f"{gene}_vs_MYCN_{group}.png"), dpi=300)
#         plt.close(fig)

# def plot_scatter_colored(df: pd.DataFrame, labels: pd.Series, gene: str, outdir: str):
#     """
#     Single figure: all samples, colored by subtype, with group-wise regressions.
#     """
#     sub = df[["MYCN", gene]].copy()
#     sub["Subtype"] = labels.values
#     fig, ax = plt.subplots(figsize=(4.5, 4.5))
#     # Points
#     sns.scatterplot(data=sub, x="MYCN", y=gene, hue="Subtype", ax=ax)
#     # Regression lines per group
#     for group in ["ETP", "nonETP"]:
#         gdf = sub[sub["Subtype"] == group]
#         if len(gdf) > 3:
#             sns.regplot(data=gdf, x="MYCN", y=gene, ax=ax, scatter=False, ci=95, truncate=True, label=f"{group} fit")
#     # Stats text
#     stats_texts = []
#     for group in ["ETP", "nonETP"]:
#         gdf = sub[sub["Subtype"] == group]
#         if len(gdf) > 3:
#             r, p = stats.pearsonr(gdf["MYCN"], gdf[gene])
#             stats_texts.append(f"{group}: R={r:.2f}, p={p:.1e}")
#     if stats_texts:
#         annotate_stats(ax, "\n".join(stats_texts), loc="upper right")
#     ax.set_title(f"{gene} vs MYCN — colored by subtype")
#     fig.tight_layout()
#     fig.savefig(os.path.join(outdir, f"{gene}_vs_MYCN_colored.png"), dpi=300)
#     plt.close(fig)

# def zscore(df: pd.DataFrame) -> pd.DataFrame:
#     return (df - df.mean(axis=1, skipna=True).values.reshape(-1,1)) / df.std(axis=1, ddof=0, skipna=True).values.reshape(-1,1)

# def heatmap_genes(df: pd.DataFrame, labels: pd.Series, genes: List[str], outpath: str, title: str):
#     """
#     df: sample x gene
#     """
#     present = [g for g in genes if g in df.columns]
#     if not present:
#         raise ValueError("None of the requested genes are present for the heatmap.")
#     mat = df[present].T  # gene x sample
#     mat_z = zscore(mat)

#     # Column colors for subtype
#     subtype_colors = labels.map({"ETP": "#d62728", "nonETP": "#1f77b4"})  # red/blue
#     col_colors = pd.DataFrame({"Subtype": subtype_colors})
#     # Cluster genes (rows) but keep samples grouped by subtype order
#     # Reorder columns: ETP first, then nonETP
#     order = list(labels[labels=="ETP"].index) + list(labels[labels=="nonETP"].index)
#     mat_z = mat_z[order]
#     col_colors = col_colors.loc[order]

#     g = sns.clustermap(
#         mat_z,
#         row_cluster=True,
#         col_cluster=False,
#         col_colors=col_colors,
#         cmap="vlag",
#         figsize=(8, max(3, 0.35*len(present))),
#         dendrogram_ratio=0.15,
#         cbar_kws={"label": "Z-score"}
#     )
#     plt.suptitle(title, y=1.02, fontsize=12)
#     plt.savefig(outpath, dpi=300, bbox_inches="tight")
#     plt.close()

# # ----------------------------
# # 2) Main orchestration
# # ----------------------------
# def main():
#     # Step 6–8: Load series matrix and metadata
#     raw_series = read_gzip_text(SERIES_MATRIX)
#     expr_df, meta_df = parse_series_matrix(raw_series)

#     # Step 8: Infer ETP vs nonETP
#     subtype = infer_etp_labels(meta_df)

#     # Step 12–15: Probe -> gene mapping via GPL
#     gpl_text = read_plain_text(GPL_ANNOT)
#     p2g = load_probe_to_gene_map(gpl_text)
#     gene_expr = collapse_probes_to_genes(expr_df, p2g)  # gene x sample

#     # Step 16–17: sample x gene and join labels
#     samp_expr = gene_expr.T  # sample x gene
#     # Keep only samples in both expr and subtype
#     common = samp_expr.index.intersection(subtype.index)
#     samp_expr = samp_expr.loc[common].copy()
#     subtype_use = subtype.loc[common].copy()

#     # Report group sizes
#     n_etp = int((subtype_use == "ETP").sum())
#     n_non = int((subtype_use == "nonETP").sum())
#     print(f"[INFO] Samples — ETP: {n_etp}, nonETP: {n_non}, total: {len(subtype_use)}")

#     # Step 18: Check genes present
#     all_needed = list(dict.fromkeys(["MYCN"] + extended_genes))
#     missing = [g for g in all_needed if g not in samp_expr.columns]
#     if missing:
#         print(f"[WARN] Missing genes (not in platform or filtered): {', '.join(missing)}")

#     # Filter to available
#     core_avail = [g for g in core_genes if g in samp_expr.columns]
#     ext_avail = [g for g in extended_genes if g in samp_expr.columns]
#     if "MYCN" not in samp_expr.columns:
#         raise RuntimeError("MYCN not found — cannot compute correlations.")

#     # Step 20–23: Correlations
#     target_for_corr = [g for g in ext_avail if g != "MYCN"]
#     corr_table = compute_correlations(samp_expr, subtype_use, target_for_corr)
#     corr_csv = os.path.join(OUT_DIR, "MYCN_PP2A_correlation_stats.csv")
#     corr_table.to_csv(corr_csv, index=False)
#     print(f"[OK] Correlation table saved: {corr_csv}")

#     # Step 24–26: Scatter plots
#     for g in target_for_corr:
#         plot_scatter_separate(samp_expr, subtype_use, g, SCATTER_DIR)
#         plot_scatter_colored(samp_expr, subtype_use, g, SCATTER_DIR)
#     print(f"[OK] Scatter plots saved to: {SCATTER_DIR}")

#     # Step 27–29: Heatmaps
#     # Subset dataframe to available genes (+ MYCN) for heatmaps
#     for name, genelist in [("core", core_avail), ("extended", ext_avail)]:
#         if genelist:
#             outpath = os.path.join(HEATMAP_DIR, f"heatmap_{name}.png")
#             title = f"Heatmap ({name}) — Z-scored genes"
#             heatmap_genes(samp_expr, subtype_use, genelist, outpath, title)
#             print(f"[OK] Heatmap saved: {outpath}")
#         else:
#             print(f"[WARN] No genes available for {name} heatmap; skipped.")

#     # Step 30: Summary
#     print("[SUMMARY]")
#     print(f"  ETP samples: {n_etp}")
#     print(f"  nonETP samples: {n_non}")
#     if missing:
#         print(f"  Missing genes: {', '.join(missing)}")
#     print(f"  Outputs: {OUT_DIR}")

# if __name__ == "__main__":
#     main()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MYCN–PP2A analysis in ETP-ALL vs non-ETP-ALL (GSE28703 + GPL13158)

What this script does (end-to-end):
1) Load GSE28703 series matrix (expression by probe) and GPL13158 platform table (probe → gene).
2) Collapse probes to gene symbols (mean across probes mapping to same gene).
3) Infer sample subtype labels (ETP vs nonETP) from sample titles/characteristics.
4) Build a sample × gene matrix and keep ONLY the slide genes (Option A set).
5) Compute Pearson correlations of MYCN vs each PP2A subunit (for reference) and
   generate one colored scatter plot per gene (MYCN on x, target gene on y).
6) Generate one heatmap for the core genes (z-scored per gene across samples).
7) Generate ETP vs nonETP boxplots for each gene (including MYCN) + a CSV of
   stats (means, fold-change, Welch t-test p-values, group sizes).

Outputs:
- outputs/scatter/*_vs_MYCN_colored.png
- outputs/heatmaps/heatmap_core.png
- outputs/boxplots/*_ETP_vs_nonETP.png
- outputs/MYCN_PP2A_correlation_stats.csv   (Pearson r & p, by all samples and by subgroup when possible)
- outputs/etp_vs_nonETP_stats.csv           (group means, fold-change, Welch t p-value)

Requirements:
    pip install pandas numpy scipy matplotlib statsmodels openpyxl
"""

import os
import io
import gzip
import math
import json
import textwrap
from typing import Tuple, List, Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


# -----------------------------
# --------- SETTINGS ----------
# -----------------------------

PROJECT_ROOT = os.path.abspath(".")
SERIES_PATH  = os.path.join(PROJECT_ROOT, "GSE28703_series_matrix.txt.gz")   # required (as downloaded from GEO)
# Use your actual annotation filename; auto-detects .gz or .txt content
GPL_ANNOT    = os.path.join(PROJECT_ROOT, "GPL13158-5065.txt")

OUT_DIR      = os.path.join(PROJECT_ROOT, "outputs")
SCATTER_DIR  = os.path.join(OUT_DIR, "scatter")
HEATMAP_DIR  = os.path.join(OUT_DIR, "heatmaps")
BOXPLOT_DIR  = os.path.join(OUT_DIR, "boxplots")
os.makedirs(SCATTER_DIR, exist_ok=True)
os.makedirs(HEATMAP_DIR, exist_ok=True)
os.makedirs(BOXPLOT_DIR, exist_ok=True)

# --- Genes to analyze (Option A: ONLY PP2A core + MYCN, from slides) ---
SLIDE_GENES = [
    "MYCN",        # ↑ in ETP (slide claim)
    "PPP2R2D",     # ↑ in ETP (slide claim)
    "PPP2R5C",     # ↓ in ETP (slide claim)
    "PPP2R5D",     # ↓ in ETP (slide claim)
    "PPP2R1B",     # ↓ in ETP (slide claim)
    "PPP2CA"       # ↓ in ETP (slide claim)
]


# -----------------------------
# ------- I/O HELPERS ---------
# -----------------------------

def read_text_auto(path: str) -> str:
    """Read a text file that may be gzipped or plain text; return as UTF-8 string."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return f.read().decode("utf-8", errors="replace")
    else:
        with open(path, "rb") as f:
            return f.read().decode("utf-8", errors="replace")


# -----------------------------
# ------ GEO PARSING ----------
# -----------------------------

def parse_series_matrix(series_text: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a GEO Series Matrix text into:
      - expr_df: probes × samples numeric expression
      - meta_df: sample metadata (columns gleaned from series text rows)

    Expression block is between:
        !series_matrix_table_begin
        !series_matrix_table_end
    """
    lower = series_text
    if "!series_matrix_table_begin" not in lower:
        raise ValueError("Series matrix table begin marker not found.")
    if "!series_matrix_table_end" not in lower:
        raise ValueError("Series matrix table end marker not found.")

    table_txt = lower.split("!series_matrix_table_begin\n", 1)[1].split("\n!series_matrix_table_end", 1)[0]
    expr_df = pd.read_csv(io.StringIO(table_txt), sep="\t", index_col=0)
    # Convert to numeric if possible
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    # Gather sample metadata (best-effort)
    meta_rows = []
    for line in series_text.splitlines():
        if line.startswith("!Sample_"):
            meta_rows.append(line)
    # Build metadata frame keyed by sample names (columns in expr_df)
    meta_dict: Dict[str, Dict[str, str]] = {s: {} for s in expr_df.columns}
    for line in meta_rows:
        # e.g., !Sample_title = GSMxxxx... or !Sample_characteristics_ch1 = key: value
        if " =" not in line:
            continue
        key, val = line.split(" =", 1)
        key = key.replace("!Sample_", "").strip()
        val = val.strip()
        # Try to align this 'val' to the right sample if it contains tab-separated fields
        # Some matrices list values in same order as samples; here we best-effort broadcast
        # If there are multiple values separated by "\t", map by index position
        if "\t" in val and len(val.split("\t")) == len(expr_df.columns):
            vals = val.split("\t")
            for s, v in zip(expr_df.columns, vals):
                meta_dict[s][key] = v
        else:
            # If a single value, just store (could be repeated lines; last wins)
            for s in expr_df.columns:
                if key not in meta_dict[s]:
                    meta_dict[s][key] = val

    meta_df = pd.DataFrame.from_dict(meta_dict, orient="index")
    meta_df.index.name = "Sample"
    return expr_df, meta_df


def parse_gpl_annotation(gpl_text: str) -> pd.DataFrame:
    """
    Parse GPL full table: return DataFrame with at least ['ID', 'Gene Symbol'] columns.
    The GPL page includes a prose header; the real table starts at line beginning with "ID\t".
    """
    lines = gpl_text.splitlines()
    header_idx = None
    for i, l in enumerate(lines):
        if l.startswith("ID\t"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find 'ID' header line in GPL annotation file.")
    tab_txt = "\n".join(lines[header_idx:])
    gpl_df = pd.read_csv(io.StringIO(tab_txt), sep="\t", dtype=str)
    # Harmonize column
    if "Gene Symbol" not in gpl_df.columns:
        # Try common alternates
        for alt in ["Gene Symbol;"]:
            if alt in gpl_df.columns:
                gpl_df["Gene Symbol"] = gpl_df[alt]
                break
    # Keep only needed columns
    keep = [c for c in ["ID", "Gene Symbol"] if c in gpl_df.columns]
    gpl_df = gpl_df[keep].copy()
    gpl_df["Gene Symbol"] = gpl_df["Gene Symbol"].fillna("").astype(str)
    return gpl_df


# -----------------------------
# ---- COLLAPSE TO GENES ------
# -----------------------------

def probes_to_genes(expr_by_probe: pd.DataFrame, gpl_df: pd.DataFrame) -> pd.DataFrame:
    """
    Map probes to gene symbols using GPL; if a probe maps to multiple symbols (e.g., 'A /// B'),
    split and keep each; then collapse to gene by taking mean across probes per gene.
    """
    # Make a mapping probe -> list of gene symbols
    annot = gpl_df.set_index("ID")["Gene Symbol"].to_dict()
    # Build mapping from gene -> list of probe rows
    gene_to_values = {}
    found = 0
    for probe, row in expr_by_probe.iterrows():
        sym = annot.get(probe, "")
        if not sym:
            continue
        # Split multi-symbol entries on common separators
        parts = [p.strip() for p in sym.replace("///", "|").replace(",", "|").split("|") if p.strip()]
        for g in parts:
            if g == "---" or g == "NA":
                continue
            found += 1
            gene_to_values.setdefault(g, []).append(row.values.astype(float))

    # Average per gene
    records = {}
    for g, lst in gene_to_values.items():
        if len(lst) == 1:
            records[g] = lst[0]
        else:
            records[g] = np.nanmean(np.vstack(lst), axis=0)
    if not records:
        raise ValueError("No probe-to-gene mapping produced any genes. Check GPL table and platform match.")
    genes_df = pd.DataFrame(records, index=expr_by_probe.columns).T  # genes × samples
    genes_df.index.name = "Gene"
    genes_df.columns.name = None
    return genes_df


# -----------------------------
# ---- SUBTYPE INFERENCE -------
# -----------------------------

def infer_etp_labels(meta_df: pd.DataFrame) -> pd.Series:
    """
    Best-effort rule: if any of the common fields contains 'ETP' (case-insensitive) -> 'ETP', else 'nonETP'.
    Uses title, source_name_ch1, characteristics fields if present.
    """
    candidates = [c for c in meta_df.columns if c.lower() in (
        "title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2",
        "characteristics_ch2", "description"
    ) or c.lower().startswith("characteristics")]
    labels = []
    for idx, row in meta_df.iterrows():
        text = " ".join([str(row[c]) for c in candidates if c in meta_df.columns and not pd.isna(row[c])]).lower()
        label = "ETP" if "etp" in text else "nonETP"
        labels.append(label)
    lab = pd.Series(labels, index=meta_df.index, name="Subtype")
    return lab


# -----------------------------
# ----- CORR & PLOTTING --------
# -----------------------------

def pearsonr_safe(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Pearson r, p with nan-safe behavior; returns (nan, nan) if insufficient variance or <3 points."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    if len(x) < 3:
        return (np.nan, np.nan)
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return (np.nan, np.nan)
    r, p = stats.pearsonr(x, y)
    return (float(r), float(p))


def scatter_colored(df: pd.DataFrame, xgene: str, ygene: str, hue: str, outfile: str):
    """
    One colored scatter plot: x=MYCN, y=target gene, colored by subtype.
    Shows Pearson r, p for ALL samples (and for subtypes if >2 samples).
    """
    use = df[[xgene, ygene, hue]].dropna()
    if use.empty:
        print(f"[WARN] No data for scatter: {xgene} vs {ygene}")
        return

    # Stats
    r_all, p_all = pearsonr_safe(use[xgene].values, use[ygene].values)
    r_e, p_e = np.nan, np.nan
    r_n, p_n = np.nan, np.nan
    if (use[hue] == "ETP").sum() >= 3:
        ue = use[use[hue] == "ETP"]
        r_e, p_e = pearsonr_safe(ue[xgene].values, ue[ygene].values)
    if (use[hue] == "nonETP").sum() >= 3:
        un = use[use[hue] == "nonETP"]
        r_n, p_n = pearsonr_safe(un[xgene].values, un[ygene].values)

    # Plot
    plt.figure(figsize=(4.4, 4.0))
    for label in ["ETP", "nonETP"]:
        d = use[use[hue] == label]
        if len(d):
            plt.scatter(d[xgene], d[ygene], alpha=0.8, label=f"{label} (n={len(d)})")
            # simple OLS fit line for each group if >= 3 points
            if len(d) >= 3 and np.nanstd(d[xgene]) > 0:
                coef = np.polyfit(d[xgene], d[ygene], 1)
                xs = np.linspace(d[xgene].min(), d[xgene].max(), 50)
                ys = coef[0] * xs + coef[1]
                plt.plot(xs, ys, alpha=0.7)

    plt.xlabel(xgene)
    plt.ylabel(ygene)
    plt.title(f"{ygene} vs {xgene}")
    txt = f"All: r={r_all:.2f}, p={p_all:.2g}"
    if not np.isnan(r_e): txt += f"\nETP: r={r_e:.2f}, p={p_e:.2g}"
    if not np.isnan(r_n): txt += f"\nnonETP: r={r_n:.2f}, p={p_n:.2g}"
    plt.gca().text(0.02, 0.98, txt, transform=plt.gca().transAxes, va="top", ha="left", fontsize=9)
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

    return {
        "gene": ygene,
        "r_all": r_all, "p_all": p_all,
        "r_ETP": r_e,  "p_ETP": p_e,
        "r_nonETP": r_n, "p_nonETP": p_n
    }


def heatmap_core(df: pd.DataFrame, genes: List[str], subtype: pd.Series, outfile: str):
    """
    Simple z-score heatmap for 'genes' across all samples, columns clustered by expression similarity.
    Adds a top color bar for subtype.
    """
    present = [g for g in genes if g in df.columns]
    if not present:
        print("[WARN] No genes present for heatmap.")
        return
    X = df[present].copy()
    # z-score rows (genes) across samples
    Z = (X - X.mean(axis=0)) / X.std(axis=0).replace(0, np.nan)

    # Reorder samples by hierarchical clustering on Z (columns)
    from scipy.cluster.hierarchy import linkage, leaves_list
    # fill NaNs with 0 for distance calc (neutral after z)
    Z2 = Z.fillna(0.0)
    # cluster on samples (columns)
    try:
        col_link = linkage(Z2.T, method="average", metric="correlation")
        order_cols = leaves_list(col_link)
    except Exception:
        order_cols = np.arange(Z2.shape[1])

    # Arrange
    Zplot = Z2.iloc[:, order_cols]
    cols = Zplot.columns

    # Subtype color bar (top)
    subtype = subtype.reindex(cols)
    col_colors = np.where(subtype == "ETP", 0.65, 0.15)  # greyscale bar: ETP lighter

    # Plot
    fig, ax = plt.subplots(figsize=(max(6, Zplot.shape[1]*0.12), 2.0 + 0.3*len(present)))
    im = ax.imshow(Zplot.T.values, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=7)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=90, fontsize=8)
    ax.set_title("Core genes (z-score across samples)")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("z-score")

    # Add subtype bar above heatmap
    ax2 = fig.add_axes([ax.get_position().x0, ax.get_position().y1 + 0.01,
                        ax.get_position().width, 0.05])
    ax2.imshow(col_colors.reshape(1, -1), aspect="auto", cmap="Greys")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("Subtype: light=ETP, dark=nonETP", fontsize=8)
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)


def etp_boxplot(expr_df: pd.DataFrame, labels: pd.Series, gene: str, outpath: str):
    """
    ETP vs nonETP boxplot with jitter and Welch's t-test p-value annotation.
    Returns a stats dict or None if insufficient data.
    """
    if gene not in expr_df.columns:
        print(f"[WARN] Missing gene for boxplot: {gene}")
        return None

    data = pd.DataFrame({"value": expr_df[gene], "Subtype": labels}).dropna()
    etp = data.loc[data["Subtype"] == "ETP", "value"]
    non = data.loc[data["Subtype"] == "nonETP", "value"]

    n_etp, n_non = len(etp), len(non)
    if n_etp < 2 or n_non < 2:
        print(f"[WARN] {gene}: too few samples (ETP={n_etp}, nonETP={n_non}). Skipping.")
        return None

    # Welch t-test
    t_stat, p_val = stats.ttest_ind(etp, non, equal_var=False, nan_policy="omit")

    # Plot
    fig, ax = plt.subplots(figsize=(4.0, 3.8))
    ax.boxplot([etp, non], labels=["ETP", "nonETP"], showfliers=False)
    # jitter
    xs = np.concatenate([np.ones(n_etp) * 1, np.ones(n_non) * 2])
    ys = pd.concat([etp.reset_index(drop=True), non.reset_index(drop=True)], ignore_index=True)
    ax.scatter(xs, ys, alpha=0.7)

    ax.set_ylabel(gene)
    ax.set_title(f"{gene}: ETP vs nonETP")

    ymax = float(np.nanmax(ys))
    ymin = float(np.nanmin(ys))
    span = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 1.0
    y_annot = ymax + span * 0.08
    ax.plot([1, 1, 2, 2], [y_annot*0.98, y_annot, y_annot, y_annot*0.98])
    ax.text(1.5, y_annot*1.02, f"p = {p_val:.3g}", ha="center", va="bottom")

    ax.set_xticklabels([f"ETP (n={n_etp})", f"nonETP (n={n_non})"])
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)

    return {
        "gene": gene,
        "mean_ETP": float(np.nanmean(etp)),
        "mean_nonETP": float(np.nanmean(non)),
        "fold_change(ETP/nonETP)": float(np.nanmean(etp) / np.nanmean(non)) if float(np.nanmean(non)) != 0 else np.nan,
        "t_pvalue": float(p_val),
        "n_ETP": int(n_etp),
        "n_nonETP": int(n_non)
    }


# -----------------------------
# ------------ MAIN -----------
# -----------------------------

def main():
    # ---- 1) Load GEO series & GPL annotation ----
    series_text = read_text_auto(SERIES_PATH)
    gpl_text    = read_text_auto(GPL_ANNOT)

    expr_by_probe, meta_df = parse_series_matrix(series_text)   # probes × samples
    gpl_df = parse_gpl_annotation(gpl_text)                     # ['ID','Gene Symbol']

    # ---- 2) Collapse probes → genes ----
    genes_by_sample = probes_to_genes(expr_by_probe, gpl_df)    # genes × samples
    # ensure numeric
    genes_by_sample = genes_by_sample.apply(pd.to_numeric, errors="coerce")

    # ---- 3) Subtype labels ----
    subtype = infer_etp_labels(meta_df)                         # index=sample, value={'ETP','nonETP'}
    # align columns/samples
    common_samples = [s for s in genes_by_sample.columns if s in subtype.index]
    if len(common_samples) < 5:
        print("[WARN] Very few overlapping samples between expression and metadata.")
    genes_by_sample = genes_by_sample[common_samples]
    subtype = subtype.loc[common_samples]

    # ---- 4) Build sample × gene matrix (keep only SLIDE_GENES) ----
    present_genes = [g for g in SLIDE_GENES if g in genes_by_sample.index]
    missing_genes = [g for g in SLIDE_GENES if g not in genes_by_sample.index]
    if missing_genes:
        print(f"[WARN] Missing genes (not on platform / not mapped): {', '.join(missing_genes)}")

    samp_expr = genes_by_sample.loc[present_genes].T   # samples × genes
    # attach subtype for convenience
    samp_expr["Subtype"] = subtype

    # ---- 5) Correlation (MYCN vs each PP2A) + one colored scatter per gene ----
    corr_rows = []
    if "MYCN" in samp_expr.columns:
        for tgt in [g for g in present_genes if g != "MYCN"]:
            out_scatter = os.path.join(SCATTER_DIR, f"{tgt}_vs_MYCN_colored.png")
            stats_row = scatter_colored(samp_expr, xgene="MYCN", ygene=tgt, hue="Subtype", outfile=out_scatter)
            if stats_row is not None:
                corr_rows.append(stats_row)
    else:
        print("[WARN] 'MYCN' not present; skipping correlation plots.")

    if corr_rows:
        corr_df = pd.DataFrame(corr_rows)
        corr_df = corr_df[["gene", "r_all", "p_all", "r_ETP", "p_ETP", "r_nonETP", "p_nonETP"]]
        corr_df.to_csv(os.path.join(OUT_DIR, "MYCN_PP2A_correlation_stats.csv"), index=False)

    # ---- 6) Heatmap (core) ----
    if len(present_genes) >= 2:
        heatmap_core(samp_expr.drop(columns=["Subtype"]), present_genes, subtype, os.path.join(HEATMAP_DIR, "heatmap_core.png"))
    else:
        print("[WARN] <2 genes present; skipping heatmap.")

    # ---- 7) ETP vs nonETP boxplots + stats (THIS IS THE NEW PART YOU ASKED FOR) ----
    rows = []
    for g in present_genes:
        outpng = os.path.join(BOXPLOT_DIR, f"{g}_ETP_vs_nonETP.png")
        row = etp_boxplot(samp_expr.drop(columns=["Subtype"]), subtype, g, outpng)
        if row is not None:
            rows.append(row)
    if rows:
        etp_stats = pd.DataFrame(rows)
        etp_stats = etp_stats.sort_values("t_pvalue")
        etp_stats.to_csv(os.path.join(OUT_DIR, "etp_vs_nonETP_stats.csv"), index=False)
        print("[OK] ETP vs nonETP boxplots saved:", BOXPLOT_DIR)
        print("[OK] ETP vs nonETP stats saved:", os.path.join(OUT_DIR, "etp_vs_nonETP_stats.csv"))
    else:
        print("[WARN] No ETP vs nonETP stats generated (likely too few samples).")

    print("\n[DONE] Outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
