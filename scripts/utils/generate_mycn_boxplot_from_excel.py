#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate MYCN boxplot from Excel file to match the style of other PP2A subunit boxplots.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# Import the boxplot function from analysis.py
from analysis import plot_boxplot_etp_vs_nonetp, welch_ttest, cohens_d

# Paths
EXCEL_FILE = os.path.join("Data", "MYCN-GDS4299_with IDs.xlsx")
OUTPUT_DIR = os.path.join("outputs", "boxplots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output file (replace the existing MYCN boxplot to keep consistency)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "MYCN_box_ETP_vs_nonETP.png")

def main():
    print("[INFO] Reading Excel file...")
    df = pd.read_excel(EXCEL_FILE)
    
    # Check column names
    print(f"[INFO] Columns: {df.columns.tolist()}")
    print(f"[INFO] Shape: {df.shape}")
    
    # Based on the structure we saw:
    # Column 0: Unnamed: 0 - GSM IDs
    # Column 1: Unnamed: 1 - Subtype labels ("ETP ALL" or "non-ETP ALL")
    # Column 2: MYCN - original expression
    # Column 3: MYCN(log2 expression) - log2 transformed values
    
    # Use column index 1 for subtype and column index 3 for expression
    subtype_col = df.columns[1]  # Unnamed: 1
    expr_col = df.columns[3]      # MYCN(log2 expression)
    
    print(f"[INFO] Using subtype column: '{subtype_col}'")
    print(f"[INFO] Using expression column: '{expr_col}'")
    
    # Extract ETP and nonETP groups
    # Clean the subtype labels
    df['Subtype_clean'] = df[subtype_col].astype(str).str.strip()
    
    # Identify ETP samples (contains "ETP" but not "non")
    etp_mask = df['Subtype_clean'].str.contains('ETP', case=False, na=False) & \
               ~df['Subtype_clean'].str.contains('non', case=False, na=False)
    non_mask = df['Subtype_clean'].str.contains('non', case=False, na=False)
    
    # Extract values
    etp_values = df.loc[etp_mask, expr_col].dropna()
    non_values = df.loc[non_mask, expr_col].dropna()
    
    print(f"[INFO] ETP samples: {len(etp_values)}")
    print(f"[INFO] nonETP samples: {len(non_values)}")
    print(f"[INFO] ETP mean: {etp_values.mean():.3f}, std: {etp_values.std():.3f}")
    print(f"[INFO] nonETP mean: {non_values.mean():.3f}, std: {non_values.std():.3f}")
    
    if len(etp_values) < 2 or len(non_values) < 2:
        print("[ERROR] Need at least 2 samples in each group.")
        return
    
    # Convert to Series for compatibility
    etp_series = pd.Series(etp_values.values, name='ETP')
    non_series = pd.Series(non_values.values, name='nonETP')
    
    # Generate boxplot using the same function as other genes
    print(f"[INFO] Generating boxplot...")
    plot_boxplot_etp_vs_nonetp(etp_series, non_series, "MYCN", OUTPUT_FILE)
    
    print(f"[OK] MYCN boxplot saved to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()

