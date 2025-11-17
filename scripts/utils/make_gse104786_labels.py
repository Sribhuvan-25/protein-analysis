# make_gse104786_labels.py
# -*- coding: utf-8 -*-
"""
Create sample labels for GSE104786 dataset.
Classification: SCPC (Small Cell Prostate Cancer = NEPC) vs Adenocarcinoma
"""
import os, re, gzip, csv

PROJECT_ROOT = "."
SERIES_MATRIX = os.path.join(PROJECT_ROOT, "Data/NEPC/GSE104786_series_matrix.txt.gz")
OUT_TEMPLATE  = os.path.join(PROJECT_ROOT, "Data/NEPC/sample_labels_gse104786_TEMPLATE.csv")
OUT_FINAL     = os.path.join(PROJECT_ROOT, "Data/NEPC/sample_labels_gse104786.csv")

def read_meta(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # pull GSM list
    gsm_ids = []
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            gsm_ids = [x.strip().strip('"') for x in ln.strip().split("\t")[1:]]
            break
    if not gsm_ids:
        raise RuntimeError("Could not find !Sample_geo_accession line")

    # gather useful sample fields
    # Note: characteristics_ch1 appears multiple times, we need to collect all
    meta = {g:{} for g in gsm_ids}
    characteristics_ch1_lines = []
    
    for ln in lines:
        if ln.startswith("!Sample_"):
            parts = ln.strip().split("\t")
            tag = parts[0][len("!Sample_"):].lower()
            vals = [v.strip('"') for v in parts[1:]]
            
            # Special handling for characteristics_ch1 (multiple lines)
            if tag == "characteristics_ch1":
                characteristics_ch1_lines.append(vals)
            else:
                # Store other fields normally
                for g, v in zip(gsm_ids, vals):
                    if tag not in meta[g] or not meta[g][tag]:  # Don't overwrite if already set
                        meta[g][tag] = v
    
    # Find the characteristics_ch1 line that contains annotations
    annotation_line = None
    for char_line in characteristics_ch1_lines:
        if any('annotation' in v.lower() for v in char_line):
            annotation_line = char_line
            break
    
    # Map annotation to samples
    if annotation_line:
        for g, ann in zip(gsm_ids, annotation_line):
            meta[g]["annotation"] = ann

    return gsm_ids, meta

def suggest_subtype(m: dict) -> tuple[str,str]:
    """
    Classification for GSE104786:
    - SCPC (Small Cell Prostate Cancer) = NEPC
    - Adenocarcinoma = High-grade adenocarcinoma
    
    Based on annotation field from characteristics_ch1.
    """
    # Check annotation field first (most reliable)
    annotation = (m.get("annotation") or "").lower()
    
    # Check for SCPC (Small Cell Prostate Cancer)
    if re.search(r"\b(small[\s-]*cell[\s-]*carcinoma|small[\s-]*cell)\b", annotation):
        return "SCPC", "matched: small cell carcinoma (from annotation)"
    
    # Check for Adenocarcinoma (but not if it also says small cell)
    if re.search(r"\b(adenocarcinoma|adeno)\b", annotation) and "small cell" not in annotation:
        return "Adenocarcinoma", "matched: adenocarcinoma (from annotation)"
    
    # Fallback: check all metadata
    text = " ".join(str(v) for v in m.values()).lower()
    if re.search(r"\b(scpc|small[\s-]*cell[\s-]*carcinoma|small[\s-]*cell)\b", text):
        return "SCPC", "matched: SCPC/small cell carcinoma (from metadata)"
    
    if re.search(r"\b(adenocarcinoma|adeno)\b", text):
        return "Adenocarcinoma", "matched: adenocarcinoma (from metadata)"
    
    return "", "WARNING: Unknown subtype - manual review needed"

def main():
    os.makedirs(os.path.dirname(OUT_TEMPLATE), exist_ok=True)
    gsm_ids, meta = read_meta(SERIES_MATRIX)

    # Build rows with auto-filled subtypes
    rows = []
    scpc_count = 0
    adeno_count = 0
    
    for g in gsm_ids:
        sug, reason = suggest_subtype(meta[g])
        if sug:
            if sug == "SCPC":
                scpc_count += 1
            elif sug == "Adenocarcinoma":
                adeno_count += 1
        
        rows.append({
            "sample_id": g,
            "subtype": sug if sug else "UNKNOWN",
            "suggest_reason": reason,
            "title": meta[g].get("title",""),
            "source": meta[g].get("source_name_ch1",""),
            "annotation": meta[g].get("annotation",""),
            "char1": meta[g].get("characteristics_ch1",""),
            "description": meta[g].get("description","")
        })

    # Write template CSV (for reference/debugging)
    with open(OUT_TEMPLATE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    
    # Write final labels file (sample_id, subtype only) - ready to use
    with open(OUT_FINAL, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "subtype"])
        w.writeheader()
        for row in rows:
            if row["subtype"] != "UNKNOWN":  # Only include classified samples
                w.writerow({"sample_id": row["sample_id"], 
                           "subtype": row["subtype"]})
    
    print(f"[OK] Template written: {OUT_TEMPLATE}")
    print(f"[OK] Labels file created: {OUT_FINAL}")
    print(f"\nClassification summary:")
    print(f"  SCPC (NEPC): {scpc_count} samples")
    print(f"  Adenocarcinoma: {adeno_count} samples")
    print(f"  Total classified: {scpc_count + adeno_count} samples")
    print(f"\nAll samples classified automatically based on metadata.")

if __name__ == "__main__":
    main()

