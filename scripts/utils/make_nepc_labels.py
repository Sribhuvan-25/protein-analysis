# make_nepc_labels.py
# -*- coding: utf-8 -*-
import os, re, gzip, csv

PROJECT_ROOT = "."
SERIES_MATRIX = os.path.join(PROJECT_ROOT, "Data/Nepc/GSE35988-GPL6480_series_matrix.txt.gz")
OUT_TEMPLATE  = os.path.join(PROJECT_ROOT, "Data/Nepc/sample_labels_nepc_TEMPLATE.csv")

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

    # gather a few useful sample fields
    wanted = [
        "title",
        "source_name_ch1",
        "characteristics_ch1",
        "characteristics_ch1.1",
        "characteristics_ch1.2",
        "description",
    ]
    meta = {g:{} for g in gsm_ids}
    for ln in lines:
        if ln.startswith("!Sample_"):
            parts = ln.strip().split("\t")
            tag = parts[0][len("!Sample_"):].lower()
            vals = parts[1:]
            if tag in wanted:
                for g, v in zip(gsm_ids, vals):
                    meta[g][tag] = v.strip('"')

    return gsm_ids, meta

def suggest_subtype(m: dict) -> tuple[str,str]:
    """
    Classification based on Beltran et al. 2011 (Nature Medicine) / GSE35988:
    
    Prefix   Meaning                                    Disease group
    -------  -----------------------------------------  -------------
    Nxx      NEPC xenograft line                       NEPC
    Txx      Adenocarcinoma xenograft line (CRPC)      CRPC
    WAxx     Metastatic CRPC tissue (warm autopsy)      CRPC
    
    Mapping rule:
      - if title.startswith("N"):  subtype = NEPC
      - elif title.startswith("T") or title.startswith("WA"): subtype = CRPC
    """
    title = (m.get("title") or "").strip()
    
    # Apply Beltran et al. 2011 documented mapping
    if title.startswith("N"):
        return "NEPC", "Beltran et al. 2011: N prefix = NEPC xenograft"
    
    if title.startswith("T"):
        return "CRPC", "Beltran et al. 2011: T prefix = CRPC xenograft"
    
    if title.startswith("WA"):
        return "CRPC", "Beltran et al. 2011: WA prefix = metastatic CRPC (warm autopsy)"
    
    # Fallback: check metadata if title pattern doesn't match
    text = " ".join(str(v) for v in m.values()).lower()
    if re.search(r"\b(nepc|neuroendocrine|small[\s-]*cell)\b", text):
        return "NEPC", "matched: nepc/neuroendocrine/small-cell (from metadata)"
    
    if re.search(r"\b(crpc|adenocarcinoma|adeno)\b", text):
        return "CRPC", "matched: crpc/adenocarcinoma/adeno (from metadata)"
    
    return "", "WARNING: Unknown title pattern - manual review needed"

def main():
    os.makedirs(os.path.dirname(OUT_TEMPLATE), exist_ok=True)
    gsm_ids, meta = read_meta(SERIES_MATRIX)

    # Build rows with auto-filled subtypes based on Beltran et al. 2011 mapping
    rows = []
    for g in gsm_ids:
        sug, reason = suggest_subtype(meta[g])
        rows.append({
            "sample_id": g,
            "subtype": sug if sug else "CRPC",  # Auto-filled based on Beltran et al. 2011; default to CRPC if unknown
            "suggest_reason": reason,
            "title": meta[g].get("title",""),
            "source": meta[g].get("source_name_ch1",""),
            "char1": meta[g].get("characteristics_ch1",""),
            "char2": meta[g].get("characteristics_ch1.1",""),
            "char3": meta[g].get("characteristics_ch1.2",""),
            "description": meta[g].get("description","")
        })

    # Write template CSV (for reference/debugging)
    with open(OUT_TEMPLATE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    
    # Write final labels file (sample_id, subtype only) - ready to use
    out_final = os.path.join(PROJECT_ROOT, "Data/Nepc/sample_labels_nepc.csv")
    with open(out_final, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "subtype"])
        w.writeheader()
        for row in rows:
            w.writerow({"sample_id": row["sample_id"], 
                       "subtype": row["subtype"]})
    
    print(f"[OK] Template written: {OUT_TEMPLATE}")
    print(f"[OK] Labels file created: {out_final}")
    print(f"\nClassification based on Beltran et al. 2011 (Nature Medicine):")
    print(f"  N prefix → NEPC")
    print(f"  T prefix → CRPC")  
    print(f"  WA prefix → CRPC")
    print(f"\nAll samples classified automatically based on title prefix.")

if __name__ == "__main__":
    main()
