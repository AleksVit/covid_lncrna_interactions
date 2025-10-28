#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert RNAInter Interactor1.Symbol -> Ensembl ID (ENSG/ENST) for Homo sapiens.
- Input: CSV with delimiter ';'
- Output: CSV with delimiter ';' and two new columns:
    Interactor1.Ensembl_ID, Interactor1.map_status
Mapping logic:
  1) If Interactor1.Symbol already looks like Ensembl (ENSG/ENST/ENSMUSG/...), keep it.
  2) Else query Ensembl REST /xrefs/symbol/homo_sapiens/{symbol}.
  3) If Category1 == miRNA and symbol is like 'hsa-miR-21-5p', normalize to 'MIR21' and retry.
  4) Fallback: query RNAcentral (returns RNAcentral:URS..., useful for rRNA/7SL/rare ncRNAs).
"""

import sys, csv, re, time, json
from typing import Optional, Tuple
from urllib.parse import quote
import urllib.request, urllib.parse

# --------- CONFIG ----------
USE_RNACENTRAL = True  # set False to disable fallback
INPUT = sys.argv[1] if len(sys.argv) > 1 else None
DELIM = ";"  # << your input CSV uses semicolons
OUTPUT_SUFFIX = "_with_ensembl.csv"
ENSEMBL_REST = "https://rest.ensembl.org"
RNACENTRAL_API = "https://rnacentral.org/api/v1/rna"
SPECIES_OK = {"Homo sapiens", "human", "HUMAN", "Homo sapiens"}  # various spellings

# Rate limiting safety (Ensembl suggests reasonable pacing)
SLEEP_BETWEEN_CALLS = 0.12  # seconds

# --------- UTILITIES ----------

def http_json(url: str, headers=None, retries: int = 3, timeout: int = 30):
    if headers is None:
        headers = {"Content-Type": "application/json", "User-Agent": "rnainter-mapper/1.0"}
    req = urllib.request.Request(url, headers=headers)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read().decode("utf-8", errors="replace")
                return json.loads(data)
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(1.5 * (attempt + 1))
    return None

def ensembl_xrefs_symbol(symbol: str):
    """Query Ensembl xrefs by symbol (human). Returns list of dicts."""
    url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{quote(symbol)}?content-type=application/json"
    time.sleep(SLEEP_BETWEEN_CALLS)
    return http_json(url)

def pick_ensembl_id_from_xrefs(xrefs) -> Optional[str]:
    """Prefer Ensembl Gene ID (ENSG...). Fallback to first ENST if no gene found."""
    if not xrefs:
        return None
    genes = [x for x in xrefs if x.get("type") == "Gene" and x.get("id","").startswith("ENSG")]
    if genes:
        genes.sort(key=lambda x: x.get("id",""))
        return genes[0]["id"]
    transcripts = [x for x in xrefs if x.get("type") == "Transcript" and x.get("id","").startswith("ENST")]
    if transcripts:
        transcripts.sort(key=lambda x: x.get("id",""))
        return transcripts[0]["id"]
    return None

def normalize_mir_symbol(sym: str) -> str:
    """
    hsa-miR-21-5p -> MIR21
    miR-155-3p   -> MIR155
    mir-16       -> MIR16
    """
    s = sym.strip()
    s = re.sub(r"(?i)^hsa[-_]", "", s)
    s = re.sub(r"(?i)^mir[-_]", "MIR", s)
    s = re.sub(r"(?i)^mir", "MIR", s)
    s = re.sub(r"[-_]?([35]p)$", "", s)
    s = re.sub(r"[-_]", "", s)
    return s

def rnacentral_search_symbol_hs(symbol: str) -> Optional[str]:
    """
    Very light fallback: text search in RNAcentral restricted to Homo sapiens.
    Returns first RNAcentral accession (URS...), or None.
    """
    q = urllib.parse.quote(f'{symbol} AND species:"Homo sapiens"')
    url = f"{RNACENTRAL_API}/?query={q}&page_size=1"
    time.sleep(SLEEP_BETWEEN_CALLS)
    data = http_json(url)
    if not data or not data.get("results"):
        return None
    res = data["results"][0]
    return res.get("rnacentral_id") or res.get("accession")

ENSEMBL_ID_RE = re.compile(r"^ENS(G|T|MUSG|MUST)\d+", re.IGNORECASE)

def map_symbol_to_id(symbol: str, category: str) -> Tuple[Optional[str], str]:
    """
    Returns (id, status):
      id  -> ENSG/ENST... or RNAcentral:URS...
      status -> mapping status string
    """
    if not symbol:
        return None, "empty_symbol"

    sym = symbol.strip()

    # If already Ensembl-like
    if ENSEMBL_ID_RE.match(sym):
        return sym, "already_ensembl"

    # 1) Try Ensembl directly
    xrefs = ensembl_xrefs_symbol(sym)
    ens = pick_ensembl_id_from_xrefs(xrefs)
    if ens:
        return ens, "mapped_via_symbol"

    # 2) Special-case miRNA normalization
    if category and category.strip().lower() == "mirna":
        norm = normalize_mir_symbol(sym)
        if norm and norm != sym:
            xrefs2 = ensembl_xrefs_symbol(norm)
            ens2 = pick_ensembl_id_from_xrefs(xrefs2)
            if ens2:
                return ens2, "mirna_normalized_symbol"

    # 3) Fallback to RNAcentral (works well for rRNA/7SL and odd ncRNAs)
    if USE_RNACENTRAL:
        rc = rnacentral_search_symbol_hs(sym)
        if not rc and category and category.strip().lower() == "mirna":
            norm = normalize_mir_symbol(sym)
            if norm:
                rc = rnacentral_search_symbol_hs(norm)
        if rc:
            return f"RNAcentral:{rc}", "rnacentral_fallback"

    return None, "not_found"

# --------- MAIN ----------

def main():
    if not INPUT:
        print("Usage: python3 map_interactor1_to_ensembl.py <rnainter.csv_with_semicolons>", file=sys.stderr)
        sys.exit(2)

    # Read CSV with semicolon delimiter
    with open(INPUT, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=DELIM)
        fieldnames = reader.fieldnames or []

        # standardize expected columns
        sym_key_candidates = ["Interactor1.Symbol", "Interactor1", "Interactor1Symbol"]
        cat_key_candidates = ["Category1", "Interactor1.Category"]
        sp_key_candidates  = ["Species1", "Interactor1.Species"]

        def pick_key(cands):
            for k in cands:
                if k in fieldnames:
                    return k
            return None

        sym_key = pick_key(sym_key_candidates)
        cat_key = pick_key(cat_key_candidates)
        sp_key  = pick_key(sp_key_candidates)

        if not sym_key or not sp_key:
            print("ERROR: Can't find columns for Interactor1.Symbol and Species1 in the header.", file=sys.stderr)
            print("Header columns:", fieldnames, file=sys.stderr)
            sys.exit(1)

        out_fields = list(fieldnames)
        if "Interactor1.Ensembl_ID" not in out_fields:
            out_fields.append("Interactor1.Ensembl_ID")
        if "Interactor1.map_status" not in out_fields:
            out_fields.append("Interactor1.map_status")

        rows_out = []
        for row in reader:
            sym = (row.get(sym_key) or "").strip()
            cat = (row.get(cat_key) or "").strip()
            sp  = (row.get(sp_key)  or "").strip()

            # Only map human
            if sp not in SPECIES_OK:
                row["Interactor1.Ensembl_ID"] = ""
                row["Interactor1.map_status"] = "non_human_or_unknown_species"
                rows_out.append(row)
                continue

            ens, status = map_symbol_to_id(sym, cat)
            row["Interactor1.Ensembl_ID"] = ens or ""
            row["Interactor1.map_status"] = status
            rows_out.append(row)

    out_name = INPUT.rsplit(".", 1)[0] + OUTPUT_SUFFIX
    with open(out_name, "w", newline="", encoding="utf-8") as fo:
        writer = csv.DictWriter(fo, fieldnames=out_fields, delimiter=DELIM)
        writer.writeheader()
        for r in rows_out:
            writer.writerow(r)

    print(f"Done: {out_name}")
    print("Statuses: already_ensembl | mapped_via_symbol | mirna_normalized_symbol | rnacentral_fallback | not_found | non_human_or_unknown_species | empty_symbol")

if __name__ == "__main__":
    main()


