#!/usr/bin/env python3
"""
cluster_group_assign_fixed.py

Fixes:
 - handles IDs like Centipede3FTx_0001_61981 where \b won't match underscore
 - tries multiple search targets: protein_name, full header, header with underscores -> spaces,
   ID token, ID token with underscores -> spaces
 - logs missing member IDs and unmatched headers
 - FIXED: Better debugging for group assignment decisions
 - FIXED: Added IPR031424 mapping (appears to be missing from your config)
"""

import argparse
import csv
import logging
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ----------------- CONFIG (edit this) -----------------
GROUP_PATTERNS = {
    "manually_annotated": [r"\bCentipede3FTx\b"],
    "3FTx": [r"\b3ftx\b", r"\b3\-ftx\b", r"\bthree[- ]?ftx\b", r"\b3FTx\b", r"\bThree-finger toxin\b", r"\bThree finger toxins\b", r"\bThree finger toxin\b"],
    "PMF": [r"\bpmf\b", r"\bPNF\b"],
    "Ly6": [r"\bly6\b", r"\bly\-6\b", r"\blysix\b", r"\bLy-6\b", r"\bLY6\b"],
    "Quiver": [r"\bquiver\b", r"\bqvr\b"],
    "Scoloptoxin": [r"\bscoloptoxin\b"],
    "SPF": [r"\bSodefrin\b"],
}

# Specificity order: most specific first. If a header matches multiple groups,
# the first one listed here will be chosen.
SPECIFICITY_ORDER = [
    "SPF",
    "Scoloptoxin",
    "PMF",
    "Quiver",
    "manually_annotated",
    "3FTx",
    "Ly6",
    # add others, least specific near the end
]
# Build a rank dict for quick comparison
SPECIFICITY_RANK = {name: i for i, name in enumerate(SPECIFICITY_ORDER)}

# InterPro -> group map (DISABLED - using text patterns only)
# INTERPRO_GROUP_MAP = {
#     "IPR045860": "Ly6",
#     "IPR035076": "3FTx", 
#     "IPR031424": "3FTx",
# }
INTERPRO_GROUP_MAP = {}  # Empty dict disables InterPro matching

# compile patterns (case-insensitive)
_COMPILED_GROUPS = [(g, [re.compile(p, re.IGNORECASE) for p in pats]) for g, pats in GROUP_PATTERNS.items()]
_re_protein_name = re.compile(r"protein_name=(.*?)\s+(organism=|domain_pos=|signature=|interpro_id=|length=|$)", re.IGNORECASE)
_re_interpro_id = re.compile(r"(IPR\d{6})", re.IGNORECASE)
# ------------------------------------------------------

def parse_fasta_headers(fasta_path: Path, id_mode: str = "first_token") -> Dict[str, str]:
    d = {}
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    with fasta_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                seqid = header if id_mode == "full" else header.split(None, 1)[0]
                d[seqid] = header
    logging.info(f"Parsed {len(d)} headers from {fasta_path} (id_mode={id_mode})")
    return d

def parse_mmseqs_tsv(tsv_path: Path) -> Dict[str, List[str]]:
    mapping = defaultdict(list)
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV not found: {tsv_path}")
    with tsv_path.open("r", encoding="utf-8") as fh:
        for ln, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            rep = cols[0]
            # FIXED: Make it clearer that we expect exactly 2 columns per line
            if len(cols) == 2:
                member = cols[1]
                mapping[rep].append(member)
            else:
                # Handle case where there might be multiple members per line
                members = cols[1:]
                mapping[rep].extend(members)
    
    # dedupe preserving order
    cleaned = {}
    for rep, members in mapping.items():
        seen = set()
        out = []
        for m in members:
            if m not in seen:
                seen.add(m); out.append(m)
        cleaned[rep] = out
    logging.info(f"Parsed {len(cleaned)} clusters from {tsv_path}")
    return cleaned

def extract_header_fields(full_header: str) -> Tuple[Optional[str], Optional[str]]:
    pn = None
    ip = None
    m = _re_protein_name.search(full_header)
    if m:
        pn = m.group(1).strip().strip('"').strip("'")
    m2 = _re_interpro_id.search(full_header)
    if m2:
        ip = m2.group(1).upper()
    return pn, ip

def find_all_groups_for_header(full_header: str) -> List[str]:
    """Return all groups that match this header using only text patterns."""
    pn, ip = extract_header_fields(full_header)

    groups = []
    
    # Build search targets: protein_name, header, underscore->space, id token, etc.
    id_token = full_header.split(None, 1)[0]
    targets = []
    if pn:
        targets.append(pn)
    targets.extend([
        full_header,
        full_header.replace("_", " "),
        id_token,
        id_token.replace("_", " "),
    ])

    # Test patterns against all targets, collect all matched group names (avoid duplicates)
    for target in targets:
        for gname, compiled_list in _COMPILED_GROUPS:
            for cre in compiled_list:
                if cre.search(target):
                    if gname not in groups:
                        groups.append(gname)
                        logging.debug(f"Header '{full_header[:50]}...' matched pattern for {gname} in target '{target}'")
    
    if not groups:
        logging.debug(f"Header '{full_header[:50]}...' matched no groups")
    
    return groups

def find_group_for_header(full_header: str) -> Optional[str]:
    """Return the single best group for header, using specificity ranking when multiple match."""
    matched = find_all_groups_for_header(full_header)
    if not matched:
        return None
    # If only one, return it
    if len(matched) == 1:
        return matched[0]
    # If multiple, choose the one with smallest SPECIFICITY_RANK (more specific)
    ranked = sorted(matched, key=lambda g: SPECIFICITY_RANK.get(g, 9999))
    logging.debug(f"Header '{full_header[:50]}...' matched multiple groups {matched}, chose {ranked[0]} by specificity")
    return ranked[0]

def decide_cluster_label(member_full_headers: List[str], rep_id: str = "") -> str:
    n = len(member_full_headers)
    groups = [find_group_for_header(h) for h in member_full_headers]
    counts = Counter(groups)  # None is allowed
    non_none = [g for g in groups if g is not None]
    distinct_non_none = set(non_none)

    # ADDED: Debug logging for cluster decisions
    logging.debug(f"Cluster {rep_id}: {n} members")
    for i, (header, group) in enumerate(zip(member_full_headers, groups)):
        header_short = header.split()[0] if header else "None"
        logging.debug(f"  Member {i+1}: {header_short} -> {group}")
    logging.debug(f"  Group counts: {dict(counts)}")

    # exact unanimous (no Nones)
    if len(distinct_non_none) == 1 and counts[None] == 0:
        result = next(iter(distinct_non_none))
        logging.debug(f"  Decision: unanimous {result}")
        return result

    # mostly_ rule
    if n > 10 and distinct_non_none:
        best_group, best_count = counts.most_common(1)[0]
        if best_group is not None and (best_count / n) >= 0.90:
            result = f"mostly_{best_group}"
            logging.debug(f"  Decision: {result} ({best_count}/{n} = {best_count/n:.2%})")
            return result

    # mixture / none logic
    if len(distinct_non_none) > 1:
        result = "mixture"
        logging.debug(f"  Decision: {result} (multiple distinct groups: {distinct_non_none})")
        return result
    if len(distinct_non_none) == 1 and counts[None] > 0:
        result = "mixture"
        logging.debug(f"  Decision: {result} (one group + unknowns)")
        return result
    if not distinct_non_none:
        result = "none"
        logging.debug(f"  Decision: {result} (no groups matched)")
        return result
    
    result = "mixture"
    logging.debug(f"  Decision: {result} (fallback)")
    return result

def build_cluster_member_headers(rep_to_members: Dict[str, List[str]],
                                 headers_by_id: Dict[str, str],
                                 rep_headers_by_id: Dict[str, str]) -> Dict[str, List[str]]:
    out = {}
    missing = set()
    for rep_id, members in rep_to_members.items():
        rep_full = rep_headers_by_id.get(rep_id) or headers_by_id.get(rep_id) or rep_id
        combined = list(members)
        if rep_id not in combined:
            combined.insert(0, rep_id)
        member_fulls = []
        for mid in combined:
            full = headers_by_id.get(mid) or rep_headers_by_id.get(mid) or mid
            if full == mid and mid not in headers_by_id and mid not in rep_headers_by_id:
                missing.add(mid)
            member_fulls.append(full)
        out[rep_full] = member_fulls
    if missing:
        logging.warning(f"{len(missing)} member IDs not matched to FASTA headers; they will be used as raw IDs.")
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug("Missing member IDs: " + ", ".join(sorted(list(missing))[:50]))
    return out

def write_output_csv(path: Path, rows: List[Dict[str, str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["identifier", "group"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    logging.info(f"Wrote {len(rows)} rows to {path}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--original-fasta", "-o", required=True, type=Path)
    p.add_argument("--rep-fasta", "-r", required=True, type=Path)
    p.add_argument("--cluster-tsv", "-c", required=True, type=Path)
    p.add_argument("--output-csv", "-O", required=True, type=Path)
    p.add_argument("--id-mode", choices=["first_token", "full"], default="first_token")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--debug-cluster", help="Debug specific cluster representative ID")
    args = p.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    headers_by_id = parse_fasta_headers(args.original_fasta, id_mode=args.id_mode)
    rep_headers_by_id = parse_fasta_headers(args.rep_fasta, id_mode=args.id_mode)
    rep_to_members = parse_mmseqs_tsv(args.cluster_tsv)

    clusters = build_cluster_member_headers(rep_to_members, headers_by_id, rep_headers_by_id)

    rows = []
    stats = Counter()
    for rep_full, member_fulls in clusters.items():
        rep_id = rep_full.split()[0]
        
        # ADDED: Debug specific cluster if requested
        if args.debug_cluster and args.debug_cluster == rep_id:
            logging.info(f"=== DEBUG CLUSTER {rep_id} ===")
            logging.getLogger().setLevel(logging.DEBUG)
        
        label = decide_cluster_label(member_fulls, rep_id)
        rows.append({"identifier": rep_id, "group": label})
        stats[label] += 1
        
        if args.debug_cluster and args.debug_cluster == rep_id:
            logging.info(f"=== END DEBUG CLUSTER {rep_id} ===")
            if not args.verbose:
                logging.getLogger().setLevel(logging.INFO)

    # reps in rep fasta but not in TSV -> singletons
    reps_in_tsv = set(rep_to_members.keys())
    reps_in_fasta = set(rep_headers_by_id.keys())
    missing_reps = reps_in_fasta - reps_in_tsv
    if missing_reps:
        logging.info(f"{len(missing_reps)} reps in rep FASTA not listed in TSV -> singletons")
    for rep_id in missing_reps:
        rep_full = rep_headers_by_id.get(rep_id, rep_id)
        label = decide_cluster_label([rep_full], rep_id)
        rows.append({"identifier": rep_full.split()[0], "group": label})
        stats[label] += 1

    write_output_csv(args.output_csv, rows)
    logging.info("Summary:")
    for k, v in stats.most_common():
        logging.info(f"  {k}: {v}")

if __name__ == "__main__":
    main()