#!/usr/bin/env python3
"""
cluster_group_assign_with_taxonomy.py

Adds last_common_ancestor column based on taxonomic classification.
Finds the most specific taxonomic level shared by all cluster members.

Fixes:
 - handles IDs like Centipede3FTx_0001_61981 where \b won't match underscore
 - tries multiple search targets: protein_name, full header, header with underscores -> spaces,
   ID token, ID token with underscores -> spaces
 - logs missing member IDs and unmatched headers
 - FIXED: Better debugging for group assignment decisions
 - FIXED: Added IPR031424 mapping (appears to be missing from your config)
 - ADDED: Monomere detection based on UniProt ID domain counts
 - NEW: Last common ancestor determination from taxonomy CSV
 - NEW: Length binning functionality with configurable number of bins
 - NEW: Configurable minimum taxonomic level for last common ancestor
"""

import argparse
import csv
import logging
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import math

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

# Taxonomic hierarchy (from most general to most specific)
TAXONOMIC_LEVELS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

# Hard-coded mapping from phylum to kingdom and domain
PHYLUM_TO_KINGDOM_DOMAIN = {
    # Animals
    "Arthropoda": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Chordata": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Mollusca": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Cnidaria": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Porifera": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Echinodermata": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Platyhelminthes": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Nematoda": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Annelida": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Rotifera": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Bryozoa": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Brachiopoda": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Ctenophora": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Hemichordata": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Onychophora": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Tardigrada": {"kingdom": "Animalia", "domain": "Eukaryota"},

    # Plants (Streptophyta can cover land plants, but I'll split for clarity)
    "Streptophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Chlorophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Tracheophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},  # vascular plants
    "Bryophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},     # mosses
    "Marchantiophyta": {"kingdom": "Plantae", "domain": "Eukaryota"}, # liverworts
    "Anthocerotophyta": {"kingdom": "Plantae", "domain": "Eukaryota"}, # hornworts
    "Cycadophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Ginkgophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Coniferophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Magnoliophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},  # flowering plants

    # Fungi
    "Ascomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Basidiomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Mucoromycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Chytridiomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Glomeromycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Blastocladiomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Neocallimastigomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},

    # Protists (optional but useful for completeness)
    "Apicomplexa": {"kingdom": "Protista", "domain": "Eukaryota"},
    "Euglenozoa": {"kingdom": "Protista", "domain": "Eukaryota"},
    "Ciliophora": {"kingdom": "Protista", "domain": "Eukaryota"},
    "Diatomea": {"kingdom": "Protista", "domain": "Eukaryota"},  # diatoms
    "Oomycota": {"kingdom": "Protista", "domain": "Eukaryota"},

    # Bacteria
    "Proteobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Firmicutes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Actinobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Bacteroidetes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Cyanobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Chlamydiae": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Spirochaetes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Acidobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Planctomycetes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Verrucomicrobia": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Deinococcus-Thermus": {"kingdom": "Bacteria", "domain": "Bacteria"},

    # Archaea
    "Euryarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Crenarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Thaumarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Korarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Nanoarchaeota": {"kingdom": "Archaea", "domain": "Archaea"}
}

# Centipede genome IDs and their assumed common taxonomy
# These 4 species are treated as having a common taxonomic ancestor
CENTIPEDE_GENOMES = [
    "GCF_000237925.1_ASM23792v2_genomic",
    "GCF_013122585.1_ASM1312258v2_genomic", 
    "GCF_016097555.1_Gae_host_genome_genomic",
    "GCF_024362695.1_ilPecGoss1.1_genomic"
]

# Assumed common taxonomy for all Centipede sequences
# TODO: Replace with actual taxonomic lookup for the 4 genome IDs above
CENTIPEDE_COMMON_TAXONOMY = {
    "domain": "Eukaryota",
    "kingdom": "Animalia",
    "phylum": "Arthropoda",
    "class": "Chilopoda", 
    "order": None,  # This might need to be adjusted based on actual species
    "family": None,  # Mixed families potentially
    "genus": None   # Mixed genera
}

# compile patterns (case-insensitive)
_COMPILED_GROUPS = [(g, [re.compile(p, re.IGNORECASE) for p in pats]) for g, pats in GROUP_PATTERNS.items()]
_re_protein_name = re.compile(r"protein_name=(.*?)\s+(organism=|domain_pos=|signature=|interpro_id=|length=|$)", re.IGNORECASE)
_re_interpro_id = re.compile(r"(IPR\d{6})", re.IGNORECASE)
_re_length = re.compile(r"length=(\d+)", re.IGNORECASE)
# ------------------------------------------------------

def parse_taxonomy_csv(csv_path: Path) -> Dict[str, Dict[str, str]]:
    """Parse taxonomy CSV file and return dict mapping protein_id to taxonomic info"""
    taxonomy_data = {}
    
    if not csv_path.exists():
        raise FileNotFoundError(f"Taxonomy CSV not found: {csv_path}")
    
    with csv_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            protein_id = row.get("protein_id", "").strip()
            if protein_id:
                # Store all taxonomic levels for this protein
                tax_info = {}
                for level in TAXONOMIC_LEVELS:
                    if level in ["domain", "kingdom"]:
                        # Skip domain and kingdom from CSV - we'll calculate them
                        continue
                    tax_info[level] = row.get(level, "").strip() or None
                
                # Add kingdom and domain based on phylum
                phylum = tax_info.get("phylum")
                if phylum and phylum in PHYLUM_TO_KINGDOM_DOMAIN:
                    mapping = PHYLUM_TO_KINGDOM_DOMAIN[phylum]
                    tax_info["kingdom"] = mapping["kingdom"]
                    tax_info["domain"] = mapping["domain"]
                else:
                    tax_info["kingdom"] = None
                    tax_info["domain"] = None
                
                taxonomy_data[protein_id] = tax_info
    
    logging.info(f"Loaded taxonomy data for {len(taxonomy_data)} proteins from {csv_path}")
    return taxonomy_data

def extract_uniprot_id(identifier: str) -> str:
    """Extract UniProt ID from identifier like A0A7J7F872_IPR042339_1_14-139
    If identifier contains 'Centipede', treat it as a unique identifier"""
    if "Centipede" in identifier:
        return identifier  # Treat the whole identifier as unique
    return identifier.split('_')[0]

def extract_length_from_header(header: str) -> Optional[int]:
    """Extract protein length from header using regex"""
    match = _re_length.search(header)
    if match:
        return int(match.group(1))
    return None

def parse_fasta_with_lengths(fasta_path: Path, id_mode: str = "first_token") -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Parse a FASTA file and return:
      - headers_by_id: mapping from seqid -> full header (without '>')
      - sequence_lengths: mapping from seqid -> sequence length (computed from sequences)
    id_mode: "first_token" (default) uses header.split()[0] as seqid, "full" uses full header as key.
    """
    headers = {}
    lengths = {}
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    with fasta_path.open("r", encoding="utf-8") as fh:
        seqid = None
        seq_lines = []
        header_text = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # store previous
                if seqid is not None:
                    seq = "".join(seq_lines)
                    lengths[seqid] = len(seq)
                    seq_lines = []
                header_text = line[1:].rstrip()
                seqid = header_text if id_mode == "full" else header_text.split(None, 1)[0]
                headers[seqid] = header_text
            else:
                seq_lines.append(line.strip())
        # last record
        if seqid is not None:
            seq = "".join(seq_lines)
            lengths[seqid] = len(seq)
    logging.info(f"Parsed {len(headers)} headers from {fasta_path} (id_mode={id_mode}), computed sequence lengths for {len(lengths)} sequences")
    return headers, lengths

def calculate_cluster_average_length(member_headers: List[str], sequence_lengths: Dict[str, int]) -> Optional[float]:
    """Calculate average length of proteins in a cluster.

    member_headers: list of full headers or IDs (the entries produced by build_cluster_member_headers)
    sequence_lengths: mapping seqid -> length (as produced by parse_fasta_with_lengths)
    """
    lengths = []
    for header in member_headers:
        # Try to extract 'length=' from header first
        length = None
        if isinstance(header, str):
            length = extract_length_from_header(header)
        if length is None:
            # try to use the first token as an ID to look up computed sequence lengths
            seqid = header.split(None, 1)[0] if isinstance(header, str) else str(header)
            # If seqid is a full header that was stored as key in sequence_lengths, it will match
            if seqid in sequence_lengths:
                length = sequence_lengths[seqid]
            else:
                # attempt replacement of spaces/underscores variations
                if seqid.replace("_", " ") in sequence_lengths:
                    length = sequence_lengths[seqid.replace("_", " ")]
        if length is not None:
            lengths.append(length)
        else:
            logging.debug(f"Could not determine length for member header/ID: '{header[:60]}'")
    if lengths:
        avg_length = sum(lengths) / len(lengths)
        logging.debug(f"Cluster average length: {avg_length:.1f} (from {len(lengths)} sequences)")
        return avg_length
    else:
        logging.debug("No length information found in cluster headers")
        return None

def create_length_bins(average_lengths: List[float], num_bins: int = 15) -> Tuple[List[Tuple[float, float]], Dict[float, str]]:
    """Create length bins and return bin ranges and mapping from average to bin label"""
    if not average_lengths:
        return [], {}
    
    min_length = min(average_lengths)
    max_length = max(average_lengths)
    
    # If all lengths are equal, create a single bin
    if min_length == max_length:
        bin_label = f"{int(min_length)}-{int(max_length)}"
        return [(min_length, max_length)], {avg: bin_label for avg in average_lengths}
    
    # Calculate bin size
    bin_size = (max_length - min_length) / num_bins
    
    # Create bin ranges
    bin_ranges = []
    length_to_bin = {}
    
    for i in range(num_bins):
        bin_start = min_length + (i * bin_size)
        bin_end = min_length + ((i + 1) * bin_size)
        
        # For the last bin, make sure it includes the maximum value
        if i == num_bins - 1:
            bin_end = max_length
        
        bin_ranges.append((bin_start, bin_end))
        bin_label = f"{int(bin_start)}-{int(bin_end)}"
        
        # Assign each average length to its corresponding bin
        for avg_len in average_lengths:
            # Use <= for upper bound on the last bin, otherwise < to avoid overlap
            if i == num_bins - 1:
                if bin_start <= avg_len <= bin_end:
                    length_to_bin[avg_len] = bin_label
            else:
                if bin_start <= avg_len < bin_end:
                    length_to_bin[avg_len] = bin_label
    
    logging.info(f"Created {num_bins} length bins from {min_length:.1f} to {max_length:.1f}")
    logging.debug(f"Bin ranges: {[(int(start), int(end)) for start, end in bin_ranges]}")
    
    return bin_ranges, length_to_bin

def find_last_common_ancestor(member_ids: List[str], taxonomy_data: Dict[str, Dict[str, str]], min_level: str = "family") -> str:
    """Find the most specific taxonomic level shared by all cluster members, with configurable minimum level
    
    min_level defines the most specific level to consider. For example:
    - If min_level="phylum", only check phylum, class, order, family, genus (skip more general levels)
    - If min_level="order", check order, family, genus (skip phylum, class and more general levels)
    """
    if not member_ids:
        return "unknown"
    
    # Validate min_level
    if min_level not in TAXONOMIC_LEVELS:
        logging.warning(f"Invalid min_level '{min_level}', using default 'family'")
        min_level = "family"
    
    # Get the index of the minimum level (most specific allowed level)
    min_level_index = TAXONOMIC_LEVELS.index(min_level)
    
    # Separate Centipede IDs from UniProt IDs
    uniprot_ids = []
    has_centipedes = False
    
    for member_id in member_ids:
        # Extract just the ID part (before any spaces)
        seq_id = member_id.split()[0]
        uniprot_id = extract_uniprot_id(seq_id)
        
        if "Centipede" in uniprot_id:
            has_centipedes = True
        else:
            uniprot_ids.append(uniprot_id)
    
    # Collect all taxonomies (including Centipede common taxonomy if present)
    taxonomies = []
    
    # Add Centipede common taxonomy if any Centipede sequences are present
    if has_centipedes:
        taxonomies.append(CENTIPEDE_COMMON_TAXONOMY)
        logging.debug("Added Centipede common taxonomy to analysis")
    
    # Add UniProt taxonomies
    missing_ids = []
    for uniprot_id in uniprot_ids:
        if uniprot_id in taxonomy_data:
            taxonomies.append(taxonomy_data[uniprot_id])
        else:
            missing_ids.append(uniprot_id)
    
    if missing_ids:
        logging.debug(f"Missing taxonomy data for {len(missing_ids)} UniProt IDs: {missing_ids[:5]}...")
    
    if not taxonomies:
        logging.debug("No taxonomy data found for any cluster members")
        return "unknown"
    
    # Create the list of levels to check, starting from min_level (most specific allowed) 
    # and going to more general levels
    levels_to_check = TAXONOMIC_LEVELS[min_level_index:]
    logging.debug(f"Checking taxonomic levels from {min_level} onwards: {levels_to_check}")
    
    # If only Centipedes in the cluster, return their common ancestor (respecting min_level)
    if len(taxonomies) == 1 and has_centipedes and not uniprot_ids:
        # Start from the most specific allowed level (min_level) and work to more general
        for level in levels_to_check:
            value = CENTIPEDE_COMMON_TAXONOMY.get(level)
            if value is not None and value != "":
                logging.debug(f"Centipede-only cluster, returning {level}: {value}")
                return value
        return "unknown"
    
    # Find the most specific level where all members agree, starting from min_level
    # and working towards more general levels
    for level in levels_to_check:
        values = set()
        all_have_value = True
        
        for tax_info in taxonomies:
            value = tax_info.get(level)
            if value is None or value == "":
                all_have_value = False
                break
            values.add(value)
        
        # If all members have a value for this level and they're all the same
        if all_have_value and len(values) == 1:
            result = values.pop()
            logging.debug(f"Found common ancestor at {level} level: {result}")
            return result
    
    logging.debug(f"No common taxonomic level found at or above {min_level}")
    return "unknown"

def count_uniprot_domains(headers_by_id: Dict[str, str]) -> Dict[str, int]:
    """Count how many domains each UniProt ID has in the dataset
    Centipede IDs are treated as unique (count = 1)"""
    uniprot_counts = Counter()
    centipede_count = 0
    
    for seq_id in headers_by_id.keys():
        uniprot_id = extract_uniprot_id(seq_id)
        if "Centipede" in uniprot_id:
            centipede_count += 1
            # Each Centipede ID is unique, so we set count to 1
            uniprot_counts[uniprot_id] = 1
        else:
            uniprot_counts[uniprot_id] += 1
    
    logging.info(f"Counted domains for {len(uniprot_counts)} UniProt IDs")
    logging.info(f"Found {centipede_count} Centipede sequences (treated as unique)")
    multimeric_proteins = sum(1 for count in uniprot_counts.values() if count > 1)
    logging.info(f"Found {multimeric_proteins} proteins with multiple domains")
    
    return dict(uniprot_counts)

def determine_monomere_status(member_ids: List[str], uniprot_counts: Dict[str, int]) -> str:
    """Determine if cluster members are all monomeric, all multimeric, or mixed
    Centipede IDs are treated as monomeric (unique)"""
    domain_types = []
    
    for member_id in member_ids:
        # Extract just the ID part (before any spaces)
        seq_id = member_id.split()[0]
        uniprot_id = extract_uniprot_id(seq_id)
        
        # Centipede IDs are always treated as unique/monomeric
        if "Centipede" in uniprot_id:
            domain_types.append("monomeric")
        else:
            domain_count = uniprot_counts.get(uniprot_id, 1)  # Default to 1 if not found
            if domain_count == 1:
                domain_types.append("monomeric")
            else:
                domain_types.append("multimeric")
    
    unique_types = set(domain_types)
    
    if len(unique_types) == 1:
        if "monomeric" in unique_types:
            return "monomere"
        else:
            return "multimere"
    else:
        return "mixture"

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
        w = csv.DictWriter(fh, fieldnames=["identifier", "group", "Oligomeric_state", "last_common_ancestor", "length"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    logging.info(f"Wrote {len(rows)} rows to {path}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--original-fasta", "-o", required=True, type=Path)
    p.add_argument("--rep-fasta", "-r", required=True, type=Path)
    p.add_argument("--cluster-tsv", "-c", required=True, type=Path)
    p.add_argument("--taxonomy-csv", "-t", required=True, type=Path, help="CSV file with taxonomy data")
    p.add_argument("--output-csv", "-O", required=True, type=Path)
    p.add_argument("--id-mode", choices=["first_token", "full"], default="first_token")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--debug-cluster", help="Debug specific cluster representative ID")
    p.add_argument("--num-bins", type=int, default=15, help="Number of length bins to create (default: 15)")
    p.add_argument("--min-tax-level", default="family", choices=TAXONOMIC_LEVELS, 
                   help=f"Minimum taxonomic level for last common ancestor (default: family)")
    args = p.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Parse FASTA files with sequence lengths
    headers_by_id, sequence_lengths = parse_fasta_with_lengths(args.original_fasta, id_mode=args.id_mode)
    rep_headers_by_id, rep_sequence_lengths = parse_fasta_with_lengths(args.rep_fasta, id_mode=args.id_mode)
    
    # Combine sequence lengths from both files
    all_sequence_lengths = {**sequence_lengths, **rep_sequence_lengths}
    
    rep_to_members = parse_mmseqs_tsv(args.cluster_tsv)
    
    # NEW: Load taxonomy data
    taxonomy_data = parse_taxonomy_csv(args.taxonomy_csv)

    # Count UniProt domain occurrences
    uniprot_counts = count_uniprot_domains(headers_by_id)

    clusters = build_cluster_member_headers(rep_to_members, headers_by_id, rep_headers_by_id)

    # NEW: First pass to collect all cluster average lengths
    cluster_avg_lengths = []
    cluster_length_data = {}  # Store rep_id -> avg_length mapping
    
    for rep_full, member_fulls in clusters.items():
        rep_id = rep_full.split()[0]
        avg_length = calculate_cluster_average_length(member_fulls, all_sequence_lengths)
        if avg_length is not None:
            cluster_avg_lengths.append(avg_length)
            cluster_length_data[rep_id] = avg_length
    
    # Process singletons for length calculation
    reps_in_tsv = set(rep_to_members.keys())
    reps_in_fasta = set(rep_headers_by_id.keys())
    missing_reps = reps_in_fasta - reps_in_tsv
    
    for rep_id in missing_reps:
        rep_full = rep_headers_by_id.get(rep_id, rep_id)
        avg_length = calculate_cluster_average_length([rep_full], all_sequence_lengths)
        if avg_length is not None:
            cluster_avg_lengths.append(avg_length)
            cluster_length_data[rep_id] = avg_length

    # NEW: Create length bins
    bin_ranges, length_to_bin = create_length_bins(cluster_avg_lengths, args.num_bins)

    rows = []
    stats = Counter()
    monomere_stats = Counter()
    taxonomy_stats = Counter()
    length_stats = Counter()
    
    for rep_full, member_fulls in clusters.items():
        rep_id = rep_full.split()[0]
        
        # ADDED: Debug specific cluster if requested
        if args.debug_cluster and args.debug_cluster == rep_id:
            logging.info(f"=== DEBUG CLUSTER {rep_id} ===")
            logging.getLogger().setLevel(logging.DEBUG)
        
        label = decide_cluster_label(member_fulls, rep_id)
        
        # Determine monomere status
        oligomeric_status = determine_monomere_status(member_fulls, uniprot_counts)
        
        # NEW: Determine last common ancestor with configurable min level
        last_common_ancestor = find_last_common_ancestor(member_fulls, taxonomy_data, args.min_tax_level)
        
        # NEW: Get length bin
        avg_length = cluster_length_data.get(rep_id)
        length_bin = length_to_bin.get(avg_length, "unknown") if avg_length is not None else "unknown"
        
        rows.append({
            "identifier": rep_id, 
            "group": label,
            "Oligomeric_state": oligomeric_status,
            "last_common_ancestor": last_common_ancestor,
            "length": length_bin
        })
        stats[label] += 1
        monomere_stats[oligomeric_status] += 1
        taxonomy_stats[last_common_ancestor] += 1
        length_stats[length_bin] += 1
        
        if args.debug_cluster and args.debug_cluster == rep_id:
            logging.info(f"=== END DEBUG CLUSTER {rep_id} ===")
            if not args.verbose:
                logging.getLogger().setLevel(logging.INFO)

    # reps in rep fasta but not in TSV -> singletons
    if missing_reps:
        logging.info(f"{len(missing_reps)} reps in rep FASTA not listed in TSV -> singletons")
    for rep_id in missing_reps:
        rep_full = rep_headers_by_id.get(rep_id, rep_id)
        label = decide_cluster_label([rep_full], rep_id)
        
        # Determine monomere status for singletons
        oligomeric_status = determine_monomere_status([rep_full], uniprot_counts)
        
        # NEW: Determine last common ancestor for singletons
        last_common_ancestor = find_last_common_ancestor([rep_full], taxonomy_data, args.min_tax_level)
        
        # NEW: Get length bin for singletons
        avg_length = cluster_length_data.get(rep_id)
        length_bin = length_to_bin.get(avg_length, "unknown") if avg_length is not None else "unknown"
        
        rows.append({
            "identifier": rep_full.split()[0], 
            "group": label,
            "Oligomeric_state": oligomeric_status,
            "last_common_ancestor": last_common_ancestor,
            "length": length_bin
        })
        stats[label] += 1
        monomere_stats[oligomeric_status] += 1
        taxonomy_stats[last_common_ancestor] += 1
        length_stats[length_bin] += 1

    write_output_csv(args.output_csv, rows)
    logging.info("Group Summary:")
    for k, v in stats.most_common():
        logging.info(f"  {k}: {v}")
    
    logging.info("Oligomeric State Summary:")
    for k, v in monomere_stats.most_common():
        logging.info(f"  {k}: {v}")
    
    logging.info("Taxonomy Summary:")
    for k, v in taxonomy_stats.most_common():
        logging.info(f"  {k}: {v}")
    
    logging.info("Length Bin Summary:")
    for k, v in length_stats.most_common():
        logging.info(f"  {k}: {v}")

if __name__ == "__main__":
    main()