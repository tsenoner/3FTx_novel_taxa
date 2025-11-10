#!/usr/bin/env python3
"""
Annotate Protein Clusters

This script annotates protein clusters with multiple characteristics:
  1. Group assignments: Pattern-based classification (3FTx, Ly6, PMF, etc.)
  2. Oligomeric state: Monomeric, multimeric, or mixture
  3. Taxonomy: Lowest common ancestor (LCA) at specified taxonomic levels
  4. Length bins: Protein length distribution

The script uses:
  - MMseqs2 cluster files for cluster membership
  - FASTA headers for protein metadata
  - taxopy for NCBI taxonomy database queries

Input:
  - Merged FASTA file with all sequences
  - Representative FASTA file with cluster representatives
  - MMseqs2 cluster TSV mapping representatives to members

Output:
  - CSV file with annotations for each cluster
  - Statistics log file with distribution summaries
"""

import argparse
import csv
import logging
import os
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import taxopy
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Centipede genome prefixes for manually annotated sequences
CENTIPEDE_PREFIXES = ("Rimm", "Lvar", "Sacu", "Cpnu")

# Pattern matching rules for protein group classification
# Note: All patterns are compiled as case-insensitive
GROUP_PATTERNS = {
    # Manually annotated sequences from centipede genomes
    "manually_annotated": [
        r"^(Rimm|Lvar|Sacu|Cpnu)_"  # Must start with centipede prefix
    ],
    # Three-finger toxins (3FTx)
    "3FTx": [
        r"\b3-?ftx\b",  # 3ftx, 3-ftx
        r"\bthree[- ]?ftx\b",  # three ftx, three-ftx
        r"\bthree[- ]?f(?:inger)?[- ]?toxins?\b",  # three finger toxin(s)
    ],
    # Pheromone/Ly6-like proteins matching Factor (PMF)
    "PMF": [
        r"\bpmf\b",  # pmf (any case)
    ],
    # Lymphocyte antigen 6 (Ly6) proteins
    "Ly6": [
        r"\bly-?6\b",  # ly6, ly-6 (any case)
        r"\blysix\b",  # lysix (spelled out)
    ],
    # Quiver proteins
    "Quiver": [
        r"\bquiver\b",  # quiver (any case)
        r"\bqvr\b",  # qvr abbreviation
    ],
    # Scoloptoxin
    "Scoloptoxin": [
        r"\bscoloptoxin\b",  # scoloptoxin (any case)
    ],
    # Sodefrin Precursor-like Factor (SPF)
    "SPF": [
        r"\bsodefrin\b",  # sodefrin (any case)
        r"\bspf\b",  # spf abbreviation
    ],
}

# Specificity ranking for resolving conflicts when multiple patterns match
# Lower index = higher priority (most specific first)
SPECIFICITY_ORDER = [
    "SPF",
    "Scoloptoxin",
    "PMF",
    "Quiver",
    "manually_annotated",
    "3FTx",
    "Ly6",
]
SPECIFICITY_RANK = {name: i for i, name in enumerate(SPECIFICITY_ORDER)}

# Taxonomic hierarchy levels (from general to specific)
TAXONOMIC_LEVELS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

# ============================================================================
# COMPILED REGEX PATTERNS
# ============================================================================

# Compile group patterns once for efficiency (case-insensitive)
_COMPILED_GROUPS = [
    (group_name, [re.compile(pattern, re.IGNORECASE) for pattern in patterns])
    for group_name, patterns in GROUP_PATTERNS.items()
]

# FASTA header field extraction patterns
_re_protein_name = re.compile(
    r"protein_name=([^=]+?)(?=\s+(?:organism=|domain_pos=|signature=|interpro_id=|length=)|$)",
    re.IGNORECASE,
)
_re_organism = re.compile(
    r"organism=([^=]+?)(?=\s+domain_pos=)",
    re.IGNORECASE,
)
_re_interpro_id = re.compile(r"\b(IPR\d{6})\b", re.IGNORECASE)
_re_length = re.compile(r"\blength=(\d+)\b", re.IGNORECASE)


# ============================================================================
# TAXONOMY DATABASE INITIALIZATION
# ============================================================================


def initialize_taxdb() -> taxopy.TaxDb:
    """
    Initialize or load the taxopy taxonomy database with intelligent caching.

    The database is cached locally and automatically refreshed if older than one week.
    Uses safe atomic file operations to prevent corruption during updates.

    Environment variable:
        PROTSPACE_TAXDB_DIR: Override default cache directory

    Returns:
        taxopy.TaxDb: Initialized taxonomy database

    Raises:
        Exception: If database cannot be downloaded or loaded
    """
    env_override = os.environ.get("PROTSPACE_TAXDB_DIR")
    db_dir = (
        Path(env_override).expanduser()
        if env_override
        else Path.home() / ".cache" / "taxopy_db"
    )
    db_dir.mkdir(parents=True, exist_ok=True)

    nodes_file = db_dir / "nodes.dmp"
    names_file = db_dir / "names.dmp"
    merged_file = db_dir / "merged.dmp"
    timestamp_file = db_dir / ".download_timestamp"

    # Determine if this is a first-time setup
    first_time_setup = not (nodes_file.exists() and names_file.exists())

    # Check if cache needs refresh based on timestamp file
    needs_refresh = False
    if timestamp_file.exists():
        try:
            with open(timestamp_file) as f:
                download_time = datetime.fromisoformat(f.read().strip())
            one_week_ago = datetime.now() - timedelta(weeks=1)
            if download_time < one_week_ago:
                logger.info(
                    "Your taxonomy dataset is more than one week old. Refreshing cache..."
                )
                needs_refresh = True
        except (ValueError, OSError) as e:
            logger.warning(f"Could not read timestamp file: {e}. Will refresh cache.")
            needs_refresh = True
    else:
        if first_time_setup:
            needs_refresh = True
        else:
            try:
                with open(timestamp_file, "w") as f:
                    f.write(datetime.now().isoformat())
            except OSError as e:
                logger.warning(
                    f"Failed to create timestamp file at first-time detection: {e}"
                )

    existing_db_present = nodes_file.exists() and names_file.exists()

    # Load or download the database with a safe refresh strategy
    if existing_db_present:
        if needs_refresh:
            logger.info(
                "Taxonomy cache is stale. Attempting safe refresh without deleting existing cache."
            )
            temp_dir_path = None
            try:
                # Download into a temporary directory first
                temp_dir_path = Path(tempfile.mkdtemp(prefix="taxopy_tmp_"))
                taxopy.TaxDb(taxdb_dir=str(temp_dir_path), keep_files=True)
                # Move refreshed files into place atomically
                for src_name, dst_path in [
                    ("nodes.dmp", nodes_file),
                    ("names.dmp", names_file),
                    ("merged.dmp", merged_file),
                ]:
                    src_path = temp_dir_path / src_name
                    if src_path.exists():
                        shutil.move(str(src_path), str(dst_path))
                # Update timestamp only after a successful refresh
                with open(timestamp_file, "w") as f:
                    f.write(datetime.now().isoformat())
            except Exception as e:
                logger.warning(
                    f"Failed to refresh taxonomy database: {e}. Falling back to existing cached database."
                )
            finally:
                if temp_dir_path and temp_dir_path.exists():
                    shutil.rmtree(temp_dir_path, ignore_errors=True)

        # Load existing (potentially refreshed) DB files
        logger.info(f"Loading taxopy database from {db_dir}")
        try:
            taxdb = taxopy.TaxDb(
                nodes_dmp=str(nodes_file),
                names_dmp=str(names_file),
                merged_dmp=str(merged_file) if merged_file.exists() else None,
            )
        except Exception as e:
            logger.error(f"Failed to load existing taxonomy database from cache: {e}")
            raise
    else:
        # First-time setup: must download
        logger.info(f"Downloading taxopy database to {db_dir}")
        try:
            taxdb = taxopy.TaxDb(taxdb_dir=str(db_dir), keep_files=True)
            # Create/update timestamp file after successful download
            with open(timestamp_file, "w") as f:
                f.write(datetime.now().isoformat())
        except Exception as e:
            logger.error(
                f"Failed to initialize taxopy database (first-time setup): {e}"
            )
            raise

    return taxdb


# ============================================================================
# FASTA HEADER PARSING
# ============================================================================


def extract_organism_name(header: str) -> Optional[str]:
    """
    Extract organism scientific name from FASTA header.

    Removes common names in parentheses and extra whitespace.

    Args:
        header: Full FASTA header line (without '>')

    Returns:
        Cleaned scientific name, or None if not found

    Example:
        'organism=Homo sapiens (Human)' -> 'Homo sapiens'
        'organism=Mus musculus domain_pos=...' -> 'Mus musculus'
    """
    match = _re_organism.search(header)
    if match:
        organism_full = match.group(1).strip()
        # Remove parenthetical common names
        organism_clean = re.sub(r"\s*\([^)]*\)\s*", "", organism_full).strip()
        return organism_clean
    return None


def get_taxonomy_from_taxopy(organism_name: str, taxdb: taxopy.TaxDb) -> Dict[str, str]:
    """Retrieve taxonomy information for an organism using taxopy.

    Returns a dict with keys: domain, kingdom, phylum, class, order, family, genus
    """
    taxonomy_info = {level: None for level in TAXONOMIC_LEVELS}

    try:
        # Get taxon ID from organism name
        taxid_list = taxopy.taxid_from_name(organism_name, taxdb)

        if not taxid_list or not taxid_list[0]:
            logger.debug(f"No taxon ID found for organism: {organism_name}")
            return taxonomy_info

        # Use the first taxon ID (handle homonyms by taking the first match)
        taxon_id = taxid_list[0]

        # Get taxon object
        taxon = taxopy.Taxon(taxon_id, taxdb)
        ranks = taxon.rank_name_dictionary

        # Map taxopy ranks to our expected levels
        taxonomy_info["domain"] = (
            ranks.get("domain", "")
            or ranks.get("realm", "")
            or ranks.get("superkingdom", "")
        )
        taxonomy_info["kingdom"] = ranks.get("kingdom", "")
        taxonomy_info["phylum"] = ranks.get("phylum", "")
        taxonomy_info["class"] = ranks.get("class", "")
        taxonomy_info["order"] = ranks.get("order", "")
        taxonomy_info["family"] = ranks.get("family", "")
        taxonomy_info["genus"] = ranks.get("genus", "")

        # Convert empty strings to None
        taxonomy_info = {k: v if v else None for k, v in taxonomy_info.items()}

        logger.debug(f"Retrieved taxonomy for {organism_name}: {taxonomy_info}")

    except Exception as e:
        logger.warning(f"Failed to get taxonomy for '{organism_name}': {e}")

    return taxonomy_info


def batch_get_taxonomy_for_organisms(
    organism_names: List[str], taxdb: taxopy.TaxDb
) -> Dict[str, Optional[taxopy.Taxon]]:
    """Batch retrieve taxonomy Taxon objects for multiple organisms.

    Returns a dict mapping organism_name -> Taxon object (or None if not found)
    """
    organism_to_taxon = {}

    # Use batch lookup for taxon IDs (more efficient than one-by-one)
    try:
        # taxopy.taxid_from_name can handle a list of names
        all_taxid_lists = taxopy.taxid_from_name(organism_names, taxdb)
    except Exception as e:
        logger.warning(
            f"Batch taxid lookup failed: {e}. Falling back to individual lookups."
        )
        all_taxid_lists = None

    # Process each organism with progress bar
    with tqdm(
        total=len(organism_names), desc="Fetching taxonomy", unit="organism"
    ) as pbar:
        for i, organism_name in enumerate(organism_names):
            try:
                # Get taxon ID from batch results or individual lookup
                if all_taxid_lists is not None:
                    taxid_list = all_taxid_lists[i] if i < len(all_taxid_lists) else []
                else:
                    # Fallback to individual lookup
                    taxid_list = taxopy.taxid_from_name(organism_name, taxdb)

                if taxid_list and taxid_list[0]:
                    # Use the first taxon ID (handle homonyms by taking the first match)
                    taxon_id = taxid_list[0]
                    taxon = taxopy.Taxon(taxon_id, taxdb)
                    organism_to_taxon[organism_name] = taxon
                    logger.debug(f"Retrieved Taxon for {organism_name}: {taxon.name}")
                else:
                    logger.debug(f"No taxon ID found for organism: {organism_name}")
                    organism_to_taxon[organism_name] = None
            except Exception as e:
                logger.warning(f"Failed to get taxonomy for '{organism_name}': {e}")
                organism_to_taxon[organism_name] = None

            pbar.update(1)

    successful = sum(1 for t in organism_to_taxon.values() if t is not None)
    logger.info(
        f"Successfully retrieved taxonomy for {successful}/{len(organism_names)} organisms"
    )
    return organism_to_taxon


def extract_all_organism_names(
    headers_by_id: Dict[str, str],
) -> Dict[str, str]:
    """Extract all organism names from headers and map seq_id -> organism_name.

    Returns a dict mapping seq_id -> organism_name
    """
    organism_by_seqid = {}
    for seq_id, header in headers_by_id.items():
        organism_name = extract_organism_name(header)
        if organism_name:
            organism_by_seqid[seq_id] = organism_name
    return organism_by_seqid


def extract_length_from_header(header: str) -> Optional[int]:
    """
    Extract protein length from FASTA header.

    Args:
        header: Full FASTA header line

    Returns:
        Length as integer, or None if not found
    """
    match = _re_length.search(header)
    if match:
        return int(match.group(1))
    return None


def is_centipede_sequence(identifier: str) -> bool:
    """
    Check if a sequence identifier is from a centipede genome annotation.

    Centipede sequences use specific prefixes defined in CENTIPEDE_PREFIXES:
      - Rimm: Strigamia maritima
      - Lvar: Lithobius variegatus
      - Sacu: Scolopendra cingulata
      - Cpnu: Cryptops punctatus

    Args:
        identifier: Sequence identifier (first token of FASTA header)

    Returns:
        True if identifier starts with a centipede prefix
    """
    return identifier.startswith(CENTIPEDE_PREFIXES)


def extract_uniprot_id(identifier: str) -> str:
    """
    Extract UniProt accession from a sequence identifier.

    For most sequences, extracts the first underscore-delimited token.
    For centipede genome sequences, returns the full identifier as they are unique.

    Args:
        identifier: Sequence identifier (e.g., 'A0A7J7F872_IPR042339_1_14-139')

    Returns:
        UniProt ID or full identifier for centipede sequences

    Examples:
        'A0A7J7F872_IPR042339_1_14-139' -> 'A0A7J7F872'
        'Rimm_gene123_domain1' -> 'Rimm_gene123_domain1' (kept as unique)
    """
    if is_centipede_sequence(identifier):
        return identifier  # Centipede IDs are treated as unique
    return identifier.split("_")[0]


# ============================================================================
# FASTA FILE PARSING
# ============================================================================


def parse_fasta_with_lengths(
    fasta_path: Path,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Parse FASTA file and extract headers and sequence lengths.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Tuple of:
          - Dict mapping sequence_id -> full header (without '>')
          - Dict mapping sequence_id -> sequence length (AA count)

    Raises:
        FileNotFoundError: If FASTA file doesn't exist
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
                seqid = header_text.split(None, 1)[0]
                headers[seqid] = header_text
            else:
                seq_lines.append(line.strip())
        # last record
        if seqid is not None:
            seq = "".join(seq_lines)
            lengths[seqid] = len(seq)
    logger.debug(
        f"Parsed {len(headers)} headers from {fasta_path}, computed sequence lengths for {len(lengths)} sequences"
    )
    return headers, lengths


# ============================================================================
# LENGTH ANALYSIS
# ============================================================================


def calculate_cluster_average_length(
    member_headers: List[str], sequence_lengths: Dict[str, int]
) -> Optional[float]:
    """
    Calculate average protein length for a cluster.

    Tries to extract length from header first, falls back to computed sequence lengths.

    Args:
        member_headers: List of FASTA headers for cluster members
        sequence_lengths: Pre-computed sequence lengths by ID

    Returns:
        Average length in amino acids, or None if no lengths found
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
            if seqid in sequence_lengths:
                length = sequence_lengths[seqid]
        if length is not None:
            lengths.append(length)
        else:
            logger.debug(
                f"Could not determine length for member header/ID: '{header[:60]}'"
            )
    if lengths:
        avg_length = sum(lengths) / len(lengths)
        logger.debug(
            f"Cluster average length: {avg_length:.1f} (from {len(lengths)} sequences)"
        )
        return avg_length
    else:
        logger.debug("No length information found in cluster headers")
        return None


def create_length_bins(
    average_lengths: List[float], num_bins: int = 15
) -> Tuple[List[Tuple[float, float]], Dict[float, str]]:
    """
    Create evenly-spaced length bins for protein length distribution.

    Args:
        average_lengths: List of average cluster lengths
        num_bins: Number of bins to create (default: 15)

    Returns:
        Tuple of:
          - List of (start, end) tuples for each bin
          - Dict mapping each average_length -> bin_label string

    Note:
        Bin labels are formatted as "start-end" (e.g., "50-72")
    """
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

    logger.debug(
        f"Created {num_bins} length bins from {min_length:.1f} to {max_length:.1f}"
    )
    logger.debug(f"Bin ranges: {[(int(start), int(end)) for start, end in bin_ranges]}")

    return bin_ranges, length_to_bin


# ============================================================================
# TAXONOMY ANALYSIS
# ============================================================================


def find_last_common_ancestor(
    member_ids: List[str],
    organism_by_seqid: Dict[str, str],
    organism_to_taxon: Dict[str, taxopy.Taxon],
    taxdb: taxopy.TaxDb,
    min_levels: List[str],
) -> Dict[str, Tuple[str, str]]:
    """
    Determine the lowest common ancestor (LCA) at specified taxonomic levels.

    For clusters with multiple organisms, uses taxopy's LCA algorithm.
    For single-organism clusters, returns that organism's taxonomy.

    Args:
        member_ids: List of sequence identifiers in the cluster
        organism_by_seqid: Mapping of sequence_id -> organism_name
        organism_to_taxon: Mapping of organism_name -> taxopy.Taxon
        taxdb: taxopy taxonomy database
        min_levels: Taxonomic levels to report (e.g., ['family', 'genus'])

    Returns:
        Dict mapping each level -> (taxon_name, actual_rank)
        Returns ("unknown", "unknown") for levels where taxonomy is unavailable

    Note:
        If the true LCA is more specific than the requested level, returns
        the ancestor at the requested level instead.
    """
    if not member_ids:
        return {level: ("unknown", "unknown") for level in min_levels}

    # Validate min_levels
    valid_levels = []
    for level in min_levels:
        if level not in TAXONOMIC_LEVELS:
            logger.warning(f"Invalid min_level '{level}', skipping")
        else:
            valid_levels.append(level)

    if not valid_levels:
        logger.warning("No valid taxonomic levels provided, using default 'family'")
        valid_levels = ["family"]

    # Collect taxon objects for all members
    taxon_objects = []
    missing_count = 0

    for member_id in member_ids:
        # Extract just the ID part (before any spaces)
        seq_id = member_id.split()[0]

        # Look up organism name and then taxon object
        organism_name = organism_by_seqid.get(seq_id)
        if organism_name and organism_name in organism_to_taxon:
            taxon_obj = organism_to_taxon[organism_name]
            if taxon_obj is not None:
                taxon_objects.append(taxon_obj)
            else:
                missing_count += 1
        else:
            missing_count += 1

    if missing_count > 0:
        logger.debug(f"Missing taxonomy data for {missing_count} cluster members")

    if not taxon_objects:
        logger.debug("No taxonomy data found for any cluster members")
        return {level: ("unknown", "unknown") for level in valid_levels}

    # Get taxonomy information for all requested levels
    result = {}

    # Handle single taxon case separately (no LCA needed)
    if len(taxon_objects) == 1:
        ranks = taxon_objects[0].rank_name_dictionary
        logger.debug(f"Single taxon cluster: {taxon_objects[0].name}")

        for level in valid_levels:
            level_index = TAXONOMIC_LEVELS.index(level)
            # Try to get the value at this level, or go up the hierarchy
            found = False
            for search_level in TAXONOMIC_LEVELS[level_index:]:
                value = ranks.get(search_level, "")
                if value:
                    logger.debug(
                        f"Single taxon cluster, returning {search_level}: {value} for level {level}"
                    )
                    result[level] = (value, search_level)
                    found = True
                    break
            if not found:
                result[level] = ("unknown", "unknown")
        return result

    # Multiple taxa: compute LCA using taxopy
    try:
        # Use taxopy's find_lca to get the true lowest common ancestor
        lca_taxon = taxopy.find_lca(taxon_objects, taxdb)
        lca_ranks = lca_taxon.rank_name_dictionary

        logger.debug(f"Taxopy LCA: {lca_taxon.name} (rank: {lca_taxon.rank})")

        # For each requested level, find the appropriate ancestor
        for level in valid_levels:
            level_index = TAXONOMIC_LEVELS.index(level)

            # Find the appropriate level to return based on min_level constraint
            # Start from min_level and go to more general levels
            found = False
            for search_level in TAXONOMIC_LEVELS[level_index:]:
                value = lca_ranks.get(search_level, "")
                if value:
                    logger.debug(
                        f"Found common ancestor at {search_level} level: {value} for level {level}"
                    )
                    result[level] = (value, search_level)
                    found = True
                    break

            if not found:
                # If no match at or above min_level, return the LCA name and rank itself
                logger.debug(
                    f"Returning LCA name directly: {lca_taxon.name} for level {level}"
                )
                result[level] = (lca_taxon.name, lca_taxon.rank)

    except Exception as e:
        logger.warning(f"Failed to compute LCA using taxopy: {e}")
        return {level: ("unknown", "unknown") for level in valid_levels}

    return result


# ============================================================================
# OLIGOMERIC STATE DETERMINATION
# ============================================================================


def count_uniprot_domains(headers_by_id: Dict[str, str]) -> Dict[str, int]:
    """
    Count how many domains each UniProt ID has in the dataset.

    Used to determine if proteins are multi-domain (multimeric) or single-domain (monomeric).
    Centipede genome sequences are always treated as unique/monomeric.

    Args:
        headers_by_id: Dict mapping sequence_id -> header

    Returns:
        Dict mapping UniProt_ID -> domain_count

    Note:
        Centipede genome IDs always get count = 1 (treated as unique)
    """
    uniprot_counts = Counter()
    centipede_count = 0

    for seq_id in headers_by_id.keys():
        uniprot_id = extract_uniprot_id(seq_id)
        if is_centipede_sequence(uniprot_id):
            centipede_count += 1
            # Each centipede genome ID is unique, so we set count to 1
            uniprot_counts[uniprot_id] = 1
        else:
            uniprot_counts[uniprot_id] += 1

    logger.debug(f"Counted domains for {len(uniprot_counts)} UniProt IDs")
    logger.debug(
        f"Found {centipede_count} centipede genome sequences (treated as unique)"
    )
    multimeric_proteins = sum(1 for count in uniprot_counts.values() if count > 1)
    logger.debug(f"Found {multimeric_proteins} proteins with multiple domains")

    return dict(uniprot_counts)


def determine_monomere_status(
    member_ids: List[str], uniprot_counts: Dict[str, int]
) -> str:
    """
    Classify cluster oligomeric state based on domain counts.

    A protein is considered:
      - Monomeric: Single domain (domain_count = 1)
      - Multimeric: Multiple domains (domain_count > 1)

    Args:
        member_ids: List of sequence identifiers in cluster
        uniprot_counts: Domain counts per UniProt ID

    Returns:
        "monomere": All members are single-domain
        "multimere": All members are multi-domain
        "mixture": Mix of single and multi-domain

    Note:
        Centipede genome sequences are always treated as monomeric
    """
    domain_types = []

    for member_id in member_ids:
        # Extract just the ID part (before any spaces)
        seq_id = member_id.split()[0]
        uniprot_id = extract_uniprot_id(seq_id)

        # Centipede genome IDs are always treated as unique/monomeric
        if is_centipede_sequence(uniprot_id):
            domain_types.append("monomeric")
        else:
            domain_count = uniprot_counts.get(
                uniprot_id, 1
            )  # Default to 1 if not found
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


# ============================================================================
# CLUSTER FILE PARSING
# ============================================================================


def parse_mmseqs_tsv(tsv_path: Path) -> Dict[str, List[str]]:
    """
    Parse MMseqs2 cluster TSV file.

    Expected format:
      - Header line (skipped)
      - Tab-separated: representative_id    member_id

    Args:
        tsv_path: Path to MMseqs2 cluster TSV file

    Returns:
        Dict mapping representative_id -> list of member_ids

    Raises:
        FileNotFoundError: If TSV file doesn't exist

    Note:
        Duplicates are removed while preserving order
    """
    mapping = defaultdict(list)
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV not found: {tsv_path}")
    with tsv_path.open("r", encoding="utf-8") as fh:
        # Skip header line
        _ = fh.readline()
        for ln, raw in enumerate(fh, 2):  # Start from line 2
            line = raw.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            rep = cols[0]
            member = cols[1]
            mapping[rep].append(member)

    # dedupe preserving order
    cleaned = {}
    for rep, members in mapping.items():
        seen = set()
        out = []
        for m in members:
            if m not in seen:
                seen.add(m)
                out.append(m)
        cleaned[rep] = out
    logger.debug(f"Parsed {len(cleaned)} clusters from {tsv_path}")
    return cleaned


# ============================================================================
# GROUP CLASSIFICATION (PATTERN MATCHING)
# ============================================================================


def extract_header_fields(full_header: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract key fields from FASTA header for pattern matching.

    Args:
        full_header: Complete FASTA header line (without '>')

    Returns:
        Tuple of (protein_name, interpro_id)
        Either value may be None if not found
    """
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
    """
    Find all group patterns that match a FASTA header.

    Tests the header against all GROUP_PATTERNS. Multiple groups may match.
    Searches in: protein_name, full_header, sequence_id (with underscores and spaces).

    Args:
        full_header: Complete FASTA header line

    Returns:
        List of matching group names (may be empty)
    """
    pn, ip = extract_header_fields(full_header)

    groups = []

    # Build search targets: protein_name, header, underscore->space, id token, etc.
    id_token = full_header.split(None, 1)[0]
    targets = []
    if pn:
        targets.append(pn)
    targets.extend(
        [
            full_header,
            full_header.replace("_", " "),
            id_token,
            id_token.replace("_", " "),
        ]
    )

    # Test patterns against all targets, collect all matched group names
    for target in targets:
        for gname, compiled_list in _COMPILED_GROUPS:
            for cre in compiled_list:
                if cre.search(target):
                    if gname not in groups:
                        groups.append(gname)
                        logger.debug(
                            f"Header '{full_header[:50]}...' matched pattern for {gname} in target '{target}'"
                        )

    if not groups:
        logger.debug(f"Header '{full_header[:50]}...' matched no groups")

    return groups


def find_group_for_header(full_header: str) -> Optional[str]:
    """
    Determine the best single group for a header when multiple patterns match.

    Uses SPECIFICITY_RANK to resolve conflicts (lower rank = higher priority).

    Args:
        full_header: Complete FASTA header line

    Returns:
        Single group name, or None if no patterns match
    """
    matched = find_all_groups_for_header(full_header)
    if not matched:
        return None
    # If only one, return it
    if len(matched) == 1:
        return matched[0]
    # If multiple, choose the one with smallest SPECIFICITY_RANK
    ranked = sorted(matched, key=lambda g: SPECIFICITY_RANK.get(g, 9999))
    logger.debug(
        f"Header '{full_header[:50]}...' matched multiple groups {matched}, chose {ranked[0]} by specificity"
    )
    return ranked[0]


def decide_cluster_label(member_full_headers: List[str], rep_id: str = "") -> str:
    """
    Assign a group label to an entire cluster based on member classifications.

    Logic:
      1. If all members match same group (unanimous): Return that group
      2. If >10 members and ≥90% match same group: Return "mostly_{group}"
      3. If multiple different groups: Return "mixture"
      4. If one group + unknowns: Return "mixture"
      5. If all unknowns: Return "none"

    Args:
        member_full_headers: List of FASTA headers for all cluster members
        rep_id: Representative ID for logging (optional)

    Returns:
        Group label string
    """
    n = len(member_full_headers)
    groups = [find_group_for_header(h) for h in member_full_headers]
    counts = Counter(groups)
    non_none = [g for g in groups if g is not None]
    distinct_non_none = set(non_none)

    logger.debug(f"Cluster {rep_id}: {n} members")
    for i, (header, group) in enumerate(zip(member_full_headers, groups)):
        header_short = header.split()[0] if header else "None"
        logger.debug(f"  Member {i + 1}: {header_short} -> {group}")
    logger.debug(f"  Group counts: {dict(counts)}")

    # exact unanimous (no Nones)
    if len(distinct_non_none) == 1 and counts[None] == 0:
        result = next(iter(distinct_non_none))
        logger.debug(f"  Decision: unanimous {result}")
        return result

    # mostly_ rule
    if n > 10 and distinct_non_none:
        best_group, best_count = counts.most_common(1)[0]
        if best_group is not None and (best_count / n) >= 0.90:
            result = f"mostly_{best_group}"
            logger.debug(
                f"  Decision: {result} ({best_count}/{n} = {best_count / n:.2%})"
            )
            return result

    # mixture / none logic
    if len(distinct_non_none) > 1:
        result = "mixture"
        logger.debug(
            f"  Decision: {result} (multiple distinct groups: {distinct_non_none})"
        )
        return result
    if len(distinct_non_none) == 1 and counts[None] > 0:
        result = "mixture"
        logger.debug(f"  Decision: {result} (one group + unknowns)")
        return result
    if not distinct_non_none:
        result = "none"
        logger.debug(f"  Decision: {result} (no groups matched)")
        return result

    result = "mixture"
    logger.debug(f"  Decision: {result} (fallback)")
    return result


def build_cluster_member_headers(
    rep_to_members: Dict[str, List[str]],
    headers_by_id: Dict[str, str],
    rep_headers_by_id: Dict[str, str],
) -> Dict[str, List[str]]:
    """
    Build complete cluster membership with full headers.

    Combines representative and member FASTA headers into complete clusters.
    Ensures the representative is included in its own cluster.

    Args:
        rep_to_members: Mapping from representative_id -> member_ids
        headers_by_id: All sequence headers (merged FASTA)
        rep_headers_by_id: Representative sequence headers (rep FASTA)

    Returns:
        Dict mapping rep_full_header -> list of member_full_headers

    Note:
        Missing IDs are used as raw strings with a warning
    """
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
            if (
                full == mid
                and mid not in headers_by_id
                and mid not in rep_headers_by_id
            ):
                missing.add(mid)
            member_fulls.append(full)
        out[rep_full] = member_fulls
    if missing:
        logger.warning(
            f"{len(missing)} member IDs not matched to FASTA headers; they will be used as raw IDs."
        )
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Missing member IDs: " + ", ".join(sorted(list(missing))[:50]))
    return out


# ============================================================================
# OUTPUT GENERATION
# ============================================================================


def write_output_csv(path: Path, rows: List[Dict[str, str]], min_tax_levels: List[str]):
    """
    Write cluster annotations to CSV file.

    Creates the output directory if it doesn't exist. The statistics file
    will be written to the same directory (see write_statistics_log).

    CSV columns:
      - identifier: Cluster representative ID
      - group: Protein family classification
      - oligomeric_state: monomere/multimere/mixture
      - lca_{level}: Taxonomy at each requested level
      - length: Length bin

    Args:
        path: Output CSV file path (parent directory will be created if needed)
        rows: List of annotation dicts
        min_tax_levels: Taxonomic levels included (for column headers)
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # Build field names with all LCA columns
    fieldnames = ["identifier", "group", "oligomeric_state"]
    fieldnames.extend([f"lca_{level}" for level in min_tax_levels])
    fieldnames.append("length")

    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    logger.info(f"Wrote {len(rows)} rows to {path}")


def write_statistics_log(
    output_path: Path,
    stats: Counter,
    monomere_stats: Counter,
    taxonomy_stats: Dict[str, Counter],
    length_stats: Counter,
    min_tax_levels: List[str],
):
    """
    Write detailed statistics summary to a text log file.

    Generates a formatted report with:
      - Total cluster count
      - Group distribution (sorted by frequency)
      - Oligomeric state distribution
      - Taxonomy distribution (separate section per level)
      - Length bin distribution (sorted numerically)

    Args:
        output_path: Output CSV path (stats file goes in same directory)
        stats: Group classification counts
        monomere_stats: Oligomeric state counts
        taxonomy_stats: Taxonomy counts per level
        length_stats: Length bin counts
        min_tax_levels: Taxonomic levels to report

    Output:
        Creates statistics.txt in the same directory as the CSV file
    """
    log_path = output_path.parent / "statistics.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as fh:
        fh.write("=" * 80 + "\n")
        fh.write("CLUSTER ANNOTATION STATISTICS\n")
        fh.write("=" * 80 + "\n\n")

        fh.write(f"Total clusters: {sum(stats.values())}\n\n")

        fh.write("-" * 80 + "\n")
        fh.write("GROUP DISTRIBUTION\n")
        fh.write("-" * 80 + "\n")
        for k, v in stats.most_common():
            fh.write(f"  {k:30s} {v:8d}\n")

        fh.write("\n" + "-" * 80 + "\n")
        fh.write("OLIGOMERIC STATE\n")
        fh.write("-" * 80 + "\n")
        for k, v in monomere_stats.most_common():
            fh.write(f"  {k:30s} {v:8d}\n")

        # Write taxonomy distribution for each requested level
        for level in min_tax_levels:
            fh.write("\n" + "-" * 80 + "\n")
            fh.write(f"TAXONOMY DISTRIBUTION (LCA at {level} level)\n")
            fh.write("-" * 80 + "\n")
            level_stats = taxonomy_stats.get(level, Counter())
            for k, v in level_stats.most_common():
                fh.write(f"  {k:60s} {v:8d}\n")

        fh.write("\n" + "-" * 80 + "\n")
        fh.write("LENGTH BINS\n")
        fh.write("-" * 80 + "\n")

        # Sort length bins by numeric value (extract start of range)
        def extract_bin_start(bin_str):
            if bin_str == "unknown":
                return float("inf")  # Put unknown at the end
            try:
                return int(bin_str.split("-")[0])
            except (ValueError, IndexError):
                return float("inf")

        sorted_bins = sorted(
            length_stats.items(), key=lambda x: extract_bin_start(x[0])
        )
        for k, v in sorted_bins:
            fh.write(f"  {k:30s} {v:8d}\n")

        fh.write("\n" + "=" * 80 + "\n")

    logger.info(f"Statistics saved to {log_path}")


# ============================================================================
# MAIN SCRIPT
# ============================================================================


def main():
    """
    Main script entry point.

    Workflow:
      1. Parse command-line arguments
      2. Initialize taxonomy database
      3. Parse FASTA files and cluster mappings
      4. Batch-retrieve taxonomy for all unique organisms
      5. Annotate each cluster with:
         - Group classification (pattern matching)
         - Oligomeric state (domain counting)
         - Taxonomy (LCA at specified levels)
         - Length bin
      6. Write CSV output and statistics log
    """
    p = argparse.ArgumentParser(
        description="Annotate protein clusters with group, oligomeric state, taxonomy, and length bins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single taxonomy level (default: family)
  %(prog)s

  # Multiple taxonomy levels
  %(prog)s --min-tax-level family genus order

  # Custom output directory with verbose output
  %(prog)s -O data/processed/my_analysis/annotations.csv -v

  # All outputs will be in: data/processed/cluster_annotations/
  #   - annotations.csv (cluster annotations)
  #   - statistics.txt (summary statistics)
        """,
    )
    p.add_argument(
        "--original-fasta",
        "-o",
        type=Path,
        default=Path("data/interm/mmseqs/merged.fasta"),
        help="Original FASTA file (default: data/interm/mmseqs/merged.fasta)",
    )
    p.add_argument(
        "--rep-fasta",
        "-r",
        type=Path,
        default=Path("data/interm/mmseqs/representatives.fasta"),
        help="Representative sequences FASTA (default: data/interm/mmseqs/representatives.fasta)",
    )
    p.add_argument(
        "--cluster-tsv",
        "-c",
        type=Path,
        default=Path("data/interm/mmseqs/clusters_with_sizes.tsv"),
        help="MMseqs2 cluster TSV file (default: data/interm/mmseqs/clusters_with_sizes.tsv)",
    )
    p.add_argument(
        "--output-csv",
        "-O",
        type=Path,
        default=Path("data/processed/cluster_annotations/annotations.csv"),
        help="Output CSV file (default: data/processed/cluster_annotations/annotations.csv)",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    p.add_argument(
        "--num-bins",
        type=int,
        default=15,
        help="Number of length bins to create (default: 15)",
    )
    p.add_argument(
        "--min-tax-level",
        nargs="+",
        default=["family"],
        choices=TAXONOMIC_LEVELS,
        help="Minimum taxonomic level(s) for last common ancestor (default: family). Can specify multiple levels.",
    )
    args = p.parse_args()

    # Ensure min_tax_level is a list
    if isinstance(args.min_tax_level, str):
        args.min_tax_level = [args.min_tax_level]

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Initialize taxopy database
    logger.info("Initializing taxonomy database...")
    taxdb = initialize_taxdb()

    # Parse FASTA files with sequence lengths
    headers_by_id, sequence_lengths = parse_fasta_with_lengths(args.original_fasta)
    rep_headers_by_id, rep_sequence_lengths = parse_fasta_with_lengths(args.rep_fasta)

    # Combine sequence lengths from both files
    all_sequence_lengths = {**sequence_lengths, **rep_sequence_lengths}

    rep_to_members = parse_mmseqs_tsv(args.cluster_tsv)

    # Count UniProt domain occurrences
    uniprot_counts = count_uniprot_domains(headers_by_id)

    clusters = build_cluster_member_headers(
        rep_to_members, headers_by_id, rep_headers_by_id
    )

    # Extract organism names from all headers
    organism_by_seqid_orig = extract_all_organism_names(headers_by_id)
    organism_by_seqid_rep = extract_all_organism_names(rep_headers_by_id)
    organism_by_seqid = {**organism_by_seqid_orig, **organism_by_seqid_rep}

    # Get unique organism names
    unique_organisms = list(set(organism_by_seqid.values()))
    logger.info(f"Retrieving taxonomy for {len(unique_organisms)} unique organisms...")

    # Batch retrieve taxonomy Taxon objects for all unique organisms
    organism_to_taxon = batch_get_taxonomy_for_organisms(unique_organisms, taxdb)

    # First pass to collect all cluster average lengths
    cluster_avg_lengths = []
    cluster_length_data = {}

    for rep_full, member_fulls in clusters.items():
        rep_id = rep_full.split()[0]
        avg_length = calculate_cluster_average_length(
            member_fulls, all_sequence_lengths
        )
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

    # Create length bins
    bin_ranges, length_to_bin = create_length_bins(cluster_avg_lengths, args.num_bins)

    rows = []
    stats = Counter()
    monomere_stats = Counter()
    taxonomy_stats = {level: Counter() for level in args.min_tax_level}
    length_stats = Counter()

    for rep_full, member_fulls in clusters.items():
        rep_id = rep_full.split()[0]

        label = decide_cluster_label(member_fulls, rep_id)

        # Determine monomere status
        oligomeric_status = determine_monomere_status(member_fulls, uniprot_counts)

        # Determine last common ancestor for all requested levels
        lca_results = find_last_common_ancestor(
            member_fulls,
            organism_by_seqid,
            organism_to_taxon,
            taxdb,
            args.min_tax_level,
        )

        # Format LCA with rank in parentheses for each level
        row = {
            "identifier": rep_id,
            "group": label,
            "oligomeric_state": oligomeric_status,
        }

        for level in args.min_tax_level:
            lca_name, lca_rank = lca_results.get(level, ("unknown", "unknown"))
            if lca_name != "unknown" and lca_rank != "unknown":
                lca_display = f"{lca_name} ({lca_rank})"
            else:
                lca_display = lca_name
            row[f"lca_{level}"] = lca_display
            taxonomy_stats[level][lca_display] += 1

        # Get length bin
        avg_length = cluster_length_data.get(rep_id)
        length_bin = (
            length_to_bin.get(avg_length, "unknown")
            if avg_length is not None
            else "unknown"
        )
        row["length"] = length_bin

        rows.append(row)
        stats[label] += 1
        monomere_stats[oligomeric_status] += 1
        length_stats[length_bin] += 1

    # Process singletons
    for rep_id in missing_reps:
        rep_full = rep_headers_by_id.get(rep_id, rep_id)
        label = decide_cluster_label([rep_full], rep_id)

        # Determine monomere status for singletons
        oligomeric_status = determine_monomere_status([rep_full], uniprot_counts)

        # Determine last common ancestor for singletons
        lca_results = find_last_common_ancestor(
            [rep_full], organism_by_seqid, organism_to_taxon, taxdb, args.min_tax_level
        )

        # Format LCA with rank in parentheses for each level
        row = {
            "identifier": rep_full.split()[0],
            "group": label,
            "oligomeric_state": oligomeric_status,
        }

        for level in args.min_tax_level:
            lca_name, lca_rank = lca_results.get(level, ("unknown", "unknown"))
            if lca_name != "unknown" and lca_rank != "unknown":
                lca_display = f"{lca_name} ({lca_rank})"
            else:
                lca_display = lca_name
            row[f"lca_{level}"] = lca_display
            taxonomy_stats[level][lca_display] += 1

        # Get length bin for singletons
        avg_length = cluster_length_data.get(rep_id)
        length_bin = (
            length_to_bin.get(avg_length, "unknown")
            if avg_length is not None
            else "unknown"
        )
        row["length"] = length_bin

        rows.append(row)
        stats[label] += 1
        monomere_stats[oligomeric_status] += 1
        length_stats[length_bin] += 1

    write_output_csv(args.output_csv, rows, args.min_tax_level)

    # Write detailed statistics to log file
    write_statistics_log(
        args.output_csv,
        stats,
        monomere_stats,
        taxonomy_stats,
        length_stats,
        args.min_tax_level,
    )

    # Print brief summary to console
    logger.info("")
    logger.info(f"Processed {len(rows)} clusters")
    logger.info(f"  - Groups: {len(stats)} types")
    for level in args.min_tax_level:
        logger.info(f"  - Taxonomies ({level}): {len(taxonomy_stats[level])} unique")
    logger.info(f"  - Length bins: {len(length_stats)}")
    logger.info("")


if __name__ == "__main__":
    main()
