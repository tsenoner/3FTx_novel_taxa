#!/usr/bin/env python3
"""
Taxonomy Utilities

Shared utility functions for taxonomy operations using taxopy.

Functions:
  - initialize_taxdb(): Initialize or load cached NCBI taxonomy database
  - extract_taxa_id(): Extract numeric taxonomy ID from FASTA headers
  - batch_get_taxonomy_for_taxids(): Batch retrieve taxonomy for multiple taxids
  - check_taxonomy_membership(): Check if taxon belongs to a taxonomic group
  - get_organism_kingdom(): Get kingdom for a taxon
"""

import logging
import os
import re
import shutil
import tempfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional

import taxopy
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Taxonomic hierarchy levels (from general to specific)
TAXONOMIC_LEVELS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

# FASTA header field extraction pattern
# Matches "taxa_id=..." (numeric ID) followed by either another field or end of line
_re_taxa_id = re.compile(
    r"taxa_id=(\d+)(?=\s+(?:protein_name=|domain_pos=|interpro_id|interpro_ids=|interpro=|signature=|length=)|$)",
    re.IGNORECASE,
)


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


def extract_taxa_id(header: str) -> Optional[int]:
    """
    Extract numeric taxonomy ID from FASTA header.

    Looks for 'taxa_id=NNNN' field format.

    Args:
        header: Full FASTA header line (without '>')

    Returns:
        Taxonomy ID as integer, or None if not found

    Example:
        'taxa_id=9606 protein_name=...' -> 9606
    """
    match = _re_taxa_id.search(header)
    if match:
        return int(match.group(1))
    return None


# ============================================================================
# TAXONOMY RETRIEVAL
# ============================================================================


def batch_get_taxonomy_for_taxids(
    taxids: List[int], taxdb: taxopy.TaxDb, show_progress: bool = True
) -> Dict[int, Optional[taxopy.Taxon]]:
    """Batch retrieve taxonomy Taxon objects for multiple taxonomy IDs.

    Args:
        taxids: List of numeric taxonomy IDs to look up
        taxdb: Initialized taxopy database
        show_progress: Whether to show progress bar (default: True)

    Returns:
        Dict mapping taxid -> Taxon object (or None if not found)
    """
    taxid_to_taxon = {}

    # Process each taxid with optional progress bar
    iterator = taxids
    if show_progress:
        iterator = tqdm(
            iterator,
            total=len(taxids),
            desc="Fetching taxonomy",
            unit="taxid",
        )

    for taxid in iterator:
        try:
            taxon = taxopy.Taxon(taxid, taxdb)
            taxid_to_taxon[taxid] = taxon
            logger.debug(f"Retrieved Taxon for taxid {taxid}: {taxon.name}")
        except Exception as e:
            logger.warning(f"Failed to get taxonomy for taxid {taxid}: {e}")
            taxid_to_taxon[taxid] = None

    successful = sum(1 for t in taxid_to_taxon.values() if t is not None)
    logger.info(
        f"Successfully retrieved taxonomy for {successful}/{len(taxids)} taxids"
    )
    return taxid_to_taxon


# ============================================================================
# TAXONOMY MEMBERSHIP CHECKING
# ============================================================================


def check_taxonomy_membership(
    taxon: taxopy.Taxon, target_taxon_name: str, taxdb: taxopy.TaxDb
) -> bool:
    """
    Check if a taxon is a member of (descends from) a target taxonomic group.

    Args:
        taxon: The taxon to check
        target_taxon_name: Name of the target taxonomic group (e.g., "Metazoa")
        taxdb: Initialized taxopy database

    Returns:
        True if taxon is a member of the target group, False otherwise

    Example:
        taxon = taxopy.Taxon(9606, taxdb)  # Homo sapiens
        check_taxonomy_membership(taxon, "Metazoa", taxdb)  # Returns True
        check_taxonomy_membership(taxon, "Viridiplantae", taxdb)  # Returns False
    """
    if taxon is None:
        return False

    try:
        # Get taxon ID for the target group
        target_taxid_list = taxopy.taxid_from_name(target_taxon_name, taxdb)
        if not target_taxid_list or not target_taxid_list[0]:
            logger.warning(
                f"Could not find taxon ID for target group: {target_taxon_name}"
            )
            return False

        target_taxid = target_taxid_list[0]
        target_taxon = taxopy.Taxon(target_taxid, taxdb)

        # Get the full lineage of our taxon
        lineage = taxon.taxid_lineage

        # Check if target taxon is in the lineage
        is_member = target_taxon.taxid in lineage

        logger.debug(
            f"{taxon.name} is {'a' if is_member else 'not a'} member of {target_taxon_name}"
        )
        return is_member

    except Exception as e:
        logger.warning(
            f"Failed to check taxonomy membership for {taxon.name} vs {target_taxon_name}: {e}"
        )
        return False


def get_organism_kingdom(taxon: taxopy.Taxon) -> Optional[str]:
    """
    Get the kingdom (or highest available rank) for a taxon.

    Args:
        taxon: The taxon to query

    Returns:
        Kingdom name, or None if not available
    """
    if taxon is None:
        return None

    try:
        ranks = taxon.rank_name_dictionary
        # Try kingdom first, then fall back to domain/superkingdom
        kingdom = (
            ranks.get("kingdom", "")
            or ranks.get("domain", "")
            or ranks.get("superkingdom", "")
        )
        return kingdom if kingdom else None
    except Exception as e:
        logger.warning(f"Failed to get kingdom for taxon: {e}")
        return None
