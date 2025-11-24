#!/usr/bin/env python3
"""
Collect and Filter InterPro Sequences by Taxonomy

Combines sequence collection from InterPro folders with taxonomy-based filtering.
Collects sequences, extracts taxonomy IDs from JSON metadata, and filters to keep
only specified taxonomic groups (default: Metazoa) in a single efficient pass.

Performance: Processes 50k+ sequences in ~10 seconds (collection + filtering).

Usage:
  uv run python src/data_prep/interpro_metazoa_sequences_collector.py \\
    data/raw/interpro data/interm/interpro --keep Metazoa
"""

import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import taxopy

from taxonomy_utils import initialize_taxdb

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
# INTERPRO DATA PARSING
# ============================================================================


def parse_fasta(fasta_file: Path) -> Dict[str, str]:
    """Parse FASTA file, return dict: seq_id -> sequence."""
    sequences = {}
    seq_id, seq_lines = None, []

    with fasta_file.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq_lines)

    return sequences


def parse_json_metadata(json_file: Path) -> Dict[str, Dict[str, str]]:
    """Extract taxonomy IDs and protein names from InterPro JSON, return dict: protein_id -> {taxid, protein_name}."""
    with json_file.open() as f:
        data = json.load(f)

    metadata = {}
    proteins_data = data.get("proteins_data", {}).get("proteins", {})

    for protein_id, protein_info in proteins_data.items():
        protein_metadata = protein_info.get("metadata", {})
        source_organism = protein_metadata.get("source_organism", {})
        protein_name = protein_metadata.get("name", "")

        if taxid := source_organism.get("taxId"):
            metadata[protein_id] = {"taxid": str(taxid), "protein_name": protein_name}

    return metadata


# ============================================================================
# TAXONOMY FILTERING
# ============================================================================


def create_taxon_objects(
    taxids: List[str], taxdb: taxopy.TaxDb
) -> Dict[str, taxopy.Taxon]:
    """Create Taxon objects for unique taxonomy IDs."""
    taxid_to_taxon = {}
    for taxid in taxids:
        try:
            taxid_to_taxon[taxid] = taxopy.Taxon(int(taxid), taxdb)
        except Exception as e:
            logger.debug(f"Failed to load taxid {taxid}: {e}")
    return taxid_to_taxon


def get_matching_taxids(
    taxid_to_taxon: Dict[str, taxopy.Taxon], target_groups: List[str]
) -> set:
    """Pre-compute which taxids belong to target taxonomic groups."""
    matching = set()
    for taxid, taxon in taxid_to_taxon.items():
        try:
            lineage_names = set(taxon.name_lineage)
            if any(group in lineage_names for group in target_groups):
                matching.add(taxid)
        except Exception as e:
            logger.debug(f"Lineage error for taxid {taxid}: {e}")
    return matching


# ============================================================================
# SEQUENCE COLLECTION AND FILTERING
# ============================================================================


def collect_and_filter_sequences(
    input_dir: Path,
    output_dir: Path,
    target_groups: List[str],
    taxdb: taxopy.TaxDb,
) -> Tuple[int, int, int]:
    """
    Collect sequences from InterPro folders and filter by taxonomy in one pass.

    Returns: (total_collected, kept, excluded)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Track sequences and metadata
    all_sequences = {}  # seq_id -> sequence
    seq_to_taxid = {}  # seq_id -> taxid
    seq_to_protein_name = {}  # seq_id -> protein_name
    seq_to_interpro = defaultdict(list)  # seq_id -> [interpro_ids]

    logger.info(f"Collecting sequences from {input_dir}...")

    # Collect sequences from all InterPro subfolders
    for subfolder in sorted(input_dir.iterdir()):
        if not subfolder.is_dir():
            continue

        interpro_id = subfolder.name

        # Find FASTA and JSON files
        fasta_file = next(
            (f for f in subfolder.iterdir() if f.suffix.lower() in [".fasta", ".fa"]),
            None,
        )
        json_file = next(
            (f for f in subfolder.iterdir() if f.suffix.lower() == ".json"), None
        )

        if not fasta_file:
            logger.debug(f"{interpro_id}: No FASTA file, skipping")
            continue

        # Load sequences
        sequences = parse_fasta(fasta_file)

        # Load taxonomy IDs and protein names from JSON
        json_metadata = {}
        if json_file:
            try:
                json_metadata = parse_json_metadata(json_file)
            except Exception as e:
                logger.warning(f"{interpro_id}: Failed to parse JSON: {e}")

        # Merge sequences and metadata
        new_count = 0
        for seq_id, sequence in sequences.items():
            if seq_id not in all_sequences:
                all_sequences[seq_id] = sequence
                new_count += 1

            # Store taxid and protein name from JSON (if available)
            if seq_id in json_metadata:
                seq_to_taxid[seq_id] = json_metadata[seq_id]["taxid"]
                if json_metadata[seq_id]["protein_name"]:
                    seq_to_protein_name[seq_id] = json_metadata[seq_id]["protein_name"]

            # Track InterPro membership
            seq_to_interpro[seq_id].append(interpro_id)

        logger.info(
            f"{interpro_id}: {len(sequences)} seqs, {new_count} new, {len(json_metadata)} with taxids"
        )

    total_collected = len(all_sequences)
    logger.info(f"\nCollected {total_collected} unique sequences")
    logger.info(f"  With taxonomy IDs: {len(seq_to_taxid)}")
    logger.info(f"  With protein names: {len(seq_to_protein_name)}")
    logger.info(f"  Missing taxonomy: {total_collected - len(seq_to_taxid)}")

    # Create Taxon objects for unique taxids
    unique_taxids = list(set(seq_to_taxid.values()))
    logger.info(f"\nCreating Taxon objects for {len(unique_taxids)} unique taxids...")
    taxid_to_taxon = create_taxon_objects(unique_taxids, taxdb)
    logger.info(f"  Loaded {len(taxid_to_taxon)}/{len(unique_taxids)} taxids")

    # Pre-compute matching taxids
    logger.info(f"\nFiltering for taxonomic groups: {', '.join(target_groups)}")
    matching_taxids = get_matching_taxids(taxid_to_taxon, target_groups)
    logger.info(f"  Matched: {len(matching_taxids)}/{len(taxid_to_taxon)} taxids")

    # Filter sequences
    logger.info(f"\nFiltering {total_collected} sequences...")
    filtered_sequences = {}
    stats = {"kept": 0, "excluded": 0, "no_taxid": 0}

    for seq_id, sequence in all_sequences.items():
        taxid = seq_to_taxid.get(seq_id)

        if not taxid:
            stats["no_taxid"] += 1
            continue

        if taxid in matching_taxids:
            protein_name = seq_to_protein_name.get(seq_id, "")
            filtered_sequences[seq_id] = (
                sequence,
                taxid,
                seq_to_interpro[seq_id],
                protein_name,
            )
            stats["kept"] += 1
        else:
            stats["excluded"] += 1

    logger.info(
        f"  Kept: {stats['kept']}, Excluded: {stats['excluded']}, No taxid: {stats['no_taxid']}"
    )

    # Write output
    output_fasta = output_dir / "interpro_metazoa_sequences.fasta"
    write_fasta(output_fasta, filtered_sequences)

    # Write report to reports subdirectory
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    write_report(
        reports_dir / "metazoa_collection_report.txt",
        total_collected,
        stats,
        target_groups,
    )

    return total_collected, stats["kept"], stats["excluded"] + stats["no_taxid"]


# ============================================================================
# OUTPUT
# ============================================================================


def write_fasta(
    output_path: Path, sequences: Dict[str, Tuple[str, str, List[str], str]]
):
    """Write sequences to FASTA with taxa_id, protein_name, and interpro fields."""
    with output_path.open("w") as f:
        for seq_id in sorted(sequences.keys()):
            sequence, taxid, interpro_ids, protein_name = sequences[seq_id]
            interpro_str = ",".join(sorted(interpro_ids))
            header = f"{seq_id} taxa_id={taxid}"
            if protein_name:
                # Replace spaces with underscores to keep header parsing simple
                protein_name_sanitized = protein_name.replace(" ", "_")
                header += f" protein_name={protein_name_sanitized}"
            header += f" interpro={interpro_str}"

            f.write(f">{header}\n")
            for i in range(0, len(sequence), 80):
                f.write(f"{sequence[i : i + 80]}\n")

    logger.info(f"\nWrote {len(sequences)} sequences to {output_path.name}")


def write_report(
    report_path: Path, total: int, stats: Dict[str, int], target_groups: List[str]
):
    """Write collection and filtering report."""
    with report_path.open("w") as f:
        f.write("=" * 80 + "\n")
        f.write("INTERPRO METAZOA SEQUENCE COLLECTION REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Target taxonomic groups: {', '.join(target_groups)}\n\n")

        f.write("-" * 80 + "\n")
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total sequences collected:  {total:8d}\n")
        f.write(
            f"Kept (in target groups):    {stats['kept']:8d} ({stats['kept'] / total * 100:.1f}%)\n"
        )
        f.write(
            f"Excluded (other taxonomy):  {stats['excluded']:8d} ({stats['excluded'] / total * 100:.1f}%)\n"
        )
        f.write(
            f"Missing taxonomy data:      {stats['no_taxid']:8d} ({stats['no_taxid'] / total * 100:.1f}%)\n"
        )
        f.write("\n" + "=" * 80 + "\n")

    logger.info(f"Report saved: {report_path.name}")


# ============================================================================
# MAIN
# ============================================================================


def main():
    """Collect InterPro sequences and filter by taxonomy in one efficient pass."""
    parser = argparse.ArgumentParser(
        description="Collect and filter InterPro sequences by taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Collect and keep only Metazoa (default)
  uv run python src/data_prep/interpro_metazoa_sequences_collector.py \\
    data/raw/interpro data/interm/interpro

  # Keep multiple taxonomic groups
  uv run python src/data_prep/interpro_metazoa_sequences_collector.py \\
    data/raw/interpro data/interm/interpro --keep Metazoa Fungi

  # With verbose logging
  uv run python src/data_prep/interpro_metazoa_sequences_collector.py \\
    data/raw/interpro data/interm/interpro -v
        """,
    )
    parser.add_argument(
        "input_dir", type=Path, help="Directory containing InterPro subfolders"
    )
    parser.add_argument(
        "output_dir", type=Path, help="Output directory for filtered sequences"
    )
    parser.add_argument(
        "--keep",
        nargs="+",
        default=["Metazoa"],
        help="Taxonomic groups to keep (default: Metazoa)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Initialize taxonomy database
    logger.info("Loading taxonomy database...")
    taxdb = initialize_taxdb()

    # Collect and filter sequences
    total, kept, excluded = collect_and_filter_sequences(
        args.input_dir, args.output_dir, args.keep, taxdb
    )

    logger.info(
        f"\nDone! Collected {total} sequences, kept {kept} ({kept / total * 100:.1f}%)"
    )


if __name__ == "__main__":
    main()
