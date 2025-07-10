#!/usr/bin/env python3
"""
InterPro Data Fetcher
Fetches protein sequences and domain boundaries from the InterPro REST API.

This script fetches protein metadata and domain locations from InterPro for a given
InterPro ID. It handles API pagination and deduplicates protein entries and their
domain locations. It can optionally fetch full protein sequences from UniProt.

Usage:
    python interpro_fetcher.py IPR000001
    python interpro_fetcher.py IPR000001 --limit 100
    python interpro_fetcher.py IPR016054 --no-sequences
"""

import argparse
import csv
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import copy

import requests
from tqdm import tqdm

# --- Constants ---
INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api"
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_PAGE_SIZE = 200


# --- Data Fetching ---
def fetch_interpro_data(
    ipr_id: str,
    max_unique_proteins: Optional[int] = None,
    existing_proteins: Optional[Dict[str, Any]] = None,
    entry_info: Optional[Dict[str, Any]] = None,
    log_filename: Optional[Path] = None,
) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """
    Fetches and processes protein data from the InterPro API for a given ID.

    This function handles pagination and appends new data to a log file incrementally.

    Args:
        ipr_id: The InterPro ID to fetch data for (e.g., "IPR000001").
        max_unique_proteins: The maximum number of unique proteins to collect.
        existing_proteins: A dict of proteins that have already been fetched.
        entry_info: Pre-fetched entry information to avoid a redundant API call.
        log_filename: The path to the JSONL file for logging new proteins.

    Returns:
        A tuple containing the entry information and the processed proteins data,
        or (None, None) if an error occurs.
    """
    print(f"Fetching data for {ipr_id}...")

    # 1. Get overall entry information
    if not entry_info:
        entry_url = f"{INTERPRO_API_URL}/entry/interpro/{ipr_id}"
        try:
            entry_response = requests.get(entry_url)
            entry_response.raise_for_status()
            entry_info = entry_response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching entry info for {ipr_id}: {e}")
            return None, None

    total_proteins = (
        entry_info.get("metadata", {}).get("counters", {}).get("proteins", 0)
    )
    print(
        f"Entry: {entry_info.get('metadata', {}).get('name', {}).get('name', 'Unknown')}"
    )
    print(f"Total proteins available in InterPro: {total_proteins:,}")

    if max_unique_proteins:
        print(f"Fetching until {max_unique_proteins:,} unique proteins are found...")
    else:
        print("Fetching ALL available proteins...")

    # 2. Paginate through protein data
    all_proteins: Dict[str, Dict[str, Any]] = existing_proteins or {}
    if existing_proteins:
        print(f"Resuming with {len(existing_proteins):,} already collected proteins.")

    total_api_entries_processed = 0
    next_page_url: Optional[str] = (
        f"{INTERPRO_API_URL}/protein/uniprot/entry/interpro/{ipr_id}/?page_size={DEFAULT_PAGE_SIZE}"
    )

    # Set up the progress bar
    progress_bar_total = max_unique_proteins or total_proteins
    with tqdm(
        total=progress_bar_total,
        desc="Fetching Proteins",
        unit="protein",
        initial=len(all_proteins),
    ) as pbar:
        while next_page_url:
            # Stop condition 1: We have collected enough unique proteins
            if max_unique_proteins and len(all_proteins) >= max_unique_proteins:
                break

            try:
                pbar.write(f"Fetching page: {next_page_url}")  # DEBUG
                response = requests.get(next_page_url, timeout=60)
                response.raise_for_status()
                data = response.json()
            except requests.exceptions.RequestException as e:
                pbar.write(f"Error fetching page: {e}")
                break

            results = data.get("results", [])
            if not results:
                pbar.write("No more results found from API.")
                break

            # Process the results from the current page
            new_proteins_on_page = []
            for protein_entry in results:
                if max_unique_proteins and len(all_proteins) >= max_unique_proteins:
                    break
                total_api_entries_processed += 1
                protein_acc = protein_entry["metadata"]["accession"]

                # Add the protein if it's new
                if protein_acc not in all_proteins:
                    # Collect new protein entry to be logged
                    new_proteins_on_page.append(protein_entry)

                    all_proteins[protein_acc] = {
                        "metadata": protein_entry["metadata"],
                        "locations": set(),  # Use a set for efficient deduplication
                    }
                    pbar.update(1)  # Update progress bar only for new unique proteins

                # Add locations for any protein we are tracking
                if protein_acc in all_proteins:
                    for interpro_entry in protein_entry.get("entries", []):
                        if interpro_entry["accession"] == ipr_id:
                            for loc in interpro_entry.get(
                                "entry_protein_locations", []
                            ):
                                for frag in loc.get("fragments", []):
                                    # Use a tuple for the location so it's hashable for the set
                                    location_tuple = (
                                        frag.get("start"),
                                        frag.get("end"),
                                        loc.get("score"),
                                    )
                                    all_proteins[protein_acc]["locations"].add(
                                        location_tuple
                                    )

            # --- Incremental Save via Logging ---
            if new_proteins_on_page and log_filename:
                with log_filename.open("a") as logfile:
                    for item in new_proteins_on_page:
                        logfile.write(json.dumps(item) + "\n")

            next_page_url = data.get("next")
            time.sleep(0.1)  # Be nice to the server

    # Final trim if the last page gave us more proteins than the limit
    if max_unique_proteins and len(all_proteins) > max_unique_proteins:
        print(
            f"Trimming collected proteins from {len(all_proteins):,} to the {max_unique_proteins:,} limit."
        )
        all_proteins = dict(list(all_proteins.items())[:max_unique_proteins])

    proteins_data = {"proteins": all_proteins}
    print(
        f"\nSuccessfully processed {total_api_entries_processed:,} API entries to store {len(all_proteins):,} unique proteins."
    )
    return entry_info, proteins_data


def fetch_protein_sequences(
    protein_ids: List[str], batch_size: int = 100
) -> Dict[str, str]:
    """
    Fetches protein sequences in FASTA format from the UniProt API.
    """
    print(f"\nFetching sequences for {len(protein_ids):,} proteins...")
    sequences: Dict[str, str] = {}
    for i in tqdm(
        range(0, len(protein_ids), batch_size), desc="Fetching Sequences", unit="batch"
    ):
        batch = protein_ids[i : i + batch_size]
        query = " OR ".join(f"accession:{pid}" for pid in batch)
        # Add size parameter to retrieve all results in one go for the batch
        url = f"{UNIPROT_API_URL}?query=({query})&format=fasta&size={len(batch)}"

        try:
            response = requests.get(url)
            response.raise_for_status()

            # Simple FASTA parsing
            fasta_content = response.text
            for entry in fasta_content.strip().split(">"):
                if not entry:
                    continue
                parts = entry.split("\n", 1)
                header = parts[0]
                seq = parts[1].replace("\n", "")
                # Assumes UniProt format like >db|accession|...
                accession = header.split("|")[1] if "|" in header else header.split()[0]
                sequences[accession] = seq
        except requests.exceptions.RequestException as e:
            tqdm.write(f"Error fetching sequence batch: {e}")
        time.sleep(0.1)
    return sequences


# --- File Saving ---
def save_to_json(data: Dict[str, Any], filename: Path):
    """Saves dictionary data to a JSON file. Creates a serializable copy to avoid side effects."""
    serializable_data = copy.deepcopy(data)

    if "fetch_args" in serializable_data:
        # Convert Path objects to strings for serialization
        if "output_dir" in serializable_data["fetch_args"] and isinstance(
            serializable_data["fetch_args"]["output_dir"], Path
        ):
            serializable_data["fetch_args"]["output_dir"] = str(
                serializable_data["fetch_args"]["output_dir"]
            )

    # Convert sets to lists for JSON serialization
    if (
        "proteins_data" in serializable_data
        and "proteins" in serializable_data["proteins_data"]
    ):
        for protein_acc in serializable_data["proteins_data"]["proteins"]:
            if (
                "locations"
                in serializable_data["proteins_data"]["proteins"][protein_acc]
            ):
                locations_set = serializable_data["proteins_data"]["proteins"][
                    protein_acc
                ]["locations"]
                if isinstance(locations_set, set):
                    serializable_data["proteins_data"]["proteins"][protein_acc][
                        "locations"
                    ] = list(locations_set)

    with filename.open("w") as f:
        json.dump(serializable_data, f, indent=2)


def save_to_csv(proteins_data: Dict[str, Any], filename: Path, ipr_id: str):
    """Saves processed protein and domain data to a CSV file."""
    print(f"Saving processed protein data to {filename}...")
    fieldnames = [
        "protein_id",
        "protein_name",
        "tax_id",
        "scientific_name",
        "sequence_length",
        "domain_start",
        "domain_end",
        "signature_name",
        "interpro_domain_name",
    ]
    with filename.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for protein_acc, protein_info in proteins_data.get("proteins", {}).items():
            meta = protein_info.get("metadata", {})
            source_organism = meta.get("source_organism", {})
            base_row = {
                "protein_id": protein_acc,
                "protein_name": meta.get("name", "Unknown"),
                "tax_id": source_organism.get("taxId"),
                "scientific_name": source_organism.get("scientificName"),
                "sequence_length": meta.get("length", 0),
                "signature_name": f"InterPro:{ipr_id}",
                "interpro_domain_name": meta.get(
                    "name"
                ),  # Re-using protein name for domain name context
            }
            locations = protein_info.get("locations")
            if locations:
                for start, end, score in locations:
                    row = base_row.copy()
                    row.update(
                        {
                            "domain_start": start,
                            "domain_end": end,
                        }
                    )
                    writer.writerow(row)
            else:
                # Still write a row for proteins with no domain locations found
                row = base_row.copy()
                row.update({"domain_start": "N/A", "domain_end": "N/A"})
                writer.writerow(row)


def save_sequences_fasta(sequences: Dict[str, str], filename: Path):
    """Saves protein sequences to a FASTA file."""
    print(f"Saving {len(sequences):,} sequences to {filename}...")
    with filename.open("w") as f:
        for protein_id, sequence in sequences.items():
            f.write(f">{protein_id}\n")
            for i in range(0, len(sequence), 80):
                f.write(f"{sequence[i : i + 80]}\n")


# --- Main Execution ---
def create_arg_parser() -> argparse.ArgumentParser:
    """Creates the argument parser for the script."""
    parser = argparse.ArgumentParser(
        description="Fetch protein data and sequences from the InterPro API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s IPR000001                          # Fetch all proteins for Kringle domain.
  %(prog)s IPR000001 --limit 500              # Limit to the first 500 unique proteins.
  %(prog)s IPR000001 --no-sequences           # Skip sequence fetching entirely.
  %(prog)s IPR002048 --output-dir custom_dir  # Save files to a custom directory.
  %(prog)s IPR000001 --limit 1000 --resume    # Resume a previous, interrupted run.
        """,
    )
    parser.add_argument(
        "interpro_id", help="InterPro ID to fetch (e.g., IPR000001, IPR002048)."
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Maximum number of *unique* proteins to fetch (default: all). More API entries may be processed to reach this limit.",
    )
    parser.add_argument(
        "--no-sequences", action="store_true", help="Skip sequence fetching entirely."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory for files (default: data/raw/interpro/IPR######/)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an interrupted job. Skips proteins already in the output file.",
    )
    return parser


def main():
    """Main function to parse arguments and run the data fetching process."""
    parser = create_arg_parser()
    args = parser.parse_args()

    ipr_id = args.interpro_id.upper()
    if not ipr_id.startswith("IPR") or len(ipr_id) != 9:
        print(f"Error: '{ipr_id}' is not a valid InterPro ID format (e.g., IPR000001).")
        return 1

    # --- Prepare for Saving ---
    date_str = datetime.now().strftime("%Y%m%d")
    output_dir = (
        args.output_dir
        if args.output_dir
        else Path("data") / "raw" / "interpro" / ipr_id
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    json_filename = output_dir / f"{ipr_id}_data_{date_str}.json"
    log_filename = output_dir / f"{ipr_id}_progress_log.jsonl"

    # --- Handle Resume ---
    existing_proteins = {}
    entry_info = None
    if args.resume and log_filename.exists():
        print(f"Resume flag set. Loading existing data from {log_filename}...")
        with log_filename.open("r") as logfile:
            for line in logfile:
                try:
                    protein_entry = json.loads(line)
                    protein_acc = protein_entry["metadata"]["accession"]
                    if protein_acc not in existing_proteins:
                        existing_proteins[protein_acc] = {
                            "metadata": protein_entry["metadata"],
                            "locations": set(),
                        }
                    # This logic assumes locations are not in the primary log
                    # but are reconstructed. For simplicity, we just add the protein.
                    # The full location data will be rebuilt on the final fetch.
                except json.JSONDecodeError:
                    tqdm.write(f"Warning: Skipping corrupted line in log: {line}")

        # To get the entry_info, we need to check the final JSON if it exists
        if json_filename.exists():
            try:
                with json_filename.open("r") as f:
                    final_data = json.load(f)
                    entry_info = final_data.get("entry_info")
            except (json.JSONDecodeError, IOError):
                entry_info = None  # Will be refetched if necessary

    # --- Fetch Data ---
    entry_info, proteins_data = fetch_interpro_data(
        ipr_id,
        max_unique_proteins=args.limit,
        existing_proteins=existing_proteins,
        entry_info=entry_info,
        log_filename=log_filename,
    )
    if not proteins_data:
        print("Failed to fetch data from InterPro. Exiting.")
        return 1

    # --- Save Final Data ---
    csv_filename = output_dir / f"{ipr_id}_proteins_{date_str}.csv"
    fasta_filename = output_dir / f"{ipr_id}_sequences_{date_str}.fasta"

    combined_data = {
        "entry_info": entry_info,
        "proteins_data": proteins_data,
        "fetch_time": datetime.now().isoformat(),
        "fetch_args": vars(args),
    }
    # Final save to ensure the file is perfectly up to date
    save_to_json(combined_data, json_filename)
    save_to_csv(proteins_data, csv_filename, ipr_id)

    sequences = {}
    protein_ids = list(proteins_data.get("proteins", {}).keys())
    if protein_ids and not args.no_sequences:
        sequences = fetch_protein_sequences(protein_ids)
        save_sequences_fasta(sequences, fasta_filename)

    # --- Final Summary ---
    print("\n----- Fetch Complete -----")
    print(
        f"Entry: {entry_info.get('metadata', {}).get('name', {}).get('name', 'Unknown')} ({ipr_id})"
    )
    print(f"Unique proteins stored: {len(protein_ids):,}")
    if not args.no_sequences:
        print(f"Sequences fetched: {len(sequences):,}")
    else:
        print("Sequence fetching skipped.")
    print("Files created:")
    print(f"  - {json_filename}")
    print(f"  - {csv_filename}")
    if not args.no_sequences and sequences:
        print(f"  - {fasta_filename}")
    print("--------------------------")

    # Clean up log file after successful completion
    if log_filename.exists():
        try:
            log_filename.unlink()
            print(f"Temporary log file {log_filename} removed.")
        except OSError as e:
            print(f"Error removing log file {log_filename}: {e}")

    return 0


if __name__ == "__main__":
    exit(main())
