#!/usr/bin/env python3
"""
UniProt to InterPro Domain Fetcher

This script takes a list of UniProt protein identifiers, or a FASTA file from
which identifiers can be parsed, and retrieves all associated InterPro domain
annotations, saving the results to a CSV file.

Usage:
  # From a list of IDs, overwriting any previous output file
  python src/interpro_retrieval/uniprot_to_interpro.py --ids P99999 Q16540 --output results.csv

  # From a FASTA file, resuming from a previous run
  python src/interpro_retrieval/uniprot_to_interpro.py --input-file proteins.fasta --output results.csv --resume
"""

import argparse
import csv
import os
import time
from typing import Any, Dict, List, Optional

import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

# --- Constants ---
INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api"


# --- API Fetching ---
def create_session() -> requests.Session:
    """Creates a requests session with retry logic."""
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    return session


def fetch_domains_for_uniprot_id(
    uniprot_id: str, session: requests.Session
) -> Optional[Dict[str, Any]]:
    """
    Fetches InterPro annotations for a single UniProt accession number.
    Tries the specific endpoint if the general one returns no entries.

    Args:
        uniprot_id: The UniProt accession (e.g., "P00750").
        session: The requests session to use for the API call.

    Returns:
        A dictionary containing the API response data, or None on failure.
    """
    # First, try the general protein endpoint
    url = f"{INTERPRO_API_URL}/protein/uniprot/{uniprot_id}"
    try:
        response = session.get(url, timeout=15)
        if response.status_code == 404:
            tqdm.write(f"Warning: No data found for UniProt ID: {uniprot_id}")
            return None
        response.raise_for_status()
        data = response.json()

        # If the primary endpoint doesn't return entries directly,
        # it's likely a metadata-only response. Use the specific endpoint for entries.
        if "entries" not in data:
            entry_url = (
                f"{INTERPRO_API_URL}/entry/interpro/protein/uniprot/{uniprot_id}"
            )
            entry_response = session.get(entry_url, timeout=15)
            # A 404 here means the protein exists but has no InterPro entries.
            if entry_response.status_code == 404:
                data["entries"] = []
            else:
                entry_response.raise_for_status()
                # The new response has the entries under the "results" key
                entry_data = entry_response.json()
                data["entries"] = entry_data.get("results", [])

        return data

    except requests.exceptions.RequestException as e:
        tqdm.write(f"Warning: Error fetching data for {uniprot_id}: {e}")
        return None


# --- FASTA Parsing ---
def parse_fasta_file(filepath: str) -> List[str]:
    """Parses UniProt IDs from a FASTA file."""
    uniprot_ids = set()  # Use a set to handle duplicates automatically
    try:
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    if "|" in header:
                        # Assumes >db|ID|... format
                        try:
                            uniprot_id = header.split("|")[1]
                            uniprot_ids.add(uniprot_id)
                        except IndexError:
                            tqdm.write(
                                f"Warning: Could not parse UniProt header: {header}"
                            )
                    else:
                        # Assumes >ID format
                        uniprot_id = header.split()[0]
                        uniprot_ids.add(uniprot_id)
    except FileNotFoundError:
        tqdm.write(f"Error: Input file not found at {filepath}")
        return []

    return list(uniprot_ids)


# --- Data Processing and Saving ---
def prepare_csv_rows(uniprot_id: str, data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Prepares a list of dictionary rows for CSV export."""
    rows = []
    protein_meta = data.get("metadata", {})

    # Extract organism info
    source_organism = protein_meta.get("source_organism", {})

    base_row_data = {
        "protein_id": uniprot_id,
        "protein_name": protein_meta.get("name"),
        "tax_id": source_organism.get("taxId"),
        "scientific_name": source_organism.get("scientificName"),
        "sequence_length": protein_meta.get("length"),
    }

    entries = data.get("entries", [])
    if not entries:
        # Still record the protein even if it has no domain entries
        row = base_row_data.copy()
        row.update(
            {
                "domain_start": "N/A",
                "domain_end": "N/A",
                "signature_name": "N/A",
                "interpro_domain_name": "N/A",
            }
        )
        rows.append(row)
        return rows

    for entry in entries:
        entry_meta = entry.get("metadata", {})
        signature_name = entry_meta.get("accession")
        interpro_domain_name = entry_meta.get("name")

        locations_found = False
        for protein_data in entry.get("proteins", []):
            if protein_data.get("accession") == uniprot_id.lower():
                for loc in protein_data.get("entry_protein_locations", []):
                    locations_found = True
                    for frag in loc.get("fragments", []):
                        row = base_row_data.copy()
                        row.update(
                            {
                                "domain_start": frag.get("start"),
                                "domain_end": frag.get("end"),
                                "signature_name": signature_name,
                                "interpro_domain_name": interpro_domain_name,
                            }
                        )
                        rows.append(row)

        if not locations_found:
            # Record the domain signature even if no locations are found for it
            row = base_row_data.copy()
            row.update(
                {
                    "domain_start": "N/A",
                    "domain_end": "N/A",
                    "signature_name": signature_name,
                    "interpro_domain_name": interpro_domain_name,
                }
            )
            rows.append(row)

    return rows


def save_to_csv(data_rows: List[Dict[str, Any]], filename: str):
    """Saves a list of data rows to a CSV file. Writes header even if no data."""
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
    with open(filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        if data_rows:
            writer.writerows(data_rows)


# --- Main Execution ---
def create_arg_parser() -> argparse.ArgumentParser:
    """Creates the argument parser for the script."""
    parser = argparse.ArgumentParser(
        description="Fetch all InterPro domain annotations for UniProt identifiers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From a list of IDs, overwriting any previous output file
  %(prog)s --ids P99999 Q16540 --output results.csv

  # From a FASTA file, resuming from a previous run
  %(prog)s --input-file proteins.fasta --output results.csv --resume
""",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--ids",
        nargs="+",
        metavar="UNIPROT_ID",
        help="One or more UniProt protein identifiers.",
    )
    group.add_argument(
        "--input-file", help="Path to a FASTA file to parse UniProt IDs from headers."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to save the results as a CSV file.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an interrupted job. Skips IDs already present in the output file.",
    )
    return parser


def main():
    """Main function to parse arguments and run the data fetching process."""
    parser = create_arg_parser()
    args = parser.parse_args()
    session = create_session()

    # 1. Get all UniProt IDs from the source
    source_ids = []
    if args.input_file:
        source_ids = parse_fasta_file(args.input_file)
        if not source_ids:
            return 1  # Exit if file not found or empty
    elif args.ids:
        source_ids = args.ids

    # 2. Handle resume logic
    processed_ids = set()
    is_resuming = args.resume and os.path.exists(args.output)
    if is_resuming:
        try:
            with open(args.output, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if "protein_id" in row:
                        processed_ids.add(row["protein_id"])
            tqdm.write(
                f"Found {len(processed_ids)} previously processed IDs in {args.output}. Resuming."
            )
        except (IOError, csv.Error) as e:
            tqdm.write(
                f"Warning: Could not read resume file {args.output}. Starting from scratch. Error: {e}"
            )
            is_resuming = False  # Force overwrite

    ids_to_process = [
        pid.upper() for pid in source_ids if pid.upper() not in processed_ids
    ]

    if not ids_to_process:
        tqdm.write("All UniProt IDs have already been processed. Nothing to do.")
        return 0

    # 3. Open file and process IDs
    file_mode = "a" if is_resuming else "w"
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

    with open(args.output, file_mode, newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not is_resuming:
            writer.writeheader()

        # The main loop with a progress bar
        for uniprot_id in tqdm(
            ids_to_process, desc="Fetching domain data", unit="protein", leave=True
        ):
            data = fetch_domains_for_uniprot_id(uniprot_id, session)
            if data:
                csv_rows = prepare_csv_rows(uniprot_id, data)
                if csv_rows:
                    writer.writerows(csv_rows)
                    csvfile.flush()  # Save progress after each ID
            time.sleep(0.1)  # Be a good API citizen

    tqdm.write(f"\nProcessing complete. Results saved to {args.output}")
    return 0


if __name__ == "__main__":
    exit(main())
