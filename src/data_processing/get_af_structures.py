#!/usr/bin/env python3
"""
Downloads available AlphaFold structures for a list of UniProt IDs from a FASTA file.
Saves PDBs to a specified directory and writes missing IDs to a text file.

Usage:
    python fetch_alphafold_structures.py \
        --fasta merged_sanitized.fasta \
        --out_dir raw/pdb \
        --missing_txt missing_ids.txt

Requirements:
    requests
    biopython
    tqdm
"""
import argparse
import os
import requests
from Bio import SeqIO
from tqdm import tqdm

ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


def check_alphafold(uniprot_id):
    """
    Check via HEAD request if the AlphaFold PDB exists for this UniProt ID.
    Returns True if available (status code 200), False otherwise.
    """
    url = ALPHAFOLD_URL_TEMPLATE.format(uniprot=uniprot_id)
    try:
        r = requests.head(url, timeout=10)
        return r.status_code == 200
    except requests.RequestException:
        return False


def download_structure(uniprot_id, out_dir):
    """
    Download the AlphaFold PDB for a given UniProt ID into out_dir.
    """
    url = ALPHAFOLD_URL_TEMPLATE.format(uniprot=uniprot_id)
    local_path = os.path.join(out_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
    try:
        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return True
    except requests.RequestException as e:
        print(f"Failed to download {uniprot_id}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download AlphaFold structures for UniProt IDs from a FASTA file."
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Input FASTA file with UniProt IDs as sequence headers"
    )
    parser.add_argument(
        "--out_dir", required=True,
        help="Directory to save downloaded PDB files"
    )
    parser.add_argument(
        "--missing_txt", required=True,
        help="Output text file to list UniProt IDs without available structures"
    )
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    # Read all records
    records = list(SeqIO.parse(args.fasta, "fasta"))
    missing = []

    # Loop with progress bar
    for rec in tqdm(records, desc="Checking AlphaFold entries", unit="entry"):
        uid = rec.id
        if check_alphafold(uid):
            success = download_structure(uid, args.out_dir)
            if not success:
                missing.append(uid)
        else:
            missing.append(uid)

    # Write missing IDs to file
    if missing:
        with open(args.missing_txt, 'w') as fh:
            for uid in missing:
                fh.write(uid + "\n")
        print(f"\nMissing AlphaFold structures for {len(missing)} IDs. See {args.missing_txt}")
    else:
        print("\nAll entries have AlphaFold structures.")

if __name__ == '__main__':
    main()
