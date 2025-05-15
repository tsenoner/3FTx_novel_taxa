#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merges two FASTA files.
- Sanitizes sequences (checks for ambiguous residues, removes duplicates).
- Cleans FASTA headers for consistency.
- Outputs a single sanitized FASTA file.
- All operations and errors are logged to 'merge_process_report.txt'
  in the output directory.

Paths are hardcoded in this script.
Output directory is derived from base Path objects.
"""

import re
import sys
from pathlib import Path
from typing import List, Set, Tuple, Optional

# --- Hardcoded Paths (using pathlib) ---
DATA_DIR = Path("data")
RAW_DATA_DIR = DATA_DIR / "raw"
INTERM_DATA_DIR = DATA_DIR / "interm"

FASTA1_PATH = RAW_DATA_DIR / "uniprot" / "3FTx_related.fasta"
CENTIPEDE_PROCESSED_DIR = INTERM_DATA_DIR / "centipede_3ftx_quiver_upar_like"
FASTA2_PATH = CENTIPEDE_PROCESSED_DIR / "centipede_3ftx_processed.fasta"
OUTPUT_DIR = INTERM_DATA_DIR / "merged"

# --- Constants ---
HEADER_CLEANER = re.compile(r"[\W\s]")
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# --- Utility Functions ---
def write_report(output_dir: Path, report_filename: str, messages: List[str]) -> None:
    """Writes report messages to a file."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        report_file_path = output_dir / report_filename
        with open(report_file_path, "w") as f:
            for msg in messages:
                f.write(msg + "\n")
    except Exception as e:
        print(
            f"CRITICAL: Failed to write report file to {output_dir / report_filename}: {e}",
            file=sys.stderr,
        )
        print("--- Accumulated Messages ---", file=sys.stderr)
        for msg in messages:
            print(msg, file=sys.stderr)


# --- FASTA Parsing ---
def custom_fasta_parser(
    file_path: Path, report_messages: List[str]
) -> List[Tuple[str, str]]:
    """
    A simple FASTA parser. Handles multi-line sequences.
    Returns a list of (header, sequence) tuples.
    Header does not include the leading '>'.
    Sequence is a single string with no newlines.
    Appends errors to report_messages.
    """
    entries: List[Tuple[str, str]] = []
    header: Optional[str] = None
    sequence_parts: List[str] = []
    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        entries.append((header, "".join(sequence_parts)))
                    header = line[1:]
                    sequence_parts = []
                elif header is not None:
                    sequence_parts.append(line)
            if header is not None:
                entries.append((header, "".join(sequence_parts)))
    except FileNotFoundError:
        report_messages.append(f"ERROR: Input FASTA file not found: {file_path}")
    return entries


# --- Sanitization Functions ---
def clean_fasta_header(header: str) -> str:
    """Clean the FASTA header.
    UniProt: extracts ACCESSION if pattern (tr|ACCESSION|... or sp|ACCESSION|...) matches.
    Otherwise: first part before whitespace, non-alphanumerics to underscores.
    """
    uniprot_match = re.match(r"^(?:tr|sp)\|([A-Za-z0-9]+)\|", header)
    if uniprot_match:
        return uniprot_match.group(1)
    else:
        base_id = header.split(None, 1)[0]
        return HEADER_CLEANER.sub("_", base_id)


def has_ambiguous_residues(sequence: str) -> Tuple[bool, str]:
    """Check for non-standard amino acid residues."""
    ambiguous_residues = set(sequence.upper()) - STANDARD_AMINO_ACIDS
    return bool(ambiguous_residues), ",".join(sorted(list(ambiguous_residues)))


def make_header_unique(
    header: str, existing_headers: Set[str], report_messages: List[str]
) -> str:
    """Ensures the cleaned header is unique, appending numbers if needed."""
    if header not in existing_headers:
        existing_headers.add(header)
        return header
    else:
        original_problem_header = header
        counter = 1
        unique_header = f"{header}_{counter}"
        while unique_header in existing_headers:
            counter += 1
            unique_header = f"{header}_{counter}"
        existing_headers.add(unique_header)
        report_messages.append(
            f"WARNING: Original header resulting in '{original_problem_header}' was not unique. Renamed to '{unique_header}'."
        )
        return unique_header


# --- Output Writing ---
def write_processed_fasta(
    file_path: Path, entries: List[Tuple[str, str]], report_messages: List[str]
) -> None:
    """Write a list of (header, sequence) tuples to a FASTA file."""
    report_messages.append(
        f"INFO: Writing {len(entries)} sanitized sequences to FASTA: {file_path}"
    )
    try:
        with file_path.open("w") as outfile:
            for header, seq in entries:
                outfile.write(f">{header}\n{seq}\n")
    except IOError as e:
        report_messages.append(f"ERROR: Could not write FASTA file to {file_path}: {e}")


# --- Main Processing Logic ---
def run_merge_process(
    fasta1_path: Path,
    fasta2_path: Path,
    output_dir_path: Path,
    report_messages: List[str],
) -> None:
    """Main function to merge and sanitize FASTA files."""

    sanitized_fasta_path = output_dir_path / "merged_sanitized.fasta"

    all_raw_entries: List[
        Tuple[str, str, str]
    ] = []  # original_header, sequence, source_tag (kept for origin tracking)

    report_messages.append(
        f"INFO: Reading sequences from {fasta1_path} (source tag: file1)..."
    )
    entries_f1 = custom_fasta_parser(fasta1_path, report_messages)
    for header, seq in entries_f1:
        all_raw_entries.append((header, seq, "file1"))  # Changed source_tag for clarity

    report_messages.append(
        f"INFO: Reading sequences from {fasta2_path} (source tag: file2)..."
    )
    entries_f2 = custom_fasta_parser(fasta2_path, report_messages)
    for header, seq in entries_f2:
        all_raw_entries.append((header, seq, "file2"))  # Changed source_tag for clarity

    report_messages.append(
        f"INFO: Read a total of {len(all_raw_entries)} raw sequences."
    )

    sanitized_entries: List[Tuple[str, str]] = []
    seen_sequences: Set[str] = set()
    unique_clean_headers: Set[str] = set()
    removed_ambiguous_count = 0
    removed_duplicate_seq_count = 0
    removed_empty_count = 0

    for original_header, sequence, source_tag in all_raw_entries:
        output_fasta_header_base = clean_fasta_header(original_header)

        if not sequence:
            report_messages.append(
                f"INFO: Removed '{original_header}' (source: {source_tag}, base header: {output_fasta_header_base}) because sequence was empty."
            )
            removed_empty_count += 1
            continue

        has_ambiguity, ambiguous_str = has_ambiguous_residues(sequence)
        if has_ambiguity:
            report_messages.append(
                f"INFO: Removed '{original_header}' (source: {source_tag}, base header: {output_fasta_header_base}) due to ambiguous residues: {ambiguous_str}."
            )
            removed_ambiguous_count += 1
            continue

        if sequence in seen_sequences:
            report_messages.append(
                f"INFO: Removed '{original_header}' (source: {source_tag}, base header: {output_fasta_header_base}) as a duplicate sequence."
            )
            removed_duplicate_seq_count += 1
            continue
        seen_sequences.add(sequence)

        final_output_header = make_header_unique(
            output_fasta_header_base, unique_clean_headers, report_messages
        )

        sanitized_entries.append((final_output_header, sequence))

    report_messages.append(
        f"INFO: FASTA sanitization complete. Kept {len(sanitized_entries)} unique, valid sequences."
    )
    if removed_ambiguous_count > 0:
        report_messages.append(
            f"INFO: Count of sequences removed due to ambiguous residues: {removed_ambiguous_count}."
        )
    if removed_empty_count > 0:
        report_messages.append(
            f"INFO: Count of sequences removed because they were empty: {removed_empty_count}."
        )
    if removed_duplicate_seq_count > 0:
        report_messages.append(
            f"INFO: Count of duplicate sequences removed: {removed_duplicate_seq_count}."
        )

    write_processed_fasta(sanitized_fasta_path, sanitized_entries, report_messages)
    num_kept_fasta = len(sanitized_entries)

    # --- Final Summary for Report ---
    num_processed_total = len(all_raw_entries)
    report_messages.extend(
        [
            "--- Summary ---",
            f"Processed {num_processed_total} sequences from input FASTA files.",
            f"Output directory: {output_dir_path}",
            f"  - Sanitized FASTA: {sanitized_fasta_path.name} ({num_kept_fasta} sequences written)",
            "Filtering/Removal Details (during FASTA sanitization):",
            f"  - Ambiguous Residues: {removed_ambiguous_count}",
            f"  - Empty Sequences: {removed_empty_count}",
            f"  - Duplicate Sequences: {removed_duplicate_seq_count}",
            "--- Done ---",
        ]
    )


def main() -> None:
    """Main execution function."""
    report_messages: List[str] = ["--- FASTA Merge and Sanitize Report ---"]
    report_file_name = "merge_process_report.txt"

    # Global Path objects (FASTA1_PATH, FASTA2_PATH, OUTPUT_DIR) are used directly.

    critical_error = False
    input_files_to_check = {
        "FASTA1": FASTA1_PATH,
        "FASTA2": FASTA2_PATH,
    }
    for name, path_obj in input_files_to_check.items():
        if not path_obj.exists():
            report_messages.append(
                f"CRITICAL ERROR: Input file for {name} not found: {path_obj}"
            )
            critical_error = True
        elif not path_obj.is_file():
            report_messages.append(
                f"CRITICAL ERROR: Input path for {name} is not a file: {path_obj}"
            )
            critical_error = True

    if not critical_error:
        try:
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            report_messages.append(
                f"INFO: Ensured output directory exists: {OUTPUT_DIR}"
            )
        except Exception as e:
            report_messages.append(
                f"CRITICAL ERROR: Could not create output directory {OUTPUT_DIR}: {e}"
            )
            critical_error = True

    if critical_error:
        write_report(OUTPUT_DIR, report_file_name, report_messages)
        print(
            f"Critical error during setup. See report: {OUTPUT_DIR / report_file_name}",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        run_merge_process(
            fasta1_path=FASTA1_PATH,
            fasta2_path=FASTA2_PATH,
            output_dir_path=OUTPUT_DIR,
            report_messages=report_messages,
        )
    except Exception as e:
        report_messages.append(
            f"CRITICAL UNEXPECTED ERROR during run_merge_process: {e}"
        )
        write_report(OUTPUT_DIR, report_file_name, report_messages)
        print(
            f"Critical unexpected error during processing. See report: {OUTPUT_DIR / report_file_name}",
            file=sys.stderr,
        )
        sys.exit(1)

    write_report(OUTPUT_DIR, report_file_name, report_messages)
    print(f"Processing complete. Sanitized FASTA and report written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
