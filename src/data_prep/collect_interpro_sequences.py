#!/usr/bin/env python3
"""
Collect Full Protein Sequences from InterPro Data

This script collects all complete protein sequences from InterPro subfolders
and merges them into a single FASTA file. Sequences are deduplicated by UniProt ID.
"""

import argparse
from pathlib import Path


def parse_fasta(fasta_file):
    """Parse FASTA file and return dict mapping sequence IDs to sequences."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(sequences, output_file, wrap_width=80):
    """Write sequences to FASTA file."""
    with open(output_file, "w") as f:
        for seq_id, sequence in sorted(sequences.items()):
            f.write(f">{seq_id}\n")
            for i in range(0, len(sequence), wrap_width):
                f.write(sequence[i : i + wrap_width] + "\n")


def collect_sequences(input_dir, output_dir):
    """Collect all sequences from InterPro subfolders."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    all_sequences = {}

    # Iterate through InterPro subfolders
    for subfolder in sorted(input_path.iterdir()):
        if not subfolder.is_dir():
            continue

        interpro_id = subfolder.name

        # Find FASTA file in subfolder
        fasta_file = next(
            (f for f in subfolder.iterdir() if f.suffix.lower() in [".fasta", ".fa"]),
            None,
        )

        if not fasta_file:
            continue

        # Load and count sequences
        sequences = parse_fasta(fasta_file)
        new_count = sum(1 for seq_id in sequences if seq_id not in all_sequences)
        all_sequences.update(sequences)

        print(f"{interpro_id}: found {len(sequences)}, added {new_count}")

    # Write merged sequences
    output_fasta = output_path / "interpro_complete_sequences.fasta"
    write_fasta(all_sequences, output_fasta)

    print(f"\nTotal: {len(all_sequences)} unique sequences -> {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Collect complete protein sequences from InterPro"
    )
    parser.add_argument("input_dir", help="Directory containing InterPro subfolders")
    parser.add_argument("output_dir", help="Output directory for merged sequences")
    args = parser.parse_args()
    collect_sequences(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
