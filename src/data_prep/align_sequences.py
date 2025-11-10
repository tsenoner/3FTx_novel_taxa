#!/usr/bin/env python3
"""
Generate a multiple sequence alignment using FAMSA2 via pyfamsa.

Usage:
    python align_sequences.py <input_fasta> <output_fasta>

Example:
    python align_sequences.py data/interm/mmseqs/representatives.fasta data/interm/mmseqs/representatives_aligned.fasta
"""

import argparse
import sys
from pathlib import Path
from pyfamsa import Aligner, Sequence


def read_fasta(fasta_path):
    """
    Read sequences from a FASTA file.

    Args:
        fasta_path: Path to input FASTA file

    Yields:
        tuple: (sequence_id, sequence) as bytes
    """
    sequence_id = None
    sequence_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence if exists
                if sequence_id is not None:
                    yield sequence_id, "".join(sequence_parts)

                # Start new sequence
                sequence_id = line[1:].encode()  # Remove '>' and convert to bytes
                sequence_parts = []
            else:
                sequence_parts.append(line)

        # Don't forget the last sequence
        if sequence_id is not None:
            yield sequence_id, "".join(sequence_parts)


def write_fasta(sequences, output_path):
    """
    Write aligned sequences to a FASTA file.

    Args:
        sequences: List of pyfamsa Sequence objects
        output_path: Path to output FASTA file
    """
    with open(output_path, "w") as f:
        for seq in sequences:
            seq_id = seq.id.decode()
            seq_str = seq.sequence.decode()
            f.write(f">{seq_id}\n")
            # Write sequence in lines of 80 characters
            for i in range(0, len(seq_str), 80):
                f.write(seq_str[i : i + 80] + "\n")


def align_sequences(input_fasta, output_fasta, guide_tree="sl", threads=0):
    """
    Align sequences using FAMSA2.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output aligned FASTA file
        guide_tree: Guide tree method ('sl', 'slink', 'upgma', 'nj')
        threads: Number of threads to use (0 = auto)
    """
    print(f"Reading sequences from {input_fasta}...")

    # Read sequences and convert to pyfamsa Sequence objects
    sequences = []
    for seq_id, seq_str in read_fasta(input_fasta):
        sequences.append(Sequence(seq_id, seq_str.encode()))

    print(f"Loaded {len(sequences)} sequences")

    # Create aligner
    print(
        f"Aligning sequences using FAMSA2 (guide tree: {guide_tree}, threads: {threads or 'auto'})..."
    )
    aligner = Aligner(guide_tree=guide_tree, threads=threads)

    # Perform alignment
    msa = aligner.align(sequences)

    # Write output
    print(f"Writing aligned sequences to {output_fasta}...")
    write_fasta(msa, output_fasta)

    print("✓ Alignment complete!")
    print(f"  Input sequences: {len(sequences)}")
    print(f"  Alignment length: {len(msa[0].sequence)}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a multiple sequence alignment using FAMSA2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default settings
  python align_sequences.py input.fasta output.fasta

  # Use UPGMA guide tree method
  python align_sequences.py input.fasta output.fasta --guide-tree upgma

  # Use 8 threads
  python align_sequences.py input.fasta output.fasta --threads 8
        """,
    )

    parser.add_argument(
        "input", type=Path, help="Input FASTA file with unaligned sequences"
    )

    parser.add_argument(
        "output", type=Path, help="Output FASTA file for aligned sequences"
    )

    parser.add_argument(
        "--guide-tree",
        choices=["nj", "sl", "slink", "upgma"],
        default="upgma",
        help="Guide tree construction method (default: upgma)",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Number of threads to use (0 = auto, default: 0)",
    )

    args = parser.parse_args()

    # Check input file exists
    if not args.input.exists():
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Run alignment
    align_sequences(args.input, args.output, args.guide_tree, args.threads)


if __name__ == "__main__":
    main()
