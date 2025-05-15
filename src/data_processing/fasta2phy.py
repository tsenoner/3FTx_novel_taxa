#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def fasta_to_relaxed_phy_aligned(
    fasta_file: str, phy_file: str, lalign: bool = False
):
    """
    Convert a FASTA file to a relaxed PHYLIP (PHY) file with aligned sequences.

    :param fasta_file: Path to the input FASTA file.
    :param phy_file: Path to the output PHYLIP file.
    :param lalign: If True, aligns the sequences by padding the sequence names with spaces.
    """
    # Initialize an empty dictionary to store the sequences
    sequences = {}

    # Define the special characters that should be replaced
    special_chars = r"[ \(\):;,\[\]/\+']"

    # Open the FASTA file
    with open(fasta_file, "r") as f:
        sequence_name = None
        for line in f:
            # Remove trailing whitespace
            line = line.strip()

            # If the line starts with '>', it's a new sequence
            if line.startswith(">"):
                # Use the whole line (except for '>') as the name
                sequence_name = line[1:]
                # Replace special characters with an underscore
                sequence_name = re.sub(special_chars, "_", sequence_name)
                sequences[sequence_name] = ""
            elif sequence_name:  # Append only if a sequence name has been found
                # Otherwise, append the line to the current sequence
                sequences[sequence_name] += line
            else:
                raise ValueError(
                    "Invalid FASTA format: sequence without a name."
                )

    # Check that all sequences are of the same length
    sequence_length = len(next(iter(sequences.values())))
    if any(len(seq) != sequence_length for seq in sequences.values()):
        raise ValueError("All sequences must have the same length.")

    # Get the number of sequences
    num_sequences = len(sequences)

    # Get the length of the longest sequence name
    max_name_length = max(len(name) for name in sequences.keys())

    # Open the output file
    with open(phy_file, "w") as f:
        # Write the first line
        f.write(f"{num_sequences} {sequence_length}\n")

        # Write each sequence
        for name, sequence in sequences.items():
            # Add spaces to the end of the name to align the sequences if lalign is True
            padded_name = name.ljust(max_name_length + 1) if lalign else name
            f.write(f"{padded_name} {sequence}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a FASTA file to a relaxed PHYLIP (PHY) file with aligned"
            " sequences."
        )
    )
    parser.add_argument(
        "fasta_file", type=str, help="Path to the input FASTA file."
    )
    parser.add_argument(
        "phy_file",
        type=str,
        nargs="?",
        help=(
            "Path to the output PHYLIP file. Defaults to the same path/name as"
            " the FASTA file with a .phy extension."
        ),
    )
    parser.add_argument(
        "--lalign",
        action="store_true",
        help="Align sequences by padding the sequence names with spaces.",
    )

    args = parser.parse_args()

    if args.phy_file is None:
        # If the phy_file argument is not provided, use the same base name as the FASTA file but with a .phy extension
        fasta_path = Path(args.fasta_file)
        args.phy_file = fasta_path.with_suffix(".phy")

    fasta_to_relaxed_phy_aligned(
        args.fasta_file, args.phy_file, lalign=args.lalign
    )
    print(f"Converted {args.fasta_file} to {args.phy_file}.")


if __name__ == "__main__":
    main()
