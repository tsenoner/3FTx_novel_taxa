#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def parse_signalp_predictions(prediction_file_path: Path, probability_cutoff: float):
    """Parses the prediction_results.txt file from SignalP.

    Args:
        prediction_file_path (Path): Path to the prediction_results.txt file.
        probability_cutoff (float): The minimum probability for SP prediction to be processed.

    Returns:
        dict: A dictionary mapping sequence IDs to their cleavage information
              (cs_end: 1-indexed end of SP, prob: probability)
              for SP predictions with probability > cutoff.
    """
    cleavage_data = {}
    # Regex to find CS position (P1 from P1-P2) and probability
    # Example: CS pos: 22-23. Pr: 0.9501
    cs_regex = re.compile(r"CS pos: (\d+)-\d+\. Pr: (\d(?:\.\d+)?)")

    try:
        with open(prediction_file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                # Expected columns: ID, Prediction, OTHER_score, SP_score, CS_Position_details
                if len(parts) < 5:
                    continue

                seq_id = parts[0]
                prediction_type = parts[1]
                cs_info_str = parts[4]

                if prediction_type == "SP" and cs_info_str:
                    match = cs_regex.search(cs_info_str)
                    if match:
                        # P1 from "CS pos: P1-P2" is the last residue of the SP (1-indexed)
                        last_sp_residue_1_indexed = int(match.group(1))
                        prob = float(match.group(2))
                        if prob > probability_cutoff:
                            cleavage_data[seq_id] = {
                                "cs_end": last_sp_residue_1_indexed,
                                "prob": prob,
                            }
    except FileNotFoundError:
        print(f"Error: Predictions file not found at {prediction_file_path}")
        return None
    except Exception as e:
        print(f"Error parsing predictions file {prediction_file_path}: {e}")
        return None
    return cleavage_data


def read_fasta_sequences(fasta_file_path: Path):
    """Reads a FASTA file and yields sequence data.

    Args:
        fasta_file_path (Path): Path to the FASTA file.

    Yields:
        tuple: (header (str, without '>'), sequence (str), sequence_id (str))
    """
    current_header = None
    current_seq_parts = []
    try:
        with open(fasta_file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_header is not None:  # Yield previous sequence if exists
                        full_header = current_header
                        seq_id = full_header.split()[0]  # ID is usually the first word
                        yield full_header, "".join(current_seq_parts), seq_id
                    current_header = line[1:]  # Store header without '>'
                    current_seq_parts = []
                elif (
                    current_header is not None
                ):  # Only append if we are inside a sequence block
                    current_seq_parts.append(line)

            if current_header is not None:  # Yield the last sequence in the file
                full_header = current_header
                seq_id = full_header.split()[0]
                yield full_header, "".join(current_seq_parts), seq_id
    except FileNotFoundError:
        print(f"Error: Original FASTA file not found at {fasta_file_path}")
        raise  # Re-raise to be caught in main
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file_path}: {e}")
        raise  # Re-raise


def main():
    parser = argparse.ArgumentParser(
        description="Processes SignalP output to create a FASTA file with signal peptides removed."
    )
    parser.add_argument(
        "output_dir",
        type=Path,  # Use Path type for argparse
        help="Path to the SignalP output directory (e.g., 'my_analysis/signalp_results'). "
        "The script expects 'prediction_results.txt' in this directory.",
    )
    parser.add_argument(
        "input_fasta",
        type=Path,
        help="Path to the original input FASTA file (e.g., 'my_data/sequences.fasta').",
    )
    parser.add_argument(
        "-p",
        "--probability_cutoff",
        type=float,
        default=0.50,
        help="Minimum probability for a signal peptide prediction to be processed (default: 0.50).",
    )
    args = parser.parse_args()

    output_dir_path = args.output_dir  # This is now a Path object
    predictions_file = output_dir_path / "prediction_results.txt"

    # Original FASTA path is derived from the output_dir_path name with .fasta extension
    # Assumes output_dir_path refers to a directory like 'toxprot_2025_subset'
    # and original fasta is 'toxprot_2025_subset.fasta' in the same parent directory as output_dir_path
    original_fasta_path = args.input_fasta  # Use the new argument

    # Output processed FASTA path, e.g., 'toxprot_2025_subset/toxprot_2025_subset_mature.fasta'
    # Save it inside the output_dir_path
    processed_fasta_path = output_dir_path / (output_dir_path.name + "_mature.fasta")

    print(f"Reading predictions from: {predictions_file}")
    print(f"Reading original FASTA from: {original_fasta_path}")
    print(f"Will write processed FASTA to: {processed_fasta_path}")

    sp_cleavage_data = parse_signalp_predictions(
        predictions_file, args.probability_cutoff
    )
    if sp_cleavage_data is None:
        return  # Error already printed by parser

    num_processed = 0
    num_total = 0
    sequences_written = 0

    try:
        with open(processed_fasta_path, "w") as outfile:
            for header, sequence, seq_id in read_fasta_sequences(original_fasta_path):
                num_total += 1
                if seq_id in sp_cleavage_data:
                    cleavage_site_1_indexed = sp_cleavage_data[seq_id]["cs_end"]
                    mature_sequence = sequence[cleavage_site_1_indexed:]

                    outfile.write(f">{header}\n")
                    outfile.write(mature_sequence + "\n")  # Write sequence as one line
                    num_processed += 1
                else:
                    outfile.write(f">{header}\n")
                    outfile.write(sequence + "\n")  # Write sequence as one line
                sequences_written += 1

        print(f"\nSuccessfully processed {num_total} sequences.")
        print(
            f"Removed signal peptides from {num_processed} sequences (where Pr > {args.probability_cutoff})."
        )
        print(f"Output written to: {processed_fasta_path}")
        if sequences_written == 0 and num_total > 0:
            print(
                "Warning: No sequences were written to the output file, though the input FASTA was read. Check parsing logic or input file contents."
            )
        elif sequences_written == 0 and num_total == 0:
            print(
                "Warning: Original FASTA file might be empty or could not be read correctly."
            )

    except FileNotFoundError:
        print(
            f"Processing aborted due to missing original FASTA file: {original_fasta_path}"
        )
    except Exception as e:
        print(f"An unexpected error occurred during FASTA processing: {e}")


if __name__ == "__main__":
    main()
