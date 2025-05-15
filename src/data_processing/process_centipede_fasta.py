import re
from pathlib import Path
import pandas as pd

# Global constants
OUTPUT_SUBDIR_NAME = "centipede_3ftx_quiver_upar_like"
INPUT_FASTA_PATH_STR = "data/raw/centipede_3ftx_quiver_upar_like/3ftx_Quiver_UPAR-like_in_Centipedes_AW_210425.fasta"
OUTPUT_PARENT_DIR_PATH_STR = "data/interm"
ASSEMBLY_ID_MAP_PATH_STR = (
    "data/raw/centipede_3ftx_quiver_upar_like/centipede_genome_tol_ids.tsv"
)


# --- Mapping File Loading ---
def get_assembly_id_map(file_path_str: str) -> dict[str, str] | None:
    """Loads the AssemblyID to OrganismId mapping."""
    file_path = Path(file_path_str)
    if not file_path.exists():
        print(f"Warning: Assembly ID mapping file not found: {file_path}.")
        return None
    try:
        df = pd.read_csv(file_path, sep="\t")
        if "AssemblyID" not in df.columns or "OrganismId" not in df.columns:
            print(
                f"Warning: Mapping file {file_path} lacks required columns 'AssemblyID' or 'OrganismId'."
            )
            return None
        mapping = dict(zip(df["AssemblyID"], df["OrganismId"]))
        return mapping
    except pd.errors.EmptyDataError:
        print(f"Warning: Mapping file {file_path} is empty.")
        return None
    except Exception as e:
        print(f"Warning: Error reading mapping file {file_path}: {e}")
        return None


# Custom FASTA parser
def custom_fasta_parser(file_path: Path):
    """
    A simple FASTA parser.
    Yields tuples of (header, sequence).
    Header does not include the leading '>'.
    Sequence is a single string with no newlines.
    """
    header = None
    sequence_parts = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence_parts)
                header = line[1:]
                sequence_parts = []
            elif header:
                sequence_parts.append(line)
    if header:
        yield header, "".join(sequence_parts)


def perform_sequence_processing(input_fasta_path: Path, assembly_id_map: dict | None):
    """
    Processes sequences from a FASTA file: filters, removes duplicates, and generates new IDs.
    New IDs are 'Centipede3FTx_NNNN', optionally with '_OrganismID' if a mapping is provided and a match is found.
    Returns a tuple: (processed_records, id_map_dict, total_initial, num_filtered, num_duplicates).
    """
    processed_records = []
    id_map_dict = {}
    filtered_out_sequences_count = 0
    duplicate_sequences_count = 0
    seen_sequences = set()
    seq_counter = 1
    total_sequences_initial_count = 0
    gca_pattern = r"(GCA_\d+\.\d+)"

    for original_id_full, sequence_str in custom_fasta_parser(input_fasta_path):
        total_sequences_initial_count += 1
        match = re.search(gca_pattern, original_id_full)
        id_for_map_lookup = match.group(1) if match else original_id_full.split(" ")[0]

        processed_sequence_str = sequence_str
        has_internal_stop = False
        if "*" in processed_sequence_str:
            if processed_sequence_str.endswith("*"):
                processed_sequence_str = processed_sequence_str[:-1]
            if "*" in processed_sequence_str:
                has_internal_stop = True

        if has_internal_stop or not processed_sequence_str:
            filtered_out_sequences_count += 1
            continue

        if processed_sequence_str in seen_sequences:
            duplicate_sequences_count += 1
            continue
        seen_sequences.add(processed_sequence_str)

        unique_numerical_id_part = f"Centipede3FTx_{seq_counter:04d}"
        final_id_candidate = unique_numerical_id_part

        if assembly_id_map:
            for map_key_assembly_id, organism_id_val in assembly_id_map.items():
                if id_for_map_lookup.startswith(
                    str(map_key_assembly_id)
                ):  # Ensure map_key is str for startswith
                    final_id_candidate = f"{unique_numerical_id_part}_{organism_id_val}"
                    break

        final_new_id = final_id_candidate
        temp_disamb_counter = 1
        while (
            final_new_id in id_map_dict
        ):  # Check against locally accumulating id_map_dict
            final_new_id = f"{final_id_candidate}_{temp_disamb_counter}"
            temp_disamb_counter += 1

        id_map_dict[final_new_id] = original_id_full
        processed_records.append((final_new_id, processed_sequence_str))
        seq_counter += 1

    return (
        processed_records,
        id_map_dict,
        total_sequences_initial_count,
        filtered_out_sequences_count,
        duplicate_sequences_count,
    )


def generate_report_and_write_outputs(
    input_fasta_path: Path,
    output_dir: Path,
    processed_records: list,
    id_map_data: dict,
    total_initial: int,
    num_filtered: int,
    num_duplicates: int,
    assembly_id_map_file_path_str: str,  # Used for report
    assembly_id_map_loaded: bool,  # Used for report
):
    """Generates report and writes output FASTA and mapping CSV files."""
    output_fasta_path = output_dir / "centipede_3ftx_processed.fasta"
    mapping_csv_path = output_dir / "sequence_mapping.csv"
    report_txt_path = output_dir / "centipede_3ftx_processing_report.txt"

    # Write FASTA
    with open(output_fasta_path, "w") as out_handle:
        for new_id, seq_str in processed_records:
            out_handle.write(f">{new_id}\n{seq_str}\n")

    # Write Mapping CSV
    if id_map_data:
        mapping_df = pd.DataFrame(
            list(id_map_data.items()), columns=["New ID", "Original Header"]
        )
        mapping_df.to_csv(mapping_csv_path, index=False)
    else:
        with open(mapping_csv_path, "w") as f:
            f.write("New ID,Original Header\n")

    # Generate Report
    num_written = len(processed_records)
    report_lines = [
        "FASTA Processing Report",
        "=======================",
        f"Input FASTA: {input_fasta_path}",
        f"Output directory: {output_dir}",
    ]
    map_status = (
        " (Loaded successfully)"
        if assembly_id_map_loaded
        else " (File not found, failed to load, or empty)"
    )
    report_lines.append(f"Assembly ID Map: {assembly_id_map_file_path_str}{map_status}")
    report_lines.extend(
        [
            "\n-------------------",
            f"Total sequences read from input: {total_initial}",
            f"Sequences written to FASTA: {num_written}",
            f"Sequences filtered (internal stop/empty): {num_filtered}",
            f"Duplicate sequences removed: {num_duplicates}",
        ]
    )

    with open(report_txt_path, "w") as report_handle:
        report_handle.write("\n".join(report_lines) + "\n")


def main():
    """Main function to orchestrate FASTA processing."""
    input_fasta_path = Path(INPUT_FASTA_PATH_STR)
    output_parent_dir = Path(OUTPUT_PARENT_DIR_PATH_STR)

    output_dir = output_parent_dir / OUTPUT_SUBDIR_NAME
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing {input_fasta_path}...")
    assembly_id_map = get_assembly_id_map(ASSEMBLY_ID_MAP_PATH_STR)
    assembly_map_loaded_successfully = assembly_id_map is not None

    (
        processed_records,
        id_map_dict,
        total_initial,
        num_filtered,
        num_duplicates,
    ) = perform_sequence_processing(input_fasta_path, assembly_id_map)

    generate_report_and_write_outputs(
        input_fasta_path=input_fasta_path,
        output_dir=output_dir,
        processed_records=processed_records,
        id_map_data=id_map_dict,
        total_initial=total_initial,
        num_filtered=num_filtered,
        num_duplicates=num_duplicates,
        assembly_id_map_file_path_str=ASSEMBLY_ID_MAP_PATH_STR,
        assembly_id_map_loaded=assembly_map_loaded_successfully,
    )

    print(f"\nProcessing complete.\nOutput files written to: {output_dir}")


if __name__ == "__main__":
    main()
