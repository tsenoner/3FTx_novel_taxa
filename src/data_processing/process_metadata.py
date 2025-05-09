import pandas as pd
import taxopy
from pathlib import Path
import sys
import logging
import json
from typing import Dict, List

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Configuration ---
# Define file paths using pathlib
DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw"
INPUT_CSV_PATH = RAW_DIR / "uniprot_transmembrane_data.csv"
INPUT_TSV_PATH = RAW_DIR / "uniprot_20250415.tsv"
INTERPRO_MAPPING_JSON = RAW_DIR / "interpro_data.json"
INTERM_DIR = DATA_DIR / "interm"
OUTPUT_CSV_PATH = INTERM_DIR / "metadata.csv"
SUMMARY_OUTPUT_CSV_PATH = INTERM_DIR / "metadata_summary.csv"

# Columns to keep initially from uniprot_transmembrane_data.csv
INITIAL_COLUMNS_TO_KEEP = [
    "uniprot_id",
    "protein_name",  # Needed temporarily for ID extraction
    "is_transmembrane",
    "is_helix",
    "is_beta",
    "has_signal",
    # "signal_start", # Removed as requested
    # "signal_end",   # Removed as requested
]

# Define InterPro ID Priority Order
INTERPRO_PRIORITY_ORDER = [
    "IPR003571",  # Snake_3FTx (Family)
    "IPR031424",  # QVR-like (Family)
    "IPR026110",  # LY6G5C (Family)
    "IPR026524",  # LY6G6d/LY6G6f (Family)
    "IPR038773",  # LY6G5B (Family)
    "IPR039457",  # LYPD6-like (Family)
    "IPR039237",  # LY6G6C (Family)
    "IPR042339",  # Ly6D (Family)
    "IPR051445",  # LY6H/LY6L_nAChR_modulators (Family)
    "IPR016054",  # LY6_UPA_recep-like (Domain)
    "IPR035076",  # Toxin/TOLIP (Domain)
    "IPR045860",  # Snake_toxin-like_sf (Homologous Superfamily)
    # Add other IDs here if needed, in descending order of priority
]

# --- Taxonomy Helper Functions (Modified) ---


def setup_taxopy_db_paths():
    """Setup and return the taxopy database paths."""
    try:
        home_dir = Path.home()
        cache_dir = home_dir / ".cache"
        db_dir = cache_dir / "taxopy_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        nodes_file = db_dir / "nodes.dmp"
        names_file = db_dir / "names.dmp"
        merged_file = db_dir / "merged.dmp"
        return db_dir, nodes_file, names_file, merged_file
    except Exception as e:
        logging.error(f"Error setting up taxopy db paths: {e}")
        raise


def initialize_taxopy_db():
    """Initialize and return the taxopy database connection."""
    db_dir, nodes_file, names_file, merged_file = setup_taxopy_db_paths()
    taxdb_params = {}
    use_dir_init = False

    if nodes_file.exists() and names_file.exists():
        logging.info(f"Found existing taxopy nodes and names files in {db_dir}")
        taxdb_params["nodes_dmp"] = str(nodes_file)
        taxdb_params["names_dmp"] = str(names_file)
        if merged_file.exists() and merged_file.is_file():
            taxdb_params["merged_dmp"] = str(merged_file)
            logging.info(f"Found existing taxopy merged file: {merged_file}")
        else:
            logging.warning(f"Merged file not found or invalid: {merged_file}")
    else:
        logging.warning(
            f"Taxonomy database files (nodes/names) not found in {db_dir}. Will attempt download/init."
        )
        use_dir_init = True

    try:
        if not use_dir_init and len(taxdb_params) >= 2:
            logging.info(
                f"Initializing taxopy TaxDb with files: {list(taxdb_params.keys())}"
            )
            taxdb = taxopy.TaxDb(**taxdb_params)
        else:
            logging.info(f"Initializing taxopy TaxDb using directory: {db_dir}")
            taxdb = taxopy.TaxDb(taxdb_dir=str(db_dir), keep_files=True)
        logging.info("Taxopy database initialized successfully.")
        return taxdb
    except Exception as e:
        logging.error(f"Failed to initialize taxopy database: {e}")
        raise


def get_taxonomy_levels(taxon_id, taxdb):
    """Get phylum, order, family, genus, and species info for a taxon ID."""
    result = {
        "phylum": pd.NA,
        "order": pd.NA,
        "family": pd.NA,
        "genus": pd.NA,
        "species": pd.NA,
    }  # Use pandas NA for missing data
    if pd.isna(taxon_id):
        return result

    try:
        taxon_id_int = int(taxon_id)
        taxon = taxopy.Taxon(taxon_id_int, taxdb)
        ranks = taxon.rank_name_dictionary
        result["phylum"] = ranks.get("phylum", pd.NA)
        result["order"] = ranks.get("order", pd.NA)
        result["family"] = ranks.get("family", pd.NA)
        result["genus"] = ranks.get("genus", pd.NA)
        result["species"] = ranks.get(
            "species", pd.NA
        )  # Or taxon.name if species rank is not explicitly found but it's the lowest known rank
    except (ValueError, taxopy.exceptions.TaxidError, TypeError) as e:
        # Log only once per unique problematic ID? For now, log lightly.
        logging.debug(f"Could not get taxonomy for ID {taxon_id}: {e}")
        pass  # Silently return NA for IDs not found or invalid
    except Exception as e:
        logging.warning(f"Unexpected error looking up taxon ID {taxon_id}: {e}")
    return result


def add_taxonomy_columns(df, taxdb, id_column="organism_id"):
    """Add phylum, order, family, genus, and species columns based on the organism ID."""
    if id_column not in df.columns:
        logging.error(f"Column '{id_column}' not found for taxonomy lookup.")
        df["phylum"] = pd.NA
        df["order"] = pd.NA
        df["family"] = pd.NA
        df["genus"] = pd.NA
        df["species"] = pd.NA
        return df

    logging.info(
        f"Fetching taxonomy (phylum, order, family, genus, species) for unique IDs in '{id_column}'."
    )
    unique_ids = df[id_column].dropna().unique()
    logging.info(f"Found {len(unique_ids)} unique organism IDs to process.")

    # Build cache for unique IDs
    taxonomy_cache = {}
    processed_count = 0
    for taxon_id in unique_ids:
        taxonomy_cache[taxon_id] = get_taxonomy_levels(taxon_id, taxdb)
        processed_count += 1
        if processed_count % 500 == 0 or processed_count == len(unique_ids):
            logging.info(
                f"  Processed {processed_count}/{len(unique_ids)} unique IDs for taxonomy cache..."
            )

    logging.info("Mapping taxonomy information to DataFrame.")
    # Map cached results
    df["phylum"] = df[id_column].map(
        lambda x: taxonomy_cache.get(x, {}).get("phylum", pd.NA)
    )
    df["order"] = df[id_column].map(
        lambda x: taxonomy_cache.get(x, {}).get("order", pd.NA)
    )
    df["family"] = df[id_column].map(
        lambda x: taxonomy_cache.get(x, {}).get("family", pd.NA)
    )
    df["genus"] = df[id_column].map(
        lambda x: taxonomy_cache.get(x, {}).get("genus", pd.NA)
    )
    df["species"] = df[id_column].map(
        lambda x: taxonomy_cache.get(x, {}).get("species", pd.NA)
    )

    logging.info("Taxonomy columns (phylum, order, family, genus, species) added.")
    return df


# --- Data Loading and Processing Functions ---


def load_transmembrane_data(file_path, columns_to_use):
    """Loads the transmembrane data CSV, keeping specified columns."""
    logging.info(f"Reading transmembrane data from: {file_path}")
    try:
        df = pd.read_csv(file_path, usecols=columns_to_use)
        logging.info(f"Successfully read {len(df)} rows from transmembrane data.")
        return df
    except FileNotFoundError:
        logging.error(f"Error: Transmembrane data file not found at {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading transmembrane data CSV: {e}")
        sys.exit(1)


def extract_uniprot_id_from_name(row):
    """Helper function to extract UniProt ID from protein_name field."""
    if pd.isna(row["uniprot_id"]):
        if isinstance(row["protein_name"], str) and "_" in row["protein_name"]:
            return row["protein_name"].split("_")[0]
    return row["uniprot_id"]


def process_transmembrane_ids(df):
    """Processes the DataFrame to fill missing uniprot_ids and clean up."""
    logging.info("Processing UniProt IDs (filling missing from protein_name)...")
    df["uniprot_id"] = df.apply(extract_uniprot_id_from_name, axis=1)

    initial_rows = len(df)
    df.dropna(subset=["uniprot_id"], inplace=True)
    dropped_rows = initial_rows - len(df)
    if dropped_rows > 0:
        logging.info(
            f"Dropped {dropped_rows} rows with missing UniProt IDs after processing."
        )

    # Drop the temporary protein_name column
    df = df.drop(columns=["protein_name"])
    logging.info("UniProt ID processing complete.")
    return df


def load_interpro_data(file_path):
    """Loads the InterPro/Pfam TSV data."""
    logging.info(f"Reading InterPro/Pfam data from: {file_path}")
    try:
        df = pd.read_csv(file_path, delimiter="\t", engine="python")
        logging.info(f"Successfully read {len(df)} rows from InterPro/Pfam data.")
        return df
    except FileNotFoundError:
        logging.error(f"Error: InterPro/Pfam data file not found at {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading InterPro/Pfam TSV file: {e}")
        sys.exit(1)


def merge_data(df_transmembrane, df_interpro):
    """Merges the two dataframes on UniProt ID."""
    logging.info("Merging transmembrane data with InterPro/Pfam data...")
    try:
        df_merged = pd.merge(
            df_transmembrane,
            df_interpro,
            left_on="uniprot_id",
            right_on="Entry",
            how="inner",
        )
        logging.info(f"Merged dataframe has {len(df_merged)} rows.")
        return df_merged
    except Exception as e:
        logging.error(f"Error during dataframe merge: {e}")
        sys.exit(1)


def load_priority_interpro_mapping(json_path: Path) -> Dict[str, str]:
    """Loads the InterPro mapping from JSON and returns a {id: short_name} dict."""
    logging.info(f"Loading InterPro priority mapping from: {json_path}")
    mapping = {}
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        for entry in data:
            if "id" in entry and "short_name" in entry:
                mapping[entry["id"]] = entry["short_name"]
        logging.info(f"Loaded {len(mapping)} priority InterPro ID mappings.")
        return mapping
    except FileNotFoundError:
        logging.error(f"Error: InterPro mapping file not found at {json_path}")
        sys.exit(1)
    except json.JSONDecodeError:
        logging.error(f"Error: InterPro mapping file {json_path} is not valid JSON.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading InterPro mapping file {json_path}: {e}")
        sys.exit(1)


def add_primary_interpro_info(
    df: pd.DataFrame, priority_order: List[str], mapping: Dict[str, str]
) -> pd.DataFrame:
    """Adds primary InterPro ID and short name based on priority order."""
    logging.info("Determining primary InterPro ID based on priority list...")

    if "interpro" not in df.columns:
        logging.warning(
            "Column 'interpro' not found. Skipping InterPro prioritization."
        )
        df["primary_interpro_id"] = pd.NA
        df["primary_interpro_short_name"] = pd.NA
        return df

    # Helper function to find the highest priority ID in a ; separated string
    def find_highest_priority_id(interpro_string: str) -> str:
        if pd.isna(interpro_string) or not isinstance(interpro_string, str):
            return pd.NA
        # Split, strip whitespace, and remove empty strings
        protein_ids = {ip.strip() for ip in interpro_string.split(";") if ip.strip()}
        for priority_id in priority_order:
            if priority_id in protein_ids:
                return priority_id
        return pd.NA  # Return NA if no priority ID is found

    df["primary_interpro_id"] = df["interpro"].apply(find_highest_priority_id)

    # Add short name using the mapping
    df["primary_interpro_short_name"] = (
        df["primary_interpro_id"].map(mapping).fillna("Unknown")
    )  # Map or mark as Unknown

    # Report how many were successfully mapped
    mapped_count = df["primary_interpro_id"].notna().sum()
    total_count = len(df)
    logging.info(
        f"Assigned primary InterPro ID to {mapped_count}/{total_count} entries."
    )
    if mapped_count < total_count:
        logging.warning(
            f"{total_count - mapped_count} entries had no match in the priority list."
        )

    # Drop the original interpro column
    if "interpro" in df.columns:
        logging.info("Dropping original 'interpro' column.")
        df = df.drop(columns=["interpro"])

    return df


def clean_and_rename_columns(df):
    """Drops redundant columns and renames others as specified."""
    logging.info("Cleaning and renaming columns...")

    # Columns to drop before final save
    columns_to_drop = [
        "Entry",
        # "organism_id", # Keep organism_id until after taxonomy lookup
    ]
    existing_cols_to_drop = [col for col in columns_to_drop if col in df.columns]
    if existing_cols_to_drop:
        df = df.drop(columns=existing_cols_to_drop)
        logging.info(f"Dropped columns: {existing_cols_to_drop}")
    # Don't warn if Entry is missing, it might have been dropped earlier if needed

    # Rename columns (headers only)
    rename_mapping = {
        "uniprot_id": "identifier",
        "Organism (ID)": "organism_id",  # Rename first, drop later
        "InterPro": "interpro",  # Rename to lowercase for subsequent processing
        "Pfam": "pfam",
    }
    # Check which columns actually exist before renaming
    actual_rename = {k: v for k, v in rename_mapping.items() if k in df.columns}
    df = df.rename(columns=actual_rename)
    missing_renames = set(rename_mapping.keys()) - set(actual_rename.keys())
    if missing_renames:
        logging.warning(f"Original columns not found for renaming: {missing_renames}")
    logging.info(f"Renamed columns: {actual_rename}")

    return df


def save_data(df, output_path):
    """Saves the final DataFrame to a CSV file."""
    logging.info(f"Saving final processed data to: {output_path}")
    try:
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        logging.info("Data saved successfully.")
    except Exception as e:
        logging.error(f"Error saving final data: {e}")
        sys.exit(1)


# --- Main Workflow ---


def main():
    """Main function to execute the data processing workflow."""
    logging.info("Starting data processing script...")

    # 0. Load InterPro Mapping Data
    interpro_map = load_priority_interpro_mapping(INTERPRO_MAPPING_JSON)

    # 1. Load Main Data
    df_transmembrane = load_transmembrane_data(INPUT_CSV_PATH, INITIAL_COLUMNS_TO_KEEP)
    df_interpro = load_interpro_data(INPUT_TSV_PATH)

    # 2. Process Transmembrane IDs
    df_transmembrane_processed = process_transmembrane_ids(df_transmembrane)

    # 3. Merge Data
    df_merged = merge_data(df_transmembrane_processed, df_interpro)

    # 4. Clean and Rename Columns (Initial pass - renames Organism ID, Pfam)
    df_cleaned = clean_and_rename_columns(df_merged)

    # 5. Add Primary InterPro Information
    df_interpro_processed = add_primary_interpro_info(
        df_cleaned, INTERPRO_PRIORITY_ORDER, interpro_map
    )

    # 6. Process Taxonomy
    try:
        taxdb = initialize_taxopy_db()
        df_with_taxonomy = add_taxonomy_columns(
            df_interpro_processed, taxdb, id_column="organism_id"
        )
    except Exception as e:
        logging.error(
            f"Taxonomy processing failed: {e}. Proceeding without taxonomy data."
        )
        # Add placeholder columns if taxonomy failed
        df_with_taxonomy = df_interpro_processed.copy()
        df_with_taxonomy["phylum"] = pd.NA
        df_with_taxonomy["order"] = pd.NA
        df_with_taxonomy["family"] = pd.NA
        df_with_taxonomy["genus"] = pd.NA
        df_with_taxonomy["species"] = pd.NA

    # 7. Save Summary Data (with sequences and full taxonomy)
    # Ensure 'Sequence' column is present if it was loaded.
    # The other new taxonomy columns ('family', 'genus', 'species') are already there.
    save_data(df_with_taxonomy, SUMMARY_OUTPUT_CSV_PATH)
    logging.info(f"Summary data saved to {SUMMARY_OUTPUT_CSV_PATH}")

    # 8. Prepare data for the original metadata.csv (protospace)
    # Drop 'Sequence' and extra taxonomy columns before saving to original path
    df_for_protospace = df_with_taxonomy.copy()
    if "Sequence" in df_for_protospace.columns:
        df_for_protospace = df_for_protospace.drop(columns=["Sequence"])
        logging.info("Dropped 'Sequence' column for protospace metadata.")

    columns_to_drop_for_protospace = ["family", "genus", "species"]
    existing_cols_to_drop_for_protospace = [
        col
        for col in columns_to_drop_for_protospace
        if col in df_for_protospace.columns
    ]
    if existing_cols_to_drop_for_protospace:
        df_for_protospace = df_for_protospace.drop(
            columns=existing_cols_to_drop_for_protospace
        )
        logging.info(
            f"Dropped extra taxonomy columns for protospace metadata: {existing_cols_to_drop_for_protospace}"
        )

    # 9. Final Column Cleanup for protospace metadata (Drop organism_id after use)
    if "organism_id" in df_for_protospace.columns:
        logging.info(
            "Dropping 'organism_id' column after taxonomy lookup for protospace metadata."
        )
        df_final_protospace = df_for_protospace.drop(columns=["organism_id"])
    else:
        logging.warning(
            "'organism_id' column not found for final drop in protospace metadata."
        )
        df_final_protospace = df_for_protospace

    # 10. Save Protospace Data
    save_data(df_final_protospace, OUTPUT_CSV_PATH)

    logging.info("Script finished successfully.")


if __name__ == "__main__":
    main()
