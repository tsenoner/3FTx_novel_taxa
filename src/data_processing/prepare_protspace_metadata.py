import pandas as pd
import taxopy
from pathlib import Path
import sys
import logging
import re  # For parsing Centipede headers
from typing import Dict, List, Any

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Configuration ---
DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw"
INTERM_DIR = DATA_DIR / "interm"

SANITIZED_FASTA_PATH = INTERM_DIR / "merged" / "merged_sanitized.fasta"
UNIPROT_TSV_PATH = RAW_DIR / "uniprot" / "3FTx_related.tsv"
TMBED_CSV_PATH = RAW_DIR / "protspace" / "tmbed_predictions.csv"
INTERPRO_MAPPING_TSV_PATH = RAW_DIR / "protspace" / "interpro_shortName_map.tsv"
MANUSCRIPT_GROUP_CSV_PATH = (
    RAW_DIR / "3FTX_manuscript" / "SD1_Dataset_and_information_updated.csv"
)

PROTOSPACE_METADATA_CSV_PATH = INTERM_DIR / "protspace" / "protspace_metadata.csv"

TMBED_COLUMNS = [
    "uniprot_id",
    "is_transmembrane",
    "is_helix",
    "is_beta",
    "has_signal",
]

INTERPRO_PRIORITY_ORDER = [
    "IPR003571",
    "IPR031424",
    "IPR026110",
    "IPR026524",
    "IPR038773",
    "IPR039457",
    "IPR039237",
    "IPR042339",
    "IPR051445",
    "IPR016054",
    "IPR035076",
    "IPR045860",
]

CENTIPEDE_HEADER_PREFIX = "Centipede3FTx"
CENTIPEDE_FUNCTIONAL_GROUP = "annotated"

SCALOPTOXIN_ACCESSIONS = {
    "I6R1R5",
    "P0DPX5",
    "P0DPX9",
    "P0DPY0",
    "P0DPY1",
    "P0DPX7",
    "P0DPX8",
    "P0DPX6",
    "P0DPU8",
}
SCALOPTOXIN_FUNCTIONAL_GROUP = "scoloptoxin"


# --- Taxonomy Helper Functions ---
def setup_taxopy_db_paths():
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
    db_dir, nodes_file, names_file, merged_file = setup_taxopy_db_paths()
    taxdb_params = {}
    use_dir_init = False
    if nodes_file.exists() and names_file.exists():
        taxdb_params["nodes_dmp"] = str(nodes_file)
        taxdb_params["names_dmp"] = str(names_file)
        if merged_file.exists() and merged_file.is_file():
            taxdb_params["merged_dmp"] = str(merged_file)
    else:
        use_dir_init = True
    try:
        if not use_dir_init and len(taxdb_params) >= 2:
            taxdb = taxopy.TaxDb(**taxdb_params)
        else:
            taxdb = taxopy.TaxDb(taxdb_dir=str(db_dir), keep_files=True)
        logging.info("Taxopy database initialized successfully.")
        return taxdb
    except Exception as e:
        logging.error(f"Failed to initialize taxopy database: {e}")
        raise


def get_taxonomy_levels(taxon_id: Any, taxdb: taxopy.TaxDb) -> Dict[str, Any]:
    result = {
        "phylum": pd.NA,
        "order": pd.NA,
        "family": pd.NA,
        "genus": pd.NA,
        "species": pd.NA,
    }
    if pd.isna(taxon_id):
        return result
    try:
        taxon_id_int = int(
            str(taxon_id).split(".")[0]
        )  # Handle potential float IDs if any
        taxon = taxopy.Taxon(taxon_id_int, taxdb)
        ranks = taxon.rank_name_dictionary
        for rank_name in result.keys():  # Iterate through phylum, order, etc.
            result[rank_name] = ranks.get(rank_name, pd.NA)
    except (ValueError, taxopy.exceptions.TaxidError, TypeError):
        logging.debug(
            f"Could not get taxonomy for ID {taxon_id}. Invalid format or not found."
        )
    except Exception as e:
        logging.warning(f"Unexpected error looking up taxon ID {taxon_id}: {e}")
    return result


def add_taxonomy_columns(
    df: pd.DataFrame, taxdb: taxopy.TaxDb, id_column: str = "organism_id_for_taxonomy"
) -> pd.DataFrame:
    taxonomy_levels_to_fetch = ["phylum", "order", "family", "genus", "species"]
    if id_column not in df.columns:
        logging.error(f"Column '{id_column}' not found for taxonomy lookup.")
        for col in taxonomy_levels_to_fetch:
            df[col] = pd.NA
        return df

    logging.info(f"Fetching taxonomy for unique IDs in '{id_column}'.")
    unique_ids = (
        df[id_column].dropna().astype(str).unique()
    )  # Ensure IDs are strings for caching keys
    logging.info(
        f"Found {len(unique_ids)} unique organism IDs to process for taxonomy."
    )

    taxonomy_cache = {}
    for uid_str in unique_ids:
        taxonomy_cache[uid_str] = get_taxonomy_levels(uid_str, taxdb)

    for col_level in taxonomy_levels_to_fetch:
        df[col_level] = (
            df[id_column]
            .astype(str)
            .map(
                lambda x_str: taxonomy_cache.get(x_str, {}).get(col_level, pd.NA)
                if x_str != "nan"
                else pd.NA
            )
        )

    logging.info(f"Taxonomy columns ({', '.join(taxonomy_levels_to_fetch)}) added.")
    return df


# --- Data Loading and Processing Functions ---
def parse_fasta_headers_and_centipede_info(fasta_path: Path) -> pd.DataFrame:
    logging.info(f"Parsing FASTA headers from: {fasta_path}")
    identifiers = []
    is_centipede_annotated_list = []
    organism_ids_from_fasta = []  # Store extracted organism IDs for centipedes

    # Regex to capture OrganismID from Centipede3FTx_NNNN_OrganismID format
    # Assumes OrganismID does not contain underscores.
    # If OrganismID can have underscores, a more specific NNNN pattern or suffix matching is needed.
    centipede_header_pattern = re.compile(rf"^{CENTIPEDE_HEADER_PREFIX}_\d+_(.+)$")

    try:
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    identifiers.append(header)

                    match = centipede_header_pattern.match(header)
                    if header.startswith(CENTIPEDE_HEADER_PREFIX) and match:
                        is_centipede_annotated_list.append(True)
                        organism_ids_from_fasta.append(
                            match.group(1)
                        )  # Captured OrganismID
                    else:
                        is_centipede_annotated_list.append(False)
                        organism_ids_from_fasta.append(pd.NA)

        df = pd.DataFrame(
            {
                "identifier": identifiers,
                "is_centipede_annotated": is_centipede_annotated_list,
                "organism_id_from_fasta": organism_ids_from_fasta,
            }
        )
        logging.info(
            f"Parsed {len(df)} headers from {fasta_path}. Found {df['is_centipede_annotated'].sum()} potential centipede entries."
        )
        return df
    except FileNotFoundError:
        logging.error(f"Sanitized FASTA file not found at {fasta_path}")
        raise FileNotFoundError(f"Sanitized FASTA file not found at {fasta_path}")
    except Exception as e:
        logging.error(f"Error parsing FASTA file {fasta_path}: {e}")
        raise RuntimeError(f"Error parsing FASTA file {fasta_path}: {e}")


def load_and_preprocess_uniprot_tsv(file_path: Path) -> pd.DataFrame:
    logging.info(f"Reading UniProt TSV data from: {file_path}")
    try:
        df = pd.read_csv(file_path, delimiter="\t", engine="python")
        id_col_original = "Entry"  # Expected UniProt ID column
        if id_col_original not in df.columns:
            logging.error(
                f"Expected UniProt ID column '{id_col_original}' not found in {file_path}."
            )
            # Allow to proceed but merges might be empty for uniprot specific data
            df["identifier"] = pd.NA
        else:
            df = df.rename(columns={id_col_original: "identifier"})

        # Standardize Organism ID and other columns, filling with NA if missing
        rename_map = {
            "Organism (ID)": "organism_id",
            "Protein names": "protein_names",
            "InterPro": "interpro",
            "Pfam": "pfam",
        }
        for old_name, new_name in rename_map.items():
            if old_name in df.columns:
                df = df.rename(columns={old_name: new_name})
            elif (
                new_name != "identifier"
            ):  # avoid re-adding identifier if it was missing
                df[new_name] = pd.NA
                logging.warning(
                    f"Column '{old_name}' (for '{new_name}') not found in {file_path}. Filled with NA."
                )

        # Ensure all target columns exist even if original ones were missing
        final_uniprot_cols = [
            "identifier",
            "protein_names",
            "organism_id",
            "interpro",
            "pfam",
        ]
        for col in final_uniprot_cols:
            if col not in df.columns:
                df[col] = pd.NA

        return df[final_uniprot_cols].copy()
    except FileNotFoundError:
        logging.error(f"UniProt TSV file not found at {file_path}")
        raise FileNotFoundError(f"UniProt TSV file not found at {file_path}")
    except Exception as e:
        logging.error(f"Error reading UniProt TSV {file_path}: {e}")
        raise RuntimeError(f"Error reading UniProt TSV {file_path}: {e}")


def load_tmbed_data(file_path: Path, columns_to_use: List[str]) -> pd.DataFrame:
    logging.info(f"Reading TMBed data from: {file_path}")
    try:
        df = pd.read_csv(file_path)
        if "uniprot_id" in df.columns:
            df = df.rename(columns={"uniprot_id": "identifier"})
        elif df.columns[0] == "protein_name" and "identifier" not in df.columns:
            df["identifier"] = df["protein_name"].apply(
                lambda x: x.split("_")[0] if isinstance(x, str) and "_" in x else x
            )
        elif "identifier" not in df.columns:
            logging.warning(
                f"No 'uniprot_id' or 'identifier' in TMBed data {file_path}. Using first column."
            )
            df = df.rename(columns={df.columns[0]: "identifier"})

        final_cols = ["identifier"]
        for col in columns_to_use:  # TMBED_COLUMNS from config
            if col == "uniprot_id":
                continue
            if col in df.columns:
                final_cols.append(col)
            else:
                df[col] = pd.NA
                final_cols.append(col)
                logging.warning(f"TMBed column '{col}' not found. Added as NA.")
        return df[final_cols].copy()
    except FileNotFoundError:
        logging.info(
            f"TMBed data file not found at {file_path}. Transmembrane info will be missing."
        )
        # Create an empty df with expected columns so merge doesn't fail
        return pd.DataFrame(
            columns=["identifier"] + [c for c in columns_to_use if c != "uniprot_id"]
        )
    except Exception as e:
        logging.error(f"Error reading TMBed CSV {file_path}: {e}")
        return pd.DataFrame(
            columns=["identifier"] + [c for c in columns_to_use if c != "uniprot_id"]
        )


def load_interpro_shortname_mapping(tsv_path: Path) -> Dict[str, str]:
    logging.info(
        f"Loading InterPro short name mapping from: {tsv_path} (expecting no header)"
    )
    mapping = {}
    try:
        df_map = pd.read_csv(tsv_path, sep="\t", header=None)
        if df_map.shape[1] >= 2:
            df_map.columns = [0, 1] + list(df_map.columns[2:])
            mapping = dict(zip(df_map[0], df_map[1]))
            logging.info(
                f"Loaded {len(mapping)} InterPro ID to short name mappings from {tsv_path}."
            )
        else:
            logging.error(
                f"InterPro mapping file {tsv_path} needs at least two columns."
            )
            # Return empty map instead of sys.exit to allow script to proceed if non-critical
            return {}
        return mapping
    except FileNotFoundError:
        logging.warning(
            f"InterPro mapping TSV file not found at {tsv_path}. Short names will be 'Unknown'."
        )
        return {}
    except pd.errors.EmptyDataError:
        logging.warning(
            f"InterPro mapping TSV file {tsv_path} is empty. Short names will be 'Unknown'."
        )
        return {}
    except Exception as e:
        logging.error(f"Error reading InterPro mapping TSV {tsv_path}: {e}")
        return {}


def load_manuscript_group_data(file_path: Path) -> pd.DataFrame:
    logging.info(f"Reading manuscript group data from: {file_path}")
    try:
        df = pd.read_csv(file_path)
        df["identifier"] = df["identifier"].str.extract(r"\|(.*?)\|")
        df["group_modified"] = df.get("group_modified", pd.NA)
        return df[["identifier", "group_modified"]].copy()
    except FileNotFoundError:
        logging.warning(
            f"Manuscript group data file not found at {file_path}. Group info will be missing."
        )
        return pd.DataFrame(columns=["identifier", "group_modified"])
    except Exception as e:
        logging.error(f"Error reading manuscript group CSV {file_path}: {e}")
        return pd.DataFrame(columns=["identifier", "group_modified"])


def add_primary_interpro_info(
    df: pd.DataFrame, priority_order: List[str], mapping: Dict[str, str]
) -> pd.DataFrame:
    logging.info("Determining primary InterPro ID based on priority list...")
    if "interpro" not in df.columns or df["interpro"].isnull().all():
        logging.warning(
            "Column 'interpro' not found or all NA. Skipping InterPro prioritization."
        )
        df["primary_interpro_id"] = pd.NA
        df["primary_interpro_short_name"] = "Unknown"
        return df

    def find_highest_priority_id(interpro_string: Any) -> Any:
        if pd.isna(interpro_string) or not isinstance(interpro_string, str):
            return pd.NA
        protein_ids = {ip.strip() for ip in interpro_string.split(";") if ip.strip()}
        for priority_id in priority_order:
            if priority_id in protein_ids:
                return priority_id
        return pd.NA

    df["primary_interpro_id"] = df["interpro"].apply(find_highest_priority_id)
    df["primary_interpro_short_name"] = (
        df["primary_interpro_id"].map(mapping).fillna("Unknown")
    )
    mapped_count = df["primary_interpro_id"].notna().sum()
    logging.info(f"Assigned primary InterPro ID to {mapped_count}/{len(df)} entries.")
    return df


def save_dataframe(df: pd.DataFrame, output_path: Path):
    logging.info(f"Saving final processed data to: {output_path}")
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        logging.info(f"Data saved successfully to {output_path}. Shape: {df.shape}")
    except Exception as e:
        logging.error(f"Error saving final data to {output_path}: {e}")
        raise RuntimeError(f"Error saving final data to {output_path}: {e}")


# --- Main Workflow ---
def main():
    logging.info("Starting Protospace metadata preparation script...")

    # 1. Load base identifiers and centipede info from sanitized FASTA
    df_base = parse_fasta_headers_and_centipede_info(SANITIZED_FASTA_PATH)
    if df_base.empty:
        logging.error("Failed to load base identifiers from FASTA. Exiting.")
        sys.exit(1)

    # 2. Load supplementary data sources
    df_uniprot = load_and_preprocess_uniprot_tsv(UNIPROT_TSV_PATH)
    df_tmbed = load_tmbed_data(TMBED_CSV_PATH, TMBED_COLUMNS)
    interpro_shortname_map = load_interpro_shortname_mapping(INTERPRO_MAPPING_TSV_PATH)
    df_manuscript_group = load_manuscript_group_data(MANUSCRIPT_GROUP_CSV_PATH)

    # 3. Merge Data: Start with base FASTA identifiers, left merge other data
    logging.info("Merging datasets with base FASTA identifiers...")
    df_merged = pd.merge(df_base, df_uniprot, on="identifier", how="left")
    df_merged = pd.merge(df_merged, df_tmbed, on="identifier", how="left")
    df_merged = pd.merge(df_merged, df_manuscript_group, on="identifier", how="left")
    logging.info(f"Shape after merges: {df_merged.shape}")

    # 4. Prepare organism ID for taxonomy lookup
    # Prioritize organism_id from FASTA for centipedes, fallback to UniProt's organism_id
    df_merged["organism_id_for_taxonomy"] = df_merged["organism_id_from_fasta"].fillna(
        df_merged["organism_id"]
    )

    # 5. Add Taxonomy (Phylum, Order, Family, Genus, Species)
    try:
        taxdb = initialize_taxopy_db()
        df_processed = add_taxonomy_columns(
            df_merged, taxdb, id_column="organism_id_for_taxonomy"
        )
    except Exception as e:
        logging.error(
            f"Taxonomy processing failed: {e}. Proceeding with NA for taxonomy columns."
        )
        df_processed = df_merged.copy()  # Ensure df_processed exists
        for col in ["phylum", "order", "family", "genus", "species"]:
            df_processed[col] = pd.NA

    # 6. Add Primary InterPro Information
    df_processed = add_primary_interpro_info(
        df_processed, INTERPRO_PRIORITY_ORDER, interpro_shortname_map
    )

    # 7. Annotate Functional Group and Protein Display Name
    logging.info("Applying functional group and protein display name annotations...")

    # Initialize functional_group from manuscript's group_modified
    df_processed["functional_group"] = df_processed["group_modified"]

    # Set general 'annotated' group for headers starting with CENTIPEDE_HEADER_PREFIX
    df_processed.loc[df_processed["is_centipede_annotated"], "functional_group"] = (
        CENTIPEDE_FUNCTIONAL_GROUP
    )

    # Override with 'scoloptoxin' for specific SCALOPTOXIN_ACCESSIONS
    df_processed.loc[
        df_processed["identifier"].isin(SCALOPTOXIN_ACCESSIONS), "functional_group"
    ] = SCALOPTOXIN_FUNCTIONAL_GROUP

    # Initialize protein_display_name from UniProt protein_names
    df_processed["protein_display_name"] = df_processed["protein_names"].apply(
        lambda x: x.split(";")[0].strip()
        if pd.notna(x) and isinstance(x, str) and ";" in x
        else (x if pd.notna(x) else pd.NA)
    )
    # For CENTIPEDE_HEADER_PREFIX entries, use their FASTA identifier as protein_display_name
    df_processed.loc[df_processed["is_centipede_annotated"], "protein_display_name"] = (
        df_processed["identifier"]
    )

    # 8. Define and Select Final Columns
    final_columns = [
        "identifier",
        "protein_display_name",
        "functional_group",
        "is_transmembrane",
        "is_helix",
        "is_beta",
        "has_signal",
        "pfam",
        "primary_interpro_id",
        "primary_interpro_short_name",
        "phylum",
        "order",
        "family",
        "genus",
        "species",
        "organism_id",  # This is the original organism_id from UniProt, kept for reference
        # organism_id_from_fasta and organism_id_for_taxonomy can be dropped if not needed in final output
    ]

    # Ensure all final columns exist, fill with NA if any are missing
    for col in final_columns:
        if col not in df_processed.columns:
            logging.warning(
                f"Final column '{col}' was not generated/retained. Filling with NA."
            )
            df_processed[col] = pd.NA

    df_final_metadata = df_processed[final_columns].copy()
    df_final_metadata.dropna(
        subset=["identifier"], inplace=True
    )  # Should not drop if FASTA was source

    # 9. Save Data
    save_dataframe(df_final_metadata, PROTOSPACE_METADATA_CSV_PATH)

    logging.info("Protospace metadata preparation script finished successfully.")


if __name__ == "__main__":
    main()
