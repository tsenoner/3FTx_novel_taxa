import pandas as pd
import plotly.express as px
from pathlib import Path
import logging

# Standard logging setup
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# --- Configuration ---
DATA_DIR = Path("data")
INTERM_DIR = DATA_DIR / "interm"
INPUT_CSV_PATH = INTERM_DIR / "metadata_summary.csv"

OUT_DIR = Path("out")
OUTPUT_HTML_PATH = OUT_DIR / "sunburst_taxonomy.html"

# Define taxonomic hierarchy and root node for the sunburst plot
BASE_TAXONOMY_LEVELS = ["phylum", "order", "family", "genus", "species"]
SUNBURST_ROOT_LABEL = "All Sequences"
SUNBURST_PATH_LEVELS = ["sunburst_root"] + BASE_TAXONOMY_LEVELS


def load_data(file_path: Path, taxonomy_levels: list[str]) -> pd.DataFrame:
    """Loads summary data and ensures specified taxonomy columns are present and filled."""
    logging.info(f"Reading summary data from: {file_path}")
    try:
        df = pd.read_csv(file_path, low_memory=False)
        logging.info(f"Successfully read {len(df)} rows from {file_path}.")
        for col in taxonomy_levels:
            if col in df.columns:
                df[col] = df[col].fillna("Unknown")
            else:
                logging.warning(
                    f"Taxonomy column '{col}' not found. Created and filled with 'Unknown'."
                )
                df[col] = "Unknown"
        return df
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Error reading CSV {file_path}: {e}")
        raise


def create_sunburst_plot(
    df: pd.DataFrame,
    base_tax_levels: list[str],
    sunburst_path_levels: list[str],
    root_label: str,
    output_path: Path,
):
    """Generates and saves an interactive sunburst plot of taxonomic distribution."""
    logging.info("Creating sunburst plot...")

    # Ensure base taxonomy columns are present (safeguard after load_data)
    for level in base_tax_levels:
        if level not in df.columns:
            df[level] = "Unknown"
        else:
            df[level] = df[level].fillna("Unknown")

    # Add the root column for the sunburst hierarchy
    df["sunburst_root"] = root_label

    # Aggregate data for sunburst plot values
    df_counts = pd.DataFrame(columns=sunburst_path_levels + ["counts"])
    if not df.empty and all(level in df.columns for level in sunburst_path_levels):
        df_counts = (
            df.groupby(sunburst_path_levels, observed=True)
            .size()
            .reset_index(name="counts")
        )

    if df_counts.empty:
        logging.warning("No data to plot. Saving an empty plot message.")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(
                "<html><body><h1>No data available for sunburst plot.</h1></body></html>"
            )
        return

    try:
        fig = px.sunburst(
            df_counts,
            path=sunburst_path_levels,
            values="counts",
            title="Taxonomic Distribution of Sequences",
            maxdepth=3,  # Initial display depth: Root -> Phylum -> Order
        )
        fig.update_traces(textinfo="label+percent entry")

        # --- Interactive Phylum Filtering Buttons ---
        phylum_col_name = base_tax_levels[0]
        unique_phyla = []
        if phylum_col_name in df.columns and not df[phylum_col_name].empty:
            unique_phyla_series = df[phylum_col_name].unique()
            unique_phyla = sorted([p for p in unique_phyla_series if p != "Unknown"])
            if "Unknown" in unique_phyla_series:
                unique_phyla.append("Unknown")

        buttons = []

        # Extract necessary data from the figure for button interactions
        sector_ids_list = []
        raw_ids = (
            fig.data[0].ids
            if fig.data and hasattr(fig.data[0], "ids") and fig.data[0].ids is not None
            else None
        )
        if raw_ids is not None:
            if hasattr(raw_ids, "tolist"):
                sector_ids_list = raw_ids.tolist()
            elif isinstance(raw_ids, (list, tuple)):
                sector_ids_list = list(raw_ids)
            elif raw_ids:
                sector_ids_list = [raw_ids]

        original_fig_values_list = []
        if sector_ids_list:  # Only proceed if there are sectors to process
            raw_values = (
                fig.data[0].values
                if fig.data
                and hasattr(fig.data[0], "values")
                and fig.data[0].values is not None
                else None
            )
            temp_values = None
            if raw_values is not None:
                if hasattr(raw_values, "tolist"):
                    temp_values = raw_values.tolist()
                elif isinstance(raw_values, (list, tuple)):
                    temp_values = list(raw_values)
                elif raw_values:
                    temp_values = [raw_values]

            if temp_values is not None and len(temp_values) == len(sector_ids_list):
                original_fig_values_list = [
                    v if isinstance(v, (int, float)) else 0 for v in temp_values
                ]
            else:
                logging.warning(
                    "Mismatch or issue with figure values vs. sector IDs. Defaulting to zeros."
                )
                original_fig_values_list = [0] * len(sector_ids_list)

        # Calculate counts for button labels
        phylum_counts = {}
        if phylum_col_name in df_counts.columns:
            phylum_counts = (
                df_counts.groupby(phylum_col_name, observed=True)["counts"]
                .sum()
                .to_dict()
            )
        grand_total_count = df_counts["counts"].sum()

        # Button for "All Phyla"
        buttons.append(
            dict(
                args=[
                    {
                        "values": [
                            original_fig_values_list.copy()
                            if original_fig_values_list
                            else []
                        ],
                        "marker.colors": [None],
                    }
                ],
                label=f"All Phyla ({grand_total_count})",
                method="restyle",
            )
        )

        # Buttons for each phylum
        for phylum_name in unique_phyla:
            new_values_for_phylum = []
            for i, sector_id in enumerate(sector_ids_list):
                val = (
                    original_fig_values_list[i]
                    if i < len(original_fig_values_list)
                    else 0
                )
                is_root = sector_id == root_label
                is_selected_phylum = sector_id == f"{root_label}/{phylum_name}"
                is_descendant = isinstance(sector_id, str) and sector_id.startswith(
                    f"{root_label}/{phylum_name}/"
                )

                if is_root or is_selected_phylum or is_descendant:
                    new_values_for_phylum.append(val)
                else:
                    new_values_for_phylum.append(0)

            count = phylum_counts.get(phylum_name, 0)
            buttons.append(
                dict(
                    args=[{"values": [new_values_for_phylum]}],
                    label=f"{str(phylum_name)} ({count})",
                    method="restyle",
                )
            )

        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="down",
                    active=0,
                    buttons=buttons,
                    pad={"r": 10, "t": 10, "b": 10},
                    showactive=True,
                    x=1.02,
                    xanchor="left",
                    y=1,
                    yanchor="top",
                )
            ],
            margin=dict(t=50, l=25, r=180, b=25),
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        logging.info(f"Saving sunburst plot to: {output_path}")
        fig.write_html(output_path)
        logging.info("Sunburst plot saved successfully.")
    except Exception as e:
        logging.error(f"Error creating or saving sunburst plot: {e}", exc_info=True)
        raise


def main():
    """Main workflow to load data and generate the sunburst plot."""
    logging.info("Starting sunburst plot generation script...")
    OUT_DIR.mkdir(
        parents=True, exist_ok=True
    )  # Ensure output dir exists before any file ops
    try:
        df_summary = load_data(INPUT_CSV_PATH, BASE_TAXONOMY_LEVELS)
        create_sunburst_plot(
            df_summary,
            BASE_TAXONOMY_LEVELS,
            SUNBURST_PATH_LEVELS,
            SUNBURST_ROOT_LABEL,
            OUTPUT_HTML_PATH,
        )
    except Exception:
        logging.error("Script failed during execution.", exc_info=True)
        return
    logging.info("Sunburst plot script finished successfully.")


if __name__ == "__main__":
    main()
