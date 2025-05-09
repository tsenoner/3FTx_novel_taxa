# Project Specification: 3FTx_novel_taxa

## 1. Overview

This document outlines the directory structure and organization of the `3FTx_novel_taxa` project. The project aims to analyse and identify possible 3FTx in other taxa outside of snakes as potentially found in centipedes and spiders.

## 2. Directory Structure

The project is organized into the following main directories:

```
.
├── .git/                   # Git version control files
├── .venv/                  # Python virtual environment
├── data/                   # Project data
│   ├── raw/                # Original, immutable input data
│   │   └── centipede_3ftx_quiver_upar_like/ # Specific dataset for centipede 3FTx Quiver UPAR-like sequences & annotations
│   ├── interm/             # Intermediate data files generated during processing
│   └── processed/          # Final, processed data ready for analysis or reporting
├── docs/                   # Project documentation
│   └── project_specification.md # This file
├── dir/                    # Currently empty directory
├── out/                    # Output files from scripts and analyses
│   ├── sunburst_taxonomy.html # Example output: Taxonomy visualization
│   └── species_uniprot_ids.txt  # Example output: List of species UniProt IDs
├── src/                    # Source code
│   ├── data_processing/    # Scripts for data cleaning, transformation, and preparation
│   │   └── process_metadata.py # Script for processing metadata
│   └── visualization/      # Scripts for generating visualizations
└── README.md               # Project overview, setup, and usage instructions
└── .gitattributes          # Git file attributes
└── .gitignore              # Specifies intentionally untracked files that Git should ignore
└── .python-version         # Specifies Python version for the project
└── pyproject.toml          # Python project configuration (PEP 518), used by uv
└── uv.lock                 # Dependency lock file for the uv resolver

```

## 3. Key Components

### 3.1. `data/` Directory

*   **`data/raw/`**: Stores all initial data. This data should be treated as read-only.
    *   **`data/raw/centipede_3ftx_quiver_upar_like/`**: Contains GFF annotation files and FASTA sequence files related to 3FTx Quiver UPAR-like proteins in centipedes.
*   **`data/interm/`**: Holds temporary files that are outputs of one processing step and inputs to another. Can be regenerated.
*   **`data/processed/`**: Contains the final datasets that are the result of the data processing pipeline. These are typically used for analysis, visualization, or reporting.

### 3.2. `src/` Directory

*   **`src/data_processing/`**: Contains Python scripts dedicated to processing and transforming data.
    *   `process_metadata.py`: A key script for handling metadata associated with the project.
*   **`src/visualization/`**: Contains scripts used to generate plots, charts, and other visual representations of the data and analysis results.

### 3.3. `out/` Directory

This directory is for storing non-data outputs, such as final reports, figures, or specific result files like `sunburst_taxonomy.html` and `species_uniprot_ids.txt`.

### 3.4. Configuration Files

*   **`pyproject.toml`**: Defines project metadata and build dependencies, used by `uv`.
*   **`uv.lock`**: Ensures reproducible Python environments by locking dependency versions with `uv`.
*   **`.python-version`**: Ensures consistent Python version usage across development environments.

## 4. Workflow (Assumed - USER: Please verify or update)

1.  **Data Collection**: Raw data is placed in `data/raw/`.
2.  **Data Processing**: Scripts in `src/data_processing/` are run to clean, transform, and prepare the data. Intermediate files are stored in `data/interm/`, and final processed data in `data/processed/`.
3.  **Analysis & Visualization**: Scripts in `src/visualization/` (and potentially other analysis scripts if added) use the processed data to generate insights and outputs, which may be saved in `out/`.

## 5. Future Considerations / TODO

*   Populate `README.md` with detailed project information, setup instructions, and how to run the code.
*   Decide on the purpose of the `dir/` directory or remove it if not needed.
*   [**USER: Add any other future plans or to-do items here**]