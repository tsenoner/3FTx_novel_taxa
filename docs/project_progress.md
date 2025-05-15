# Project Progress: 3FTx_novel_taxa

## Introduction

This document outlines the key steps undertaken in the `3FTx_novel_taxa` project. The primary goal of this research is to investigate the evolution and occurrence of Three-Finger Toxin (3FTx)-like proteins beyond snakes, with a particular focus on invertebrates like centipedes.

The project was spurred by the intriguing discovery of toxin protein sequences in centipedes that occupy a similar structural space to snake 3FTx in ProtT5 embeddings. This is significant because venomous snake 3FTx are understood to have evolved from the Ly6 protein family after the divergence of snakes and lizards. Centipedes, being invertebrate Arthropods, present a fascinating case: did their 3FTx-like proteins also evolve from Ly6 in a separate specialization event, or is there another evolutionary pathway? This project aims to explore these questions, seeking to understand when and how these proteins evolved in centipedes and whether similar sequences exist in other arthropods.

## Steps Completed

Alongside the specific data processing steps detailed below, several Python scripts within the `src/data_processing/` directory underwent significant refactoring. These efforts focused on simplifying the codebase, improving readability, and standardizing operations. Key changes included hardcoding file paths to streamline execution, removing command-line argument parsing (e.g., `argparse`) in favor of fixed configurations, simplifying console output and reporting mechanisms, and restructuring code for better maintainability (e.g., separating core logic from I/O operations). These refinements apply to scripts such as `process_centipede_fasta.py` and `merge_data.py`.

### 1. Sequence Gathering from UniProt

To compile a comprehensive dataset of known and potentially related sequences, the following actions were taken:

*   **InterPro ID Selection**: We identified and selected all InterPro IDs associated with the 3FTx family (e.g., `IPR003571`, `IPR045860`) and the related LY6/UPAR domain (`IPR016054`, `IPR035076`) and its associated families. These IDs are linked to the characteristic three-finger fold.
*   **Inclusion of Novel Centipede Proteins**: The newly discovered centipede protein sequences, which prompted this investigation, were incorporated into the search strategy.
*   **UniProt Query**: A detailed query was constructed to search the UniProt database. The query, as documented in `data/raw/uniprot/README.md`, is:
    ```
    (((xref:interpro-IPR016054 OR xref:interpro-IPR026110 OR xref:interpro-IPR026524 OR xref:interpro-IPR031424 OR xref:interpro-IPR038773 OR xref:interpro-IPR039457 OR xref:interpro-IPR039237 OR xref:interpro-IPR042339 OR xref:interpro-IPR051445 OR xref:interpro-IPR035076 OR xref:interpro-IPR003571 OR xref:interpro-IPR045860) AND fragment:false) OR (accession:I6R1R5 OR accession:P0DPX5 OR accession:P0DPX9 OR accession:P0DPY0 OR accession:P0DPY1 OR accession:P0DPX7 OR accession:P0DPX8 OR accession:P0DPX6 OR accession:P0DPU8)) AND (length:[51 TO 199]) AND (taxonomy_id:33208)
    ```
*   **Data Download**: The search results were downloaded from UniProt in several formats:
    *   FASTA file (for sequences)
    *   TSV file (with headers: `Entry`, `Organism (ID)`, `InterPro`, `Pfam`)
    *   ProtT5 per-protein embeddings
    *   JSON file (as a backup for additional information if needed later)

### 2. Centipede Genome Search and Embedding Generation

To validate the identified centipede sequences and discover any other similar genes within the centipede lineage:

*   **Genome Search**: The protein sequences found via the UniProt query were searched against available centipede genome data. This step aimed to confirm their authenticity and to extract all genes encoding similar protein sequences.
*   **ProtT5 Embeddings**: For the sequences extracted directly from the centipede genome analysis, ProtT5 embeddings were generated.

### 3. Sequence Merging and Cleaning

To create a refined dataset for downstream analysis:

*   **Data Consolidation and Cleaning**: Manually annotated centipede sequences (from the centipede genome search) were merged with the sequences obtained from the UniProt download using the `src/data_processing/merge_data.py` script. This script also performed quality control on the merged FASTA dataset by:
    *   Removing duplicate sequences.
    *   Eliminating sequences with ambiguous amino acid characters or internal stop codons (trailing stop codons were trimmed).
    *   Cleaning FASTA headers for consistency.

### 4. Transmembrane Domain Prediction

To identify potential transmembrane regions within the protein sequences, which can be relevant for understanding their function and localization:

*   **TMBed Prediction**: The TMBed tool was run on the curated sequences from `data/interm/merged/merged_sanitized.fasta` to predict transmembrane helices. The results are stored in `data/raw/protspace/tmbed_predictions.csv`.

### 5. Comprehensive Metadata Compilation for Protospace

To consolidate all relevant information for each protein sequence into a single, analyzable table for Protospace visualization and further research, a dedicated script `src/data_processing/prepare_protspace_metadata.py` was developed. This script performs the following key operations:

*   **Primary Input**: Uses `data/interm/merged/merged_sanitized.fasta` as the primary source of sequence identifiers.
*   **Data Integration**: It parses and merges data from multiple sources:
    *   FASTA headers from `merged_sanitized.fasta`, with special parsing for sequences prefixed with `Centipede3FTx` to extract `OrganismID` and flag them as centipede-annotated.
    *   UniProt metadata from `data/raw/uniprot/3FTx_related.tsv` (providing protein names, UniProt organism ID, InterPro annotations, Pfam annotations).
    *   Transmembrane predictions from `data/raw/protspace/tmbed_predictions.csv`.
    *   InterPro short name mappings from `data/raw/protspace/interpro_shortName_map.tsv` to get descriptive names for InterPro IDs.
    *   Functional group information from the manuscript's supplementary data (`data/raw/3FTX_manuscript/SD1_Dataset_and_information_updated.csv`).
*   **Taxonomy Assignment**: Adds detailed taxonomic information (phylum, order, family, genus, species) for each sequence by querying the NCBI taxonomy database (via `taxopy`) using the `OrganismID`. For `Centipede3FTx` sequences, the `OrganismID` extracted from the FASTA header is prioritized; otherwise, the UniProt organism ID is used.
*   **InterPro Prioritization**: Determines a `primary_interpro_id` and `primary_interpro_short_name` for each sequence by selecting the highest-priority InterPro ID from a predefined list if multiple are present.
*   **Functional Group Annotation**: Assigns a `functional_group` to each protein:
    *   Initially populated from the manuscript's `group_modified` column.
    *   Sequences with FASTA headers starting with `Centipede3FTx` are assigned the `functional_group` "annotated".
    *   A specific set of proteins, identified by their accessions (defined as `SCALOPTOXIN_ACCESSIONS` in the script), are assigned the `functional_group` "scoloptoxin".
*   **Protein Display Name Generation**: Creates a `protein_display_name` for each entry:
    *   Typically derived from UniProt protein names.
    *   For sequences with headers starting with `Centipede3FTx`, their full FASTA identifier is used as the display name.
*   **Output**: The script generates a comprehensive metadata table saved as `data/interm/protspace/protspace_metadata.csv`.

### 6. ProtT5 Embedding Consolidation and Dimensionality Reduction

With the curated sequences and comprehensive metadata, the next step involves working with protein embeddings for structural and functional analysis:

*   **Embedding Association**: The pre-computed ProtT5 per-protein embeddings (gathered during initial UniProt download or generated for new centipede sequences) are associated with the entries in `protspace_metadata.csv`.
*   **Dimensionality Reduction**: To visualize and analyze the high-dimensional ProtT5 embeddings, various dimensionality reduction techniques are applied. These include:
    *   UMAP (Uniform Manifold Approximation and Projection)
    *   PCA (Principal Component Analysis)
    *   PacMAP (Pairwise Controlled Manifold Approximation)
    These techniques help in projecting the embeddings into lower-dimensional spaces, facilitating the identification of clusters and relationships between sequences.

### 8. Interactive Visualization with Protospace

To interactively explore the sequence landscape based on their ProtT5 embeddings and associated metadata, the Protospace visualization tool is utilized. This involves two main steps:

*   **Generating Protospace JSON Data**: The `protspace-json` command-line tool is used to process the ProtT5 embeddings (from `data/interm/merged/merged_sanitized.h5`) and the comprehensive metadata file (`data/interm/protspace/protspace_metadata.csv`). It performs dimensionality reduction using specified methods (PCA, UMAP, PacMAP with 2 components each) and integrates the metadata to create a single JSON file (`data/interm/protspace/protspace.json`) required by the Protospace viewer. Key parameters for UMAP/PacMAP like `n_neighbors` (100) and `min_dist` (0.5) are specified to control the layout of the manifold.

    ```bash
    protspace-json -i data/interm/merged/merged_sanitized.h5 -m data/interm/protspace/protspace_metadata.csv -o data/interm/protspace/protspace.json --methods pca2 umap2 pacmap2 --n_neighbors 100 --min_dist 0.5
    ```

*   **Applying Feature Styling**: To enhance the visualization with custom color schemes and styles for different metadata features, the `protspace-feature-colors` tool is employed. This tool takes a predefined styling configuration JSON (`data/interm/protspace/styling.json`) and applies it to the `protspace.json` file, producing a styled JSON file (`data/processed/protspace_styled.json`) that can then be loaded into the Protospace web viewer for interactive analysis.

    ```bash
    protspace-feature-colors --feature_styles data/interm/protspace/styling.json data/interm/protspace/protspace.json data/processed/protspace_styled.json
    ```

This structured approach aims to build a robust dataset and analytical framework to investigate the evolutionary origins and characteristics of 3FTx-like proteins in novel taxa.

### 9. Data Preparation for Phylogenetic Analysis (ExaBayes)

To prepare the curated sequences for phylogenetic tree construction using ExaBayes, the following steps were performed:

*   **Multiple Sequence Alignment (MSA)**: The merged and sanitized sequences were aligned using MAFFT. This step is crucial for identifying homologous regions across different sequences, which is a prerequisite for phylogenetic inference. The following command was used:
    ```bash
    mafft --auto --dpparttree --thread 10 data/interm/merged/merged_sanitized.fasta > data/interm/exabayes/merged_aligned_auto.fasta
    ```
    This command automatically selects the appropriate MAFFT strategy (`--auto`), uses a progressive method with a distance-based partition tree (`--dpparttree`), and utilizes 10 threads for parallel processing (`--thread 10`). The output is an aligned FASTA file stored at `data/interm/exabayes/merged_aligned_auto.fasta`.

*   **Conversion to PHYLIP Format**: The aligned FASTA file was then converted to the PHYLIP format, which is required by many phylogenetic software packages, including ExaBayes. This conversion was performed using the `src/data_processing/fasta2phy.py` script. This script takes the aligned FASTA file (`data/interm/exabayes/merged_aligned_auto.fasta`) as input and produces a corresponding `.phy` file.

This structured approach aims to build a robust dataset and analytical framework to investigate the evolutionary origins and characteristics of 3FTx-like proteins in novel taxa.