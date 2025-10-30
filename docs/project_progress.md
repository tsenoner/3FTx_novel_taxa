# Project Progress: 3FTx_novel_taxa

## Introduction

This document outlines the key steps undertaken in the `3FTx_novel_taxa` project. The primary goal of this research is to investigate the evolution and occurrence of Three-Finger Toxin (3FTx)-like proteins beyond snakes, with a particular focus on invertebrates like centipedes.

The project was spurred by the intriguing discovery of toxin protein sequences in centipedes that occupy a similar structural space to snake 3FTx in ProtT5 embeddings. This is significant because venomous snake 3FTx are understood to have evolved from the Ly6 protein family after the divergence of snakes and lizards. Centipedes, being invertebrate Arthropods, present a fascinating case: did their 3FTx-like proteins also evolve from Ly6 in a separate specialization event, or is there another evolutionary pathway? This project aims to explore these questions, seeking to understand when and how these proteins evolved in centipedes and whether similar sequences exist in other arthropods.

## Starting Sequences

The project was initiated based on the discovery of the following centipede protein sequences that showed 3FTx-like characteristics:

- I6R1R5
- P0DPX5
- P0DPX9
- P0DPY0
- P0DPY1
- P0DPX7
- P0DPX8
- P0DPX6
- P0DPU8

accession:I6R1R5 OR accession:P0DPX5 OR accession:P0DPX9 OR accession:P0DPY0 OR accession:P0DPY1 OR accession:P0DPX7 OR accession:P0DPX8 OR accession:P0DPX6 OR accession:P0DPU8

## Workflow Overview

### 1. InterPro Data Retrieval

- Used `interpro_fetcher.py` to programmatically retrieve proteins from InterPro database
- Targeted the following InterPro IDs associated with 3FTx and Ly6/UPAR domains:
  - IPR016054 (Ly6/UPAR domain)
  - IPR026110 (LY6G5C)
  - IPR026524 (LY6G6d/LY6G6f)
  - IPR031424 (QVR-like)
  - IPR038773 (LY6G5B)
  - IPR039457 (LYPD6-like)
  - IPR039237 (LY6G6C)
  - IPR042339 (Ly6D)
  - IPR051445 (LY6H/LY6L_nAChR_modulators)
  - IPR035076 (Toxin/TOLIP)
  - IPR003571 (Three-finger toxin)
  - IPR045860 (Snake_toxin-like_sf)
- Retrieved protein sequences, metadata, and domain annotations:
  - Ran `uv run python src/interpro_retrieval/interpro_fetcher.py <INTERPRO_ID>` for:
    IPR016054, IPR026110, IPR026524, IPR031424, IPR038773, IPR039457, IPR039237, IPR042339, IPR051445, IPR035076, IPR003571, IPR045860
- Data stored in `data/raw/interpro/`
- all from the original 9 sequences are present, but one (P0DPX8)

### 2. Sequence Collection

```bash
uv run python src/data_prep/collect_interpro_sequences.py data/raw/interpro data/interm/interpro
```

- Collected complete protein sequences from all InterPro subfolders
- Deduplicated sequences by UniProt ID (52,828 unique sequences from 12 InterPro IDs)
- Output: `data/interm/interpro/interpro_complete_sequences.fasta`

### 3. Domain Extraction

```bash
uv run python src/data_prep/domain_extracter.py data/raw/interpro data/interm/interpro
```

- Extracted optimal non-overlapping protein domains from InterPro data
- Resolved overlapping domain annotations using optimization algorithms
- Output: `data/interm/interpro/extracted_domains.fasta`

### 4. Centipede Genome Extraction

**Manual annotation step:**

- Identified and annotated 3FTx-like gene regions in centipede genomes using sequence similarity searches
- Created GFF annotation files with gene coordinates and exon structures
- Files stored in `data/raw/centipede_genome/annotations/` (44 GFF files for 4 species)

**Automated extraction:**

```bash
uv run python src/data_prep/process_centipede_genomes.py
```

- Downloads reference genomes from NCBI (Cylindrodesmus punicus, Rhysida immarginata, Strigamia acuminata, Lithobius variegatus)
- Parses GFF annotations and deduplicates genes based on genomic coordinates (34 duplicates removed)
- Extracts and translates multi-exon genes, truncating at first stop codon
- Filters sequences < 50 amino acids
- **Output**: 151 protein sequences → `data/interm/centipede_genome/centipede_3ftx_proteins.fasta`

### 5. Data Merging and Clustering

- Merged all domain sequences and synteny data into a single dataset
- Used MMseqs2 to create representative sequences at 70% sequence similarity
- Resulting dataset: X representative sequences

### 6. Phylogenetic Analysis

- Generated multiple sequence alignments using FAMSA2
- Constructed phylogenetic trees using IQ-TREE3 with parameters: X
- Performed X iterations with X bootstrap replicates

### 7. Tree Annotation and Visualization

- Parsed metadata information for all sequences
- Generated iTOL-compatible annotation files
- Visualized phylogenetic trees with taxonomic and functional annotations in iTOL
- Created interactive tree visualizations for further investigation

Analyses:

- purity
- taxonomy LCA
- cluster size
- domain length

## Key Results

- Total sequences processed: X
- Representative sequences at 70% identity: X
- Final phylogenetic tree contains: X sequences
- Identified X distinct clades/groups of 3FTx-like proteins
- Confirmed presence of 3FTx-like proteins in X arthropod species
