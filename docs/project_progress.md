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

## Workflow Overview

### 1. InterPro Data Retrieval

- Used `interpro_fetcher.py` to programmatically retrieve proteins from InterPro database
- Targeted the following InterPro IDs associated with 3FTx and Ly6/UPAR domains:
  - IPR016054 (Ly6/UPAR domain)
  - IPR026110
  - IPR026524
  - IPR031424
  - IPR038773
  - IPR039457
  - IPR039237
  - IPR042339
  - IPR051445
  - IPR035076
  - IPR003571 (Three-finger toxin)
  - IPR045860
- Retrieved protein sequences, metadata, and domain annotations
- Data stored in `data/raw/data_for_tree/interpro/`

### 2. Domain Extraction

- Used `domain_extracter.py` to extract optimal non-overlapping protein domains
- Processed all InterPro data to identify and extract domain sequences
- Resolved overlapping domain annotations using optimization algorithms
- Generated FASTA files with extracted domain sequences

### 3. Centipede Synteny Analysis

- Performed synteny analysis to identify additional 3FTx-like sequences in centipede genomes
- Extracted similar sequences from centipede genomic data
- Results stored in `data/raw/centipede_3ftx_quiver_upar_like/`

### 4. Data Merging and Clustering

- Merged all domain sequences and synteny data into a single dataset
- Used MMseqs2 to create representative sequences at 70% sequence similarity
- Resulting dataset: X representative sequences

### 5. Phylogenetic Analysis

- Generated multiple sequence alignments using FAMSA2
- Constructed phylogenetic trees using IQ-TREE3 with parameters: X
- Performed X iterations with X bootstrap replicates

### 6. Tree Annotation and Visualization

- Parsed metadata information for all sequences
- Generated iTOL-compatible annotation files
- Visualized phylogenetic trees with taxonomic and functional annotations in iTOL
- Created interactive tree visualizations for further investigation

## Key Results

- Total sequences processed: X
- Representative sequences at 70% identity: X
- Final phylogenetic tree contains: X sequences
- Identified X distinct clades/groups of 3FTx-like proteins
- Confirmed presence of 3FTx-like proteins in X arthropod species
