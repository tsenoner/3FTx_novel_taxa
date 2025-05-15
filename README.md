# 3FTx_novel_taxa

Investigates Three-Finger Toxins (3FTx) evolution beyond snakes.
Full details: `docs/project_specification.md`.

## Setup

1.  Requires `uv` ([https://github.com/astral-sh/uv](https://github.com/astral-sh/uv)).
2.  Clone repo & `cd` to root.
3.  Install dependencies:
    ```bash
    uv venv
    uv sync
    ```

## Example: Process FASTA

Clean raw centipede FASTA (remove stops/duplicates, new IDs) and output to default `data/interm/centipede_3ftx_quiver_upar_like/`:
```bash
uv run python src/data_processing/process_centipede_fasta.py -i data/raw/centipede_3ftx_quiver_upar_like/3ftx_Quiver_UPAR-like_in_Centipedes_AW_210425.fasta
```
For more, see script help (`-h`) or project docs.

## Data Sources

Raw data is located in `data/raw/`:
*   `uniprot/`: Sequences queried from UniProt (see `README.md` within for details).
*   `centipede_3ftx_quiver_upar_like/`: Manually annotated centipede FASTA/GFF files.
