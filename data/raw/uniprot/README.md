# UniProt Data for 3FTx Analysis

This directory contains data retrieved from the UniProt database.

## Retrieval Details

-   **Source:** UniProt (https://www.uniprot.org/)
-   **Date of Retrieval:** 12-05-2025
-   **Query Used:**
    ```
    (((xref:interpro-IPR016054 OR xref:interpro-IPR026110 OR xref:interpro-IPR026524 OR xref:interpro-IPR031424 OR xref:interpro-IPR038773 OR xref:interpro-IPR039457 OR xref:interpro-IPR039237 OR xref:interpro-IPR042339 OR xref:interpro-IPR051445 OR xref:interpro-IPR035076 OR xref:interpro-IPR003571 OR xref:interpro-IPR045860) AND fragment:false) OR (accession:I6R1R5 OR accession:P0DPX5 OR accession:P0DPX9 OR accession:P0DPY0 OR accession:P0DPY1 OR accession:P0DPX7 OR accession:P0DPX8 OR accession:P0DPX6 OR accession:P0DPU8)) AND (length:[51 TO 199]) AND (taxonomy_id:33208)
    ```

This query was designed to retrieve specific protein sequences relevant to the project:

*   **InterPro IDs (`xref:interpro-...`)**: These IDs were selected because they are all related to the Three-Finger Toxin (3FTx) family (`IPR003571`, `IPR045860`) or the related LY6/UPAR domain (`IPR016054`, `IPR035076`) and its associated families, which are considered potential ancestors or relatives of 3FTx. Including these aims to capture a broad set of sequences for correlation analysis. The `fragment:false` condition ensures complete sequences are retrieved where possible.
*   **Accession Numbers (`accession:...`)**: These specific accessions (I6R1R5, P0DPX5, etc.) correspond to Scoloptoxins found in centipedes. Preliminary analysis using ProtSpace suggested these might be structurally related to 3FTx, forming the initial motivation for this project.
*   **Length Constraint (`length:[51 TO 199]`)**: This filter restricts the results to proteins within a typical size range for single-domain 3FTx/LY6 proteins, aiming to exclude multi-domain proteins or very short fragments.
*   **Taxonomy ID (`taxonomy_id:33208`)**: This restricts the search to Metazoa (animals), excluding sequences from bacteria or other kingdoms, as the focus is on the evolution and correlation within the animal kingdom.