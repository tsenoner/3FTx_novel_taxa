#!/usr/bin/env python3
"""
Process Centipede Genome Annotations

Extracts and translates protein sequences from manually annotated centipede genomes.

Workflow:
  1. Parse GFF annotation files with gene coordinates and exon structures
  2. Download reference genomes from NCBI (auto-skips if present)
  3. Deduplicate genes with identical genomic coordinates
  4. Extract and concatenate exon sequences
  5. Translate to proteins and truncate at first stop codon
  6. Filter sequences by minimum length (default: 50 AA)
  7. Output FASTA with standardized headers matching InterPro domain format

Output Format:
  Header: >ID taxa_id=NNNN protein_name=NAME domain_pos=1-LEN domain_length=LEN
  Compatible with annotate_clusters.py and merge_and_cluster.py
"""

import argparse
import csv
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import requests
from pyfaidx import Fasta


# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Species abbreviation mapping for GFF filenames
SPECIES_ABBREV_MAP = {
    "Cylindrodesmus": "Cpun",
    "Rhysida": "Rimm",
    "Strigamia": "Sacu",
    "Lithobius": "Lvar",
}

# Genetic code (standard codon table)
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class FastaSequence:
    """
    Represents a protein sequence with standardized metadata.

    Attributes:
        id: Unique sequence identifier
        taxa_id: NCBI taxonomy ID
        protein_name: Descriptive protein name
        sequence: Amino acid sequence
        domain_pos: Sequence position range (default: 1-length for full proteins)
        domain_length: Sequence length in amino acids
    """

    def __init__(
        self,
        seq_id: str,
        taxa_id: int,
        protein_name: str,
        sequence: str,
        domain_pos: str = None,
        domain_length: int = None,
    ):
        self.id = seq_id
        self.taxa_id = taxa_id
        self.protein_name = protein_name
        self.sequence = sequence.upper()
        self.domain_pos = domain_pos or f"1-{len(sequence)}"
        self.domain_length = domain_length or len(sequence)

    def __len__(self):
        return len(self.sequence)


class GFFFeature:
    """Represents a single GFF feature (exon)."""

    def __init__(
        self, seqid: str, start: int, end: int, strand: str, name: str, note: str
    ):
        self.seqid = seqid
        self.start = start  # 1-based, inclusive
        self.end = end  # 1-based, inclusive
        self.strand = strand
        self.name = name
        self.note = note

    def __repr__(self):
        return (
            f"GFFFeature({self.seqid}:{self.start}-{self.end}:"
            f"{self.strand}, {self.name})"
        )


class Gene:
    """Represents a gene with potentially multiple exons."""

    def __init__(self, name: str, species: str):
        self.name = name
        self.species = species
        self.exons: List[GFFFeature] = []

    def add_exon(self, exon: GFFFeature):
        """Add an exon to the gene."""
        self.exons.append(exon)

    def sort_exons(self):
        """Sort exons by genomic position."""
        # For negative strand, exons are typically in reverse order in GFF
        # Sort by position (ascending)
        self.exons.sort(key=lambda x: x.start)

    def get_strand(self) -> str:
        """Get the strand of the gene."""
        if self.exons:
            return self.exons[0].strand
        return "+"

    def get_seqid(self) -> str:
        """Get the sequence ID (chromosome/contig)."""
        if self.exons:
            return self.exons[0].seqid
        return ""


def parse_fasta(fasta_path: Path) -> Fasta:
    """
    Parse a FASTA file using pyfaidx.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Fasta object from pyfaidx
    """
    return Fasta(str(fasta_path))


def write_fasta(sequences: List[FastaSequence], output_path: Path, width: int = 80):
    """
    Write sequences to FASTA file with standardized header format.

    Header format matches interpro_domain_extractor.py output for consistency
    across all FASTA files in the pipeline.

    Args:
        sequences: List of FastaSequence objects
        output_path: Output file path
        width: Line width for sequence wrapping (default: 80)
    """
    with open(output_path, "w") as f:
        for seq in sequences:
            header = (
                f">{seq.id} taxa_id={seq.taxa_id} protein_name={seq.protein_name} "
                f"domain_pos={seq.domain_pos} domain_length={seq.domain_length}"
            )
            f.write(header + "\n")
            # Write sequence with line wrapping
            for i in range(0, len(seq.sequence), width):
                f.write(seq.sequence[i : i + width] + "\n")


def reverse_complement(dna_seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        dna_seq: DNA sequence string

    Returns:
        Reverse complement sequence
    """
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "N": "N",
        "n": "n",
    }

    rev_comp = "".join(complement.get(base, "N") for base in reversed(dna_seq))
    return rev_comp


def translate_dna(dna_seq: str) -> str:
    """
    Translate DNA sequence to protein.

    Args:
        dna_seq: DNA sequence string (must be uppercase)

    Returns:
        Protein sequence string
    """
    protein = []

    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i : i + 3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, "X")
            protein.append(aa)

    return "".join(protein)


def parse_gff_file(gff_path: Path) -> Tuple[List[Gene], str]:
    """
    Parse GFF file and extract gene annotations.

    Groups exons by gene name and extracts species abbreviation from filename.
    Normalizes gene names by removing '_exon_N' suffixes.

    Args:
        gff_path: Path to GFF annotation file

    Returns:
        Tuple of (list of Gene objects, species abbreviation)

    Example:
        Filename: "3ftx_UPAR-like_in_Cpun_OZ222540_annotations.gff"
        Returns: ([Gene(...), ...], "Cpun")
    """
    # Extract species abbreviation from filename pattern: "_in_SPECIES_"
    match = re.search(r"_in_([A-Za-z]+)_", gff_path.stem)
    species = match.group(1) if match else "Unknown"

    genes_dict = defaultdict(lambda: Gene(name="", species=species))

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            seqid, start, end, strand, attributes = (
                parts[0],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                parts[8],
            )

            # Parse GFF attributes (name and note)
            name_match = re.search(r"name=([^;]+)", attributes)
            if not name_match:
                continue

            gene_name = name_match.group(1)
            note_match = re.search(r"note=([^;]+)", attributes)
            note = note_match.group(1) if note_match else ""

            # Normalize: "gene_exon_1" -> "gene"
            normalized_name = re.sub(r"_exon_\d+$", "", gene_name)

            if normalized_name not in genes_dict:
                genes_dict[normalized_name] = Gene(normalized_name, species)

            exon = GFFFeature(seqid, start, end, strand, gene_name, note)
            genes_dict[normalized_name].add_exon(exon)

    # Sort exons by genomic position
    genes = list(genes_dict.values())
    for gene in genes:
        gene.sort_exons()

    return genes, species


def download_genome(assembly_id: str, output_dir: Path) -> Path:
    """
    Download genome from NCBI using Assembly ID.

    Args:
        assembly_id: NCBI Assembly ID (e.g., GCA_965125795.1)
        output_dir: Directory to save genome

    Returns:
        Path to downloaded genome file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{assembly_id}_genomic.fna"

    # Check if already downloaded
    if output_file.exists():
        logger.info(f"Genome {assembly_id} already exists, skipping download")
        return output_file

    logger.info(f"Downloading genome {assembly_id}...")

    # NCBI Datasets API URL
    base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
    url = f"{base_url}/{assembly_id}/download"

    try:
        # Using NCBI Datasets API
        params = {
            "include_annotation_type": ["GENOME_FASTA"],
            "filename": f"{assembly_id}.zip",
        }

        response = requests.get(url, params=params, stream=True, timeout=300)
        response.raise_for_status()

        # Save as zip and extract
        import zipfile
        from io import BytesIO

        zip_data = BytesIO(response.content)
        with zipfile.ZipFile(zip_data) as zip_ref:
            # Find the genomic.fna file in the zip
            fna_files = [
                f
                for f in zip_ref.namelist()
                if f.endswith("_genomic.fna") or f.endswith(".fna")
            ]

            if fna_files:
                # Extract the first .fna file
                with zip_ref.open(fna_files[0]) as source:
                    with open(output_file, "wb") as target:
                        target.write(source.read())
                logger.info(f"Successfully downloaded {assembly_id}")
            else:
                logger.error(f"No .fna file found in zip for {assembly_id}")
                raise Exception("No .fna file in zip")

    except Exception as e:
        logger.error(f"Failed to download {assembly_id}: {e}")
        logger.warning(f"Please download manually from NCBI and place at {output_file}")
        logger.warning(
            f"URL: https://www.ncbi.nlm.nih.gov/datasets/genome/{assembly_id}/"
        )
        raise

    return output_file


def extract_sequence(gene: Gene, genome_seqs: Fasta) -> str:
    """
    Extract and concatenate exon sequences for a gene.

    Handles both positive and negative strand genes. For negative strand,
    exons are reversed and the sequence is reverse-complemented.

    Args:
        gene: Gene object with exon annotations
        genome_seqs: Indexed genome sequences (pyfaidx Fasta object)

    Returns:
        Concatenated DNA sequence (empty string on error)
    """
    if not gene.exons:
        return ""

    seqid = gene.get_seqid()
    if seqid not in genome_seqs.keys():
        logger.error(f"Sequence {seqid} not found in genome")
        return ""

    # Extract all exon sequences (GFF: 1-based, pyfaidx slicing: 0-based)
    exon_sequences = []
    for exon in gene.exons:
        exon_seq = str(genome_seqs[seqid][exon.start - 1 : exon.end])
        exon_sequences.append(exon_seq)

    # Concatenate and handle strand orientation
    strand = gene.get_strand()
    if strand == "-":
        # Negative strand: reverse order and reverse complement
        concatenated = "".join(reversed(exon_sequences))
        return reverse_complement(concatenated)
    else:
        return "".join(exon_sequences)


def translate_and_process(dna_seq: str, gene_name: str) -> Tuple[str, bool]:
    """
    Translate DNA to protein, truncate at first stop, and check length.

    Args:
        dna_seq: DNA sequence
        gene_name: Name of the gene (for logging)

    Returns:
        Tuple of (protein sequence string, is_valid boolean)
    """
    if len(dna_seq) == 0:
        return "", False

    if len(dna_seq) % 3 != 0:
        logger.warning(
            f"{gene_name}: Length {len(dna_seq)} not divisible by 3, "
            f"truncating to {len(dna_seq) - len(dna_seq) % 3}"
        )
        dna_seq = dna_seq[: len(dna_seq) - len(dna_seq) % 3]

    # Translate
    try:
        protein = translate_dna(dna_seq.upper())
    except Exception as e:
        logger.error(f"{gene_name}: Translation failed: {e}")
        return "", False

    # Find first stop codon and truncate
    stop_pos = protein.find("*")
    if stop_pos != -1:
        protein = protein[:stop_pos]

    # Check length
    if len(protein) < 50:
        logger.debug(f"{gene_name}: Sequence too short ({len(protein)} aa), skipping")
        return protein, False

    return protein, True


def deduplicate_genes(genes: List[Gene]) -> Tuple[List[Gene], int]:
    """
    Remove duplicate genes based on genomic coordinates.

    Two genes are considered duplicates if they have:
    - Same scaffold (seqid)
    - Same strand
    - Identical exon coordinates (all start/end positions)

    Args:
        genes: List of Gene objects

    Returns:
        Tuple of (deduplicated genes list, number of duplicates removed)
    """
    unique_genes = []
    seen_signatures = set()
    duplicates_removed = 0

    for gene in genes:
        # Create a signature based on genomic coordinates
        if not gene.exons:
            continue

        seqid = gene.get_seqid()
        strand = gene.get_strand()

        # Sort exons by position to ensure consistent signature
        sorted_exons = sorted(gene.exons, key=lambda e: (e.start, e.end))
        exon_coords = tuple((e.start, e.end) for e in sorted_exons)

        signature = (seqid, strand, exon_coords)

        if signature not in seen_signatures:
            seen_signatures.add(signature)
            unique_genes.append(gene)
        else:
            duplicates_removed += 1
            logger.debug(
                f"Removing duplicate: {gene.name} at {seqid}:{strand} "
                f"(same coordinates as existing gene)"
            )

    return unique_genes, duplicates_removed


def load_assembly_mapping(tsv_path: Path) -> Dict[str, Dict]:
    """
    Load assembly and taxonomy information from TSV file.

    Maps species abbreviations (used in GFF filenames) to their assembly IDs,
    full species names, and NCBI taxonomy IDs.

    Args:
        tsv_path: Path to TSV with columns: AssemblyID, ToLID, OrganismId, SpeciesName

    Returns:
        Dict mapping species_abbrev -> {assembly_id, species_name, tol_id, taxa_id}

    Example:
        {"Cpun": {"assembly_id": "GCA_965125795.1", "taxa_id": 61981, ...}}
    """
    mapping = {}

    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species_name = row["SpeciesName"]
            taxa_id = int(row["OrganismId"])

            # Determine species abbreviation using genus name
            abbrev = None
            for genus, abbrev_code in SPECIES_ABBREV_MAP.items():
                if genus in species_name:
                    abbrev = abbrev_code
                    break

            if abbrev:
                mapping[abbrev] = {
                    "assembly_id": row["AssemblyID"],
                    "species_name": species_name,
                    "tol_id": row["ToLID"],
                    "taxa_id": taxa_id,
                }

    return mapping


def main():
    """
    Main entry point for centipede genome processing.

    Orchestrates the complete workflow from GFF parsing to FASTA output.
    """
    parser = argparse.ArgumentParser(
        description="Process centipede genome annotations and extract proteins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default processing
  %(prog)s

  # Custom paths with verbose logging
  %(prog)s --annotations-dir data/raw/annotations --verbose

  # Custom minimum length threshold
  %(prog)s --min-length 40
        """,
    )
    parser.add_argument(
        "--annotations-dir",
        type=Path,
        default=Path("data/raw/centipede_genome/annotations"),
        help="Directory containing GFF annotation files",
    )
    parser.add_argument(
        "--genomes-dir",
        type=Path,
        default=Path("data/raw/centipede_genome/genomes"),
        help="Directory to store downloaded genomes",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/interm/centipede_genome"),
        help="Directory for output FASTA file",
    )
    parser.add_argument(
        "--assembly-tsv",
        type=Path,
        default=Path("data/raw/centipede_genome/centipede_genome_tol_ids.tsv"),
        help="TSV file with assembly ID mappings",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=50,
        help="Minimum protein length in amino acids",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load assembly mapping
    logger.info(f"Loading assembly mapping from {args.assembly_tsv}")
    assembly_mapping = load_assembly_mapping(args.assembly_tsv)
    logger.info(f"Loaded mappings for {len(assembly_mapping)} species")

    # Download genomes if needed (automatically skips if already present)
    logger.info("Checking/downloading genomes...")
    for species_abbrev, info in assembly_mapping.items():
        assembly_id = info["assembly_id"]
        try:
            download_genome(assembly_id, args.genomes_dir)
        except Exception as e:
            logger.warning(f"Could not download {assembly_id}: {e}")
            logger.warning("Continuing with available genomes...")

    # Process each GFF file
    logger.info(f"Processing GFF files from {args.annotations_dir}")
    gff_files = sorted(args.annotations_dir.glob("*.gff"))

    all_proteins = []
    stats = {
        "total_genes": 0,
        "valid_proteins": 0,
        "too_short": 0,
        "failed_extraction": 0,
    }

    # Group GFF files by species to avoid reloading genomes
    gff_by_species = defaultdict(list)
    for gff_file in gff_files:
        genes, species = parse_gff_file(gff_file)
        gff_by_species[species].append((gff_file, genes))

    # Collect all genes per species for deduplication
    genes_by_species = defaultdict(list)
    for species, gff_list in gff_by_species.items():
        for gff_file, genes in gff_list:
            genes_by_species[species].extend(genes)

    # Deduplicate genes within each species based on genomic coordinates
    logger.info("\nDeduplicating genes based on genomic coordinates...")
    total_duplicates = 0
    for species in genes_by_species:
        before_count = len(genes_by_species[species])
        genes_by_species[species], dup_count = deduplicate_genes(
            genes_by_species[species]
        )
        after_count = len(genes_by_species[species])
        if dup_count > 0:
            logger.info(
                f"  {species}: Removed {dup_count} duplicates "
                f"({before_count} → {after_count} genes)"
            )
        total_duplicates += dup_count

    if total_duplicates > 0:
        logger.info(f"Total duplicates removed: {total_duplicates}")
    else:
        logger.info("No duplicates found")

    # Process each species
    for species in genes_by_species.keys():
        logger.info(f"\nProcessing species: {species}")

        # Find genome file
        if species not in assembly_mapping:
            logger.warning(f"No assembly mapping found for {species}, skipping")
            continue

        assembly_id = assembly_mapping[species]["assembly_id"]
        taxa_id = assembly_mapping[species]["taxa_id"]
        genome_file = args.genomes_dir / f"{assembly_id}_genomic.fna"

        if not genome_file.exists():
            logger.warning(f"Genome file not found: {genome_file}, skipping {species}")
            continue

        # Load genome
        logger.info(f"Loading genome from {genome_file}")
        try:
            genome_seqs = parse_fasta(genome_file)
            logger.info(f"Loaded {len(genome_seqs.keys())} sequences")
        except Exception as e:
            logger.error(f"Failed to load genome: {e}")
            continue

        # Process deduplicated genes for this species
        deduplicated_genes = genes_by_species[species]
        logger.info(f"Processing {len(deduplicated_genes)} unique genes")

        for gene in deduplicated_genes:
            stats["total_genes"] += 1

            # Extract sequence
            dna_seq = extract_sequence(gene, genome_seqs)
            if len(dna_seq) == 0:
                stats["failed_extraction"] += 1
                continue

            # Translate and process
            protein_seq, is_valid = translate_and_process(dna_seq, gene.name)

            if is_valid:
                # Create sequence record with standardized header
                record_id = f"{species}_{gene.name}_{gene.get_seqid()}"
                protein_name = gene.name.replace(" ", "_")  # FASTA header compatible
                record = FastaSequence(
                    seq_id=record_id,
                    taxa_id=taxa_id,
                    protein_name=protein_name,
                    sequence=protein_seq,
                )
                all_proteins.append(record)
                stats["valid_proteins"] += 1
            else:
                if len(protein_seq) > 0:
                    stats["too_short"] += 1
                else:
                    stats["failed_extraction"] += 1

    # Write output
    output_file = args.output_dir / "centipede_3ftx_proteins.fasta"
    logger.info(f"\nWriting {len(all_proteins)} proteins to {output_file}")
    write_fasta(all_proteins, output_file)

    # Print statistics
    logger.info("\n" + "=" * 50)
    logger.info("Processing Statistics:")
    logger.info(f"Total genes processed: {stats['total_genes']}")
    logger.info(f"Valid proteins (>= {args.min_length} aa): {stats['valid_proteins']}")
    logger.info(f"Too short (< {args.min_length} aa): {stats['too_short']}")
    logger.info(f"Failed extraction/translation: {stats['failed_extraction']}")
    logger.info("=" * 50)

    logger.info(f"\nDone! Output saved to {output_file}")


if __name__ == "__main__":
    main()
