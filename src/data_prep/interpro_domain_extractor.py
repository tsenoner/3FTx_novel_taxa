#!/usr/bin/env python3
"""
Extract optimal non-overlapping protein domains from InterPro annotations.

This script processes InterPro domain annotations and protein sequences to:
1. Filter domains by minimum length
2. Resolve overlapping domains by finding optimal non-overlapping sets
3. Extract domain sequences with metadata
4. Generate comprehensive reports

The optimization uses either enumeration (for ≤10 domains) or dynamic programming
with weighted interval scheduling (for larger sets) to maximize coverage.
"""

import csv
from pathlib import Path
from collections import defaultdict
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# ============================================================================
# FASTA I/O FUNCTIONS
# ============================================================================


def parse_fasta(fasta_file):
    """
    Parse FASTA file with metadata in headers.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        dict: {seq_id: (sequence, metadata_dict)}
    """
    sequences = {}
    current_id = None
    current_metadata = {}
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = ("".join(current_seq), current_metadata)

                parts = line[1:].split()
                current_id = parts[0]
                current_metadata = {}

                for part in parts[1:]:
                    if "=" in part:
                        key, value = part.split("=", 1)
                        current_metadata[key] = value

                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ("".join(current_seq), current_metadata)

    return sequences


def write_fasta(records, output_file, wrap_width=80):
    """
    Write domain records to FASTA file.

    Args:
        records: List of dicts with 'id', 'sequence', and optional 'description'
        output_file: Output path
        wrap_width: Characters per line (default: 80)
    """
    with open(output_file, "w") as f:
        for record in records:
            header = f">{record['id']}"
            if description := record.get("description"):
                header += f" {description}"
            f.write(f"{header}\n")

            sequence = record["sequence"]
            for i in range(0, len(sequence), wrap_width):
                f.write(f"{sequence[i : i + wrap_width]}\n")


# ============================================================================
# DOMAIN PROCESSING FUNCTIONS
# ============================================================================


def process_single_protein(
    protein_id, domains, full_sequence, protein_metadata, min_domain_length
):
    """
    Extract optimal non-overlapping domains from a single protein.

    This function is designed for parallel execution.

    Args:
        protein_id: Protein identifier
        domains: List of domain annotations
        full_sequence: Full protein sequence
        protein_metadata: Dict with taxa_id, protein_name, interpro
        min_domain_length: Minimum domain length threshold

    Returns:
        tuple: (extracted_domains, overlap_report, stats, domain_counts)
    """
    extracted_domains = []
    overlap_report = []
    stats = {"total_domains": 0, "filtered_out": 0, "kept_domains": 0}
    domain_counts = defaultdict(int)

    stats["total_domains"] = len(domains)

    # Filter by minimum length
    filtered_domains = [
        d
        for d in domains
        if (d["domain_end"] - d["domain_start"] + 1) >= min_domain_length
    ]
    stats["filtered_out"] = len(domains) - len(filtered_domains)
    stats["kept_domains"] = len(filtered_domains)

    if not filtered_domains:
        return extracted_domains, overlap_report, stats, domain_counts

    # Identify overlaps
    for i, d1 in enumerate(filtered_domains):
        for d2 in filtered_domains[i + 1 :]:
            if not (
                d1["domain_end"] < d2["domain_start"]
                or d2["domain_end"] < d1["domain_start"]
            ):
                overlap_report.append(
                    {
                        "protein_id": protein_id,
                        "domain1": f"{d1['interpro_id']}:{d1['domain_start']}-{d1['domain_end']}",
                        "domain2": f"{d2['interpro_id']}:{d2['domain_start']}-{d2['domain_end']}",
                        "domain1_length": d1["domain_end"] - d1["domain_start"] + 1,
                        "domain2_length": d2["domain_end"] - d2["domain_start"] + 1,
                    }
                )

    # Find optimal non-overlapping set
    extractor = DomainExtractor.__new__(DomainExtractor)
    optimal_domains = extractor.find_optimal_domain_set(filtered_domains)

    # Extract domain sequences
    for i, domain in enumerate(optimal_domains):
        domain_start = domain["domain_start"]
        domain_end = domain["domain_end"]
        domain_seq = full_sequence[domain_start - 1 : domain_end]
        domain_length = len(domain_seq)

        domain_id = (
            f"{protein_id}_{domain['interpro_id']}_{i + 1}_{domain_start}-{domain_end}"
        )

        description_parts = []
        if "taxa_id" in protein_metadata:
            description_parts.append(f"taxa_id={protein_metadata['taxa_id']}")
        if "protein_name" in protein_metadata:
            description_parts.append(f"protein_name={protein_metadata['protein_name']}")
        if "interpro" in protein_metadata:
            description_parts.append(f"interpro={protein_metadata['interpro']}")
        description_parts.extend(
            [
                f"domain_pos={domain_start}-{domain_end}",
                f"domain_length={domain_length}",
            ]
        )

        extracted_domains.append(
            {
                "id": domain_id,
                "sequence": domain_seq,
                "description": " ".join(description_parts),
            }
        )
        domain_counts[domain["interpro_id"]] += 1

    return extracted_domains, overlap_report, stats, domain_counts


# ============================================================================
# DOMAIN EXTRACTOR CLASS
# ============================================================================


class DomainExtractor:
    """Extract optimal non-overlapping protein domains from InterPro data."""

    def __init__(self, input_dir, output_dir, filtered_fasta=None):
        """
        Initialize the domain extractor.

        Args:
            input_dir: Directory containing InterPro subfolders with CSV files
            output_dir: Output directory for results
            filtered_fasta: Optional path to pre-filtered FASTA sequences
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.reports_dir = self.output_dir / "reports"
        self.reports_dir.mkdir(exist_ok=True)

        self.all_sequences = {}
        self.protein_metadata = {}
        self.protein_annotations = defaultdict(list)
        self.extraction_stats = defaultdict(int)

        self.filtered_fasta = Path(filtered_fasta) if filtered_fasta else None

    def load_sequences(self, fasta_file):
        """Load sequences from FASTA file."""
        try:
            return parse_fasta(fasta_file)
        except Exception as e:
            logger.error(f"Error loading sequences from {fasta_file}: {e}")
            return {}

    def load_domain_annotations(self, csv_file, interpro_id):
        """Load domain annotations from CSV file."""
        annotations = []
        try:
            with open(csv_file) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    annotations.append(
                        {
                            "protein_id": row["protein_id"],
                            "interpro_id": interpro_id,
                            "domain_start": int(row["domain_start"]),
                            "domain_end": int(row["domain_end"]),
                        }
                    )
        except Exception as e:
            logger.error(f"Error loading annotations from {csv_file}: {e}")
        return annotations

    def collect_all_data(self):
        """Collect sequences and annotations from InterPro data."""
        if self.filtered_fasta and self.filtered_fasta.exists():
            logger.info(f"Loading sequences from: {self.filtered_fasta}")
            sequences = self.load_sequences(self.filtered_fasta)
            for seq_id, (sequence, metadata) in sequences.items():
                self.all_sequences[seq_id] = sequence
                self.protein_metadata[seq_id] = metadata
            logger.info(f"Loaded {len(self.all_sequences):,} sequences")
        else:
            logger.info("Collecting sequences from InterPro subfolders...")
            for subfolder in self.input_dir.iterdir():
                if not subfolder.is_dir():
                    continue

                fasta_file = next(
                    (
                        f
                        for f in subfolder.iterdir()
                        if f.suffix.lower() in [".fasta", ".fa"]
                    ),
                    None,
                )
                if fasta_file:
                    sequences = self.load_sequences(fasta_file)
                    for seq_id, (sequence, metadata) in sequences.items():
                        if seq_id not in self.all_sequences:
                            self.all_sequences[seq_id] = sequence
                            self.protein_metadata[seq_id] = metadata

        logger.info("Loading domain annotations...")
        for subfolder in self.input_dir.iterdir():
            if not subfolder.is_dir():
                continue

            interpro_id = subfolder.name
            csv_file = next(subfolder.glob("*.csv"), None)
            if csv_file:
                annotations = self.load_domain_annotations(csv_file, interpro_id)
                for ann in annotations:
                    self.protein_annotations[ann["protein_id"]].append(ann)

    def domains_overlap(self, domain1, domain2):
        """Check if two domains overlap."""
        return not (
            domain1["domain_end"] < domain2["domain_start"]
            or domain2["domain_end"] < domain1["domain_start"]
        )

    def get_coverage(self, domains):
        """Calculate total non-overlapping coverage of domain set."""
        if not domains:
            return 0

        sorted_domains = sorted(domains, key=lambda x: x["domain_start"])
        coverage = 0
        last_end = 0

        for domain in sorted_domains:
            start = max(domain["domain_start"], last_end + 1)
            end = domain["domain_end"]
            if start <= end:
                coverage += end - start + 1
                last_end = end

        return coverage

    def find_optimal_domain_set(self, domains):
        """
        Find optimal non-overlapping domain set that maximizes coverage.

        Uses enumeration for small sets (≤10 domains) or dynamic programming
        with weighted interval scheduling for larger sets.
        """
        if not domains:
            return []
        if len(domains) == 1:
            return domains

        if len(domains) <= 10:
            return self._find_optimal_by_enumeration(domains)
        else:
            return self._find_optimal_by_interval_scheduling(domains)

    def _find_optimal_by_enumeration(self, domains):
        """Find optimal set by trying all combinations (for ≤10 domains)."""
        from itertools import combinations

        n = len(domains)
        overlap_matrix = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if self.domains_overlap(domains[i], domains[j]):
                    overlap_matrix[i][j] = overlap_matrix[j][i] = True

        best_set = []
        best_coverage = 0

        for r in range(1, n + 1):
            for combo_indices in combinations(range(n), r):
                has_overlap = any(
                    overlap_matrix[i][j]
                    for idx, i in enumerate(combo_indices)
                    for j in combo_indices[idx + 1 :]
                )

                if not has_overlap:
                    combo_domains = [domains[i] for i in combo_indices]
                    coverage = self.get_coverage(combo_domains)
                    if coverage > best_coverage:
                        best_coverage = coverage
                        best_set = combo_domains

        return best_set

    def _find_optimal_by_interval_scheduling(self, domains):
        """Weighted interval scheduling with DP (O(n log n))."""
        if not domains:
            return []

        sorted_domains = sorted(domains, key=lambda x: x["domain_end"])
        n = len(sorted_domains)
        weights = [d["domain_end"] - d["domain_start"] + 1 for d in sorted_domains]
        end_positions = [d["domain_end"] for d in sorted_domains]

        def latest_non_overlapping(i):
            """Binary search for latest compatible domain."""
            start = sorted_domains[i]["domain_start"]
            left, right = 0, i - 1
            result = -1

            while left <= right:
                mid = (left + right) // 2
                if end_positions[mid] < start:
                    result = mid
                    left = mid + 1
                else:
                    right = mid - 1

            return result

        dp = [0] * n
        include_domain = [False] * n
        latest_compatible = [latest_non_overlapping(i) for i in range(n)]

        dp[0] = weights[0]
        include_domain[0] = True

        for i in range(1, n):
            exclude_weight = dp[i - 1]
            include_weight = weights[i]
            if (latest := latest_compatible[i]) >= 0:
                include_weight += dp[latest]

            if include_weight > exclude_weight:
                dp[i] = include_weight
                include_domain[i] = True
            else:
                dp[i] = exclude_weight

        result = []
        i = n - 1
        while i >= 0:
            if include_domain[i]:
                result.append(sorted_domains[i])
                i = latest_compatible[i]
            else:
                i -= 1

        return result[::-1]

    def extract_domains(self, min_domain_length=0, n_workers=None):
        """
        Extract optimal domains for all proteins in parallel.

        Args:
            min_domain_length: Minimum domain length threshold
            n_workers: Number of parallel workers (None = all CPUs)

        Returns:
            tuple: (extracted_domains, overlap_report, filtered_stats)
        """
        logger.info(f"Finding optimal domain sets (min length: {min_domain_length} AA)")

        protein_tasks = []
        skipped_proteins = 0

        for protein_id, domains in self.protein_annotations.items():
            if protein_id not in self.all_sequences:
                skipped_proteins += 1
                continue
            metadata = self.protein_metadata.get(protein_id, {})
            protein_tasks.append(
                (protein_id, domains, self.all_sequences[protein_id], metadata)
            )

        if skipped_proteins > 0:
            logger.info(
                f"Skipped {skipped_proteins:,} proteins not in filtered sequences"
            )

        logger.info(f"Processing {len(protein_tasks):,} proteins...")

        all_extracted_domains = []
        overlap_report = []
        filtered_stats = {"total_domains": 0, "filtered_out": 0, "kept_domains": 0}

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_protein = {
                executor.submit(
                    process_single_protein,
                    protein_id,
                    domains,
                    sequence,
                    metadata,
                    min_domain_length,
                ): protein_id
                for protein_id, domains, sequence, metadata in protein_tasks
            }

            for future in tqdm(
                as_completed(future_to_protein),
                total=len(protein_tasks),
                desc="Processing proteins",
                unit="protein",
            ):
                try:
                    domains, overlaps, stats, domain_counts = future.result()
                    all_extracted_domains.extend(domains)
                    overlap_report.extend(overlaps)

                    filtered_stats["total_domains"] += stats["total_domains"]
                    filtered_stats["filtered_out"] += stats["filtered_out"]
                    filtered_stats["kept_domains"] += stats["kept_domains"]

                    for interpro_id, count in domain_counts.items():
                        self.extraction_stats[interpro_id] += count

                except Exception as e:
                    protein_id = future_to_protein[future]
                    logger.error(f"Error processing {protein_id}: {e}")

        logger.info(f"Extracted {len(all_extracted_domains):,} optimal domains")
        logger.info(f"  Total domains: {filtered_stats['total_domains']:,}")
        logger.info(f"  Filtered out: {filtered_stats['filtered_out']:,}")
        logger.info(f"  Kept: {filtered_stats['kept_domains']:,}")

        return all_extracted_domains, overlap_report, filtered_stats

    def generate_reports(self, overlap_report, extracted_domains, filtered_stats):
        """Generate summary reports."""
        if overlap_report:
            overlap_file = self.reports_dir / "domain_overlap.csv"
            with open(overlap_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=overlap_report[0].keys())
                writer.writeheader()
                writer.writerows(overlap_report)
            logger.info(f"Overlap report: {overlap_file}")

        domain_lengths = [len(d["sequence"]) for d in extracted_domains]
        avg_length = sum(domain_lengths) / len(domain_lengths) if domain_lengths else 0
        min_length = min(domain_lengths) if domain_lengths else 0
        max_length = max(domain_lengths) if domain_lengths else 0

        summary_file = self.reports_dir / "domain_extraction.txt"
        with open(summary_file, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("DOMAIN EXTRACTION SUMMARY\n")
            f.write("=" * 80 + "\n\n")

            f.write("OVERVIEW\n")
            f.write("-" * 80 + "\n")
            f.write(
                f"Total proteins in annotations:      {len(self.protein_annotations):>10,}\n"
            )
            f.write(
                f"Total sequences loaded:              {len(self.all_sequences):>10,}\n"
            )
            f.write(
                f"Total domains extracted:             {sum(self.extraction_stats.values()):>10,}\n"
            )
            f.write(
                f"Proteins with overlapping domains:   {len(set(o['protein_id'] for o in overlap_report)):>10,}\n"
            )
            f.write(
                f"Total overlap cases resolved:        {len(overlap_report):>10,}\n\n"
            )

            f.write("DOMAIN FILTERING\n")
            f.write("-" * 80 + "\n")
            f.write(
                f"Total domains found:                 {filtered_stats['total_domains']:>10,}\n"
            )
            f.write(
                f"Domains filtered out (too short):    {filtered_stats['filtered_out']:>10,}\n"
            )
            f.write(
                f"Domains kept for processing:         {filtered_stats['kept_domains']:>10,}\n"
            )
            f.write(
                f"Final domains after optimization:    {len(extracted_domains):>10,}\n\n"
            )

            f.write("DOMAIN LENGTH STATISTICS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Minimum domain length:               {min_length:>10} AA\n")
            f.write(f"Maximum domain length:               {max_length:>10} AA\n")
            f.write(f"Average domain length:               {avg_length:>10.1f} AA\n\n")

            f.write("DOMAINS PER INTERPRO ID\n")
            f.write("-" * 80 + "\n")
            total = sum(self.extraction_stats.values())
            for interpro_id, count in sorted(
                self.extraction_stats.items(), key=lambda x: x[1], reverse=True
            ):
                percentage = (count / total) * 100
                f.write(
                    f"{interpro_id:>12}:  {count:>6,} domains  ({percentage:>5.1f}%)\n"
                )

            f.write("\n" + "=" * 80 + "\n")

        logger.info(f"Summary report: {summary_file}")

    def run(self, min_domain_length=50, n_workers=None):
        """Run the complete domain extraction pipeline."""
        logger.info("=" * 80)
        logger.info("DOMAIN EXTRACTION PIPELINE")
        logger.info("=" * 80)
        logger.info(f"Minimum domain length: {min_domain_length} AA")
        logger.info(f"Parallel workers: {n_workers or 'all available CPUs'}")
        logger.info("")

        self.collect_all_data()
        logger.info(f"Loaded {len(self.all_sequences):,} sequences")
        logger.info(f"Processing {len(self.protein_annotations):,} proteins\n")

        extracted_domains, overlap_report, filtered_stats = self.extract_domains(
            min_domain_length, n_workers
        )

        output_fasta = self.output_dir / "extracted_domains.fasta"
        write_fasta(extracted_domains, output_fasta)

        self.generate_reports(overlap_report, extracted_domains, filtered_stats)

        logger.info("")
        logger.info("=" * 80)
        logger.info("EXTRACTION COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Extracted domains: {output_fasta}")
        logger.info(f"Total domains: {len(extracted_domains):,}")


def main():
    """Main entry point for domain extraction."""
    parser = argparse.ArgumentParser(
        description="Extract optimal non-overlapping protein domains",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Extract domains from filtered sequences
  python interpro_domain_extractor.py \\
    data/raw/interpro \\
    data/interm/interpro \\
    --filtered-fasta data/interm/interpro/interpro_metazoa_sequences.fasta \\
    --min-domain-length 50

  # Use all InterPro sequences (not recommended)
  python interpro_domain_extractor.py \\
    data/raw/interpro \\
    data/interm/interpro \\
    --min-domain-length 50
        """,
    )
    parser.add_argument("input_dir", help="Directory with InterPro subfolders")
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument(
        "--filtered-fasta",
        help="Pre-filtered FASTA file (e.g., Metazoa-only sequences)",
    )
    parser.add_argument(
        "--min-domain-length",
        type=int,
        default=50,
        help="Minimum domain length in amino acids (default: 50)",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        help="Number of parallel workers (default: all CPUs)",
    )

    args = parser.parse_args()

    extractor = DomainExtractor(args.input_dir, args.output_dir, args.filtered_fasta)
    extractor.run(args.min_domain_length, args.n_workers)


if __name__ == "__main__":
    main()
