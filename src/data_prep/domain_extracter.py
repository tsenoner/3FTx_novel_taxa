#!/usr/bin/env python3

import csv
import json
from pathlib import Path
from collections import defaultdict
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequences.

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        dict: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)

                # Start new sequence
                current_id = line[1:].split()[0]  # Get ID (first word after >)
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(records, output_file, wrap_width=80):
    """
    Write sequence records to a FASTA file.

    Args:
        records: List of dictionaries with 'id', 'sequence', and optional 'description'
        output_file: Path to output FASTA file
        wrap_width: Number of characters per line for sequence (default: 80)
    """
    with open(output_file, "w") as f:
        for record in records:
            seq_id = record["id"]
            sequence = record["sequence"]
            description = record.get("description", "")

            # Write header
            if description:
                f.write(f">{seq_id} {description}\n")
            else:
                f.write(f">{seq_id}\n")

            # Write sequence with wrapping
            for i in range(0, len(sequence), wrap_width):
                f.write(sequence[i : i + wrap_width] + "\n")


def process_single_protein(protein_id, domains, full_sequence, min_domain_length):
    """
    Process a single protein to extract optimal domains.
    This function is designed to be called in parallel.

    Returns:
        tuple: (extracted_domains, overlap_report, stats)
    """
    extracted_domains = []
    overlap_report = []
    stats = {"total_domains": 0, "filtered_out": 0, "kept_domains": 0}
    domain_counts = defaultdict(int)

    stats["total_domains"] = len(domains)

    # Filter domains by minimum length
    filtered_domains = []
    for domain in domains:
        domain_length = domain["domain_end"] - domain["domain_start"] + 1
        if domain_length >= min_domain_length:
            filtered_domains.append(domain)
        else:
            stats["filtered_out"] += 1

    stats["kept_domains"] = len(filtered_domains)

    if not filtered_domains:
        return extracted_domains, overlap_report, stats, domain_counts

    # Find overlaps between filtered domains
    n = len(filtered_domains)
    for i in range(n):
        for j in range(i + 1, n):
            d1, d2 = filtered_domains[i], filtered_domains[j]
            # Quick overlap check
            if not (
                d1["domain_end"] < d2["domain_start"]
                or d2["domain_end"] < d1["domain_start"]
            ):
                overlap_report.append(
                    {
                        "protein_id": protein_id,
                        "domain1": f"{d1['interpro_id']}:{d1['domain_start']}-{d1['domain_end']}",
                        "domain2": f"{d2['interpro_id']}:{d2['domain_start']}-{d2['domain_end']}",
                        "signature1": d1["signature_name"],
                        "signature2": d2["signature_name"],
                        "domain1_length": d1["domain_end"] - d1["domain_start"] + 1,
                        "domain2_length": d2["domain_end"] - d2["domain_start"] + 1,
                    }
                )

    # Find optimal domain set
    extractor = DomainExtractor.__new__(DomainExtractor)
    optimal_domains = extractor.find_optimal_domain_set(filtered_domains)

    # Extract sequences for optimal domains
    for i, domain in enumerate(optimal_domains):
        domain_seq = full_sequence[domain["domain_start"] - 1 : domain["domain_end"]]
        domain_id = f"{protein_id}_{domain['interpro_id']}_{i + 1}_{domain['domain_start']}-{domain['domain_end']}"

        description = (
            f"protein_name={domain['protein_name']} "
            f"organism={domain['organism']} "
            f"domain_pos={domain['domain_start']}-{domain['domain_end']} "
            f"signature={domain['signature_name']} "
            f"interpro_id={domain['interpro_id']} "
            f"length={len(domain_seq)}"
        )

        if "tax_id" in domain and domain["tax_id"]:
            description += f" tax_id={domain['tax_id']}"
        if "interpro_domain_name" in domain and domain["interpro_domain_name"]:
            description += f" interpro_domain_name={domain['interpro_domain_name']}"

        domain_record = {
            "id": domain_id,
            "sequence": domain_seq,
            "description": description,
        }

        extracted_domains.append(domain_record)
        domain_counts[domain["interpro_id"]] += 1

    return extracted_domains, overlap_report, stats, domain_counts


class DomainExtractor:
    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Create output directories
        self.reports_dir = self.output_dir / "reports"
        self.reports_dir.mkdir(exist_ok=True)

        # Store all protein sequences and annotations
        self.all_sequences = {}
        self.protein_annotations = defaultdict(list)  # protein_id -> list of domains
        self.extraction_stats = defaultdict(int)

    def load_sequences(self, fasta_file):
        """Load sequences from FASTA file"""
        try:
            sequences = parse_fasta(fasta_file)
        except Exception as e:
            logger.error(f"Error loading sequences from {fasta_file}: {e}")
            sequences = {}
        return sequences

    def load_domain_annotations(self, csv_file, interpro_id):
        """Load domain annotations from CSV file"""
        annotations = []
        try:
            with open(csv_file, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Handle both CSV formats
                    if "organism" in row:
                        annotation = {
                            "protein_id": row["protein_id"],
                            "interpro_id": interpro_id,
                            "protein_name": row["protein_name"],
                            "organism": row["organism"],
                            "sequence_length": int(row["sequence_length"]),
                            "domain_start": int(row["domain_start"]),
                            "domain_end": int(row["domain_end"]),
                            "domain_score": row.get("domain_score", ""),
                            "signature_name": row["signature_name"],
                        }
                    else:
                        annotation = {
                            "protein_id": row["protein_id"],
                            "interpro_id": interpro_id,
                            "protein_name": row["protein_name"],
                            "organism": row["scientific_name"],
                            "tax_id": row.get("tax_id", ""),
                            "sequence_length": int(row["sequence_length"]),
                            "domain_start": int(row["domain_start"]),
                            "domain_end": int(row["domain_end"]),
                            "domain_score": "",
                            "signature_name": row["signature_name"],
                            "interpro_domain_name": row.get("interpro_domain_name", ""),
                        }
                    annotations.append(annotation)
        except Exception as e:
            logger.error(f"Error loading annotations from {csv_file}: {e}")
        return annotations

    def collect_all_data(self):
        """Collect all sequences and annotations from all subfolders"""
        logger.info("Collecting all data from subfolders...")

        for subfolder in self.input_dir.iterdir():
            if not subfolder.is_dir():
                continue

            interpro_id = subfolder.name
            logger.info(f"Processing {interpro_id}")

            # Find files
            csv_file = next(subfolder.glob("*.csv"), None)
            fasta_file = next(
                (
                    f
                    for f in subfolder.iterdir()
                    if f.suffix.lower() in [".fasta", ".fa"]
                ),
                None,
            )

            if not csv_file or not fasta_file:
                logger.warning(f"Missing files in {interpro_id}")
                continue

            # Load sequences (merge with existing)
            sequences = self.load_sequences(fasta_file)
            self.all_sequences.update(sequences)

            # Load annotations and group by protein_id
            annotations = self.load_domain_annotations(csv_file, interpro_id)
            for ann in annotations:
                self.protein_annotations[ann["protein_id"]].append(ann)

    def domains_overlap(self, domain1, domain2):
        """Check if two domains overlap"""
        return not (
            domain1["domain_end"] < domain2["domain_start"]
            or domain2["domain_end"] < domain1["domain_start"]
        )

    def get_coverage(self, domains):
        """Calculate total coverage (non-overlapping) of a set of domains"""
        if not domains:
            return 0

        # Sort by start position
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
        Find the set of non-overlapping domains that maximizes coverage.
        Uses dynamic programming approach for optimal solution.
        """
        if not domains:
            return []

        if len(domains) == 1:
            return domains

        # For small numbers of domains, try all combinations for optimal solution
        # For larger sets, use faster interval scheduling with binary search
        if len(domains) <= 10:  # 2^10 = 1024 combinations, very fast
            return self._find_optimal_by_enumeration(domains)
        else:
            # For larger sets, use interval scheduling approach (much faster)
            return self._find_optimal_by_interval_scheduling(domains)

    def _find_optimal_by_enumeration(self, domains):
        """Find optimal set by trying all possible combinations (for small sets)"""
        from itertools import combinations

        # Pre-compute overlap matrix for faster checking
        n = len(domains)
        overlap_matrix = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if self.domains_overlap(domains[i], domains[j]):
                    overlap_matrix[i][j] = True
                    overlap_matrix[j][i] = True

        best_set = []
        best_coverage = 0

        # Try all possible combinations
        for r in range(1, n + 1):
            for combo_indices in combinations(range(n), r):
                # Check if this combination has overlaps using pre-computed matrix
                has_overlap = False
                for i_idx, i in enumerate(combo_indices):
                    for j in combo_indices[i_idx + 1 :]:
                        if overlap_matrix[i][j]:
                            has_overlap = True
                            break
                    if has_overlap:
                        break

                if not has_overlap:
                    combo_domains = [domains[i] for i in combo_indices]
                    coverage = self.get_coverage(combo_domains)
                    if coverage > best_coverage:
                        best_coverage = coverage
                        best_set = combo_domains

        return best_set

    def _find_optimal_by_interval_scheduling(self, domains):
        """
        For larger sets, use weighted interval scheduling approach with binary search.
        Sort by end position and use dynamic programming - O(n log n) complexity.
        """
        if not domains:
            return []

        # Sort domains by end position
        sorted_domains = sorted(domains, key=lambda x: x["domain_end"])
        n = len(sorted_domains)

        # Calculate weights (length of each domain)
        weights = [d["domain_end"] - d["domain_start"] + 1 for d in sorted_domains]

        # Pre-compute end positions for binary search
        end_positions = [d["domain_end"] for d in sorted_domains]

        # Find the latest non-overlapping domain using binary search - O(log n)
        def latest_non_overlapping(i):
            """Binary search for latest domain that doesn't overlap with domain i"""
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

        # DP arrays
        dp = [0] * n  # Maximum weight ending at or before position i
        include_domain = [False] * n  # Whether domain i is included

        # Pre-compute latest non-overlapping for all domains
        latest_compatible = [latest_non_overlapping(i) for i in range(n)]

        # Fill DP table
        dp[0] = weights[0]
        include_domain[0] = True

        for i in range(1, n):
            # Option 1: Don't include current domain
            exclude_weight = dp[i - 1]

            # Option 2: Include current domain
            include_weight = weights[i]
            latest = latest_compatible[i]
            if latest >= 0:
                include_weight += dp[latest]

            if include_weight > exclude_weight:
                dp[i] = include_weight
                include_domain[i] = True
            else:
                dp[i] = exclude_weight
                include_domain[i] = False

        # Reconstruct solution efficiently
        result = []
        i = n - 1
        while i >= 0:
            if include_domain[i]:
                result.append(sorted_domains[i])
                i = latest_compatible[i]
            else:
                i -= 1

        return result[::-1]  # Reverse to get domains in original order

    def extract_domains(self, min_domain_length=0, n_workers=None):
        """
        Extract optimal domain sets for each protein using parallel processing.

        Args:
            min_domain_length: Minimum domain length to keep
            n_workers: Number of parallel workers (None = use all CPUs)
        """
        logger.info(
            f"Finding optimal domain sets for each protein (min length: {min_domain_length})"
        )

        # Prepare data for parallel processing
        protein_tasks = []
        for protein_id, domains in self.protein_annotations.items():
            if protein_id not in self.all_sequences:
                logger.warning(f"No sequence found for {protein_id}")
                continue
            protein_tasks.append((protein_id, domains, self.all_sequences[protein_id]))

        logger.info(f"Processing {len(protein_tasks)} proteins in parallel...")

        all_extracted_domains = []
        overlap_report = []
        filtered_stats = {"total_domains": 0, "filtered_out": 0, "kept_domains": 0}

        # Process proteins in parallel
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Submit all tasks
            future_to_protein = {
                executor.submit(
                    process_single_protein,
                    protein_id,
                    domains,
                    sequence,
                    min_domain_length,
                ): protein_id
                for protein_id, domains, sequence in protein_tasks
            }

            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_protein):
                protein_id = future_to_protein[future]
                try:
                    domains, overlaps, stats, domain_counts = future.result()

                    all_extracted_domains.extend(domains)
                    overlap_report.extend(overlaps)

                    # Update stats
                    filtered_stats["total_domains"] += stats["total_domains"]
                    filtered_stats["filtered_out"] += stats["filtered_out"]
                    filtered_stats["kept_domains"] += stats["kept_domains"]

                    # Update extraction stats
                    for interpro_id, count in domain_counts.items():
                        self.extraction_stats[interpro_id] += count

                    completed += 1
                    if completed % 100 == 0:
                        logger.info(
                            f"Processed {completed}/{len(protein_tasks)} proteins"
                        )

                except Exception as e:
                    logger.error(f"Error processing protein {protein_id}: {e}")

        # Log filtering statistics
        logger.info("Domain filtering stats:")
        logger.info(f"  Total domains found: {filtered_stats['total_domains']}")
        logger.info(
            f"  Domains filtered out (< {min_domain_length} AA): {filtered_stats['filtered_out']}"
        )
        logger.info(f"  Domains kept for processing: {filtered_stats['kept_domains']}")

        return all_extracted_domains, overlap_report

    def generate_reports(self, overlap_report):
        """Generate summary reports"""

        # Overlap report
        if overlap_report:
            overlap_file = self.reports_dir / "overlap_report.csv"
            with open(overlap_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=overlap_report[0].keys())
                writer.writeheader()
                writer.writerows(overlap_report)
            logger.info(f"Overlap report saved to {overlap_file}")

        # Extraction statistics
        stats_data = {
            "extraction_stats": dict(self.extraction_stats),
            "total_domains_extracted": sum(self.extraction_stats.values()),
            "total_proteins_processed": len(self.protein_annotations),
            "total_sequences_loaded": len(self.all_sequences),
            "overlapping_cases": len(overlap_report),
        }

        stats_file = self.reports_dir / "extraction_stats.json"
        with open(stats_file, "w") as f:
            json.dump(stats_data, f, indent=2)

        # Summary report
        summary_file = self.reports_dir / "summary_report.txt"
        with open(summary_file, "w") as f:
            f.write("Domain Extraction Summary Report\n")
            f.write("=" * 40 + "\n\n")

            f.write(f"Total proteins processed: {len(self.protein_annotations)}\n")
            f.write(f"Total sequences loaded: {len(self.all_sequences)}\n")
            f.write(f"Total domains extracted: {sum(self.extraction_stats.values())}\n")
            f.write(
                f"Overlapping domain cases resolved: {len(set(o['protein_id'] for o in overlap_report))}\n\n"
            )

            f.write("Domains extracted per InterPro ID:\n")
            for interpro_id, count in sorted(self.extraction_stats.items()):
                f.write(f"  {interpro_id}: {count} domains\n")

        logger.info(f"Summary report saved to {summary_file}")

    def run(self, min_domain_length=0, n_workers=None):
        """Run the domain extraction process"""
        logger.info("Starting domain extraction...")
        logger.info(f"Minimum domain length: {min_domain_length} amino acids")
        if n_workers:
            logger.info(f"Using {n_workers} parallel workers")
        else:
            logger.info("Using all available CPU cores for parallel processing")

        # Collect all data first
        self.collect_all_data()

        logger.info(f"Loaded {len(self.all_sequences)} sequences")
        logger.info(f"Processing {len(self.protein_annotations)} proteins")

        # Extract optimal domains with length filtering and parallel processing
        extracted_domains, overlap_report = self.extract_domains(
            min_domain_length, n_workers
        )

        # Write output FASTA
        output_fasta = self.output_dir / "extracted_domains.fasta"
        if extracted_domains:
            write_fasta(extracted_domains, output_fasta)
            logger.info(f"Extracted {len(extracted_domains)} optimal domains")

        # Generate reports
        self.generate_reports(overlap_report)

        logger.info("Domain extraction completed!")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info(f"Extracted domains: {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract optimal non-overlapping protein domains with parallel processing"
    )
    parser.add_argument("input_dir", help="Directory containing InterPro subfolders")
    parser.add_argument("output_dir", help="Output directory for extracted domains")
    parser.add_argument(
        "--min-domain-length",
        type=int,
        default=50,
        help="Minimum domain length in amino acids (default: 50)",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: use all CPU cores)",
    )

    args = parser.parse_args()

    extractor = DomainExtractor(args.input_dir, args.output_dir)
    extractor.run(args.min_domain_length, args.n_workers)


if __name__ == "__main__":
    main()
