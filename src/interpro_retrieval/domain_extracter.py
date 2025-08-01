#!/usr/bin/env python3

import os
import csv
import json
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DomainExtractor:
    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create output directories
        self.domains_dir = self.output_dir / "extracted_domains"
        self.reports_dir = self.output_dir / "reports"
        self.domains_dir.mkdir(exist_ok=True)
        self.reports_dir.mkdir(exist_ok=True)
        
        # Store all protein sequences and annotations
        self.all_sequences = {}
        self.protein_annotations = defaultdict(list)  # protein_id -> list of domains
        self.extraction_stats = defaultdict(int)

    def load_sequences(self, fasta_file):
        """Load sequences from FASTA file"""
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences[record.id] = str(record.seq)
        except Exception as e:
            logger.error(f"Error loading sequences from {fasta_file}: {e}")
        return sequences

    def load_domain_annotations(self, csv_file, interpro_id):
        """Load domain annotations from CSV file"""
        annotations = []
        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Handle both CSV formats
                    if 'organism' in row:
                        annotation = {
                            'protein_id': row['protein_id'],
                            'interpro_id': interpro_id,
                            'protein_name': row['protein_name'],
                            'organism': row['organism'],
                            'sequence_length': int(row['sequence_length']),
                            'domain_start': int(row['domain_start']),
                            'domain_end': int(row['domain_end']),
                            'domain_score': row.get('domain_score', ''),
                            'signature_name': row['signature_name']
                        }
                    else:
                        annotation = {
                            'protein_id': row['protein_id'],
                            'interpro_id': interpro_id,
                            'protein_name': row['protein_name'],
                            'organism': row['scientific_name'],
                            'tax_id': row.get('tax_id', ''),
                            'sequence_length': int(row['sequence_length']),
                            'domain_start': int(row['domain_start']),
                            'domain_end': int(row['domain_end']),
                            'domain_score': '',
                            'signature_name': row['signature_name'],
                            'interpro_domain_name': row.get('interpro_domain_name', '')
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
            fasta_file = next((f for f in subfolder.iterdir() 
                             if f.suffix.lower() in ['.fasta', '.fa']), None)
            
            if not csv_file or not fasta_file:
                logger.warning(f"Missing files in {interpro_id}")
                continue
            
            # Load sequences (merge with existing)
            sequences = self.load_sequences(fasta_file)
            self.all_sequences.update(sequences)
            
            # Load annotations and group by protein_id
            annotations = self.load_domain_annotations(csv_file, interpro_id)
            for ann in annotations:
                self.protein_annotations[ann['protein_id']].append(ann)

    def domains_overlap(self, domain1, domain2):
        """Check if two domains overlap"""
        return not (domain1['domain_end'] < domain2['domain_start'] or 
                   domain2['domain_end'] < domain1['domain_start'])

    def get_coverage(self, domains):
        """Calculate total coverage (non-overlapping) of a set of domains"""
        if not domains:
            return 0
        
        # Sort by start position
        sorted_domains = sorted(domains, key=lambda x: x['domain_start'])
        coverage = 0
        last_end = 0
        
        for domain in sorted_domains:
            start = max(domain['domain_start'], last_end + 1)
            end = domain['domain_end']
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
        
####################################################################################

#       IF THIS TAKES TOO LONG; PUT BACK IN THE CODE; IF A PROTEIN HAS MORE THAN 15 DOMAINS; IT WILL USE A WORSE BUT FASTER ALGORITHM

        # For small numbers of domains, try all combinations
        #if len(domains) <= 15:  # 2^15 = 32K combinations, manageable
        return self._find_optimal_by_enumeration(domains)
        #else:
            # For larger sets, use interval scheduling approach
        #    return self._find_optimal_by_interval_scheduling(domains)
    
#####################################################################################

    def _find_optimal_by_enumeration(self, domains):
        """Find optimal set by trying all possible combinations (for small sets)"""
        from itertools import combinations
        
        best_set = []
        best_coverage = 0
        
        # Try all possible combinations
        for r in range(1, len(domains) + 1):
            for combo in combinations(domains, r):
                # Check if this combination has overlaps
                has_overlap = False
                for i, d1 in enumerate(combo):
                    for j, d2 in enumerate(combo[i+1:], i+1):
                        if self.domains_overlap(d1, d2):
                            has_overlap = True
                            break
                    if has_overlap:
                        break
                
                if not has_overlap:
                    coverage = self.get_coverage(list(combo))
                    if coverage > best_coverage:
                        best_coverage = coverage
                        best_set = list(combo)
        
        return best_set
    
    def _find_optimal_by_interval_scheduling(self, domains):
        """
        For larger sets, use weighted interval scheduling approach.
        Sort by end position and use dynamic programming.
        """
        if not domains:
            return []
        
        # Sort domains by end position
        sorted_domains = sorted(domains, key=lambda x: x['domain_end'])
        n = len(sorted_domains)
        
        # Calculate weights (length of each domain)
        weights = [d['domain_end'] - d['domain_start'] + 1 for d in sorted_domains]
        
        # Find the latest non-overlapping domain for each domain
        def latest_non_overlapping(i):
            for j in range(i-1, -1, -1):
                if sorted_domains[j]['domain_end'] < sorted_domains[i]['domain_start']:
                    return j
            return -1
        
        # DP arrays
        dp = [0] * n  # Maximum weight ending at or before position i
        parent = [-1] * n  # To reconstruct solution
        
        # Fill DP table
        dp[0] = weights[0]
        
        for i in range(1, n):
            # Option 1: Don't include current domain
            exclude_weight = dp[i-1]
            
            # Option 2: Include current domain
            include_weight = weights[i]
            latest = latest_non_overlapping(i)
            if latest >= 0:
                include_weight += dp[latest]
            
            if include_weight > exclude_weight:
                dp[i] = include_weight
                parent[i] = latest
            else:
                dp[i] = exclude_weight
                parent[i] = -2  # Indicates we didn't include domain i
        
        # Reconstruct solution
        result = []
        i = n - 1
        while i >= 0:
            if parent[i] != -2:  # Domain i was included
                result.append(sorted_domains[i])
                i = parent[i]
            else:
                i -= 1
        
        return result

    def extract_domains(self, min_domain_length=0):
        """Extract optimal domain sets for each protein"""
        logger.info(f"Finding optimal domain sets for each protein (min length: {min_domain_length})")
        
        all_extracted_domains = []
        overlap_report = []
        filtered_stats = {'total_domains': 0, 'filtered_out': 0, 'kept_domains': 0}
        
        for protein_id, domains in self.protein_annotations.items():
            if protein_id not in self.all_sequences:
                logger.warning(f"No sequence found for {protein_id}")
                continue
            
            filtered_stats['total_domains'] += len(domains)
            
            # Filter domains by minimum length
            filtered_domains = []
            for domain in domains:
                domain_length = domain['domain_end'] - domain['domain_start'] + 1
                if domain_length >= min_domain_length:
                    filtered_domains.append(domain)
                else:
                    filtered_stats['filtered_out'] += 1
            
            filtered_stats['kept_domains'] += len(filtered_domains)
            
            if not filtered_domains:
                logger.debug(f"No domains pass length filter for {protein_id}")
                continue
            
            # Report overlaps before resolution (only for domains that passed filter)
            overlaps = []
            for i, d1 in enumerate(filtered_domains):
                for j, d2 in enumerate(filtered_domains[i+1:], i+1):
                    if self.domains_overlap(d1, d2):
                        overlaps.append({
                            'protein_id': protein_id,
                            'domain1': f"{d1['interpro_id']}:{d1['domain_start']}-{d1['domain_end']}",
                            'domain2': f"{d2['interpro_id']}:{d2['domain_start']}-{d2['domain_end']}",
                            'signature1': d1['signature_name'],
                            'signature2': d2['signature_name'],
                            'domain1_length': d1['domain_end'] - d1['domain_start'] + 1,
                            'domain2_length': d2['domain_end'] - d2['domain_start'] + 1
                        })
            overlap_report.extend(overlaps)
            
            # Find optimal domain set from filtered domains
            optimal_domains = self.find_optimal_domain_set(filtered_domains)
            
            # Extract sequences for optimal domains
            full_sequence = self.all_sequences[protein_id]
            
            for i, domain in enumerate(optimal_domains):
                # Extract domain sequence (1-based to 0-based)
                domain_seq = full_sequence[domain['domain_start']-1:domain['domain_end']]
                
                # Create domain ID
                domain_id = f"{protein_id}_{domain['interpro_id']}_{i+1}_{domain['domain_start']}-{domain['domain_end']}"
                
                # Create description
                description = (f"protein_name={domain['protein_name']} "
                             f"organism={domain['organism']} "
                             f"domain_pos={domain['domain_start']}-{domain['domain_end']} "
                             f"signature={domain['signature_name']} "
                             f"interpro_id={domain['interpro_id']} "
                             f"length={len(domain_seq)}")
                
                # Add optional fields
                if 'tax_id' in domain and domain['tax_id']:
                    description += f" tax_id={domain['tax_id']}"
                if 'interpro_domain_name' in domain and domain['interpro_domain_name']:
                    description += f" interpro_domain_name={domain['interpro_domain_name']}"
                
                domain_record = SeqRecord(
                    Seq(domain_seq),
                    id=domain_id,
                    description=description
                )
                
                all_extracted_domains.append(domain_record)
                self.extraction_stats[domain['interpro_id']] += 1
        
        # Log filtering statistics
        logger.info(f"Domain filtering stats:")
        logger.info(f"  Total domains found: {filtered_stats['total_domains']}")
        logger.info(f"  Domains filtered out (< {min_domain_length} AA): {filtered_stats['filtered_out']}")
        logger.info(f"  Domains kept for processing: {filtered_stats['kept_domains']}")
        
        return all_extracted_domains, overlap_report

    def generate_reports(self, overlap_report):
        """Generate summary reports"""
        
        # Overlap report
        if overlap_report:
            overlap_file = self.reports_dir / "overlap_report.csv"
            with open(overlap_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=overlap_report[0].keys())
                writer.writeheader()
                writer.writerows(overlap_report)
            logger.info(f"Overlap report saved to {overlap_file}")
        
        # Extraction statistics
        stats_data = {
            'extraction_stats': dict(self.extraction_stats),
            'total_domains_extracted': sum(self.extraction_stats.values()),
            'total_proteins_processed': len(self.protein_annotations),
            'total_sequences_loaded': len(self.all_sequences),
            'overlapping_cases': len(overlap_report)
        }
        
        stats_file = self.reports_dir / "extraction_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats_data, f, indent=2)
        
        # Summary report
        summary_file = self.reports_dir / "summary_report.txt"
        with open(summary_file, 'w') as f:
            f.write("Domain Extraction Summary Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Total proteins processed: {len(self.protein_annotations)}\n")
            f.write(f"Total sequences loaded: {len(self.all_sequences)}\n")
            f.write(f"Total domains extracted: {sum(self.extraction_stats.values())}\n")
            f.write(f"Overlapping domain cases resolved: {len(set(o['protein_id'] for o in overlap_report))}\n\n")
            
            f.write("Domains extracted per InterPro ID:\n")
            for interpro_id, count in sorted(self.extraction_stats.items()):
                f.write(f"  {interpro_id}: {count} domains\n")
        
        logger.info(f"Summary report saved to {summary_file}")

    def run(self, min_domain_length=0):
        """Run the domain extraction process"""
        logger.info("Starting domain extraction...")
        logger.info(f"Minimum domain length: {min_domain_length} amino acids")
        
        # Collect all data first
        self.collect_all_data()
        
        logger.info(f"Loaded {len(self.all_sequences)} sequences")
        logger.info(f"Processing {len(self.protein_annotations)} proteins")
        
        # Extract optimal domains with length filtering
        extracted_domains, overlap_report = self.extract_domains(min_domain_length)
        
        # Write output FASTA
        output_fasta = self.domains_dir / "all_optimal_domains.fasta"
        if extracted_domains:
            SeqIO.write(extracted_domains, output_fasta, "fasta")
            logger.info(f"Extracted {len(extracted_domains)} optimal domains")
        
        # Generate reports
        self.generate_reports(overlap_report)
        
        logger.info("Domain extraction completed!")
        logger.info(f"Results saved to: {self.output_dir}")
        logger.info(f"Optimal domains: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Extract optimal non-overlapping protein domains")
    parser.add_argument("input_dir", help="Directory containing InterPro subfolders")
    parser.add_argument("output_dir", help="Output directory for extracted domains")
    parser.add_argument("--min-domain-length", type=int, default=0,
                       help="Minimum domain length in amino acids (default: 0, no filtering)")
    
    args = parser.parse_args()
    
    extractor = DomainExtractor(args.input_dir, args.output_dir)
    extractor.run(args.min_domain_length)

if __name__ == "__main__":
    main()