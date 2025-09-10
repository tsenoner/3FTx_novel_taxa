#!/usr/bin/env python3
"""
cluster_group_assign_with_taxonomy.py

Updated: normalization + robust taxonomy lookup + domain->parent handling.
Behavior: find_last_common_ancestor checks the requested min level first,
and if no consensus or missing, steps up to more general levels until domain.
"""
import argparse
import csv
import logging
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import math

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ----------------- CONFIGURATION -----------------

# Group patterns for protein classification
GROUP_PATTERNS = {
    "manually_annotated": [r"\bCentipede3FTx\b"],
    "3FTx": [r"\b3ftx\b", r"\b3\-ftx\b", r"\bthree[- ]?ftx\b", r"\b3FTx\b",
             r"\bThree-finger toxin\b", r"\bThree finger toxins\b", r"\bThree finger toxin\b"],
    "PMF": [r"\bpmf\b", r"\bPNF\b"],
    "Ly6": [r"\bly6\b", r"\bly\-6\b", r"\blysix\b", r"\bLy-6\b", r"\bLY6\b", r"\bUPAR/Ly6\b"],
    "Quiver": [r"\bquiver\b", r"\bqvr\b"],
    "Scoloptoxin": [r"\bscoloptoxin\b"],
    "SPF": [r"\bSodefrin\b"],
}

# Specificity order: most specific first
SPECIFICITY_ORDER = [
    "SPF", "Scoloptoxin", "PMF", "Quiver",
    "manually_annotated", "3FTx", "Ly6"
]

# InterPro domain to group mapping
INTERPRO_GROUP_MAP = {
    # fill as needed
}

# Taxonomic hierarchy (general → specific)
TAXONOMIC_LEVELS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

# Phylum to kingdom/domain mapping
PHYLUM_TO_HIGHER_TAXA = {
    # Animals
    "Arthropoda": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Chordata": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Mollusca": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Cnidaria": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Nematoda": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Annelida": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Echinodermata": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Platyhelminthes": {"kingdom": "Animalia", "domain": "Eukaryota"},
    "Porifera": {"kingdom": "Animalia", "domain": "Eukaryota"},

    # Plants
    "Streptophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Chlorophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Tracheophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},
    "Magnoliophyta": {"kingdom": "Plantae", "domain": "Eukaryota"},

    # Fungi
    "Ascomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Basidiomycota": {"kingdom": "Fungi", "domain": "Eukaryota"},
    "Mucoromycota": {"kingdom": "Fungi", "domain": "Eukaryota"},

    # Bacteria
    "Proteobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Firmicutes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Actinobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Bacteroidetes": {"kingdom": "Bacteria", "domain": "Bacteria"},
    "Cyanobacteria": {"kingdom": "Bacteria", "domain": "Bacteria"},

    # Archaea
    "Euryarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Crenarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
    "Thaumarchaeota": {"kingdom": "Archaea", "domain": "Archaea"},
}

# Centipede common taxonomy
CENTIPEDE_COMMON_TAXONOMY = {
    "domain": "Eukaryota",
    "kingdom": "Animalia",
    "phylum": "Arthropoda",
    "class": "Chilopoda",
    "order": None,
    "family": None,
    "genus": None
}

# Compiled regex patterns
SPECIFICITY_RANK = {name: i for i, name in enumerate(SPECIFICITY_ORDER)}
_COMPILED_GROUPS = [(g, [re.compile(p, re.IGNORECASE) for p in pats])
                    for g, pats in GROUP_PATTERNS.items()]
_re_protein_name = re.compile(
    r"protein_name=(.*?)\s+(?:organism=|domain_pos=|signature=|interpro_id=|length=|$)", re.IGNORECASE
)
_re_interpro_id = re.compile(r"(IPR\d{6})", re.IGNORECASE)
_re_length = re.compile(r"length=(\d+)", re.IGNORECASE)

# ----------------- ID normalization helpers -----------------
_DOMAIN_SUFFIX_RE = re.compile(r'(_IPR\d{6}(?:_[0-9]+)?(?:_[0-9\-]+)?)$')
_PIPE_RE = re.compile(r'^[ps]p\|([^|]+)\|')
_ISOFORM_RE = re.compile(r'^([A-NR-Z0-9]+)(?:-\d+)$')
_VERSION_RE = re.compile(r'^([A-NR-Z0-9]+)(?:\.\d+)$')

def normalize_protein_id(raw_id: str) -> str:
    if not raw_id:
        return raw_id
    s = raw_id.split(None, 1)[0].strip()

    if "Centipede" in s:
        return s

    m = _PIPE_RE.match(s)
    if m:
        return m.group(1)

    s = _DOMAIN_SUFFIX_RE.sub('', s)

    m = _ISOFORM_RE.match(s)
    if m:
        s = m.group(1)
    m = _VERSION_RE.match(s)
    if m:
        s = m.group(1)

    if "_" in s and re.match(r'^[A-Z0-9]{5,}', s):
        s = s.split("_", 1)[0]

    return s.strip()

# ----------------- CORE FUNCTIONS -----------------

class TaxonomyParser:
    @staticmethod
    def parse_csv(csv_path: Path) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
        taxonomy_data = {}

        if not csv_path.exists():
            raise FileNotFoundError(f"Taxonomy CSV not found: {csv_path}")

        with csv_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            csv_fieldnames = [c for c in (reader.fieldnames or [])]
            csv_cols_lower = {c.lower(): c for c in csv_fieldnames}

            for row in reader:
                protein_id = (row.get("protein_id") or row.get("Protein_id") or row.get("Protein_ID") or "").strip()
                if not protein_id:
                    continue

                tax_info = {}
                for level in TAXONOMIC_LEVELS:
                    if level in ["domain", "kingdom"]:
                        continue
                    col = csv_cols_lower.get(level)
                    if col:
                        tax_info[level] = (row.get(col) or "").strip() or None
                    else:
                        alt = None
                        for candidate in (level, level.capitalize(), level + "_name", level + "Name"):
                            colalt = csv_cols_lower.get(candidate.lower())
                            if colalt:
                                alt = colalt
                                break
                        tax_info[level] = (row.get(alt) or "").strip() or None if alt else None

                phylum = (tax_info.get("phylum") or "").strip()
                if phylum:
                    mapping = PHYLUM_TO_HIGHER_TAXA.get(phylum) or PHYLUM_TO_HIGHER_TAXA.get(phylum.title())
                    if mapping:
                        tax_info.update(mapping)
                    else:
                        tax_info["kingdom"] = None
                        tax_info["domain"] = None
                else:
                    tax_info["kingdom"] = None
                    tax_info["domain"] = None

                taxonomy_data[protein_id] = tax_info

        taxonomy_by_norm: Dict[str, Dict[str, str]] = {}
        collisions = 0
        for pid, tinfo in taxonomy_data.items():
            n = normalize_protein_id(pid)
            if n in taxonomy_by_norm and taxonomy_by_norm[n] != tinfo:
                collisions += 1
                logging.debug(f"Normalized taxonomy key collision for {n}: existing != new for {pid}")
            else:
                taxonomy_by_norm[n] = tinfo

        logging.info(f"Loaded taxonomy data for {len(taxonomy_data)} proteins ({len(taxonomy_by_norm)} normalized keys, {collisions} collisions)")
        return taxonomy_data, taxonomy_by_norm

    @staticmethod
    def find_last_common_ancestor(member_ids: List[str],
                                  taxonomy_data: Dict[str, Dict[str, str]],
                                  taxonomy_by_norm: Optional[Dict[str, Dict[str, str]]] = None,
                                  min_level: str = "family") -> str:
        if not member_ids:
            return "unknown"

        if min_level not in TAXONOMIC_LEVELS:
            logging.warning(f"Invalid min_level '{min_level}', using 'family'")
            min_level = "family"

        min_index = TAXONOMIC_LEVELS.index(min_level)

        # Build levels_to_check starting at min_level and moving to more general
        # e.g. min_level='genus' -> ['genus','family','order','class','phylum','kingdom','domain']
        levels_to_check = [TAXONOMIC_LEVELS[i] for i in range(min_index, -1, -1)]

        # Normalize member IDs and separate Centipede
        uniprot_ids = []
        has_centipedes = False
        for member_id in member_ids:
            seq_id = member_id.split()[0]
            norm_id = normalize_protein_id(seq_id)
            if "Centipede" in norm_id:
                has_centipedes = True
            else:
                uniprot_ids.append(norm_id)

        taxonomies = []
        if has_centipedes:
            taxonomies.append(CENTIPEDE_COMMON_TAXONOMY)
            logging.debug("Added Centipede taxonomy")

        missing_count = 0

        if taxonomy_by_norm is None:
            taxonomy_by_norm = {}
            for pid, tinfo in taxonomy_data.items():
                taxonomy_by_norm.setdefault(normalize_protein_id(pid), tinfo)

        seen_parents = set()
        for uniprot_id in uniprot_ids:
            if uniprot_id in seen_parents:
                continue
            seen_parents.add(uniprot_id)

            tinfo = taxonomy_by_norm.get(uniprot_id)
            if tinfo is None:
                alt = re.sub(r'(-\d+)$', '', uniprot_id)
                alt = normalize_protein_id(alt)
                tinfo = taxonomy_by_norm.get(alt)
            if tinfo is None:
                missing_count += 1
                logging.debug(f"Missing taxonomy for uid (cluster lookup): {uniprot_id}")
            else:
                taxonomies.append(tinfo)

        if missing_count > 0:
            logging.info(f"Missing taxonomy for {missing_count} UniProt IDs in cluster")

        if not taxonomies:
            return "unknown"

        # Check levels in order from most specific (min_level) to most general (domain)
        for level in levels_to_check:
            values = set()
            all_have_value = True

            for tax_info in taxonomies:
                value = tax_info.get(level)
                if not value:
                    all_have_value = False
                    break
                values.add(value.strip())

            if all_have_value and len(values) == 1:
                result = values.pop()
                logging.debug(f"LCA at {level}: {result}")
                return f"{level}:{result}"

        return "unknown"


class ProteinAnalyzer:
    @staticmethod
    def extract_uniprot_id(identifier: str) -> str:
        return normalize_protein_id(identifier)

    @staticmethod
    def extract_length_from_header(header: str) -> Optional[int]:
        match = _re_length.search(header)
        return int(match.group(1)) if match else None

    @staticmethod
    def extract_header_fields(full_header: str) -> Tuple[Optional[str], Optional[str]]:
        protein_name = None
        interpro_id = None

        pn_match = _re_protein_name.search(full_header)
        if pn_match:
            protein_name = pn_match.group(1).strip().strip('"').strip("'")

        ip_match = _re_interpro_id.search(full_header)
        if ip_match:
            interpro_id = ip_match.group(1).upper()

        return protein_name, interpro_id

    @staticmethod
    def count_uniprot_domains(headers_by_id: Dict[str, str]) -> Dict[str, int]:
        uniprot_counts = Counter()

        for seq_id in headers_by_id.keys():
            parent = normalize_protein_id(seq_id)
            if "Centipede" not in parent:
                uniprot_counts[parent] += 1

        for seq_id in headers_by_id.keys():
            parent = normalize_protein_id(seq_id)
            if "Centipede" in parent:
                uniprot_counts[parent] = 1

        multimeric_count = sum(1 for count in uniprot_counts.values() if count > 1)
        logging.info(f"Found {multimeric_count} proteins with multiple domains (normalized)")

        return dict(uniprot_counts)

    @staticmethod
    def determine_oligomeric_state(member_ids: List[str], uniprot_counts: Dict[str, int]) -> str:
        domain_types = []
        for member_id in member_ids:
            seq_id = member_id.split()[0]
            parent = normalize_protein_id(seq_id)
            domain_count = uniprot_counts.get(parent, 1)
            domain_types.append("monomeric" if domain_count == 1 else "multimeric")

        unique_types = set(domain_types)
        if len(unique_types) == 1:
            return "monomer" if "monomeric" in unique_types else "multimer"
        else:
            return "mixture"


class GroupClassifier:
    @staticmethod
    def find_all_groups_for_header(full_header: str) -> List[str]:
        protein_name, interpro_id = ProteinAnalyzer.extract_header_fields(full_header)
        groups = []

        if interpro_id and interpro_id in INTERPRO_GROUP_MAP:
            group = INTERPRO_GROUP_MAP[interpro_id]
            groups.append(group)
            logging.debug(f"InterPro {interpro_id} -> {group}")

        id_token = full_header.split(None, 1)[0]
        targets = []
        if protein_name:
            targets.append(protein_name)
        targets.extend([
            full_header,
            full_header.replace("_", " "),
            id_token,
            id_token.replace("_", " "),
        ])

        for target in targets:
            for group_name, compiled_patterns in _COMPILED_GROUPS:
                for pattern in compiled_patterns:
                    if pattern.search(target) and group_name not in groups:
                        groups.append(group_name)
                        logging.debug(f"Pattern match: {group_name} in '{target[:50]}'")

        return groups

    @staticmethod
    def find_best_group_for_header(full_header: str) -> Optional[str]:
        matched_groups = GroupClassifier.find_all_groups_for_header(full_header)
        if not matched_groups:
            return None
        if len(matched_groups) == 1:
            return matched_groups[0]
        best_group = min(matched_groups, key=lambda g: SPECIFICITY_RANK.get(g, 9999))
        logging.debug(f"Multiple matches {matched_groups}, chose {best_group}")
        return best_group

    @staticmethod
    def decide_cluster_label(member_headers: List[str], rep_id: str = "") -> str:
        n = len(member_headers)
        groups = [GroupClassifier.find_best_group_for_header(h) for h in member_headers]
        counts = Counter(groups)
        non_none_groups = [g for g in groups if g is not None]
        distinct_groups = set(non_none_groups)

        logging.debug(f"Cluster {rep_id}: {dict(counts)}")
        if len(distinct_groups) == 1 and counts[None] == 0:
            return next(iter(distinct_groups))
        if n > 10 and distinct_groups:
            best_group, best_count = counts.most_common(1)[0]
            if best_group is not None and (best_count / n) >= 0.90:
                return f"mostly_{best_group}"
        if len(distinct_groups) > 1 or (len(distinct_groups) == 1 and counts[None] > 0):
            return "mixture"
        return "none"


class LengthBinner:
    def __init__(self, num_bins: int = 15):
        self.num_bins = num_bins
        self.bin_ranges = []

    def create_bins(self, average_lengths: List[float]) -> None:
        if not average_lengths:
            self.bin_ranges = []
            return
        min_length = min(average_lengths)
        max_length = max(average_lengths)
        if min_length == max_length:
            self.bin_ranges = [(min_length, max_length)]
            return
        bin_size = (max_length - min_length) / self.num_bins
        self.bin_ranges = []
        for i in range(self.num_bins):
            start = min_length + i * bin_size
            end = min_length + (i + 1) * bin_size if i < self.num_bins - 1 else max_length
            self.bin_ranges.append((start, end))
        logging.info(f"Created {self.num_bins} bins: {min_length:.1f}-{max_length:.1f}")

    def get_bin_label(self, avg_length: Optional[float]) -> str:
        if avg_length is None or not self.bin_ranges:
            return "unknown"
        for start, end in self.bin_ranges:
            if start <= avg_length <= end:
                return f"{int(start)}-{int(end)}"
        return "unknown"


class ClusterProcessor:
    @staticmethod
    def parse_fasta_with_lengths(fasta_path: Path, id_mode: str = "first_token") -> Tuple[Dict[str, str], Dict[str, int]]:
        headers = {}
        lengths = {}
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA not found: {fasta_path}")
        with fasta_path.open("r", encoding="utf-8") as fh:
            seqid = None
            seq_lines = []
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seqid is not None:
                        seq = "".join(seq_lines)
                        lengths[seqid] = len(seq)
                        seq_lines = []
                    header_text = line[1:].rstrip()
                    seqid = header_text if id_mode == "full" else header_text.split(None, 1)[0]
                    headers[seqid] = header_text
                else:
                    seq_lines.append(line.strip())
            if seqid is not None:
                seq = "".join(seq_lines)
                lengths[seqid] = len(seq)
        logging.info(f"Parsed {len(headers)} sequences from {fasta_path}")
        return headers, lengths

    @staticmethod
    def parse_mmseqs_tsv(tsv_path: Path) -> Dict[str, List[str]]:
        mapping = defaultdict(list)
        if not tsv_path.exists():
            raise FileNotFoundError(f"TSV not found: {tsv_path}")
        with tsv_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split("\t")
                if len(cols) >= 2:
                    rep = cols[0]
                    member = cols[1]
                    mapping[rep].append(member)
        cleaned = {}
        for rep, members in mapping.items():
            seen = set()
            unique_members = []
            for m in members:
                if m not in seen:
                    seen.add(m)
                    unique_members.append(m)
            cleaned[rep] = unique_members
        logging.info(f"Parsed {len(cleaned)} clusters from TSV")
        return cleaned

    @staticmethod
    def calculate_cluster_average_length(member_headers: List[str], sequence_lengths: Dict[str, int]) -> Optional[float]:
        lengths = []
        for header in member_headers:
            length = ProteinAnalyzer.extract_length_from_header(header)
            if length is None:
                seqid = header.split(None, 1)[0]
                length = sequence_lengths.get(seqid) or sequence_lengths.get(seqid.replace("_", " "))
            if length is not None:
                lengths.append(length)
        return sum(lengths) / len(lengths) if lengths else None

    @staticmethod
    def build_cluster_member_headers(rep_to_members: Dict[str, List[str]],
                                     headers_by_id: Dict[str, str],
                                     rep_headers_by_id: Dict[str, str]) -> Dict[str, List[str]]:
        result = {}
        missing_count = 0
        for rep_id, members in rep_to_members.items():
            rep_header = rep_headers_by_id.get(rep_id) or headers_by_id.get(rep_id) or rep_id
            combined_ids = [rep_id] + [m for m in members if m != rep_id]
            member_headers = []
            for member_id in combined_ids:
                header = headers_by_id.get(member_id) or rep_headers_by_id.get(member_id) or member_id
                if header == member_id and member_id not in headers_by_id and member_id not in rep_headers_by_id:
                    missing_count += 1
                member_headers.append(header)
            result[rep_header] = member_headers
        if missing_count > 0:
            logging.warning(f"{missing_count} member IDs not found in FASTA files")
        return result


def write_output_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["identifier", "group", "oligomeric_state", "last_common_ancestor", "length"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logging.info(f"Wrote {len(rows)} rows to {path}")


def main():
    parser = argparse.ArgumentParser(description="Assign taxonomic groups to protein clusters")
    parser.add_argument("--original-fasta", "-o", required=True, type=Path,
                        help="Original FASTA file with all sequences")
    parser.add_argument("--rep-fasta", "-r", required=True, type=Path,
                        help="Representative sequences FASTA file")
    parser.add_argument("--cluster-tsv", "-c", required=True, type=Path,
                        help="MMseqs2 clustering TSV file")
    parser.add_argument("--taxonomy-csv", "-t", required=True, type=Path,
                        help="Taxonomy CSV file")
    parser.add_argument("--output-csv", "-O", required=True, type=Path,
                        help="Output CSV file")
    parser.add_argument("--id-mode", choices=["first_token", "full"], default="first_token",
                        help="How to extract sequence IDs from FASTA headers")
    parser.add_argument("--num-bins", type=int, default=15,
                        help="Number of length bins")
    parser.add_argument("--min-tax-level", default="family", choices=TAXONOMIC_LEVELS,
                        help="Minimum taxonomic level for LCA")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose logging")
    parser.add_argument("--debug-cluster",
                        help="Enable detailed debugging for specific cluster")
    parser.add_argument("--debug-taxonomy", action="store_true",
                        help="Enable detailed logging of taxonomy matching (for debugging)")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        headers_by_id, sequence_lengths = ClusterProcessor.parse_fasta_with_lengths(
            args.original_fasta, args.id_mode)
        rep_headers_by_id, rep_sequence_lengths = ClusterProcessor.parse_fasta_with_lengths(
            args.rep_fasta, args.id_mode)

        all_sequence_lengths = {**sequence_lengths, **rep_sequence_lengths}

        rep_to_members = ClusterProcessor.parse_mmseqs_tsv(args.cluster_tsv)
        taxonomy_data, taxonomy_by_norm = TaxonomyParser.parse_csv(args.taxonomy_csv)

        if getattr(args, "debug_taxonomy", False) or args.verbose:
            sample_headers = list(headers_by_id.values())
            present = 0
            missing = 0
            missing_examples = []
            for h in sample_headers[:200]:
                uid = ProteinAnalyzer.extract_uniprot_id(h.split()[0])
                if uid in taxonomy_by_norm:
                    present += 1
                else:
                    missing += 1
                    if len(missing_examples) < 10:
                        missing_examples.append((h.split()[0], uid))
            logging.debug(f"Taxonomy debug: {present} present, {missing} missing (checked up to 200 headers)")
            if missing_examples:
                logging.debug(f"Missing examples (raw,normalized): {missing_examples}")

        uniprot_counts = ProteinAnalyzer.count_uniprot_domains(headers_by_id)
        clusters = ClusterProcessor.build_cluster_member_headers(
            rep_to_members, headers_by_id, rep_headers_by_id)

        cluster_avg_lengths = []
        cluster_length_data = {}

        for rep_header, member_headers in clusters.items():
            rep_id = rep_header.split()[0]
            avg_length = ClusterProcessor.calculate_cluster_average_length(
                member_headers, all_sequence_lengths)
            if avg_length is not None:
                cluster_avg_lengths.append(avg_length)
                cluster_length_data[rep_id] = avg_length

        reps_in_tsv = set(rep_to_members.keys())
        reps_in_fasta = set(rep_headers_by_id.keys())
        singletons = reps_in_fasta - reps_in_tsv

        for rep_id in singletons:
            rep_header = rep_headers_by_id.get(rep_id, rep_id)
            avg_length = ClusterProcessor.calculate_cluster_average_length(
                [rep_header], all_sequence_lengths)
            if avg_length is not None:
                cluster_avg_lengths.append(avg_length)
                cluster_length_data[rep_id] = avg_length

        length_binner = LengthBinner(args.num_bins)
        length_binner.create_bins(cluster_avg_lengths)

        results = []
        stats = Counter()
        oligomeric_stats = Counter()
        taxonomy_stats = Counter()
        length_stats = Counter()

        logging.info(f"Processing {len(clusters)} clusters and {len(singletons)} singletons")

        for rep_header, member_headers in clusters.items():
            rep_id = rep_header.split()[0]
            if args.debug_cluster and args.debug_cluster == rep_id:
                logging.info(f"=== DEBUG CLUSTER {rep_id} ===")
                original_level = logging.getLogger().level
                logging.getLogger().setLevel(logging.DEBUG)
                logging.debug(f"Cluster has {len(member_headers)} members:")
                for i, header in enumerate(member_headers):
                    logging.debug(f"  {i+1}. {header[:60]}...")

            group_label = GroupClassifier.decide_cluster_label(member_headers, rep_id)
            oligomeric_state = ProteinAnalyzer.determine_oligomeric_state(member_headers, uniprot_counts)
            lca = TaxonomyParser.find_last_common_ancestor(member_headers, taxonomy_data, taxonomy_by_norm, args.min_tax_level)

            avg_length = cluster_length_data.get(rep_id)
            length_bin = length_binner.get_bin_label(avg_length)

            results.append({
                "identifier": rep_id,
                "group": group_label,
                "oligomeric_state": oligomeric_state,
                "last_common_ancestor": lca,
                "length": length_bin
            })

            stats[group_label] += 1
            oligomeric_stats[oligomeric_state] += 1
            taxonomy_stats[lca] += 1
            length_stats[length_bin] += 1

            if args.debug_cluster and args.debug_cluster == rep_id:
                logging.info(f"Final classification: group={group_label}, oligomeric={oligomeric_state}, LCA={lca}")
                logging.info(f"=== END DEBUG CLUSTER {rep_id} ===")
                logging.getLogger().setLevel(original_level)

        if singletons:
            logging.info(f"Processing {len(singletons)} singleton sequences")

        for rep_id in singletons:
            rep_header = rep_headers_by_id.get(rep_id, rep_id)
            group_label = GroupClassifier.decide_cluster_label([rep_header], rep_id)
            oligomeric_state = ProteinAnalyzer.determine_oligomeric_state([rep_header], uniprot_counts)
            lca = TaxonomyParser.find_last_common_ancestor([rep_header], taxonomy_data, taxonomy_by_norm, args.min_tax_level)

            avg_length = cluster_length_data.get(rep_id)
            length_bin = length_binner.get_bin_label(avg_length)

            results.append({
                "identifier": rep_id,
                "group": group_label,
                "oligomeric_state": oligomeric_state,
                "last_common_ancestor": lca,
                "length": length_bin
            })

            stats[group_label] += 1
            oligomeric_stats[oligomeric_state] += 1
            taxonomy_stats[lca] += 1
            length_stats[length_bin] += 1

        write_output_csv(args.output_csv, results)

        logging.info("\n=== FINAL STATISTICS ===")
        logging.info("Group Classification:")
        for group, count in stats.most_common():
            logging.info(f"  {group}: {count}")

        logging.info("Oligomeric State:")
        for state, count in oligomeric_stats.most_common():
            logging.info(f"  {state}: {count}")

        logging.info("Taxonomic Classification (top 20):")
        for taxon, count in taxonomy_stats.most_common()[:20]:
            logging.info(f"  {taxon}: {count}")

        logging.info("Length Distribution:")
        for length_bin, count in length_stats.most_common():
            logging.info(f"  {length_bin}: {count}")

        logging.info(f"\nProcessed {len(results)} total entries successfully")

        if getattr(args, "debug_taxonomy", False) or args.verbose:
            missing_uids = Counter()
            for header in list(headers_by_id.values())[:5000]:
                seqid = header.split()[0]
                uid_norm = normalize_protein_id(seqid)
                if uid_norm not in taxonomy_by_norm:
                    missing_uids[uid_norm] += 1
            logging.debug(f"Top missing normalized taxonomy IDs (sample 50): {missing_uids.most_common(50)[:20]}")

    except Exception as e:
        logging.error(f"Error during processing: {e}")
        raise


if __name__ == "__main__":
    main()
