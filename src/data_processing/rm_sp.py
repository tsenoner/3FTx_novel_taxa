#!/usr/bin/env python3
"""
Removes signal peptides from a directory of PDB files based on SignalP6 predictions,
outputs mature FASTA and per-chain mature PDBs with unique filenames.

Usage:
    python rm_sp.py \
        --pdb_dir input_pdb_directory \
        --signalp signalp_results.txt \
        --fasta input_full_length.fasta \
        --out_pdb_dir mature_pdb_output_dir \
        --out_fasta mature_sequences.fasta \
        [--threshold 0.8]

SignalP6 results (none format) must have whitespace-separated:
    sequence_id  cleavage_site  prediction_score

Requires:
    biopython (pip install biopython)
"""
import argparse
import os
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class ChainMatureSelect(Select):
    def __init__(self, chain_id, cleavage_map):
        self.chain_id = chain_id
        self.cleavage_map = cleavage_map

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        if self.chain_id in self.cleavage_map:
            if residue.id[1] <= self.cleavage_map[self.chain_id]:
                return False
        return True


def parse_signalp_results_gff3(results_file, threshold):
    """
    Parses SignalP 6.0 GFF3 output to extract cleavage sites with scores over the threshold.
    Returns a dict: {sequence_id: cleavage_site}
    """
    cleavages = {}
    with open(results_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            seq_id = parts[0]
            feature_type = parts[2]
            end = int(parts[4])  # SP ends at this position (inclusive, 1-based)
            score = float(parts[5])
            if feature_type == "signal_peptide" and score >= threshold:
                cleavages[seq_id] = end
    return cleavages


def write_mature_fasta(input_fasta, cleavage_map, output_fasta):
    records = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        seqid = rec.id
        if seqid in cleavage_map:
            cut = cleavage_map[seqid]
            seq = rec.seq[cut:]
        else:
            seq = rec.seq
        records.append(SeqRecord(seq, id=rec.id, description=rec.description))
    SeqIO.write(records, output_fasta, 'fasta')
    print(f"Written mature FASTA to {output_fasta}")


def write_mature_pdbs_for_file(input_pdb, cleavage_map, output_dir):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(os.path.basename(input_pdb), input_pdb)
    io = PDBIO()
    io.set_structure(struct)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    base = os.path.splitext(os.path.basename(input_pdb))[0]
    for chain in struct.get_chains():
        cid = chain.id
        select = ChainMatureSelect(cid, cleavage_map)
        out_file = os.path.join(output_dir, f"{base}_{cid}_mature.pdb")
        io.save(out_file, select=select)
        print(f"Written mature PDB: {out_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Remove signal peptides from a directory of PDBs and produce mature outputs."
    )
    parser.add_argument(
        "--pdb_dir", required=True,
        help="Directory containing input PDB files (*.pdb)"
    )
    parser.add_argument(
        "--signalp", required=True,
        help="SignalP6 results file (format none)"
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Input FASTA file with full-length sequences"
    )
    parser.add_argument(
        "--out_pdb_dir", required=True,
        help="Directory to write mature PDB files"
    )
    parser.add_argument(
        "--out_fasta", required=True,
        help="Output FASTA file for mature sequences"
    )
    parser.add_argument(
        "--threshold", type=float, default=0.8,
        help="Prediction score threshold (default: 0.8)"
    )
    args = parser.parse_args()

    cleavages = parse_signalp_results_gff3(args.signalp, args.threshold)
    if not cleavages:
        print("No signal peptides detected over threshold; outputs will be unchanged.")

    write_mature_fasta(args.fasta, cleavages, args.out_fasta)

    pdb_files = [os.path.join(args.pdb_dir, f)
                 for f in os.listdir(args.pdb_dir)
                 if f.endswith('.pdb')]
    for pdb_path in pdb_files:
        write_mature_pdbs_for_file(pdb_path, cleavages, args.out_pdb_dir)

def extract_fasta_by_uniprot_ids(input_fasta, id_file, output_fasta):
    """
    Extracts sequences from a FASTA file that match UniProt IDs listed in a text file.

    Parameters:
    - input_fasta (str): Path to the source FASTA file.
    - id_file (str): Path to the text file containing UniProt IDs (one per line).
    - output_fasta (str): Path to the output FASTA file to write matched sequences.
    """
    # Load UniProt IDs from file
    with open(id_file) as f:
        missing_ids = set(line.strip() for line in f if line.strip())

    matched = 0
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Extract ID from FASTA header (e.g. >sp|P12345|... or >P12345 ...)
            header_id = record.id.split('|')[1] if '|' in record.id else record.id
            if header_id in missing_ids:
                SeqIO.write(record, out_f, "fasta")
                matched += 1

    print(f"✓ Extracted {matched} sequences to {output_fasta}")


def count_fasta_sequences(filepath):
    """
    Counts the number of sequences in a FASTA file.

    Parameters:
    - filepath (str): Path to the FASTA file.

    Returns:
    - int: Number of sequences.
    """
    count = sum(1 for _ in SeqIO.parse(filepath, "fasta"))
    print(f"{filepath}: {count} sequences")
    return count


if __name__ == '__main__':
    # usage:
    # python rm_sp.py
    # --pdb_dir ../../data/raw/pdb/
    # --signalp ../../data/interm/sp6/output.gff3
    # --fasta ../../data/interm/backup/merged_sanitized.fasta
    # --out_pdb_dir ../../data/interm/pdb
    # --out_fasta ../../data/interm/pdb
    # --threshold 0.8
    #main()

    extract_fasta_by_uniprot_ids(
        input_fasta="../../data/interm/pdb/mature_sequences.fasta",
        id_file="../../data/raw/pdb/missing.txt",
        output_fasta="../../data/interm/pdb/mature_missing_sequences.fasta"
    )

    # Compare counts
    # processed_count = count_fasta_sequences("../../data/interm/backup/merged_sanitized.fasta")
    #
    # processed_count = count_fasta_sequences("../../data/interm/pdb/mature_missing_sequences.fasta")
    # merged_count = count_fasta_sequences("../../data/interm/pdb/mature_sequences.fasta")
    #
    # missing = merged_count - processed_count
    # print(f"\nMissing sequences: {missing}")
