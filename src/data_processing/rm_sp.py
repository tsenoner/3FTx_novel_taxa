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


def parse_signalp_results(results_file, threshold):
    cleavages = {}
    with open(results_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            seq_id, site_str, score_str = parts[:3]
            try:
                site = int(site_str)
                score = float(score_str)
            except ValueError:
                continue
            if score >= threshold:
                cleavages[seq_id] = site
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

    cleavages = parse_signalp_results(args.signalp, args.threshold)
    if not cleavages:
        print("No signal peptides detected over threshold; outputs will be unchanged.")

    write_mature_fasta(args.fasta, cleavages, args.out_fasta)

    pdb_files = [os.path.join(args.pdb_dir, f)
                 for f in os.listdir(args.pdb_dir)
                 if f.endswith('.pdb')]
    for pdb_path in pdb_files:
        write_mature_pdbs_for_file(pdb_path, cleavages, args.out_pdb_dir)

if __name__ == '__main__':
    main()
