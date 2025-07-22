#!/usr/bin/env python3
"""
Use case: Extract domains from monomer and multimer proteins (to compare monomeric repeats within multimers), then reduce redundancy in those domain sequences via aggressive MMseqs2 clustering for downstream phylogenetic analysis (ExaBayes tree).
"""
import argparse
import subprocess
import sys
import shutil
from pathlib import Path


def check_mmseqs():
    if not shutil.which("mmseqs"):
        sys.exit("ERROR: mmseqs2 binary not found in $PATH. Please install mmseqs2 first.")


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.check_call(cmd)


def main():
    parser = argparse.ArgumentParser(
        description="MMseqs2 aggressive clustering for domain redundancy reduction"
    )
    parser.add_argument("input_fasta", type=Path,
                        help="FASTA of extracted domains from InterPro (non-overlapping)")
    parser.add_argument("output_dir", type=Path,
                        help="Directory for clustering outputs and final representatives")
    parser.add_argument("--min-seq-id", type=float, default=0.30,
                        help="Min sequence identity for clustering (default: 0.30)")
    parser.add_argument("--cov", type=float, default=0.80,
                        help="Coverage threshold fraction (default: 0.80)")
    parser.add_argument("--cov-mode", type=int, choices=[0,1,2], default=2,
                        help="Coverage mode: 0=both,1=target,2=query (default: 2)")
    parser.add_argument("--cluster-mode", type=int, choices=[0,1,2,3], default=1,
                        help="Cluster mode: 0=set-cover,1=connected,2/3=length-based (default: 1 for connected components)")
    parser.add_argument("--reassign", action="store_true",
                        help="Enable --cluster-reassign to correct cascaded clustering errors")
    parser.add_argument("-s", dest="sensitivity", type=float, default=9,
                        help="Prefilter sensitivity (1=faster,7.5=sensitive; default: 9)")
    parser.add_argument("--threads", type=int, default=6,
                        help="Number of threads to use (default: 6)")
    args = parser.parse_args()

    check_mmseqs()

    # Prepare directories
    db_dir = args.output_dir / "db"
    tmp_dir = args.output_dir / "tmp"
    rep_fa = args.output_dir / "representatives.fasta"
    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    db_dir.mkdir(parents=True)
    tmp_dir.mkdir(parents=True)

    # 1) Create DB
    db_in = db_dir / "inputDB"
    run(["mmseqs", "createdb", str(args.input_fasta), str(db_in)])

    # 2) Aggressive clustering
    db_clu = db_dir / "clusters"
    cmd = [
        "mmseqs", "cluster",
        str(db_in), str(db_clu), str(tmp_dir),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.cov),
        "--cov-mode", str(args.cov_mode),
        "--cluster-mode", str(args.cluster_mode),
        "-s", str(args.sensitivity),
        "--split-mode", "2",
        "--threads", str(args.threads)
    ]
    if args.reassign:
        cmd.append("--cluster-reassign")
    run(cmd)

    # 3) Extract representatives
    db_rep = db_dir / "repDB"
    run(["mmseqs", "result2repseq", str(db_in), str(db_clu), str(db_rep)])

    # 4) Convert to FASTA
    run(["mmseqs", "convert2fasta", str(db_rep), str(rep_fa)])

    print(f"\n✓ DONE — final representatives at {rep_fa}")

if __name__ == "__main__":
    main()