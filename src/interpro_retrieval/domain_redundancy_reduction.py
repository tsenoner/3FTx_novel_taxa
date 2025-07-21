#!/usr/bin/env python3
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
        description="MMseqs2-based hierarchical clustering and representative extraction"
    )
    parser.add_argument("input_fasta", type=Path,
                        help="Input FASTA of domain sequences")
    parser.add_argument("output_dir", type=Path,
                        help="Directory to write MMseqs2 DBs & final reps")
    parser.add_argument("--min-seq-id", type=float, default=0.40,
                        help="Minimum sequence identity for clustering (0–1.0, default: 0.40)")
    parser.add_argument("--cov", type=float, default=0.40,
                        help="Coverage fraction for clustering (0–1.0, default: 0.40)")
    parser.add_argument("--cov-mode", type=int, choices=[0,1,2], default=2,
                        help="Coverage mode (0: both, 1: target, 2: query; default: 2)")
    parser.add_argument("--cluster-mode", type=int, choices=[0,1,2,3], default=2,
                        help="Clustering algorithm (0: set-cover, 1: connected, 2/3: length-based; default: 2)")
    parser.add_argument("--kmer-per-seq", type=int, default=15,
                        help="Number of k-mers per sequence for clustering sensitivity (default: 15)")
    parser.add_argument("--threads", type=int, default=6,
                        help="Number of threads to use")
    args = parser.parse_args()

    check_mmseqs()

    # Prepare directories
    db_dir  = args.output_dir / "db"
    tmp_dir = args.output_dir / "tmp"
    rep_fa  = args.output_dir / "representatives.fasta"
    # Clean
    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Create MMseqs2 DB
    db_in = db_dir / "inputDB"
    run(["mmseqs", "createdb", str(args.input_fasta), str(db_in)])

    # Hierarchical clustering (single stage)
    db_clu = db_dir / "clusters"
    run([
        "mmseqs", "cluster",
        str(db_in), str(db_clu), str(tmp_dir),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.cov),
        "--cov-mode", str(args.cov_mode),
        "--cluster-mode", str(args.cluster_mode),
        "--kmer-per-seq", str(args.kmer_per_seq),
        "--threads", str(args.threads)
    ])

    # Extract one representative per cluster
    db_rep = db_dir / "repDB"
    run([
        "mmseqs", "result2repseq",
        str(db_in), str(db_clu), str(db_rep)
    ])

    # Convert repDB back to FASTA
    run([
        "mmseqs", "convert2fasta", str(db_rep), str(rep_fa)
    ])

    print(f"\n✓ DONE — representatives written to: {rep_fa}")

if __name__ == "__main__":
    main()