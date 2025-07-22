#!/usr/bin/env python3
"""
Scan MMseqs2 clustering over multiple parameter combinations and plot the elbow curves (#representatives vs %ID) for each.
"""
import argparse
import subprocess
import shutil
import sys
from pathlib import Path
import matplotlib.pyplot as plt


def check_mmseqs():
    if not shutil.which("mmseqs"):
        sys.exit("ERROR: mmseqs2 not found in $PATH")


def run(cmd):
    print(">>>", " ".join(cmd))
    subprocess.check_call(cmd)


def count_fasta_seqs(fasta_path):
    return sum(1 for line in open(fasta_path) if line.startswith(">"))


def main():
    p = argparse.ArgumentParser(
        description="Scan MMseqs2 cluster thresholds and plot #representatives for various parameter sets"
    )
    p.add_argument("--input", required=True, type=Path,
                   help="Input FASTA of domains")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Base directory for outputs")
    p.add_argument("--ids", nargs="+", type=float, required=True,
                   help="List of min-seq-id thresholds (e.g. 0.30 0.35 0.40 …)")
    p.add_argument("--cov", type=float, default=0.80,
                   help="Coverage fraction for all runs")
    p.add_argument("--cov-mode", type=int, default=2, choices=[0,1,2],
                   help="Coverage mode: 0=both,1=target,2=query")
    p.add_argument("--cluster-modes", nargs="+", type=int, choices=[0,1,2,3], default=[1],
                   help="List of cluster-mode values to test")
    p.add_argument("--reassign-options", nargs="+", type=int, choices=[0,1], default=[0],
                   help="Whether to include --cluster-reassign: 1=on, 0=off")
    p.add_argument("--sensitivities", nargs="+", type=float, default=[9.0],
                   help="List of MMseqs prefilter sensitivity values to test")
    p.add_argument("--threads", type=int, default=6,
                   help="Threads to use")
    args = p.parse_args()

    check_mmseqs()

    # Prepare for plotting
    curves = {}

    # Iterate over combinations of cluster-mode, reassign, sensitivity
    for cluster_mode in args.cluster_modes:
        for reassign_flag in args.reassign_options:
            for sens in args.sensitivities:
                label = f"mode{cluster_mode}_reassign{reassign_flag}_s{sens}"  
                curves[label] = []
                print(f"\n=== Running config: {label} ===")

                for id_val in args.ids:
                    work = args.outdir / label / f"{int(id_val*100)}"
                    db_in   = work / "db_in"
                    db_clu  = work / "db_clu"
                    tmp     = work / "tmp"
                    rep_fa  = work / "reps.fasta"

                    # fresh directories
                    if work.exists():
                        shutil.rmtree(work)
                    work.mkdir(parents=True)
                    tmp.mkdir()

                    # run mmseqs
                    run(["mmseqs", "createdb", str(args.input), str(db_in)])
                    cmd = [
                        "mmseqs", "cluster",
                        str(db_in), str(db_clu), str(tmp),
                        "--min-seq-id", str(id_val),
                        "-c", str(args.cov),
                        "--cov-mode", str(args.cov_mode),
                        "--cluster-mode", str(cluster_mode),
                        "-s", str(sens),
                        "--split-mode", "2",
                        "--threads", str(args.threads)
                    ]
                    if reassign_flag == 1:
                        cmd.append("--cluster-reassign")
                    run(cmd)

                    run(["mmseqs", "result2repseq", str(db_in), str(db_clu), str(work/"db_rep")])
                    run(["mmseqs", "convert2fasta", str(work/"db_rep"), str(rep_fa)])

                    n = count_fasta_seqs(rep_fa)
                    print(f"  -> {n} reps at {id_val*100:.0f}% ID")
                    curves[label].append((id_val*100, n))

    # Plot all curves
    plt.figure(figsize=(8,6))
    for label, data in curves.items():
        xs, ys = zip(*data)
        plt.plot(xs, ys, marker='o', label=label)
    plt.xlabel("Sequence identity cutoff (%)")
    plt.ylabel("Number of representative sequences")
    plt.title("MMseqs2 clustering elbow across parameter sets")
    plt.legend()
    plt.tight_layout()
    out_png = args.outdir / "elbow_multi.png"
    plt.savefig(out_png)
    print(f"\nElbow plot with multiple curves written to {out_png}")

if __name__ == "__main__":
    main()
