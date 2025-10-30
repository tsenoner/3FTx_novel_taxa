#!/usr/bin/env python3
"""
Merge domain sequences and centipede proteins, then cluster using MMseqs2.

This script:
1. Merges extracted_domains.fasta and centipede_3ftx_proteins.fasta
2. Deduplicates sequences by ID
3. Runs MMseqs2 clustering at specified identity threshold
4. Outputs representative sequences and cluster information
"""

import argparse
import sys
import subprocess
import shutil
from pathlib import Path
from collections import OrderedDict


def read_fasta(fasta_path):
    """
    Read a FASTA file and return an OrderedDict of {id: sequence}.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        OrderedDict mapping sequence IDs to sequences
    """
    sequences = OrderedDict()
    current_id = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                # Extract sequence ID (everything after '>')
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(sequences, output_path):
    """
    Write sequences to a FASTA file.

    Args:
        sequences: Dict mapping sequence IDs to sequences
        output_path: Path to output FASTA file
    """
    with open(output_path, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def merge_fasta_files(fasta_paths, output_path):
    """
    Merge multiple FASTA files, deduplicating by sequence ID.

    Args:
        fasta_paths: List of paths to input FASTA files
        output_path: Path to output merged FASTA file

    Returns:
        Number of unique sequences in merged file
    """
    all_sequences = OrderedDict()
    duplicates = 0

    print("\nMerging FASTA files...")
    for fasta_path in fasta_paths:
        print(f"  Reading: {fasta_path}")
        sequences = read_fasta(fasta_path)
        print(f"    Found {len(sequences):,} sequences")

        for seq_id, seq in sequences.items():
            if seq_id in all_sequences:
                duplicates += 1
            else:
                all_sequences[seq_id] = seq

    print(f"\nTotal unique sequences: {len(all_sequences):,}")
    if duplicates > 0:
        print(f"Duplicates removed: {duplicates:,}")

    print(f"Writing merged file: {output_path}")
    write_fasta(all_sequences, output_path)

    return len(all_sequences)


def check_mmseqs():
    """Check if MMseqs2 is installed and available."""
    if not shutil.which("mmseqs"):
        print("ERROR: mmseqs2 not found in $PATH", file=sys.stderr)
        print(
            "Please install MMseqs2: https://github.com/soedinglab/MMseqs2",
            file=sys.stderr,
        )
        sys.exit(1)


def run_command(cmd, description=None):
    """Run a shell command and handle errors."""
    if description:
        print(f"  {description}...")
    print(f"    Running: {' '.join(str(c) for c in cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=False)
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Error running command: {e}", file=sys.stderr)
        raise


def parse_mmseqs_clusters(cluster_tsv_path):
    """Parse MMseqs2 cluster TSV to get cluster statistics."""
    cluster_map = {}  # rep_id -> list of member_ids

    with open(cluster_tsv_path, "r") as f:
        for line in f:
            rep_id, member_id = line.strip().split("\t")
            if rep_id not in cluster_map:
                cluster_map[rep_id] = []
            cluster_map[rep_id].append(member_id)

    return cluster_map


def run_mmseqs_clustering(
    input_fasta,
    output_dir,
    min_seq_id=0.70,
    coverage=0.80,
    cov_mode=0,
    cluster_mode=0,
    sensitivity=9.0,
    cluster_reassign=False,
    threads=14,
):
    """
    Run MMseqs2 clustering using native mmseqs command-line tool.

    Args:
        input_fasta: Path to input FASTA file
        output_dir: Directory for output files
        min_seq_id: Minimum sequence identity (0.0-1.0)
        coverage: Coverage threshold (0.0-1.0)
        cov_mode: Coverage mode (0=both, 1=target, 2=query)
        cluster_mode: Clustering mode (0=set cover, 1=connected component, 2=greedy, 3=greedy incremental)
        sensitivity: MMseqs sensitivity parameter
        cluster_reassign: Whether to reassign clusters
        threads: Number of threads

    Returns:
        Tuple of (representatives_fasta_path, tsv_path)
    """
    # Check MMseqs2 is installed
    check_mmseqs()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Temporary directory for all intermediate files
    tmp_dir = output_dir / "tmp"

    # Recreate tmp directory (removes all old temporary files)
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    # Database paths (inside tmp directory)
    db_in = tmp_dir / "db_in"
    db_clu = tmp_dir / "db_clu"
    db_rep = tmp_dir / "db_rep"

    # Output paths (in main output directory)
    rep_fasta = output_dir / "representatives.fasta"
    cluster_tsv = output_dir / "clusters.tsv"

    # Clean up old output files if they exist
    for path in [rep_fasta, cluster_tsv]:
        if path.exists():
            path.unlink()

    print("\n" + "=" * 60)
    print("Running MMseqs2 Clustering")
    print("=" * 60)
    print(f"Input: {input_fasta}")
    print(f"Output directory: {output_dir}")
    print("\nParameters:")
    print(f"  Min sequence identity: {min_seq_id * 100:.1f}%")
    print(f"  Coverage: {coverage * 100:.1f}%")
    print(
        f"  Coverage mode: {cov_mode} ({'both' if cov_mode == 0 else 'target' if cov_mode == 1 else 'query'})"
    )
    print(
        f"  Cluster mode: {cluster_mode} ({'set cover' if cluster_mode == 0 else 'connected component' if cluster_mode == 1 else 'greedy' if cluster_mode == 2 else 'greedy incremental'})"
    )
    print(f"  Sensitivity: {sensitivity}")
    print(f"  Cluster reassign: {cluster_reassign}")
    print(f"  Threads: {threads}")
    print()

    # Count input sequences
    print("Counting input sequences...")
    num_input = sum(1 for line in open(input_fasta) if line.startswith(">"))
    print(f"  Total sequences: {num_input:,}")

    # Run clustering
    print("\nClustering sequences (this may take a while)...")
    try:
        # Step 1: Create database
        run_command(
            ["mmseqs", "createdb", str(input_fasta), str(db_in)],
            "Creating MMseqs2 database",
        )

        # Step 2: Cluster
        cluster_cmd = [
            "mmseqs",
            "cluster",
            str(db_in),
            str(db_clu),
            str(tmp_dir),
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
            "--cluster-mode",
            str(cluster_mode),
            "-s",
            str(sensitivity),
            "-a",  # Add backtrace for cluster assignment
            "--threads",
            str(threads),
        ]
        if cluster_reassign:
            cluster_cmd.append("--cluster-reassign")

        run_command(cluster_cmd, "Running clustering")

        # Step 3: Extract representatives
        run_command(
            ["mmseqs", "result2repseq", str(db_in), str(db_clu), str(db_rep)],
            "Extracting representative sequences",
        )

        # Step 4: Convert representatives to FASTA
        run_command(
            ["mmseqs", "convert2fasta", str(db_rep), str(rep_fasta)],
            "Converting representatives to FASTA",
        )

        # Step 5: Create TSV cluster mapping
        run_command(
            [
                "mmseqs",
                "createtsv",
                str(db_in),
                str(db_in),
                str(db_clu),
                str(cluster_tsv),
            ],
            "Creating cluster mapping",
        )

        print("\n✓ Clustering complete!")

        # Count representatives
        num_reps = sum(1 for line in open(rep_fasta) if line.startswith(">"))
        print("\nResults:")
        print(f"  Representatives: {num_reps:,}")
        print(
            f"  Reduction: {num_input - num_reps:,} sequences ({(1 - num_reps / num_input) * 100:.1f}%)"
        )

        # Parse cluster statistics
        print("\nParsing cluster statistics...")
        cluster_map = parse_mmseqs_clusters(cluster_tsv)
        cluster_sizes = [len(members) for members in cluster_map.values()]

        print(f"  Total clusters: {len(cluster_map):,}")
        print(f"  Min cluster size: {min(cluster_sizes)}")
        print(f"  Max cluster size: {max(cluster_sizes)}")
        print(f"  Mean cluster size: {sum(cluster_sizes) / len(cluster_sizes):.2f}")

        # Create enhanced TSV with cluster sizes
        enhanced_tsv = output_dir / "clusters_with_sizes.tsv"
        print(f"\nCreating enhanced cluster mapping: {enhanced_tsv}")
        with open(enhanced_tsv, "w") as out_f:
            out_f.write("representative\tmember\tcluster_size\n")
            with open(cluster_tsv, "r") as in_f:
                for line in in_f:
                    rep_id, member_id = line.strip().split("\t")
                    cluster_size = len(cluster_map[rep_id])
                    out_f.write(f"{rep_id}\t{member_id}\t{cluster_size}\n")

        print("✓ Complete!")

        return rep_fasta, enhanced_tsv

    except Exception as e:
        print(f"\n✗ Error during clustering: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Merge FASTA files and cluster sequences using MMseqs2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/output arguments
    parser.add_argument(
        "--interpro-fasta",
        type=Path,
        default=Path("data/interm/interpro/extracted_domains.fasta"),
        help="Path to InterPro extracted domains FASTA",
    )
    parser.add_argument(
        "--centipede-fasta",
        type=Path,
        default=Path("data/interm/centipede_genome/centipede_3ftx_proteins.fasta"),
        help="Path to centipede proteins FASTA",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/interm/mmseqs"),
        help="Output directory for merged and clustered sequences",
    )

    # MMseqs2 parameters
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.70,
        help="Minimum sequence identity threshold (0.0-1.0)",
    )
    parser.add_argument(
        "--coverage", type=float, default=0.80, help="Coverage threshold (0.0-1.0)"
    )
    parser.add_argument(
        "--cov-mode",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Coverage mode: 0=both, 1=target, 2=query",
    )
    parser.add_argument(
        "--cluster-mode",
        type=int,
        choices=[0, 1, 2, 3],
        default=0,
        help="Cluster mode: 0=set cover, 1=connected component, 2=greedy, 3=greedy incremental",
    )
    parser.add_argument(
        "--sensitivity",
        type=float,
        default=9.0,
        help="MMseqs sensitivity parameter (higher = more sensitive, slower)",
    )
    parser.add_argument(
        "--cluster-reassign",
        action="store_true",
        help="Enable cluster reassignment",
    )
    parser.add_argument(
        "--threads", type=int, default=14, help="Number of threads to use"
    )
    parser.add_argument(
        "--skip-merge",
        action="store_true",
        help="Skip merging step and use existing merged.fasta",
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Merge FASTA files
    merged_fasta = args.output_dir / "merged.fasta"

    if args.skip_merge:
        if not merged_fasta.exists():
            print(
                f"ERROR: --skip-merge specified but {merged_fasta} does not exist",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Skipping merge step, using existing: {merged_fasta}")
    else:
        # Check input files exist
        if not args.interpro_fasta.exists():
            print(
                f"ERROR: InterPro FASTA not found: {args.interpro_fasta}",
                file=sys.stderr,
            )
            sys.exit(1)
        if not args.centipede_fasta.exists():
            print(
                f"ERROR: Centipede FASTA not found: {args.centipede_fasta}",
                file=sys.stderr,
            )
            sys.exit(1)

        merge_fasta_files([args.interpro_fasta, args.centipede_fasta], merged_fasta)

    # Step 2: Run MMseqs2 clustering
    rep_fasta, cluster_tsv = run_mmseqs_clustering(
        input_fasta=merged_fasta,
        output_dir=args.output_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        cov_mode=args.cov_mode,
        cluster_mode=args.cluster_mode,
        sensitivity=args.sensitivity,
        cluster_reassign=args.cluster_reassign,
        threads=args.threads,
    )

    print("\n" + "=" * 60)
    print("All steps completed successfully!")
    print("=" * 60)
    print(f"Merged sequences: {merged_fasta}")
    print(f"Representatives: {rep_fasta}")
    print(f"Cluster mapping: {cluster_tsv}")
    print()


if __name__ == "__main__":
    main()
