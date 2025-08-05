#!/usr/bin/env bash
# ------------------------------------------------------------------
# mmseqs_one_vs_all.sh — near-exhaustive one-vs-all search with header
#
# Usage:
#   ./mmseqs_one_vs_all.sh [-f] <sequences.fasta> <QUERY_ID> <out_dir> <tmp_dir>
#
# Options:
#   -f | --force   Remove any existing MMseqs databases before rebuilding them.
#
# Output:
#   <out_dir>/<QUERY_ID>_vs_all.tsv   —  TSV with header
# ------------------------------------------------------------------

set -euo pipefail

# ---- Parse optional -f / --force flag ---------------------------------------
FORCE=0
if [[ ${1:-} == "-f" || ${1:-} == "--force" ]]; then
    FORCE=1
    shift
fi

# ---- Check positional arguments --------------------------------------------
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 [-f] <sequences.fasta> <QUERY_ID> <out_dir> <tmp_dir>" >&2
    exit 1
fi

FASTA="$1"
QUERY_ID="$2"
OUT_DIR="$3"
TMP_DIR="$4"

mkdir -p "$OUT_DIR" "$TMP_DIR"

ALL_DB="$TMP_DIR/allDB"
QUERY_FASTA="$TMP_DIR/query.fasta"
QUERY_DB="$TMP_DIR/queryDB"
RESULT_DB="$TMP_DIR/resultDB"
RESULT_TSV="$OUT_DIR/${QUERY_ID}_vs_all.tsv"

# ---- Optional clean-up ------------------------------------------------------
if [[ $FORCE -eq 1 ]]; then
    echo "[0/4] Removing existing MMseqs databases …"
    rm -f "${ALL_DB}"* "${QUERY_DB}"* "${RESULT_DB}"*
fi

# ---- 1) Build the database with all sequences ------------------------------
echo "[1/4] Building MMseqs database from ${FASTA} …"
mmseqs createdb "$FASTA" "$ALL_DB"

# ---- 2) Extract the query sequence -----------------------------------------
echo "[2/4] Extracting sequence with ID '${QUERY_ID}' …"
awk -v id="$QUERY_ID" '
    BEGIN        { hdr = ">" id; keep = 0 }
    /^>/         { keep = ($0 == hdr) }
    keep         { print }
' "$FASTA" >"$QUERY_FASTA"

if [[ ! -s "$QUERY_FASTA" ]]; then
    echo "❌  Query ID '${QUERY_ID}' not found in the FASTA file." >&2
    exit 1
fi

mmseqs createdb "$QUERY_FASTA" "$QUERY_DB"

# ---- 3) Run the most permissive search MMseqs allows -----------------------
echo "[3/4] Running near-exhaustive MMseqs2 search …"
mmseqs search \
    "$QUERY_DB" "$ALL_DB" "$RESULT_DB" "$TMP_DIR" \
    --max-seqs 1000000 \
    -s 7.5 \
    -e inf \
    -c 0.0 --cov-mode 0 \
    --min-ungapped-score 0 \
    --prefilter-mode 2 \
    --exhaustive-search-filter 1 \
    -a 1 \
    --num-iterations 4

# ---- 4) Convert to TSV with header -----------------------------------------
echo "[4/4] Writing TSV (with header) to ${RESULT_TSV} …"
mmseqs convertalis \
    "$QUERY_DB" "$ALL_DB" "$RESULT_DB" "$RESULT_TSV" \
    --format-mode 4 \
    --format-output "query,target,fident,alnlen,nident,mismatch,gapopen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits"

echo "✓ Finished. Results saved to: ${RESULT_TSV}"
