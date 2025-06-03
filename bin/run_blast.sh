#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_blast.sh data/refs.fasta data/proteome_db/custom_proteomes results
#
# Arguments:
#   refs.fasta   – FASTA file of query sequences
#   db_prefix    – prefix of BLAST DB (e.g. custom_proteomes)
#   out_dir      – directory to write intermediate and final results
#
# This script will:
#  1. Run blastp of refs.fasta against the specified DB
#  2. Filter hits by: percent identity ≥ 30%, alignment coverage ≥ 70%, e‐value ≤ 1e‐5
#  3. For each query, keep only the single best hit (highest bitscore)
#  4. Write out:
#     • blast_raw.tsv           (all raw hits)
#     • blast_filtered.tsv      (after applying pident/coverage/e‐value filters)
#     • blast_best_per_query.tsv (best hit per query)

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 refs.fasta db_prefix out_dir"
  exit 1
fi

REFS_FA="$1"
DB_PREFIX="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

RAW_OUT="$OUT_DIR/blast_raw.tsv"
FILTERED_OUT="$OUT_DIR/blast_filtered.tsv"
BEST_OUT="$OUT_DIR/blast_best_per_query.tsv"

echo "1) Running blastp..."
blastp \
  -query "$REFS_FA" \
  -db "$DB_PREFIX" \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -evalue 1e-5 \
  -num_threads 4 \
  -out "$RAW_OUT"

echo "2) Filtering hits (pident>=30%, coverage>=70%, evalue<=1e-5)..."
awk '{
  qid   = $1
  sid   = $2
  pident = $3
  alnlen = $4
  qlen   = $5
  evalue = $7
  # compute percent coverage of query: alnlen/qlen
  if (pident >= 30 && (alnlen/qlen) >= 0.70 && evalue <= 1e-5) {
    print $0
  }
}' "$RAW_OUT" > "$FILTERED_OUT"

echo "3) Selecting best hit per query (highest bitscore)..."
# sort by query (field 1) and descending bitscore (field 8), then pick first per query
sort -k1,1 -k8,8nr "$FILTERED_OUT" \
  | awk '!seen[$1]++' \
  > "$BEST_OUT"

echo "Done."
echo "  • Raw BLAST output:        $RAW_OUT"
echo "  • Filtered hits:          $FILTERED_OUT"
echo "  • Best hit per query:     $BEST_OUT"
