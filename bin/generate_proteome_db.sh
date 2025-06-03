#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./generate_proteome_db.sh data/species_taxids.tsv data/proteome_db
#
# Arguments:
#   species_taxids.txt  – whitespace-separated file with columns: species (can have spaces) and taxid (last field)
#   out_dir             – directory to hold downloaded proteomes, BLAST DB, and missing list
#
# Requirements:
#   • NCBI Datasets CLI installed and on PATH (`datasets` command)
#   • blast+ installed (`makeblastdb` command)
#   • unzip

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 species_taxids.txt out_dir"
  exit 1
fi

SPECIES_FILE="$1"
OUT_DIR="$2"
ZIP_DIR="$OUT_DIR/zips"
EXTRACT_DIR="$OUT_DIR/extracted"
COMBINED_FASTA="$OUT_DIR/all_proteomes.faa"
DB_PREFIX="$OUT_DIR/custom_proteomes"
MISSING_FILE="$OUT_DIR/missing_taxids.txt"

mkdir -p "$ZIP_DIR" "$EXTRACT_DIR" "$OUT_DIR"
> "$MISSING_FILE"

# Step 1: Download reference proteome ZIP for each taxID from species_taxids.txt
# We assume each line’s last whitespace-separated field is the taxID.
tail -n +2 "$SPECIES_FILE" | while read -r LINE; do
  # Extract taxid = last field
  taxid="$(echo "$LINE" | awk '{print $NF}')"
  # Extract species = everything before the last field
  species="$(echo "$LINE" | sed "s/[[:space:]]\+$taxid\$//")"
  species="$(echo "$species" | xargs)"  # trim leading/trailing spaces

  if [[ -z "$taxid" ]]; then
    continue
  fi

  # Sanitize species for filenames
  safe_name="$(echo "$species" | sed 's/[[:space:]/]/_/g')"
  ZIP_PATH="$ZIP_DIR/${safe_name}_${taxid}_proteome.zip"
  echo "Downloading proteome for '$species' (taxID $taxid) → $ZIP_PATH"

  # Try strict --reference
  if ! datasets download genome taxon "$taxid" \
       --reference \
       --include protein \
       --filename "$ZIP_PATH"; then

    echo "  ● No 'reference' found for $species ($taxid); falling back to complete/chromosome"
    if ! datasets download genome taxon "$taxid" \
         --assembly-level complete,chromosome \
         --include protein \
         --filename "$ZIP_PATH"; then
      echo "  ✖ Failed to download any proteome for $species ($taxid)"
      echo "$species	$taxid" >> "$MISSING_FILE"
      rm -f "$ZIP_PATH"
      continue
    fi
  fi
done

# Step 2: Unzip each downloaded ZIP into its own subdirectory
for ZIP in "$ZIP_DIR"/*_proteome.zip; do
  [ -f "$ZIP" ] || continue
  BASENAME="$(basename "$ZIP" .zip)"
  DEST="$EXTRACT_DIR/$BASENAME"
  echo "Extracting $ZIP → $DEST"
  mkdir -p "$DEST"
  unzip -qq "$ZIP" -d "$DEST"
done

# Step 3: Concatenate all protein.faa into one combined FASTA
echo "Concatenating protein.faa files into $COMBINED_FASTA"
rm -f "$COMBINED_FASTA"
FOUND_ANY=false
for FAA in "$EXTRACT_DIR"/*/ncbi_dataset/data/*/protein.faa; do
  if [[ -f "$FAA" ]]; then
    cat "$FAA" >> "$COMBINED_FASTA"
    FOUND_ANY=true
  fi
done

if ! $FOUND_ANY; then
  echo "No proteomes were extracted; all taxa may have failed. See $MISSING_FILE"
  exit 1
fi

# Step 4: Build BLAST database
echo "Building BLAST DB with makeblastdb"
makeblastdb -in "$COMBINED_FASTA" -dbtype prot -out "$DB_PREFIX"

echo "Done. BLAST DB files are at:"
echo "  ${DB_PREFIX}.phr"
echo "  ${DB_PREFIX}.pin"
echo "  ${DB_PREFIX}.psq"
echo "Missing entries are listed in $MISSING_FILE"
