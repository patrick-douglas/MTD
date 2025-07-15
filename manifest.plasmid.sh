#!/bin/bash
set -euo pipefail

###############################################################################
# BASIC CONFIGURATION
###############################################################################
FTP_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_DIR="genomes/refseq/plasmid/"
LOCAL_DIR="/media/me/18TB_BACKUP_LBN/lbn_workspace/RNA-Seq-LBN/viral-rna-seq/MTD/Compressed/MTD/Kraken2DB_micro/library/plasmid"

# Temporary files and failed download list
TMP_DIR="$HOME"
mkdir -p "$TMP_DIR"
FAILED_DL="$TMP_DIR/failed_downloads.txt"

TMP_REMOTE_LIST="$TMP_DIR/ncbi_remote_list.txt"
TMP_REMOTE_NAMES="$TMP_DIR/ncbi_remote_names.txt"
TMP_LOCAL_LIST="$TMP_DIR/local_list.txt"

###############################################################################
# CREATE / CHECK DESTINATION FOLDER
###############################################################################
mkdir -p "$LOCAL_DIR"

###############################################################################
# FETCH REMOTE & LOCAL FILE LISTS
###############################################################################
ftp -n "$FTP_SERVER" <<END_SCRIPT >"$TMP_REMOTE_LIST"
quote USER anonymous
quote PASS anonymous
cd $FTP_DIR
ls
quit
END_SCRIPT
awk '{print $NF}' "$TMP_REMOTE_LIST" >"$TMP_REMOTE_NAMES"

ls "$LOCAL_DIR" >"$TMP_LOCAL_LIST"

###############################################################################
# DEFINE FILES TO BE DOWNLOADED
###############################################################################
mapfile -t ALL_REMOTE <"$TMP_REMOTE_NAMES"

MISSING=()
for f in "${ALL_REMOTE[@]}"; do
  [[ -f "$LOCAL_DIR/$f" ]] || MISSING+=("$f")
done

TOTAL=${#ALL_REMOTE[@]}
NEEDED=${#MISSING[@]}
LOCAL_NOW=$(ls "$LOCAL_DIR" | wc -l)

echo -e "\n────────  SUMMARY  ────────
Total on FTP     : $TOTAL
Already local    : $LOCAL_NOW
To be downloaded : $NEEDED
───────────────────────────"

[[ $NEEDED -eq 0 ]] && { echo "🎉 Nothing new to download."; rm -f "$TMP_REMOTE_LIST" "$TMP_REMOTE_NAMES" "$TMP_LOCAL_LIST"; exit 0; }

###############################################################################
# FUNCTIONS: PROGRESS + DOWNLOAD
###############################################################################
draw_bar() {
  local done=$1 need=$2 width=40
  local pct=$((100*done/need))
  local fill=$((pct*width/100))
  printf "\rProgress: [%s%s] %d/%d (%d%%)" \
         "$(printf '=%.0s' $(seq 1 $fill))" \
         "$(printf ' %.0s' $(seq 1 $((width-fill))))" \
         "$done" "$need" "$pct"
}

download_one() {
  local file="$1"
  aria2c --continue --auto-file-renaming=false -x16 -s16 \
         -d "$LOCAL_DIR" -o "$file" \
         "ftp://$FTP_SERVER/$FTP_DIR/$file" \
         > /dev/null 2>&1 || {
    echo "❌ $file" >>"$FAILED_DL"
    return 1
  }
}

export -f download_one

###############################################################################
# LOOP / PARALLEL WITH PROGRESS BAR
###############################################################################
> "$FAILED_DL"
DL_DONE=0

if command -v parallel &>/dev/null; then
  printf "%s\n" "${MISSING[@]}" | parallel --no-notice --halt now,fail=1 -j4 --bar \
    'aria2c --continue --auto-file-renaming=false -x16 -s16 \
     -d "'"$LOCAL_DIR"'" -o {} "ftp://'"$FTP_SERVER"'/'"$FTP_DIR"'{}" > /dev/null 2>&1 || \
     { echo "❌ {}" >> "'"$FAILED_DL"'"; rm -f "'"$LOCAL_DIR"'/{}"; }'
else
  # sequential fallback with custom progress bar
  for file in "${MISSING[@]}"; do
    download_one "$file" || true
    ((DL_DONE++))
    draw_bar "$DL_DONE" "$NEEDED"
  done
  echo
fi
###############################################################################
# VERIFICAÇÃO DE INTEGRIDADE E RE-DOWNLOAD DE ARQUIVOS CORROMPIDOS
###############################################################################
echo -e "\n🔎 Verifying integrity of downloaded .gz files..."

CORRUPTED_LIST="$TMP_DIR/corrupted_files.txt"
> "$CORRUPTED_LIST"

for file in "$LOCAL_DIR"/*.gz; do
  gzip -t "$file" 2>/dev/null || echo "$(basename "$file")" >> "$CORRUPTED_LIST"
done

CORRUPTED_COUNT=$(wc -l < "$CORRUPTED_LIST")
if [[ $CORRUPTED_COUNT -gt 0 ]]; then
  echo -e "⚠️  Found $CORRUPTED_COUNT corrupted files. Attempting re-download with wget (3 attempts max)...\n"
  while read -r file; do
    echo "🔁 Re-downloading $file..."

    for attempt in {1..3}; do
      echo "  Attempt $attempt..."
      wget -q -O "$LOCAL_DIR/$file" "ftp://$FTP_SERVER/$FTP_DIR/$file"
      if gzip -t "$LOCAL_DIR/$file" 2>/dev/null; then
        echo "  ✅ Integrity check passed on attempt $attempt"
        break
      else
        echo "  ❌ Integrity check failed on attempt $attempt"
        (( attempt == 3 )) && {
          echo "  ⚠️  Final attempt failed — removing $file and logging as failed."
          echo "$file" >> "$FAILED_DL"
          rm -f "$LOCAL_DIR/$file"
        }
      fi
    done

  done < "$CORRUPTED_LIST"
else
  echo "✅ All .gz files passed integrity check."
fi

rm -f "$CORRUPTED_LIST"

###############################################################################
# FINAL REPORT
###############################################################################
if [[ -s "$FAILED_DL" ]]; then
  echo -e "\n⚠️  Some downloads failed:"
  cat "$FAILED_DL"
else
  echo -e "\n🎉 All downloads completed successfully!"
fi

# Cleanup
rm -f "$TMP_REMOTE_LIST" "$TMP_REMOTE_NAMES" "$TMP_LOCAL_LIST"
