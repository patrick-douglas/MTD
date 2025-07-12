#!/usr/bin/env bash
set -euo pipefail

offline_files_folder="/media/me/4TB_BACKUP_LBN/Compressed/MTD"
new_download_dir="$offline_files_folder/Kraken2DB_micro/library/bacteria/all"
home_download_dir="$HOME/MTD/kraken2DB_micro/library/bacteria/all"

mkdir -p "$new_download_dir" "$home_download_dir"

assembly_summary_file="$offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary_bacteria.txt"
manifest_list="$offline_files_folder/Kraken2DB_micro/library/bacteria/manifest.list.txt"
failed_downloads="$offline_files_folder/failed_downloads.txt"
rm -rf $assembly_summary_file $manifest_list


###############################################################################
# 1. Baixar assembly_summary.txt e gerar manifest.list.txt
###############################################################################
echo "🔄 Downloading latest assembly_summary.txt (bacteria)..."
curl -fsSL -o "$assembly_summary_file" \
  ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

awk -F '\t' '
    BEGIN { OFS="/" }
    /^#/ { next }
    ($12 == "Complete Genome" || $12 == "Chromosome") && $20 != "na" {
        ftp_path=$20
        sub(/^ftp:/,"https:",ftp_path)
        split(ftp_path,a,"/")
        asm=a[length(a)]
        print ftp_path, asm "_genomic.fna.gz"
    }
' "$assembly_summary_file" > "$manifest_list"

echo "✅ $(wc -l < "$manifest_list") genomes available at NCBI servers"

###############################################################################
# 3. Verificar arquivos já existentes
###############################################################################
echo "🔄 Checking for new available bacterial genomes, this may take a while, please wait..."
mapfile -t all_urls < "$manifest_list"

total_files=${#all_urls[@]}

# Criar lista dos arquivos faltantes (que não existem no new_download_dir)
missing_urls=()
for url in "${all_urls[@]}"; do
  filename=$(basename "$url")
  if [[ ! -f "$new_download_dir/$filename" ]]; then
    missing_urls+=("$url")
  fi
done

missing_count=${#missing_urls[@]}

echo "Files to download: $missing_count / $total_files"

###############################################################################
# 4. Download missing files with adaptive = progress bar
###############################################################################
cd "$new_download_dir"
> "$failed_downloads"
download_count=0

# Define a frequência adaptativa de atualização da barra
if (( missing_count <= 100 )); then
    update_freq=1
elif (( missing_count <= 1000 )); then
    update_freq=$(( missing_count / 100 ))  # 1%
elif (( missing_count <= 10000 )); then
    update_freq=$(( missing_count / 200 ))  # 0.5%
else
    update_freq=500  # fixo
fi

draw_progress_bar() {
    local current=$1
    local total=$2
    local bar_width=40
    local percent=$(( 100 * current / total ))
    local filled=$(( (percent * bar_width) / 100 ))
    local empty=$(( bar_width - filled ))
    local bar=$(printf "%0.s=" $(seq 1 $filled))
    local space=$(printf "%0.s " $(seq 1 $empty))
    printf "\rProgress: [${bar}${space}] %d/%d (%d%%)" "$current" "$total" "$percent"
}

download_one () {
  local url="$1"
  local file
  file="$(basename "$url")"

  aria2c --continue --auto-file-renaming=false -x16 -s16 -o "$file" "$url" \
    > /dev/null 2>&1

  local exit_code=$?

  if [[ $exit_code -ne 0 ]]; then
      echo "❌ Failed: $file" >> "$failed_downloads"
      rm -f "$file"
  else
      [[ -f "$home_download_dir/$file" ]] || cp -p "$file" "$home_download_dir/"
  fi

  ((download_count++))

  if (( download_count % update_freq == 0 || download_count == missing_count )); then
    draw_progress_bar "$download_count" "$missing_count"
  fi
}

export -f download_one
export -f draw_progress_bar
export failed_downloads home_download_dir missing_count download_count update_freq

if [[ $missing_count -gt 0 ]]; then
  echo ""
  if command -v parallel &>/dev/null; then
    printf "%s\n" "${missing_urls[@]}" | parallel --halt now,fail=1 -j 4 download_one
  else
    for url in "${missing_urls[@]}"; do
      download_one "$url"
    done
  fi
  draw_progress_bar "$missing_count" "$missing_count"
  echo -e "\n✅ All downloads completed!"
else
  echo "No files to download."
fi
echo -e "\n"

###############################################################################
# 4. Relatório final
###############################################################################
if [[ -s "$failed_downloads" ]]; then
  echo "⚠️ Some downloads failed:"
  cat "$failed_downloads"
else
  echo "🎉 All downloads completed successfully!"
fi
