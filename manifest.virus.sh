#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0. Diretórios e nomes de arquivos
###############################################################################
offline_files_folder="/media/me/18TB_BACKUP_LBN/lbn_workspace/RNA-Seq-LBN/viral-rna-seq/MTD/Compressed/MTD/"

# pasta definitiva para TODOS os .gz e para o FASTA combinado
new_download_dir="$offline_files_folder/Kraken2DB_micro/library/viral/all"

mkdir -p "$new_download_dir"

assembly_summary_file="$offline_files_folder/Kraken2DB_micro/library/viral/assembly_summary_viral.txt"
manifest_list="$offline_files_folder/Kraken2DB_micro/library/viral/manifest_viral.list.txt"
failed_downloads="$offline_files_folder/failed_downloads.txt"
rm -rf "$assembly_summary_file" "$manifest_list" "$failed_downloads"

###############################################################################
# 1. Baixar assembly_summary.txt (viral) e gerar manifest
###############################################################################
echo "🔄 Downloading latest assembly_summary.txt (viral)…"
curl -fsSL -o "$assembly_summary_file" \
  https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt

awk -F '\t' '
    BEGIN { OFS="/" }
    /^#/ { next }
    $20 != "na" {
        ftp_path=$20
        sub(/^ftp:/,"https:",ftp_path)
        split(ftp_path,a,"/")
        asm=a[length(a)]
        print ftp_path, asm "_genomic.fna.gz"
    }
' "$assembly_summary_file" > "$manifest_list"

echo "✅ $(wc -l < "$manifest_list") viral genomes listed on NCBI servers."

###############################################################################
# 2. Verificar quais arquivos ainda não estão localmente
###############################################################################
echo "🔄 Checking local vs. remote…"
mapfile -t all_urls < "$manifest_list"
total_files=${#all_urls[@]}

missing_urls=()
for url in "${all_urls[@]}"; do
  fname=$(basename "$url")
  [[ -f "$new_download_dir/$fname" ]] || missing_urls+=("$url")
done

missing_count=${#missing_urls[@]}
available_local=$(ls -1 "$new_download_dir" 2>/dev/null | wc -l)

echo "
Total remote   : $total_files
Already local  : $available_local
To download    : $missing_count
"

###############################################################################
# 3. Baixar arquivos faltantes
###############################################################################
cd "$new_download_dir"
> "$failed_downloads"

download_one() {
    local url="$1"
    local file; file="$(basename "$url")"

    for attempt in {1..3}; do
        aria2c --continue --auto-file-renaming=false -x16 -s16 -o "$file" "$url" > /dev/null 2>&1
        if [[ -f "$file" ]] && gzip -t "$file" 2>/dev/null; then
            return 0
        else
            echo "❌ Attempt $attempt failed for $file"
            rm -f "$file"
            sleep 1
        fi
    done

    echo "❌ $file" >> "$failed_downloads"
}
export -f download_one
export failed_downloads

if (( missing_count == 0 )); then
    echo "🎉 Nothing new to download."
else
    echo "📥 Downloading $missing_count genome(s)…"
    if command -v parallel &>/dev/null && command -v aria2c &>/dev/null; then
        printf "%s\n" "${missing_urls[@]}" \
            | parallel --bar -j 4 download_one
    else
        echo "⚠️  aria2c/parallel not found – falling back to serial wget."
        for url in "${missing_urls[@]}"; do
          file=$(basename "$url")
          wget -q -c -O "$file" "$url" || echo "❌ Failed: $file" >> "$failed_downloads"
        done
    fi
    echo "✅ Download stage completed."
fi

###############################################################################
# 4. Checar integridade e re-baixar (até 3 tentativas) usando wget
###############################################################################
echo -e "\n🔎 Verifying integrity of .gz files…"
corrupted_list="$offline_files_folder/corrupted_viral.txt"
> "$corrupted_list"

for gz in "$new_download_dir"/*.gz; do
  gzip -t "$gz" 2>/dev/null || echo "$(basename "$gz")" >> "$corrupted_list"
done

corrupted_count=$(wc -l < "$corrupted_list")
if (( corrupted_count > 0 )); then
  echo "⚠️  $corrupted_count corrupted file(s) found. Re-downloading with wget (max 3 attempts)…"
  while read -r fname; do
    for attempt in {1..3}; do
      echo "  🔁 $fname (attempt $attempt)…"
      wget -q -O "$new_download_dir/$fname" \
        "https://ftp.ncbi.nih.gov/genomes/refseq/viral/*/*/*/$fname" || true
      if gzip -t "$new_download_dir/$fname" 2>/dev/null; then
        echo "    ✅ Integrity OK after attempt $attempt"
        break
      elif [[ $attempt -eq 3 ]]; then
        echo "    ❌ Still corrupted – removing and logging"
        echo "❌ Failed after 3 attempts: $fname" >> "$failed_downloads"
        rm -f "$new_download_dir/$fname"
      fi
    done
  done < "$corrupted_list"
else
  echo "✅ All .gz files passed integrity check."
fi
rm -f "$corrupted_list"

###############################################################################
# 5. Descompactar e concatenar em um único FASTA
###############################################################################
combined_fasta="$new_download_dir/all_viral_genomes.fna"
echo -e "\n🧬 Building combined FASTA → $(basename "$combined_fasta")"
> "$combined_fasta"
for gz in "$new_download_dir"/*.gz; do
  zcat "$gz" >> "$combined_fasta"
done
mv $combined_fasta $offline_files_folder/Kraken2DB_micro/library/viral/
echo "✅ Combined FASTA created: $combined_fasta"

###############################################################################
# 6. Relatório final
###############################################################################
if [[ -s "$failed_downloads" ]]; then
  echo -e "\n⚠️  The following genomes could not be retrieved:"
  cat "$failed_downloads"
else
  echo -e "\n🎉 All viral genomes downloaded, verified and concatenated successfully!"
fi
