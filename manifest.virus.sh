#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0. Diretórios e nomes de arquivos
###############################################################################
offline_files_folder="/media/me/18TB_BACKUP_LBN/lbn_workspace/RNA-Seq-LBN/viral-rna-seq/MTD/Compressed/MTD"

# pasta definitiva para TODOS os .gz e para o FASTA combinado
new_download_dir="$offline_files_folder/Kraken2DB_micro/library/viral/all"

mkdir -p "$new_download_dir"
mkdir -p "$offline_files_folder/Kraken2DB_micro/library/viral"

assembly_summary_file="$offline_files_folder/Kraken2DB_micro/library/viral/assembly_summary_viral.txt"
manifest_list="$offline_files_folder/Kraken2DB_micro/library/viral/manifest_viral.list.txt"
failed_downloads="$offline_files_folder/failed_downloads.txt"
corrupted_list="$offline_files_folder/corrupted_viral.txt"
to_download_list="$offline_files_folder/to_download_viral.txt"
obsolete_list="$offline_files_folder/obsolete_local_viral.txt"

rm -f "$assembly_summary_file" "$manifest_list" "$failed_downloads" \
      "$corrupted_list" "$to_download_list" "$obsolete_list"

###############################################################################
# 1. Baixar assembly_summary.txt (viral) e gerar manifest
###############################################################################
echo "🔄 Downloading latest assembly_summary.txt (viral)…"
curl -4 --retry 5 --retry-delay 2 --connect-timeout 20 -fsSL \
  -o "$assembly_summary_file" \
  "https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt"

awk -F '\t' '
    /^#/ { next }
    $20 != "na" {
        ftp_path = $20
        sub(/^ftp:/, "https:", ftp_path)
        gsub(/\/+$/, "", ftp_path)

        n = split(ftp_path, a, "/")
        asm = a[n]

        if (asm != "") {
            print ftp_path "/" asm "_genomic.fna.gz"
        }
    }
' "$assembly_summary_file" > "$manifest_list"

echo "✅ $(wc -l < "$manifest_list") viral genomes listed on NCBI servers."

###############################################################################
# 2. Sincronizar pasta local com o NCBI
#    - remover obsoletos
#    - detectar faltantes
#    - detectar alterados por tamanho
###############################################################################
echo "🔄 Syncing local folder against NCBI manifest…"

mapfile -t remote_urls < "$manifest_list"
total_remote=${#remote_urls[@]}

declare -A remote_map
declare -A remote_size_map

for url in "${remote_urls[@]}"; do
    fname="$(basename "$url")"
    remote_map["$fname"]="$url"
done

shopt -s nullglob
local_paths=("$new_download_dir"/*.gz)
shopt -u nullglob

: > "$obsolete_list"
: > "$to_download_list"
: > "$failed_downloads"

# 2A. Remover arquivos locais que não existem mais no NCBI
obsolete_count=0
for path in "${local_paths[@]}"; do
    fname="$(basename "$path")"
    if [[ -z "${remote_map[$fname]:-}" ]]; then
        echo "❌ Removing obsolete local file: $fname"
        echo "$fname" >> "$obsolete_list"
        rm -f "$path"
        ((obsolete_count+=1))
    fi
done

# Atualizar lista local após remoções
shopt -s nullglob
local_paths=("$new_download_dir"/*.gz)
shopt -u nullglob

declare -A local_map
for path in "${local_paths[@]}"; do
    fname="$(basename "$path")"
    local_map["$fname"]="$path"
done

# Função para obter tamanho remoto
get_remote_size() {
    local url="$1"
    curl -4 -fsSI --connect-timeout 20 "$url" \
      | awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {gsub("\r","",$2); print $2; exit}'
}

# 2B. Detectar faltantes e alterados
missing_count=0
changed_count=0

for fname in "${!remote_map[@]}"; do
    url="${remote_map[$fname]}"

    if [[ -z "${local_map[$fname]:-}" ]]; then
        echo "$url" >> "$to_download_list"
        ((missing_count+=1))
        continue
    fi

    # compara tamanho remoto x local para detectar alteração/renomeação parcial/substituição
    local_size="$(stat -c%s "${local_map[$fname]}" 2>/dev/null || echo 0)"
    remote_size="$(get_remote_size "$url" || echo "")"

    if [[ -n "$remote_size" && "$remote_size" != "$local_size" ]]; then
        echo "🔁 Remote file changed, re-downloading: $fname"
        rm -f "${local_map[$fname]}"
        echo "$url" >> "$to_download_list"
        ((changed_count+=1))
    fi
done

available_local=$(find "$new_download_dir" -maxdepth 1 -type f -name "*.gz" | wc -l)
to_download_count=$(wc -l < "$to_download_list" 2>/dev/null || echo 0)

echo "
Total remote      : $total_remote
Already local     : $available_local
Obsolete removed  : $obsolete_count
Missing files     : $missing_count
Changed files     : $changed_count
To download now   : $to_download_count
"

###############################################################################
# 3. Baixar arquivos faltantes/alterados
###############################################################################
cd "$new_download_dir"

download_one() {
    local url="$1"
    local file
    file="$(basename "$url")"

    for attempt in {1..3}; do
        aria2c --disable-ipv6=true \
               --continue=true \
               --auto-file-renaming=false \
               -x16 -s16 \
               -o "$file" \
               "$url" > /dev/null 2>&1

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

if (( to_download_count == 0 )); then
    echo "🎉 Local folder is already synchronized with NCBI."
else
    echo "📥 Downloading $to_download_count genome(s)…"
    mapfile -t download_urls < "$to_download_list"

    if command -v parallel >/dev/null 2>&1 && command -v aria2c >/dev/null 2>&1; then
        printf "%s\n" "${download_urls[@]}" | parallel --bar -j 4 download_one
    else
        echo "⚠️  aria2c/parallel not found – falling back to serial wget."
        for url in "${download_urls[@]}"; do
            file="$(basename "$url")"
            wget -4 -q -c -O "$file" "$url" || echo "❌ Failed: $file" >> "$failed_downloads"
        done
    fi
    echo "✅ Download stage completed."
fi

###############################################################################
# 4. Checar integridade e re-baixar (até 3 tentativas) usando a URL do manifest
###############################################################################
echo -e "\n🔎 Verifying integrity of .gz files…"
: > "$corrupted_list"

shopt -s nullglob
gz_files=("$new_download_dir"/*.gz)
shopt -u nullglob

if (( ${#gz_files[@]} > 0 )); then
    for gz in "${gz_files[@]}"; do
        gzip -t "$gz" 2>/dev/null || echo "$(basename "$gz")" >> "$corrupted_list"
    done
else
    echo "⚠️  No .gz files found yet in $new_download_dir"
fi

corrupted_count=$(wc -l < "$corrupted_list")

if (( corrupted_count > 0 )); then
    echo "⚠️  $corrupted_count corrupted file(s) found. Re-downloading with wget (max 3 attempts)…"

    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue

        url="${remote_map[$fname]:-}"

        if [[ -z "$url" ]]; then
            echo "❌ URL not found in manifest for $fname" >> "$failed_downloads"
            rm -f "$new_download_dir/$fname"
            continue
        fi

        for attempt in {1..3}; do
            echo "  🔁 $fname (attempt $attempt)…"
            wget -4 -q -O "$new_download_dir/$fname" "$url" || true

            if gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                echo "    ✅ Integrity OK after attempt $attempt"
                break
            elif [[ $attempt -eq 3 ]]; then
                echo "    ❌ Still corrupted – removing and logging"
                echo "❌ Failed after 3 attempts: $fname" >> "$failed_downloads"
                rm -f "$new_download_dir/$fname"
            else
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
final_fasta="$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna"

echo -e "\n🧬 Building combined FASTA → $(basename "$final_fasta")"
: > "$combined_fasta"

shopt -s nullglob
gz_files=("$new_download_dir"/*.gz)
shopt -u nullglob

if (( ${#gz_files[@]} == 0 )); then
    echo "❌ No .gz files available to concatenate."
    exit 1
fi

for gz in "${gz_files[@]}"; do
    zcat "$gz" >> "$combined_fasta"
done

mv "$combined_fasta" "$final_fasta"
echo "✅ Combined FASTA created: $final_fasta"

###############################################################################
# 6. Relatório final
###############################################################################
echo -e "\n📋 Final report"
echo "Remote genomes listed : $total_remote"
echo "Obsolete removed      : $obsolete_count"
echo "Missing detected      : $missing_count"
echo "Changed detected      : $changed_count"
echo "Downloaded this run   : $to_download_count"

if [[ -s "$failed_downloads" ]]; then
    echo -e "\n⚠️  The following genomes could not be retrieved:"
    cat "$failed_downloads"
else
    echo -e "\n🎉 All viral genomes synchronized, verified and concatenated successfully!"
fi
