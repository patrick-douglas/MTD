#!/bin/bash

# List of files to check
files=(
    "gene2pubmed.gz"
    "gene2accession.gz"
    "gene2refseq.gz"
    "gene_info.gz"
    "gene2go.gz"
)

# Base URL for NCBI FTP
ftp_base="https://ftp.ncbi.nlm.nih.gov/gene/DATA"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --MakeOrgpkgDir) offline_files_folder="$2"; shift ;;
        *) echo "❌ Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

# Check if folder was provided
if [ -z "$offline_files_folder" ]; then
    echo "❌ Error: You must provide the folder using --MakeOrgpkgDir"
    exit 1
fi

# Check if folder exists
if [ ! -d "$offline_files_folder" ]; then
    echo "❌ The folder $offline_files_folder does not exist!"
    exit 1
fi

# Function to get remote file size
get_remote_size() {
    curl -sI "$1" | grep -i "^Content-Length" | awk '{print $2}' | tr -dc '0-9'
}

# Check and download each file
for file in "${files[@]}"; do
    echo "🔍 Checking $file..."

    local_file="${offline_files_folder}/${file}"
    remote_file_url="${ftp_base}/${file}"
    remote_size=$(get_remote_size "$remote_file_url")

    if [ -z "$remote_size" ]; then
        echo "⚠️  Could not retrieve remote size for $file. Skipping..."
        continue
    fi

    if [ ! -f "$local_file" ]; then
        echo "📥 $file not found locally. Downloading..."
        curl -# -o "$local_file" "$remote_file_url"
        continue
    fi

    local_size=$(stat -c %s "$local_file")

    if [ "$local_size" != "$remote_size" ]; then
        echo "🔄 $file is outdated (local: $local_size bytes vs remote: $remote_size bytes). Updating..."
        curl -# -o "$local_file" "$remote_file_url"
    else
        echo "✅ $file is up to date."
    fi
done

echo "🟢 All files checked."
