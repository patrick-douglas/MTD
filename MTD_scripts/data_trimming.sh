max_jobs=$(nproc)
max_fastp_cores=16

if [ "$threads" -gt "$max_fastp_cores" ]; then
    fastp_threads=$max_fastp_cores
else
    fastp_threads=$threads
fi

# Processar cada amostra
for i in $lsn; do
    # Encontre o arquivo fastq correspondente (suporta .fq.gz, .fastq.gz, .fq, ou .fastq)
    fq=$(find $inputdr -name "${i}*.fq.gz" -o -name "${i}*.fastq.gz" -o -name "${i}*.fq" -o -name "${i}*.fastq" -type f)

    if [ -z "$no_trimm" ]; then
        # Se no_trimm não for definido, use o fastp para limpeza dos dados
        echo 'Trimming fastq files with fastp'
        fastp --trim_poly_x \
              --length_required $length \
              --thread $fastp_threads \
              -i $fq \
              -o $outputdr/temp/Trimmed_${i}.fq.gz
    fi
done

# Compressão paralela e cópia se no_trimm for definido
if [ -n "$no_trimm" ]; then
    echo 'Compressing fastq files to .gz'
    find $inputdr -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" | xargs -I {} -P $max_jobs sh -c '
        input_file="$1"
        base_name=$(basename "${input_file%.*}")
        
        if [ "${input_file##*.}" = "gz" ]; then
            # Se o arquivo já estiver comprimido, apenas processe o nome
            base_name_no_suffix=$(echo "$base_name" | sed "s/_R[0-9]$//")
            output_file="$2/Trimmed_${base_name_no_suffix}.fq.gz"
            # Copie o arquivo comprimido para o diretório de saída com o novo nome
            cp "$input_file" "$output_file"
        else
            # Se o arquivo não estiver comprimido, remova o sufixo e comprima
            base_name_no_suffix=$(echo "$base_name" | sed "s/_R[0-9]$//")
            output_file="$2/Trimmed_${base_name_no_suffix}.fq.gz"
            gzip --fast -c "$input_file" > "$output_file"
        fi
    ' _ {} "$outputdr/temp"
fi
