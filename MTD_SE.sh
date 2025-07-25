#!/bin/bash
# Defining colors
w=$(tput sgr0) 
r=$(tput setaf 1)
g=$(tput setaf 2) 
y=$(tput setaf 3) 
p=$(tput setaf 5) 
echo "${w}"
# default settings
pdm="spearman" # method in HALLA
length=35 # read length trimming by fastp
read_len=75 # the read length in bracken
threads=`nproc`
blast="hisat"  # Define "hisat" como padrão
while getopts i:o:h:m:p:l:r:bt option
do
    case "${option}" in
        i) inputdr=${OPTARG};;
        o) outputdr=${OPTARG};;
        h) hostid=${OPTARG};;
        m) metadata=${OPTARG};;
        p) pdm=${OPTARG};;
        l) length=${OPTARG};;
        r) read_len=${OPTARG};;
        b) blast="blast";;  # Se -b for usado, define blast como "blast"
        t) no_trimm=1;;
    esac
done

#Example of the latest cmd
#bash ~/MTD/MTD_SE.sh -i $input_dir/samplesheet.csv -o ~/MTD/A.macularius_vs_C.pusilla/ -h 8839 -b -t

# inputdr=~/RNAseq_raw_data/samplesheet.csv # select input directory; must store singe-end .fq.gz (eg. DJ01.fq.gz) of each sample in the same folder as the samplesheet.csv
# outputdr=~/MTD_Results/test1 # select outputdr directory
# hostid=9544 # Enter host species taxonomy ID; initally supporting 9544 (rhesus monkey), 9606 (human), and 10090 (mouse).
# threads=20 # CPU threads; suggest >=16, eg. 20
# pdm= spearman or pearson or mi or nmi or xicor or dcor # pairwise distance metrics refer to HALLA mannual
if [ -d "$outputdr" ]; then
    echo "The output directory '$outputdr' already exists."
    read -p "Do you want to delete it and overwrite the files? (y/n): " answer
    if [[ "$answer" =~ ^[Yy]$ ]]; then
        rm -rf "$outputdr"
        echo "Directory deleted."
    else
        echo "Operation cancelled by the user. Exiting."
        exit 1
    fi
fi
# get MTD.sh script file path (in the MTD folder)
MTDIR=$(dirname $(readlink -f $0))
#parentname="$(dirname "$MTDIR")"
echo "MTD directory is $MTDIR"

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh
conda deactivate # aviod multiple conda environment
conda activate MTD

inputdr=$(dirname $inputdr)
mkdir -p $outputdr
mkdir -p $outputdr/temp
cd $outputdr/temp

# Step 0: Host database auto selection
if [[ $hostid == 9606 ]]; then
    DB_host=$MTDIR/kraken2DB_human # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_human/genome_tran #for hisat2
    DB_blast=$MTDIR/human_blastdb/human_blastdb # for blast
    gtf=$MTDIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz # for featureCounts

elif [[ $hostid == 9544 ]]; then
    DB_host=$MTDIR/kraken2DB_rhesus # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_rhesus/genome_tran #for hisat2
    DB_blast=$MTDIR/rhesus_blastdb/rhesus_blastdb # for blast
    gtf=$MTDIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz # for featureCounts

elif [[ $hostid == 10090 ]]; then
    DB_host=$MTDIR/kraken2DB_mice
t    DB_hisat2=$MTDIR/hisat2_index_mouse/genome_tran
    DB_blast=$MTDIR/mouse_blastdb/mouse_blastdb # for blast
    gtf=$MTDIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz # for featureCounts

elif [[ -d "$MTDIR/kraken2DB_${hostid}" ]]; then # test if customized host species exist
    DB_host=$MTDIR/kraken2DB_${hostid}
    DB_hisat2=$MTDIR/hisat2_index_${hostid}/genome_tran
    DB_blast=$MTDIR/blastdb_${hostid}/blastdb_${hostid} # for blast
    gtf=$MTDIR/ref_${hostid}/ref_${hostid}.gtf.gz

else
    echo "Host species is not supported. You can use bash Customized_host.sh for building."
    exit 1
fi
DB_micro=$MTDIR/kraken2DB_micro # customized kraken database for microbiome

# for SRR input samples in the samplesheet.csv; download SRR samples
cd $inputdr
if [ ! -z "$(cat samplesheet.csv | cut -f 1 -d ','| grep ^SRR)" ]; then
    for s in $(cat samplesheet.csv | cut -f 1 -d ','| grep ^SRR); do
        # check if fastq files of SRR sample exists
        if [[ ! -f ${s}_1.fastq || ! -f ${s}_2.fastq ]]; then
            echo "File ${s} fastq files NOT exists. Start downloading..."
            echo 'download SRA files...'
            prefetch -X 999G ${s}
            echo 'split SRA files to fastq files...'
            fasterq-dump -p --split-files ${s}
            rm -rf ${s}
        fi
    done
fi
cd $outputdr/temp

# To extract sample names from input fastq files (support .fq.gz, .fastq.gz, .fq, or .fastq)
files=$(find $inputdr -name "*.fq.gz" -or -name "*.fastq.gz" -or -name "*.fq" -or -name "*.fastq" -type f)
#b=$(basename -a $inputdr) # store basenames of input directories into variable b to make a list of input sample names (eg. DJ01 EM77...)
for i in $files; do
    fn=$(basename $i) #Extract file name, eg. DJ01_1.fq.gz
    sn=$(echo $fn | awk -F '_R1' '{print $(NF-1)}') #Extract sample name, eg. DJ01
    lsn=$lsn" "$sn #Make a list of sample names; store basenames of input directories into variable lsn to make a list of input sample names (eg. DJ01 EM77...)
done

# check if input files match the samplesheet.csv
fastq_files=$(echo $lsn | tr " " "\n" | sort | tr "\n" " ")
SamplesInSheet=$(cat $inputdr/samplesheet.csv | cut -f 1 -d ',' | tail -n +2 | sort | tr "\n" " ")
if [[ "$fastq_files" != "$SamplesInSheet" ]]; then
    echo "The samples' fastq files in the input folder do not match with your samplesheet.csv"
    echo "Please double-check with the samplesheet.csv and input files. Please ensure no other fastq files are under the input folder and its subfolders. You can refer to the user guide on https://github.com/FEI38750/MTD."
    exit 1
fi
species_name=$(awk -F, -v taxid="$hostid" '$1 == taxid {print $3}' "$MTDIR/HostSpecies.csv")

# Verifique se o nome da espécie foi encontrado
if [ -z "$species_name" ]; then
  echo "Error: species name not found for Taxon_ID $hostid."
  exit 1
fi

echo "${g}"
echo "============================================"
echo -e "Selected host species:${w}\e[3m $species_name\e[0m${g}"
echo "Taxon ID:${w} $hostid ${g}"
echo ''
echo "============================================"
echo "Main study design:${w}"
awk -F',' 'NR>1 {groups[$2]++; if ($5=="vs") comparisons[$2" vs "$6]++} END { 
  for (g in groups) printf "Group: %s - Number of samples: %d\n", g, groups[g]; 
}' $inputdr/samplesheet.csv
echo "${g}============================================"
echo "${w}"

# Assuming $metadata is the path to the metadata file passed via -m flag
if [ ! -z "$metadata" ]; then
    echo "============================================"

    # Reading the header to dynamically identify columns
    header=$(head -n 1 "$metadata")
    IFS=',' read -ra columns <<< "$header"

    # Ignoring the first two columns
    for ((i=3; i<=${#columns[@]}; i++)); do
        col="${columns[$i-1]}"
        echo "Metadata column: $col,"
        echo "Meta-groups:"
        awk -v col_index="$i" -F',' '
        NR > 1 {
            values[$col_index]++;
        }
        END {
            for (value in values) {
                printf "  %s: %d\n", value, values[value];
            }
        }' "$metadata"
    done

    echo "============================================"
    echo ''
fi

echo "${g}MTD running  progress:"
echo ">>                  [10%]"

echo "Raw reads trimming${w}"
choice="execute"
#choice="skip" Just copy pre compressed files path is required
#choice="execute" Perform the filtering or not based if the parameter -t is declared or not
case $choice in
  execute)
max_jobs=$(( $(nproc) / 2 ))
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
    echo "${g}Compressing fastq files to .gz"
    echo "${y}WARNING: As the parameter -t was declared the data will not be trimmed/filtered with fastp${g}"
    echo 'Skipping trimming step...'
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
    ;;
  skip)
CUSTOM_PATH=/media/me/4TB_BACKUP_LBN/temp/A.macularius
    echo "WARNING: USING UNTRIMMED DATA FROM $CUSTOM_PATH
Skipping trimming with fastp step..."
    cp $CUSTOM_PATH/* .
    ;;
esac

#$MTDIR/MTD_scripts/data_trimming.sh 

echo "${g}MTD running  progress:"
echo '>>>>                [20%]'

echo "Reads classification by kraken2; 1st step for host ${w}"
for i in $lsn; do # store input sample name in i; eg. DJ01
    kraken2 --db $DB_host --use-names \
        --report Report_host_$i.txt \
        --threads $threads \
        --gzip-compressed \
        --classified-out ${i}_host.fq \
        --unclassified-out ${i}_non-host_raw.fq \
        Trimmed_${i}.fq.gz \
        > Report_host_$i.kraken
done

echo "${g}MTD running  progress:"
echo '>>>>>               [25%]'

echo "Reads classification by kraken2; 2nd step for non-host reads ${w}"
for i in $lsn; do # store input sample name in i; eg. DJ01
    kraken2 --db $DB_micro --use-names \
        --report Report_non-host.raw_$i.txt \
        --threads $threads \
        --classified-out ${i}_raw_cseqs.fq \
        --unclassified-out ${i}_raw_ucseqs.fq \
        ${i}_non-host_raw.fq \
        > Report_non-host_raw_$i.kraken
done

echo "${g}MTD running  progress:"
echo '>>>>>>              [30%]'

echo "Decontamination step${w}"
conta_file=$MTDIR/conta_ls.txt
if test -f "$conta_file"; then
    tls=$(awk -F '\t' '{print $2}' $conta_file)
    conta_ls="${tls//$'\r\n'/ }"
    for i in $lsn; do
        python $MTDIR/Tools/KrakenTools/extract_kraken_reads.py \
            -k Report_non-host_raw_${i}.kraken \
            -s1 ${i}_non-host_raw.fq \
            -o ${i}_non-host.fq \
            -r Report_non-host.raw_${i}.txt \
            --taxid $conta_ls --exclude --include-children
    done

    echo "${g}MTD running  progress:"
    echo '>>>>>>>             [35%]'

echo "Reads classification by kraken2; 3rd step for decontaminated non-host reads to get reports ${w}"
    for i in $lsn; do
        kraken2 --db $DB_micro --use-names \
            --report Report_non-host_$i.txt \
            --threads $threads \
            --classified-out ${i}_cseqs.fq \
            --unclassified-out ${i}_ucseqs.fq \
            ${i}_non-host.fq \
            > Report_non-host_$i.kraken
    done
fi

echo "${g}MTD running  progress:"
echo '>>>>>>>>            [40%]'

echo "Bracken analysis ${w}"
for i in $lsn; do # store input sample name in i; eg. DJ01
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.phylum.bracken -r $read_len -l P -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.genus.bracken -r $read_len -l G -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.species.bracken -r $read_len -l S -t $threads
done

echo "${g}MTD running  progress:"
echo '>>>>>>>>>           [45%]'
echo "combined .bracken files (table like) into a single outputdr for Deseq2 ${w}"
python $MTDIR/Tools/combine_bracken_outputs.py --files *.phylum.bracken -o $outputdr/bracken_phylum_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.genus.bracken -o $outputdr/bracken_genus_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.species.bracken -o $outputdr/bracken_species_all

echo "${g}Move _bracken report files (tree like) to a separate folder${w}"
mkdir -p Report_non-host_bracken_species_normalized
mv *_bracken_species.txt Report_non-host_bracken_species_normalized
cd Report_non-host_bracken_species_normalized

echo "${g}Trim the name of _bracken report files (tree like) to the sample name (eg. DJ01) ${w}"
for i in $lsn; do
    mv *${i}_* $i
done

echo "${g}Converted original _bracken report files (tree like) into .biom file for ANCOMBC and diversity analysis in phyloseq (R) etc. in DEG_Anno_Plot.R ${w}"
kraken-biom * -o $outputdr/temp/bracken_species_all0.biom --fmt json

echo "${g}Adjust bracken file (tree like) by normalizated reads counts; for additional visualization (.biom, .mpa, .krona) ${w}"
conda deactivate
conda activate R412
Rscript $MTDIR/Normalization_afbr.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $outputdr/temp/Report_non-host_bracken_species_normalized $metadata

conda deactivate
conda activate MTD

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>          [50%]'

echo "Converted adjusted _bracken report files (tree like) into .biom file for graph visualization: graphlan, MPA, krona ${w}"
kraken-biom * -o $outputdr/bracken_species_all.biom --fmt json
#Converted original _bracken report files (tree like) into .biom file
#kraken-biom * -o $outputdr/temp/bracken_species_all0.biom --fmt json
#kraken-biom *_bracken_phylum -o bracken_phylum_all.biom --fmt json
#kraken-biom *_bracken_genus -o bracken_genus_all.biom --fmt json

echo "${g}Remove "sp. " in the .biom file; correct improper format before run export2graphlan.py"
sed -i 's/sp. //g' $outputdr/bracken_species_all.biom

echo "Go to temp folder${w}"
cd ../

mkdir -p ../graphlan
cd ../graphlan

#source $condapath/etc/profile.d/conda.sh
conda deactivate
conda activate py2

python $MTDIR/Tools/export2graphlan/export2graphlan.py \
    -i ../bracken_species_all.biom \
    -a annot.txt -t tree.txt \
    --discard_otus --most_abundant 50 \
    --annotations 2,3,4,5,6 \
    --external_annotations 7 --internal_levels --max_clade_size 300

conda deactivate
conda activate MTD

cd ../temp

echo "${g}DEG & Annotation & Plots & Diversity & Preprocess for Microbiome ${w}"
conda deactivate
conda activate R412
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $hostid $MTDIR/HostSpecies.csv $metadata
conda deactivate
conda activate MTD

cd $outputdr/temp
mkdir -p bracken_raw_results # save the raw output from bracken (table like)
mv ../bracken_*_all bracken_raw_results

cd ../graphlan
echo "${g}Applying a fix for both tree.txt and Annot.txt${w}"
python $MTDIR/Tools/graphlan/verify_and_correct_annotations.py tree.txt annot.txt corrected_annot.txt
mv annot.txt annot_original.txt
mv corrected_annot.txt annot.txt

python $MTDIR/Tools/graphlan/graphlan_annotate.py --annot annot.txt tree.txt outtree.txt # attach annotation to the tree
python $MTDIR/Tools/graphlan/graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png # generate the graphlan png
python $MTDIR/Tools/graphlan/graphlan.py outtree.txt outimg.pdf # generate the graphlan pdf

cd ../temp

echo "${g}Visualization preprocess"
echo "For krona${w}"
mkdir -p ../krona
for i in $lsn; do # store input sample name in i; eg. DJ01
    python $MTDIR/Tools/KrakenTools/kreport2krona.py \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ../krona/${i}-bracken.krona
done

echo "${g}To make MPA style file${w}"
for i in $lsn; do # store input sample name in i; eg. DJ01
    python $MTDIR/Tools/KrakenTools/kreport2mpa.py \
        --display-header \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ${i}-bracken.mpa.txt
done
echo "${g}Combine MPA files${w}"
python $MTDIR/Tools/KrakenTools/combine_mpa.py \
    -i *.mpa.txt \
    -o ../Combined.mpa

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>         [55%]'

echo "HUMAnN3${w}"
mkdir -p HUMAnN_output

for n1 in *\_non-host.fq; do
    cp $n1 HUMAnN_output/$n1
done

cd HUMAnN_output

for file in *; do #trim the file name
    mv $file ${file/_non-host/}
done

echo "${g}Run HUMAnN3${w}"
for i in *.fq; do
    humann --input $i \
        --output hmn_output \
        --threads $threads \
        --verbose
done
echo "${g}>>>>>>>>>>>>        [60%]"

echo "Join all gene family and pathway abudance files${w}"
humann_join_tables -i hmn_output/ -o humann_pathabundance.tsv --file_name pathabundance
humann_join_tables -i hmn_output/ -o humann_genefamilies.tsv --file_name genefamilies

# #Normalizing RPKs to CPM
# humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_cpm.tsv --units cpm --update-snames
# humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_cpm.tsv --units cpm --update-snames

echo "${g}Normalizing RPKs to "relab" (relative abundance)${w}"
humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_relab.tsv --units relab --update-snames
humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_relab.tsv --units relab --update-snames

echo "${g}Generate stratified tables; This utility will split a table into two files (one stratified and one unstratified). ${w}"
humann_split_stratified_table --input humann_pathabundance_relab.tsv --output ./
humann_split_stratified_table --input humann_genefamilies_relab.tsv --output ./
    echo "${g}Stratify unnormalized table (for Deseq2)${w}"
    humann_split_stratified_table --input humann_pathabundance.tsv --output ./
    humann_split_stratified_table --input humann_genefamilies.tsv --output ./

echo "${g}Regroup gene familites table into KEGG orthologs and GO terms${w}"
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_ko --output humann_genefamilies_relAbundance_kegg.tsv
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_go --output humann_genefamilies_relAbundance_go.tsv
    echo "${g}Regroup unnormalized table (for Deseq2${w}"
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_ko --output humann_genefamilies_Abundance_kegg.tsv
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_go --output humann_genefamilies_Abundance_go.tsv

echo "${g}Translate KEGG and GO ID to human readable terms${w}"
conda deactivate
conda activate R412
#Rscript $MTDIR/humann_ID_translation.R \
Rscript $MTDIR/humann_ID_translation_adjusted.R $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_kegg.tsv $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_go.tsv $MTDIR
    # Tranlate unnormalized table (for Deseq2)
#    Rscript $MTDIR/humann_ID_translation.R \
Rscript $MTDIR/humann_ID_translation_adjusted.R $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_kegg.tsv $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_go.tsv $MTDIR
conda deactivate
conda activate MTD

#Cleaning up file structure
mkdir -p $outputdr/hmn_pathway_abundance_files
mkdir -p $outputdr/hmn_genefamily_abundance_files
mv *pathabundance* $outputdr/hmn_pathway_abundance_files/
mv *genefamilies* $outputdr/hmn_genefamily_abundance_files/

# #Translate KEGG and GO ID to human readable terms
# Rscript $MTDIR/humann_ID_translation.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg.tsv \
#     $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go.tsv

echo "${g}DEG & Annotation & Plots & Diversity & Preprocess${w}"
cd $outputdr/hmn_genefamily_abundance_files
conda deactivate
conda activate R412
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg_translated.tsv $inputdr/samplesheet.csv
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go_translated.tsv $inputdr/samplesheet.csv
conda deactivate && conda activate MTD

#humann_barplot
# humann_barplot --input $outputdr/hmn_pathway_abundance_files/humann_pathabundance_cpm_stratified.tsv \
#     --focal-metadatum Group --last-metadatum Group \
#     --focal-feature PWY-3781 \
#     --output $outputdr/hmn_pathway_abundance_files/humann_pathabundance_barplot.png
# humann_barplot --input $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_cpm_stratified.tsv \
#     --output $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_barplot.png

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>       [65%]'

echo "Starting to process the host reads...${w}"
## continue to process the host reads
cd $outputdr/temp
if [[ $blast == blast ]]; then
    echo "${g}Magic-BLAST${w}"
    for i in $lsn; do # store input sample name in i; eg. DJ01
        magicblast -query ${i}_host.fq \
        -db $DB_blast \
        -infmt fastq \
        -out $i.sam \
        -num_threads 8 #$threads
    done
#for i in $lsn; do magicblast -query ${i}_host.fq -db $DB_blast -infmt fastq -out $i.sam -num_threads 8; done
else
    echo "${g}HISAT2 alignment${w}"
    for i in $lsn; do # store input sample name in i; eg. DJ01
        hisat2 -p $threads -q \
            -x $DB_hisat2 \
            --summary-file ${i}_hisat2_summary.txt \
            -U Trimmed_${i}.fq.gz \
            -S $i.sam
    done
fi


echo "${g}featureCounts${w}"
featureCounts -T $threads -a $gtf -o $outputdr/host_counts.txt *.sam


for i in $lsn; do
    samtools view -bS $i.sam > $i.bam -@ $threads
    samtools sort $i.bam -o $i.sorted.bam -@ $threads
    samtools index $i.sorted.bam -@ $threads
done
#Comando abaixo [e o mesmo acima, mas em uma unica linha
#for i in $lsn; do samtools view -bS $i.sam > $i.bam -@ $threads && samtools sort $i.bam -o $i.sorted.bam -@ $threads && samtools index $i.sorted.bam -@ $threads; done

mkdir -p BAM
mv *.sorted.bam *.sorted.bam.bai BAM/

cd $outputdr
# trim the featureCounts output(host_counts.txt) for downstream analysis
echo "${g}Delete the first line/row of a file then trim the sample name${w}"
sed '1d; 2 s/\.sam//g' host_counts.txt > tmpfile; mv tmpfile host_counts.txt

echo "${g}DEG & Annotation & Plots & preprocess for host${w}"
conda deactivate
conda activate R412
cd $outputdr
echo "${r}"
echo "before DEG_Anno_Plot.R "
read -p "PRESS ENTER"
echo "${w}"
echo $MTDIR
echo $outputdr
echo $inputdr
echo $hostid 
echo $metadata

Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/host_counts.txt $inputdr/samplesheet.csv $hostid $MTDIR/HostSpecies.csv $metadata
#Aqui o arquivo definido pela variavel $metadata pode causar erros na analise DE, principalemnte se tiver grupos com apenas 1 fator, melhor rodar sem o $metadata e usar apenas do samplessheet.csv

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>     [75%]'

echo "ssGSEA${w}"
Rscript $MTDIR/gct_making.R $outputdr/Host_DEG/host_counts_TPM.csv $inputdr/samplesheet.csv

Rscript $MTDIR/Tools/ssGSEA2.0/ssgsea-cli.R -i $outputdr/ssGSEA/host.gct -o $outputdr/ssGSEA/ssgsea-results -d $MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt -y $MTDIR/Tools/ssGSEA2.0/config.yaml -u $threads

Rscript $MTDIR/for_halla.R $outputdr/ssGSEA/ssgsea-results-scores.gct $inputdr/samplesheet.csv $metadata

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>>    [80%]'
echo "MTD DEG analyses are done. Starting microbiome x host association analyses..."

echo "halla: association analysis${w}"
#mkdir -p $outputdr/Associations
conda deactivate
conda activate halla0820
echo "${g}Analyzing microbiome x host_genes associations...${w}"
#mkdir -p $outputdr/halla/host_gene # need to create a new directory for output to avoid "exists; deleting..." issue by halla
halla -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/host_gene --x_dataset_label Microbiomes --y_dataset_label Host_gene --diagnostic_plot -m ${pdm}

#O script abaixo parece ser inutil, pois so gera um heatmap com as medias do hostgene e micromiomas, esse script nao e original foi em quem fiz
#python $MTDIR/generate_halla_heatmap.py -m $outputdr/halla/Microbiomes.txt -g $outputdr/halla/Host_gene.txt  -o $outputdr/halla/host_gene/hallagram_all.pdf

#Abaixo uma abordagem diferente para verificar as associaçoes entre host gene vs microbiomas
python $MTDIR/pls_da_analysis.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/pls_da_results.pdf

#usando k-means
python $MTDIR/kmeans_clustering.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/kmeans_results.pdf -k 3

#Abaixo as tres opcoes de correlacoes para ver se alguma da um valor significativo pois a default de pearson nao deu nada
halla -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/pearson --x_dataset_label Microbiomes --y_dataset_label Host_gene --diagnostic_plot -m pearson --num_threads 12
halla -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/spearman --x_dataset_label Microbiomes --y_dataset_label Host_gene --diagnostic_plot -m spearman --num_threads 12
 
   # show all clusters
    if [[ $pdm == "spearman" ]]; then
        pdm_name='Pairwise Spearman'
    elif [[ $pdm == "pearson" ]]; then
        pdm_name='Pairwise Pearson'
    elif [[ $pdm == "mi" ]]; then
        pdm_name='mi'
    elif [[ $pdm == "nmi" ]]; then
        pdm_name='nmi'
    elif [[ $pdm == "xicor" ]]; then
        pdm_name='xicor'
    elif [[ $pdm == "dcor" ]]; then
        pdm_name='dcor'
    fi

#    hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_all.png --block_num -1
    hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_all.pdf --block_num -1
        # if hallagram_all.png not exist, show top 300 blocks
        if [[ ! -f $outputdr/halla/host_gene/hallagram_all.pdf ]]; then
hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_Top5.pdf --block_num 5
hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_Top10.pdf --block_num 10
hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_Top25.pdf --block_num 25
hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_Top50.pdf --block_num 50
hallagram -i $outputdr/halla/host_gene --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_gene --output $outputdr/halla/host_gene/hallagram_Top300.pdf --block_num 300
        fi
echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'

echo 'Analyzing microbiome x host_pathways associations...'
echo "${w}"
# for microbiome x host_pathways(ssGSEA)
#mkdir -p $outputdr/halla/pathway
halla -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_score.txt -o $outputdr/halla/pathway --x_dataset_label Microbiomes --y_dataset_label Host_pathway --diagnostic_plot -m ${pdm}

# show all clusters
hallagram -i $outputdr/halla/pathway --cbar_label "${pdm_name[@]}" --x_dataset_label Microbiomes --y_dataset_label Host_pathway --output $outputdr/halla/pathway_hallagram_all.pdf --block_num -1
echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "MTD running is finished"
echo -e "${w}"
