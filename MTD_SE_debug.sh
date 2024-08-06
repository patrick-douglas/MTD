#!/bin/bash

# default settings
pdm="spearman" # method in HALLA
length=35 # read length trimming by fastp
read_len=75 # the read length in bracken
threads=`nproc`
while getopts i:o:h:m:p:l:r:b:t option
do
    case "${option}" in
        i) inputdr=${OPTARG};;
        o) outputdr=${OPTARG};;
        h) hostid=${OPTARG};;
        m) metadata=${OPTARG};;
        p) pdm=${OPTARG};;
        l) length=${OPTARG};;
        r) read_len=${OPTARG};;
        b) blast=${OPTARG};;
        t) no_trimm=1;;
    esac
done

# inputdr=~/RNAseq_raw_data/samplesheet.csv # select input directory; must store singe-end .fq.gz (eg. DJ01.fq.gz) of each sample in the same folder as the samplesheet.csv
# outputdr=~/MTD_Results/test1 # select outputdr directory
# hostid=9544 # Enter host species taxonomy ID; initally supporting 9544 (rhesus monkey), 9606 (human), and 10090 (mouse).
# threads=20 # CPU threads; suggest >=16, eg. 20
# pdm= spearman or pearson or mi or nmi or xicor or dcor # pairwise distance metrics refer to HALLA mannual

# get MTD.sh script file path (in the MTD folder)
MTDIR=$(dirname $(readlink -f $0))
#parentname="$(dirname "$MTDIR")"
echo "MTD directory is $MTDIR"

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
echo "Conda path is $condapath"

# activate MTD conda environment
echo "Activating MTD conda environment..."
source $condapath/etc/profile.d/conda.sh
conda deactivate # avoid multiple conda environment
conda activate MTD
echo "Conda environment MTD activated."

inputdr=$(dirname $inputdr)
mkdir -p $outputdr
mkdir -p $outputdr/temp
cd $outputdr/temp
echo "Directories created and moved to $outputdr/temp."

# Step 0: Host database auto selection
echo "Step 0: Auto selecting host database based on hostid $hostid..."
if [[ $hostid == 9606 ]]; then
    DB_host=$MTDIR/kraken2DB_human # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_human/genome_tran #for hisat2
    DB_blast=$MTDIR/human_blastdb/human_blastdb # for blast
    gtf=$MTDIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz # for featureCounts
    echo "Selected human database."
elif [[ $hostid == 9544 ]]; then
    DB_host=$MTDIR/kraken2DB_rhesus # for kraken2
    DB_hisat2=$MTDIR/hisat2_index_rhesus/genome_tran #for hisat2
    DB_blast=$MTDIR/rhesus_blastdb/rhesus_blastdb # for blast
    gtf=$MTDIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz # for featureCounts
    echo "Selected rhesus monkey database."
elif [[ $hostid == 10090 ]]; then
    DB_host=$MTDIR/kraken2DB_mice
    DB_hisat2=$MTDIR/hisat2_index_mouse/genome_tran
    DB_blast=$MTDIR/mouse_blastdb/mouse_blastdb # for blast
    gtf=$MTDIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz # for featureCounts
    echo "Selected mouse database."
elif [[ -d "$MTDIR/kraken2DB_${hostid}" ]]; then # test if customized host species exist
    DB_host=$MTDIR/kraken2DB_${hostid}
    DB_hisat2=$MTDIR/hisat2_index_${hostid}/genome_tran
    DB_blast=$MTDIR/blastdb_${hostid}/blastdb_${hostid} # for blast
    gtf=$MTDIR/ref_${hostid}/ref_${hostid}.gtf.gz
    echo "Selected customized host database."
else
    echo "Host species is not supported. You can use bash Customized_host.sh for building."
    exit 1
fi
DB_micro=$MTDIR/kraken2DB_micro # customized kraken database for microbiome
echo "Database paths set up."

# for SRR input samples in the samplesheet.csv; download SRR samples
echo "Checking for SRR samples in $inputdr..."
cd $inputdr
if [ ! -z "$(cat samplesheet.csv | cut -f 1 -d ','| grep ^SRR)" ]; then
    for s in $(cat samplesheet.csv | cut -f 1 -d ','| grep ^SRR); do
        # check if fastq files of SRR sample exists
        if [[ ! -f ${s}_1.fastq || ! -f ${s}_2.fastq ]]; then
            echo "File ${s} fastq files NOT exists. Start downloading..."
            echo 'Downloading SRA files...'
            prefetch -X 999G ${s}
            echo 'Splitting SRA files to fastq files...'
            fasterq-dump -p --split-files ${s}
            rm -rf ${s}
            echo "Fastq files for ${s} downloaded and split."
        else
            echo "Fastq files for ${s} already exist."
        fi
    done
else
    echo "No SRR samples found in samplesheet.csv."
fi
cd $outputdr/temp
echo "Finished checking and processing SRR samples."

# To extract sample names from input fastq files (support .fq.gz, .fastq.gz, .fq, or .fastq)
echo "Extracting sample names from input fastq files..."
files=$(find $inputdr -name "*.fq.gz" -or -name "*.fastq.gz" -or -name "*.fq" -or -name "*.fastq" -type f)
for i in $files; do
    fn=$(basename $i) # Extract file name, eg. DJ01_1.fq.gz
    sn=$(echo $fn | awk -F '_R1' '{print $(NF-1)}') # Extract sample name, eg. DJ01
    lsn=$lsn" "$sn # Make a list of sample names
done
echo "Sample names extracted."

# check if input files match the samplesheet.csv
echo "Checking if input files match the samplesheet.csv..."
fastq_files=$(echo $lsn | tr " " "\n" | sort | tr "\n" " ")
SamplesInSheet=$(cat $inputdr/samplesheet.csv | cut -f 1 -d ',' | tail -n +2 | sort | tr "\n" " ")
if [[ "$fastq_files" != "$SamplesInSheet" ]]; then
    echo "The samples' fastq files in the input folder do not match with your samplesheet.csv"
    echo "Please double-check with the samplesheet.csv and input files. Ensure no other fastq files are under the input folder and its subfolders. Refer to the user guide on https://github.com/FEI38750/MTD."
    exit 1
fi
echo "Input files match with samplesheet.csv."
echo 'MTD running progress:'
echo '>>                  [10%] - Raw reads trimming'

# Raw reads trimming
max_jobs=$(nproc)
max_fastp_cores=16

if [ "$threads" -gt "$max_fastp_cores" ]; then
    fastp_threads=$max_fastp_cores
else
    fastp_threads=$threads
fi

# Processar cada amostra
for i in $lsn; do
    echo "Processing sample $i..."
    
    # Encontre o arquivo fastq correspondente (suporta .fq.gz, .fastq.gz, .fq, ou .fastq)
    fq=$(find $inputdr -name "${i}*.fq.gz" -o -name "${i}*.fastq.gz" -o -name "${i}*.fq" -o -name "${i}*.fastq" -type f)

    if [ -z "$no_trimm" ]; then
        # Se no_trimm não for definido, use o fastp para limpeza dos dados
        echo "Trimming fastq files for sample $i with fastp"
        fastp --trim_poly_x \
              --length_required $length \
              --thread $fastp_threads \
              -i $fq \
              -o $outputdr/temp/Trimmed_${i}.fq.gz
        echo "Sample $i trimmed and saved as Trimmed_${i}.fq.gz"
    else
        echo "No trimming required for sample $i."
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
            echo "Compressed file $input_file copied to $output_file"
        else
            # Se o arquivo não estiver comprimido, remova o sufixo e comprima
            base_name_no_suffix=$(echo "$base_name" | sed "s/_R[0-9]$//")
            output_file="$2/Trimmed_${base_name_no_suffix}.fq.gz"
            gzip -c "$input_file" > "$output_file"
            echo "File $input_file compressed to $output_file"
        fi
    ' _ {} "$outputdr/temp"
    echo 'Compression and copying completed.'
fi

echo 'MTD running progress:'
echo '>>>>                [20%] - Reads classification by kraken2; 1st step for host'

# Reads classification by kraken2; 1st step for host
for i in $lsn; do
    echo "Classifying reads for sample $i with kraken2 (host database)"
    kraken2 --db $DB_host --use-names \
        --report Report_host_$i.txt \
        --threads $threads \
        --gzip-compressed \
        --classified-out ${i}_host.fq \
        --unclassified-out ${i}_non-host_raw.fq \
        Trimmed_${i}.fq.gz \
        > Report_host_$i.kraken
    echo "Sample $i classified; reports saved as Report_host_$i.kraken and Report_host_$i.txt"
done

echo 'MTD running progress:'
echo '>>>>>               [25%] - Reads classification by kraken2; 2nd step for non-host reads'

# Reads classification by kraken2; 2nd step for non-host reads
for i in $lsn; do
    echo "Classifying non-host reads for sample $i with kraken2 (microbiome database)"
    kraken2 --db $DB_micro --use-names \
        --report Report_non-host.raw_$i.txt \
        --threads $threads \
        --classified-out ${i}_raw_cseqs.fq \
        --unclassified-out ${i}_raw_ucseqs.fq \
        ${i}_non-host_raw.fq \
        > Report_non-host_raw_$i.kraken
    echo "Non-host reads for sample $i classified; reports saved as Report_non-host_raw_$i.kraken and Report_non-host.raw_$i.txt"
done

echo 'MTD running progress:'
echo '>>>>>>              [30%] - Process completed'
# Decontamination step
conta_file=$MTDIR/conta_ls.txt
if test -f "$conta_file"; then
    echo "Decontamination file found. Processing decontamination..."

    tls=$(awk -F '\t' '{print $2}' $conta_file)
    conta_ls="${tls//$'\r\n'/ }"
    
    for i in $lsn; do
        echo "Extracting reads for sample $i using decontamination list"
        python $MTDIR/Tools/KrakenTools/extract_kraken_reads.py \
            -k Report_non-host_raw_${i}.kraken \
            -s1 ${i}_non-host_raw.fq \
            -o ${i}_non-host.fq \
            -r Report_non-host.raw_${i}.txt \
            --taxid $conta_ls --exclude --include-children
        echo "Decontaminated reads for sample $i saved as ${i}_non-host.fq"
    done

    echo 'MTD running progress:'
    echo '>>>>>>>             [35%] - Decontamination complete, proceeding with classification'

    # Reads classification by kraken2; 3rd step for decontaminated non-host reads to get reports
    for i in $lsn; do
        echo "Classifying decontaminated reads for sample $i with kraken2 (microbiome database)"
        kraken2 --db $DB_micro --use-names \
            --report Report_non-host_$i.txt \
            --threads $threads \
            --classified-out ${i}_cseqs.fq \
            --unclassified-out ${i}_ucseqs.fq \
            ${i}_non-host.fq \
            > Report_non-host_$i.kraken
        echo "Decontaminated reads for sample $i classified; reports saved as Report_non-host_$i.kraken and Report_non-host_$i.txt"
    done
fi

echo 'MTD running progress:'
echo '>>>>>>>>            [40%] - Bracken analysis'

# Bracken analysis
for i in $lsn; do
    echo "Running Bracken analysis for sample $i"
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.phylum.bracken -r $read_len -l P -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.genus.bracken -r $read_len -l G -t $threads
    bracken -d $DB_micro -i Report_non-host_${i}.txt -o Report_$i.species.bracken -r $read_len -l S -t $threads
    echo "Bracken analysis for sample $i completed; reports saved as Report_$i.phylum.bracken, Report_$i.genus.bracken, and Report_$i.species.bracken"
done

echo 'MTD running progress:'
echo '>>>>>>>>>           [45%] - Combining Bracken files and normalization'

# Combined .bracken files (table-like) into a single output directory for DESeq2
echo "Combining Bracken output files into single files"
python $MTDIR/Tools/combine_bracken_outputs.py --files *.phylum.bracken -o $outputdr/bracken_phylum_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.genus.bracken -o $outputdr/bracken_genus_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.species.bracken -o $outputdr/bracken_species_all
echo "Bracken files combined; results saved in $outputdr"

# Move _bracken report files (tree-like) to a separate folder
echo "Organizing Bracken report files"
mkdir -p Report_non-host_bracken_species_normalized
mv *_bracken_species.txt Report_non-host_bracken_species_normalized
cd Report_non-host_bracken_species_normalized

# Trim the name of _bracken report files (tree-like) to the sample name (e.g., DJ01)
for i in $lsn; do
    echo "Renaming Bracken report files for sample $i"
    mv *${i}_* $i
done

# Convert original _bracken report files (tree-like) into .biom file for ANCOMBC and diversity analysis in phyloseq (R) etc. in DEG_Anno_Plot.R
echo "Converting Bracken report files to .biom format"
kraken-biom * -o $outputdr/temp/bracken_species_all0.biom --fmt json

# Adjust bracken file (tree-like) by normalizing read counts; for additional visualization (.biom, .mpa, .krona)
echo "Normalizing Bracken data"
conda deactivate
conda activate R412
Rscript $MTDIR/Normalization_afbr.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $outputdr/temp/Report_non-host_bracken_species_normalized $metadata
conda deactivate
conda activate MTD

echo 'MTD running progress:'
echo '>>>>>>>>>>          [50%] - Normalization complete and process finalized'
# Converted adjusted _bracken report files (tree-like) into .biom file for graph visualization: graphlan, MPA, krona
echo "Converting adjusted _bracken report files to .biom format for visualization"
kraken-biom * -o $outputdr/bracken_species_all.biom --fmt json

# Remove "sp. " in the .biom file; correct improper format before running export2graphlan.py
echo "Removing 'sp. ' from .biom file to correct format"
sed -i 's/sp. //g' $outputdr/bracken_species_all.biom

# Go to temp folder
echo "Navigating to temp folder and preparing for graph visualization"
cd ../
mkdir -p ../graphlan
cd ../graphlan

# Activate conda environment for graphlan
echo "Activating conda environment for graphlan"
conda deactivate
conda activate py2

# Export to graphlan format
echo "Exporting .biom file to graphlan format"
python $MTDIR/Tools/export2graphlan/export2graphlan.py \
    -i ../bracken_species_all.biom \
    -a annot.txt -t tree.txt \
    --discard_otus --most_abundant 50 \
    --annotations 2,3,4,5,6 \
    --external_annotations 7 --internal_levels --max_clade_size 300

echo "Graphlan format export complete"
conda deactivate
conda activate MTD

# Move to temp directory
cd ../temp

# DEG & Annotation & Plots & Diversity & Preprocess for Microbiome
echo "Running DEG & Annotation & Plotting & Diversity Analysis"
conda deactivate
conda activate R412
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/bracken_species_all $inputdr/samplesheet.csv $hostid $MTDIR/HostSpecies.csv $metadata
echo "DEG & Annotation & Plotting complete"
conda deactivate
conda activate MTD

# Organize and move HUMAnN3 results
echo "Organizing HUMAnN3 results"
cd $outputdr/temp
mkdir -p bracken_raw_results # Save the raw output from bracken (table-like)
mv ../bracken_*_all bracken_raw_results

cd ../graphlan
echo "Annotating and generating Graphlan image"
python $MTDIR/Tools/graphlan/graphlan_annotate.py --annot annot.txt tree.txt outtree.txt # Attach annotation to the tree
python $MTDIR/Tools/graphlan/graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png # Generate the Graphlan image
echo "Graphlan image generated"

# Visualization preprocess for Krona
echo "Preprocessing for Krona visualization"
mkdir -p ../krona
for i in $lsn; do # Store input sample name in i; e.g., DJ01
    python $MTDIR/Tools/KrakenTools/kreport2krona.py \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ../krona/${i}-bracken.krona
    echo "Krona report for sample $i generated"
done

# To make MPA style file
echo "Converting Bracken reports to MPA format"
for i in $lsn; do # Store input sample name in i; e.g., DJ01
    python $MTDIR/Tools/KrakenTools/kreport2mpa.py \
        --display-header \
        -r Report_non-host_bracken_species_normalized/${i} \
        -o ${i}-bracken.mpa.txt
    echo "MPA file for sample $i created"
done

# Combine MPA files
echo "Combining MPA files"
python $MTDIR/Tools/KrakenTools/combine_mpa.py \
    -i *.mpa.txt \
    -o ../Combined.mpa

echo 'MTD running progress:'
echo '>>>>>>>>>>>>        [55%] - HUMAnN3 and MPA files processing complete'

# HUMAnN3 processing
echo "Running HUMAnN3 analysis"
mkdir -p HUMAnN_output

for n1 in *_non-host.fq; do
    cp $n1 HUMAnN_output/$n1
done

cd HUMAnN_output

for file in *; do # Trim the file name
    mv $file ${file/_non-host/}
done

# Run HUMAnN3
for i in *.fq; do
    echo "Running HUMAnN3 for file $i"
    humann --input $i \
        --output hmn_output \
        --threads $threads \
        --verbose
    echo "HUMAnN3 analysis for file $i complete"
done

echo 'MTD running progress:'
echo '>>>>>>>>>>>>        [60%] - HUMAnN3 analysis complete'

# Join all gene family and pathway abundance files
echo "Joining gene family and pathway abundance files"
humann_join_tables -i hmn_output/ -o humann_pathabundance.tsv --file_name pathabundance
humann_join_tables -i hmn_output/ -o humann_genefamilies.tsv --file_name genefamilies

# Normalize RPKs to "relab" (relative abundance)
echo "Normalizing RPKs to relative abundance"
humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_relab.tsv --units relab --update-snames
humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_relab.tsv --units relab --update-snames

# Generate stratified tables
echo "Generating stratified tables"
humann_split_stratified_table --input humann_pathabundance_relab.tsv --output ./
humann_split_stratified_table --input humann_genefamilies_relab.tsv --output ./
humann_split_stratified_table --input humann_pathabundance.tsv --output ./
humann_split_stratified_table --input humann_genefamilies.tsv --output ./

# Regroup gene families table into KEGG orthologs and GO terms
echo "Regrouping gene families into KEGG orthologs and GO terms"
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_ko \
    --output humann_genefamilies_relAbundance_kegg.tsv
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_go \
    --output humann_genefamilies_relAbundance_go.tsv
humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_ko \
    --output humann_genefamilies_Abundance_kegg.tsv
humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_go \
    --output humann_genefamilies_Abundance_go.tsv

# Translate KEGG and GO ID to human readable terms
echo "Translating KEGG and GO IDs to human readable terms"
conda deactivate
conda activate R412
Rscript $MTDIR/humann_ID_translation.R \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_kegg.tsv \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_go.tsv \
    $MTDIR
Rscript $MTDIR/humann_ID_translation.R \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_kegg.tsv \
    $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_go.tsv \
    $MTDIR
conda deactivate
conda activate MTD

# Cleaning up file structure
echo "Cleaning up file structure"
mkdir $outputdr/hmn_pathway_abundance_files
mkdir $outputdr/hmn_genefamily_abundance_files
mv *pathabundance* $outputdr/hmn_pathway_abundance_files/
mv *genefamilies* $outputdr/hmn_genefamily_abundance_files/

# DEG & Annotation & Plots & Diversity & Preprocess
echo "Running DEG & Annotation & Plotting for HUMAnN3 results"
cd $outputdr/hmn_genefamily_abundance_files
conda deactivate
conda activate R412
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg_translated.tsv $inputdr/samplesheet.csv
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go_translated.tsv $inputdr/samplesheet.csv
conda deactivate
conda activate MTD

# Optionally generate barplots (commented out)
# echo "Generating HUMAnN3 barplots"
# humann_barplot --input $outputdr/hmn_pathway_abundance_files/humann_pathabundance_cpm_stratified.tsv \
#     --focal-metadatum Group --last-metadatum Group \
#     --focal-feature PWY-3781 \
#     --output $outputdr/hmn_pathway_abundance_files/humann_pathabundance_barplot.png
# humann_barplot --input $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_cpm_stratified.tsv \
#     --output $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_barplot.png

echo 'Starting to process the host reads...'
echo '>>>>>>>>>>>>        [65%] - Processing host reads'

## Continue to process the host reads
cd $outputdr/temp

if [[ $blast == blast ]]; then
    # Magic-BLAST
    echo "Running Magic-BLAST on host reads"
    for i in $lsn; do # store input sample name in i; eg. DJ01
        magicblast -query ${i}_host.fq \
        -db $DB_blast \
        -infmt fastq \
        -out $i.sam \
        -num_threads $threads
        echo "Magic-BLAST for sample $i completed"
    done
else
    # HISAT2 alignment
    echo "Running HISAT2 alignment on host reads"
    for i in $lsn; do # store input sample name in i; eg. DJ01
        hisat2 -p $threads -q \
            -x $DB_hisat2 \
            --summary-file ${i}_hisat2_summary.txt \
            -U Trimmed_${i}.fq.gz \
            -S $i.sam
        echo "HISAT2 alignment for sample $i completed"
    done
fi

# featureCounts
echo "Running featureCounts on SAM files"
featureCounts -T $threads \
   -a $gtf \
   -o $outputdr/host_counts.txt \
   *.sam

# Convert SAM to BAM, sort, and index
for i in $lsn; do
    echo "Processing SAM to BAM for sample $i"
    samtools view -bS $i.sam > $i.bam -@ $threads
    samtools sort $i.bam -o $i.sorted.bam -@ $threads
    samtools index $i.sorted.bam -@ $threads
    echo "SAM to BAM conversion and sorting for sample $i completed"
done

mkdir -p BAM
mv *.sorted.bam *.sorted.bam.bai BAM/

cd $outputdr
# Trim the featureCounts output for downstream analysis
echo "Trimming the featureCounts output for downstream analysis"
sed '1d; 2 s/\.sam//g' host_counts.txt > tmpfile
mv tmpfile host_counts.txt

# DEG & Annotation & Plots & preprocess for host
echo "Running DEG & Annotation & Plotting for host counts"
conda deactivate
conda activate R412
Rscript $MTDIR/DEG_Anno_Plot.R $outputdr/host_counts.txt $inputdr/samplesheet.csv $hostid $MTDIR/HostSpecies.csv $metadata
conda deactivate
conda activate MTD

echo 'MTD running progress:'
echo '>>>>>>>>>>>>>>>>    [75%] - DEG & Annotation & Plotting complete'

# ssGSEA
echo "Running ssGSEA analysis"
Rscript $MTDIR/gct_making.R $outputdr/Host_DEG/host_counts_TPM.csv $inputdr/samplesheet.csv

Rscript $MTDIR/Tools/ssGSEA2.0/ssgsea-cli.R \
    -i $outputdr/ssGSEA/host.gct \
    -o $outputdr/ssGSEA/ssgsea-results \
    -d $MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt \
    -y $MTDIR/Tools/ssGSEA2.0/config.yaml \
    -u $threads
echo "ssGSEA analysis completed"

Rscript $MTDIR/for_halla.R $outputdr/ssGSEA/ssgsea-results-scores.gct $inputdr/samplesheet.csv $metadata

echo 'MTD running progress:'
echo '>>>>>>>>>>>>>>>>    [80%] - ssGSEA analysis complete'

echo "MTD DEG analyses are done. Starting microbiome x host association analyses..."

# halla: association analysis
echo "Analyzing microbiome x host gene associations"
conda deactivate
conda activate halla0820
halla -x $outputdr/halla/Microbiomes.txt \
    -y $outputdr/halla/Host_gene.txt \
    -o $outputdr/halla/host_gene \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_gene \
    --diagnostic_plot -m ${pdm}

    # Show all clusters
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

    hallagram \
        -i $outputdr/halla/host_gene \
        --cbar_label "${pdm_name[@]}" \
        --x_dataset_label Microbiomes \
        --y_dataset_label Host_gene \
        --output $outputdr/halla/host_gene/hallagram_all.png \
        --block_num -1

    # If hallagram_all.png not exist, show top 300 blocks
    if [[ ! -f $outputdr/halla/host_gene/hallagram_all.png ]]; then
        hallagram \
            -i $outputdr/halla/host_gene \
            --cbar_label "${pdm_name[@]}" \
            --x_dataset_label Microbiomes \
            --y_dataset_label Host_gene \
            --output $outputdr/halla/host_gene/hallagram_Top300.png \
            --block_num 300
    fi

echo 'MTD running progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%] - Microbiome x host gene associations complete'

echo 'Analyzing microbiome x host pathways associations...'
# for microbiome x host pathways (ssGSEA)
halla -x $outputdr/halla/Microbiomes.txt \
    -y $outputdr/halla/Host_score.txt \
    -o $outputdr/halla/pathway \
    --x_dataset_label Microbiomes \
    --y_dataset_label Host_pathway \
    --diagnostic_plot -m ${pdm}

    # Show all clusters
    hallagram \
        -i $outputdr/halla/pathway \
        --cbar_label "${pdm_name[@]}" \
        --x_dataset_label Microbiomes \
        --y_dataset_label Host_pathway \
        --output $outputdr/halla/pathway_hallagram_all.png \
        --block_num -1

echo 'MTD running progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%] - Microbiome x host pathways associations complete'
echo "MTD running is finished"

