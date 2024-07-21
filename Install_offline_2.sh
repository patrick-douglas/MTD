#!/bin/bash

# Defining colors
w=$(tput sgr0)
r=$(tput setaf 1)
g=$(tput setaf 2)
y=$(tput setaf 3)
p=$(tput setaf 5)
echo "${w}"

# Setting default values
kmer="" # --kmer-len in kraken2-build
min_l="" # --minimizer-len in kraken2-build
min_s="" # --minimizer-spaces in kraken2-build
read_len=75 # the read length in bracken-build
threads=`nproc`
condapath=~/miniconda3
offline_files_folder=""
logfile="installation_log.txt"

# Function to display usage message
usage() {
    echo "Usage: $0 -p <condapath> -o <offline_files_folder> [-k <kmer>] [-m <minimizer_length>] [-s <minimizer_spaces>] [-r <read_length>]"
    exit 1
}

# Checking if the required parameters are provided
if [ $# -lt 4 ]; then
    usage
fi

# Processing arguments
while getopts ":p:o:k:m:s:r:" option; do
    case "${option}" in
        p)
            condapath=${OPTARG}
            ;;
        o)
            offline_files_folder=${OPTARG}
            ;;
        k)
            kmer=${OPTARG}
            ;;
        m)
            min_l=${OPTARG}
            ;;
        s)
            min_s=${OPTARG}
            ;;
        r)
            read_len=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# Verifying if the required parameters are provided
if [ -z "$condapath" ] || [ -z "$offline_files_folder" ]; then
    usage
fi

# Initialize log file
echo "Installation Log - $(date)" > $logfile

# Function to check if the last command was successful
check_status() {
    if [ $? -ne 0 ]; then
        echo "$1 failed!" | tee -a $logfile
        exit 1
    else
        echo "$1 succeeded!" | tee -a $logfile
    fi
}

# Get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place
touch condaPath
echo "$condapath" > $dir/condaPath

source $condapath/etc/profile.d/conda.sh

# Update and install Conda environments
echo 'Installing Conda environments...' | tee -a $logfile
sudo apt-get update
check_status "APT update"

conda deactivate
conda env create -f Installation/MTD.yml
check_status "Creating MTD environment"
conda env create -f Installation/py2.yml
check_status "Creating py2 environment"
conda env create -f Installation/halla0820.yml
check_status "Creating halla0820 environment"
conda activate halla0820
pip install --upgrade setuptools pip
pip install -r Installation/pip.requirements
pip install jenkspy matplotlib numpy pandas PyYAML scipy seaborn
pip install --no-deps halla==0.8.20
check_status "Installing Python packages"
conda deactivate
conda env create -f Installation/R412.yml
check_status "Creating R412 environment"

echo 'MTD installation progress: [15%]' | tee -a $logfile

# Check R package installations
conda activate halla0820
conda install -n halla0820 -y -c conda-forge pkg-config
R -e "install.packages('lattice',repos = 'http://cran.us.r-project.org', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages(c('XICOR','mclust','BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e 'BiocManager::install("preprocessCore", ask = FALSE)'
R -e "install.packages('eva', INSTALL_opts = '--no-lock', repos='http://cran.us.r-project.org', Ncpus=$threads)"
check_status "Installing R packages"

# Check R package versions
R_ver=`R --version | grep version | grep R | awk '{print $3}'`
echo "${g}*************************************************"
echo "${w}R $R_ver packages versions${g}"
echo "*************************************************${g}"
for package in lattice MASS Matrix XICOR mclust BiocManager eva; do
    version=$(R --no-restore -e "packageVersion('$package')" | grep packageVersion -A 1 | grep '[1]' | awk '{print $2}' | sed -r 's/^.{1}//' | sed 's/.$//')
    echo -n "${g}$package: ${w}"; echo $version
done
echo "${g}*************************************************${w}" | tee -a $logfile

conda deactivate

echo 'Conda environments installed' | tee -a $logfile

echo 'MTD installation progress: [20%]' | tee -a $logfile
echo 'Downloading virome database...' | tee -a $logfile
conda activate MTD
sudo apt-get update
sudo apt-get install rsync -y
check_status "Installing rsync"
conda deactivate
conda activate MTD
conda install -n MTD -y -c conda-forge pkg-config
check_status "Installing pkg-config in MTD environment"

# Download and prepare virome database
cp -f $offline_files_folder/virushostdb.genomic.fna.gz .
unpigz -f virushostdb.genomic.fna.gz
cat Installation/M33262_SIVMM239.fa virushostdb.genomic.fna > viruses4kraken.fa

# Update rsync scripts for kraken2-build
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/download_genomic_library.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library.sh $condapath/envs/MTD/libexec/download_genomic_library.sh

echo 'MTD installation progress: [30%]' | tee -a $logfile
echo 'Preparing microbiome (virus, bacteria, archaea, protozoa, fungi, plasmid, UniVec_Core) database...' | tee -a $logfile
# Kraken2 database building - Microbiome

DBNAME=kraken2DB_micro
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 taxonomy for microbiome"

# Substituição de FTP por wget para o download da biblioteca de plasmídeos
echo 'Downloading plasmid files from FTP...' | tee -a $logfile
wget -r -nd -np -A 'plasmid.*' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/ -P $DBNAME/library/plasmid/
check_status "Downloading Kraken2 plasmid library"

kraken2-build --download-library archaea --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 archaea library"
kraken2-build --download-library bacteria --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 bacteria library"
kraken2-build --download-library protozoa --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 protozoa library"
kraken2-build --download-library fungi --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 fungi library"

kraken2-build --download-library UniVec_Core --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Downloading Kraken2 UniVec_Core library"
kraken2-build --add-to-library viruses4kraken.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Adding virome library to Kraken2 microbiome database"
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s
check_status "Building Kraken2 microbiome database"

# Prepare host databases
for host in rhesus human mouse; do
    DBNAME=kraken2DB_$host
    mkdir -p $DBNAME
    cd $DBNAME || exit
    cp $offline_files_folder/GCF_000001635.27_GRCm39_genomic.fna.gz .
    unpigz GCF_000001635.27_GRCm39_genomic.fna.gz
    mv GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.fa
    kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
    check_status "Downloading Kraken2 taxonomy for host $host"
    kraken2-build --add-to-library $DBNAME/GCF_000001635.27_GRCm39_genomic.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
    check_status "Adding host $host library to Kraken2 database"
    kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s
    check_status "Building Kraken2 host $host database"
    cd ..
done

# Prepare Bracken database
echo 'Preparing Bracken database...' | tee -a $logfile
if [[ -z $kmer ]]; then
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len
else
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len -k $kmer
fi
check_status "Building Bracken database"

# Install HUMAnN3 databases
echo 'Installing HUMAnN3 databases...' | tee -a $logfile
mkdir -p $dir/HUMAnN/ref_database/
cd $dir/HUMAnN/ref_database/
for file in full_chocophlan.v201901_v31.tar.gz uniref90_annotated_v201901b_full.tar.gz full_mapping_v201901b.tar.gz; do
    cp $offline_files_folder/$file .
    tar xzvf $file -C $(basename $file .tar.gz)/
    check_status "Extracting $file"
done
humann_config --update database_folders nucleotide $dir/HUMAnN/ref_database/chocophlan
humann_config --update database_folders protein $dir/HUMAnN/ref_database/uniref
humann_config --update database_folders utility_mapping $dir/HUMAnN/ref_database/utility_mapping
check_status "Updating HUMAnN3 database configuration"

# Download host GTF files
echo 'Downloading host references...' | tee -a $logfile
for ref in rhesus human mouse; do
    mkdir -p ref_$ref
    cp $offline_files_folder/Macaca_mulatta.Mmul_10.104.gtf.gz ref_rhesus
    cp $offline_files_folder/Homo_sapiens.GRCh38.104.gtf.gz ref_human
    cp $offline_files_folder/Mus_musculus.GRCm39.104.gtf.gz ref_mouse
done
check_status "Downloading GTF files"

# Build HISAT2 indexes
echo 'Building HISAT2 indexes...' | tee -a $logfile
for host in rhesus mouse; do
    mkdir -p hisat2_index_$host
    cd hisat2_index_$host || exit
    cp ../ref_$host/Macaca_mulatta.Mmul_10.104.gtf.gz .
    gzip -d Macaca_mulatta.Mmul_10.104.gtf.gz
    mv Macaca_mulatta.Mmul_10.104.gtf genome.gtf
    python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
    python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
    cp $offline_files_folder/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz .
    gzip -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
    mv Macaca_mulatta.Mmul_10.dna.toplevel.fa genome.fa
    hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
    check_status "Building HISAT2 index for $host"
    cd ..
done

echo 'Installation complete!' | tee -a $logfile

