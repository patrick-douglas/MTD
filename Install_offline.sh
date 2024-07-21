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
#threads=$(($(nproc) - 2))
condapath=~/miniconda3
offline_files_folder=""

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

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place
touch condaPath
echo "$condapath" > $dir/condaPath

source $condapath/etc/profile.d/conda.sh
sudo apt-get update
conda deactivate
echo 'installing conda environments...'
conda env create -f Installation/MTD.yml
conda env create -f Installation/py2.yml
conda env create -f Installation/halla0820.yml
conda activate halla0820
pip install --upgrade setuptools pip
pip install -r Installation/pip.requirements
pip install jenkspy matplotlib numpy pandas PyYAML scipy seaborn
pip install --no-deps halla==0.8.20
conda deactivate
conda env create -f Installation/R412.yml

echo 'MTD installation progress:'
echo '>>                  [10%]'

# conda activate py2 # install dependencies of py2 in case pip does work in conda yml
# pip install backports-functools-lru-cache==1.6.1 biom-format==2.0.1 cycler==0.10.0 h5py==2.10.0 hclust2==1.0.0 kiwisolver==1.1.0 matplotlib==2.2.5 numpy==1.16.6 pandas==0.24.2 pyparsing==2.4.7 pyqi==0.3.2 python-dateutil==2.8.1 pytz==2021.1 scipy==1.2.3 six==1.15.0 subprocess32==3.5.4
# conda deactivate

conda activate halla0820 # install dependencies of halla
#halla0820
conda install -n halla0820 -y -c conda-forge pkg-config
R -e "install.packages('lattice',repos = 'http://cran.us.r-project.org', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages(c('XICOR','mclust','BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e 'BiocManager::install("preprocessCore", ask = FALSE)'
R -e "install.packages('eva', INSTALL_opts = '--no-lock', repos='http://cran.us.r-project.org', Ncpus=$threads)"

#Check Dependencies installation
R_ver=`R --version | grep version | grep R | awk '{print $3}'`
echo "${g}*************************************************"
echo "${w}R $R_ver packages versions${g}"
echo "*************************************************${g}"
echo -n "${g}lattice:                      	        ${w}" ; R --no-restore -e 'packageVersion("lattice")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}MASS: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("MASS")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Matrix:                       	        ${w}" ; R --no-restore -e 'packageVersion("Matrix")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}XICOR: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("XICOR")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}mclust:                       	        ${w}" ; R --no-restore -e 'packageVersion("mclust")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocManager:                   	        ${w}" ; R --no-restore -e 'packageVersion("BiocManager")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}eva: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("eva")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo "${g}*************************************************${w}"

conda deactivate

echo 'conda environments installed'

echo 'MTD installation progress:'
echo '>>>                 [15%]'
echo 'downloading virome database...'
conda activate MTD
sudo apt-get update
sudo apt-get install rsync -y
conda deactivate
##conda install -y python=3.10 
#conda install -n MTD -y -c bioconda metaphlan=4.0.6=pyhca03a8a_0 #Instalar no env MTD
conda activate MTD
conda install -n MTD -y -c conda-forge pkg-config

#Check if the file exists and have the same size before download
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/virushostdb.genomic.fna.gz
cp -f $offline_files_folder/virushostdb.genomic.fna.gz .
#wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.genomic.fna.gz
unpigz -f virushostdb.genomic.fna.gz
cat Installation/M33262_SIVMM239.fa virushostdb.genomic.fna > viruses4kraken.fa

# debug rsync error of kraken2-build
# debug rsync error of kraken2-build
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/download_genomic_library.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library.sh $condapath/envs/MTD/libexec/download_genomic_library.sh

echo 'MTD installation progress:'
echo '>>>>                [20%]'
echo 'Preparing microbiome (virus, bacteria, archaea, protozoa, fungi, plasmid, UniVec_Core) database...'
# Kraken2 database building - Microbiome

DBNAME=kraken2DB_micro
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
#Use local files for bacteria
kraken2-build --download-library archaea --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/Installation/rsync_from_ncbi_bacteria.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
kraken2-build --download-library bacteria --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
kraken2-build --download-library protozoa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library fungi --threads $threads --db $DBNAME $kmer $min_l $min_s

#Use local files for plasmid
cp -f $dir/Installation/download_genomic_library_plasmid.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library_plasmid.sh $condapath/envs/MTD/libexec/download_genomic_library.sh
kraken2-build --download-library plasmid --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/Installation/download_genomic_library.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library.sh $condapath/envs/MTD/libexec/download_genomic_library.sh

kraken2-build --download-library UniVec_Core --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library viruses4kraken.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo 'MTD installation progress:'
echo '>>>>>>              [30%]'
echo 'Preparing host (human) database...'
# Kraken2 database building - Human
DBNAME=kraken2DB_human
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --download-library human --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo 'MTD installation progress:'
echo '>>>>>>>             [35%]'
echo 'Preparing host (mouse) database...'
# Kraken2 database building - Mouse
DBNAME=kraken2DB_mice
mkdir -p $DBNAME
cd $DBNAME
#wget -T 300 -t 5 -N --no-if-modified-since https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/GCF_000001635.27_GRCm39_genomic.fna.gz .
cp $offline_files_folder/GCF_000001635.27_GRCm39_genomic.fna.gz .
unpigz GCF_000001635.27_GRCm39_genomic.fna.gz
mv GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.fa
cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library $DBNAME/GCF_000001635.27_GRCm39_genomic.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo 'MTD installation progress:'
echo '>>>>>>>>            [40%]'
echo 'Preparing host (rhesus monkey) database...'
# Kraken2 database building - Rhesus macaque
DBNAME=kraken2DB_rhesus
mkdir -p $DBNAME
cd $DBNAME
#wget -T 300 -t 5 -N --no-if-modified-since https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/GCF_003339765.1_Mmul_10_genomic.fna.gz .
cp $offline_files_folder/GCF_003339765.1_Mmul_10_genomic.fna.gz .
unpigz GCF_003339765.1_Mmul_10_genomic.fna.gz
mv GCF_003339765.1_Mmul_10_genomic.fna GCF_003339765.1_Mmul_10_genomic.fa
cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library $DBNAME/GCF_003339765.1_Mmul_10_genomic.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo 'MTD installation progress:'
echo '>>>>>>>>>           [45%]'
echo 'Bracken database building...'
# Bracken database building
if [[ $kmer == "" ]]; then
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len
else
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len -k $kmer
fi

echo 'MTD installation progress:'
echo '>>>>>>>>>>>         [55%]'
echo 'installing HUMAnN3 databases...'
# install HUMAnN3 databases
mkdir -p $dir/HUMAnN/ref_database/
cd $dir/HUMAnN/ref_database/
#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v296_201901.tar.gz
#Link working but slow
#wget -T 300 -t 5 -N --no-if-modified-since http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_chocophlan.v296_201901.tar.gz
#Temporary cp solution
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/full_chocophlan.v296_201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/full_chocophlan.v201901_v31.tar.gz .
cp $offline_files_folder/full_chocophlan.v201901_v31.tar.gz .
#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz
#Link working but slow
#wget -T 300 -t 5 -N --no-if-modified-since http://cmprod1.cibio.unitn.it/databases/HUMAnN/uniref90_annotated_v201901.tar.gz
#Temporary cp solution
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/uniref90_annotated_v201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/uniref90_annotated_v201901b_full.tar.gz .
cp $offline_files_folder/uniref90_annotated_v201901b_full.tar.gz .

#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_v201901.tar.gz
#Link working but slow
#wget -c http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_mapping_v201901.tar.gz
#Link source forge patrick-douglas
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/HUMAnN/ref_database/full_mapping_v201901.tar.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/full_mapping_v201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/full_mapping_v201901b.tar.gz .
cp $offline_files_folder/full_mapping_v201901b.tar.gz .
mkdir -p $dir/HUMAnN/ref_database/chocophlan
#tar xzvf full_chocophlan.v296_201901.tar.gz -C chocophlan/
tar xzvf full_chocophlan.v201901_v31.tar.gz -C chocophlan/
#mkdir -p $dir/HUMAnN/ref_database/full_UniRef90
mkdir -p $dir/HUMAnN/ref_database/uniref
#tar xzvf uniref90_annotated_v201901.tar.gz -C full_UniRef90/
tar xzvf uniref90_annotated_v201901b_full.tar.gz -C uniref/
mkdir -p $dir/HUMAnN/ref_database/utility_mapping
#tar xzvf full_mapping_v201901.tar.gz -C utility_mapping/full_mapping_v201901b.tar.gz
tar xzvf full_mapping_v201901b.tar.gz -C utility_mapping/
cd $dir

humann_config --update database_folders nucleotide $dir/HUMAnN/ref_database/chocophlan
humann_config --update database_folders protein $dir/HUMAnN/ref_database/uniref
humann_config --update database_folders utility_mapping $dir/HUMAnN/ref_database/utility_mapping

echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>      [70%]'
echo 'Downloading host (default: rhesus, human, mouse) references...'
# install host references
# download host GTF
    # download rhesus macaque GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz -P ref_rhesus
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz -P ref_rhesus
mkdir -p ref_rhesus && cp $offline_files_folder/Macaca_mulatta.Mmul_10.104.gtf.gz ref_rhesus

    # download human GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz -P ref_human
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_human/Homo_sapiens.GRCh38.104.gtf.gz -P ref_human
mkdir -p ref_human && cp $offline_files_folder/Homo_sapiens.GRCh38.104.gtf.gz ref_human

    # download mouse GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz -P ref_mouse
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz -P ref_mouse
mkdir -p ref_mouse && cp $offline_files_folder/Mus_musculus.GRCm39.104.gtf.gz ref_mouse

# Building indexes for hisat2
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>     [75%]'
echo 'Building host indexes (rhesus monkey) for hisat2...'
# rhesus macaques
mkdir -p hisat2_index_rhesus
cd hisat2_index_rhesus
cp ../ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz .
gzip -d Macaca_mulatta.Mmul_10.104.gtf.gz
mv Macaca_mulatta.Mmul_10.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_rhesus/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
cp $offline_files_folder/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz .
gzip -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
mv Macaca_mulatta.Mmul_10.dna.toplevel.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>    [80%]'
echo 'Building host indexes (mouse) for hisat2...'
# mouse
mkdir -p hisat2_index_mouse
cd hisat2_index_mouse
cp ../ref_mouse/Mus_musculus.GRCm39.104.gtf.gz .
gzip -d Mus_musculus.GRCm39.104.gtf.gz
mv Mus_musculus.GRCm39.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
cp $offline_files_folder/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz .
gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>   [85%]'
echo 'Building host indexes (human) for hisat2...'
# human
mkdir -p hisat2_index_human
cd hisat2_index_human
cp ../ref_human/Homo_sapiens.GRCh38.104.gtf.gz .
gzip -d Homo_sapiens.GRCh38.104.gtf.gz
mv Homo_sapiens.GRCh38.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
cp $offline_files_folder/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

# # download preduild index for hisat2 # prebuild is from NCBI, may be not compatiable with featureCount
#     # H. sapiens
#     mkdir -p hisat2_index_human
#     cd hisat2_index_human
#     wget https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz
#     pigz -dc grch38_tran.tar.gz | tar xf -
#     cd ..

# Create a BLAST database for Magic-BLAST
makeblastdb -in $dir/hisat2_index_human/genome.fa -dbtype nucl -parse_seqids -out $dir/human_blastdb/human_blastdb
makeblastdb -in $dir/hisat2_index_mouse/genome.fa -dbtype nucl -parse_seqids -out $dir/mouse_blastdb/mouse_blastdb
makeblastdb -in $dir/hisat2_index_rhesus/genome.fa -dbtype nucl -parse_seqids -out $dir/rhesus_blastdb/rhesus_blastdb

echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'
echo 'installing R packages...'
# install R packages
conda deactivate
conda install -n py2 -y -c conda-forge pkg-config
conda activate R412
#$dir/update_fix/update_conda_pkgs.sh
conda install -n R412 -y -c conda-forge pkg-config
~/MTD/update_fix/Install.R.packages.R412.sh

# debug in case libcurl cannot be located in the conda R environment
wget -T 300 -t 5 -N --no-if-modified-since https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.2.tar.gz
# if /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc exists, use it
if [ -f /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc ]; then
    locate_lib=/usr/lib/x86_64-linux-gnu/pkgconfig
    else 
    locate_lib=$(dirname $(locate libcurl | grep '\.pc'))
fi
R CMD INSTALL --configure-vars='LIB_DIR='"$locate_lib" curl_4.3.2.tar.gz

Rscript $dir/Installation/R_packages_installation.R
~/MTD/update_fix/check_R_pkg.R412.sh
#~/miniconda3/envs/MTD/opt/krona/updateTaxonomy.sh
conda deactivate
conda activate MTD
#~/miniconda3/envs/MTD/opt/krona/updateTaxonomy.sh
chmod +x MTD.sh
cp $dir/update_fix/hclust2.py ~/miniconda3/envs/py2/lib/python2.7/site-packages/hclust2.py
echo "*********************************"
echo "R packages version for conda envs"
echo "*********************************"
conda run -n R412 ~/MTD/update_fix/check_R_pkg.R412.sh
conda run -n halla0820 ~/MTD/update_fix/check_R_pkg_halla0820.sh
echo "*********************************"
echo ""
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "MTD installation is finished"
