#!/bin/bash
MTDIR=~/MTD
condapath=~/miniconda3
threads=`nproc`

# Função para mostrar instruções de uso
usage() {
    echo ""
    echo "Usage:"
    echo "  bash $0 --genome <genome.fa.gz> --gtf-file <annotations.gtf.gz> --customized <TaxonID> --offline-folder <path>"
    echo ""
    echo "Required options:"
    echo "  -d, --genome           Path to host genome FASTA (.fa.gz)"
    echo "  -g, --gtf-file         Path to GTF annotation file (.gtf.gz)"
    echo "  -c, --ncbi-taxon-id    NCBI Taxon ID of the host species"
    echo "  -o, --offline-folder   Folder where offline files will be stored"
    echo "      --help             Show this help message"
    echo ""
    echo "Example:"
    echo "  bash $0 --genome ~/genomes/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa.gz --gtf-file ~/annotations/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa.gz --ncbi-taxon-id 8839 --offline-folder ~/makeOrgPackageFromNCBI/"
    echo ""
    exit 1
}

# Define as opções curtas e longas
OPTIONS=c:d:g:o:
LONGOPTIONS=ncbi-taxon-id:,genome:,gtf-file:,offline-folder:,help

# Processa os argumentos
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    usage
fi

eval set -- "$PARSED"

# Variáveis
while true; do
    case "$1" in
        -c|--ncbi-taxon-id)
            customized="$2"
            shift 2
            ;;
        -d|--genome)
            download="$2"
            shift 2
            ;;
        -g|--gtf-file)
            gtf="$2"
            shift 2
            ;;
        -o|--offline-folder)
            offline_files_folder="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: invalid arguments: $1"
            usage
            ;;
    esac
done

# Verificação básica de argumentos obrigatórios
if [[ -z "$download" || -z "$gtf" || -z "$customized" || -z "$offline_files_folder" ]]; then
    echo "Error: Missing arguments"
    usage
fi

# Remove any previously created directories or files from previous runs
rm -rf $MTDIR/kraken2DB_${customized} $MTDIR/ref_${customized} $MTDIR/hisat2_index_${customized} $MTDIR/blastdb_${customized}

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh


conda activate MTD

# Kraken2 database building - Customized
DBNAME=kraken2DB_${customized}
rm -rf $DBNAME 
mkdir -p $DBNAME
cd $DBNAME
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Calidris_pugnax.ASM143184v1.dna.toplevel.fa.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz .
#wget -c $download #http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
cp -u $download .
#bash ~/MTD/Customized_host.sh -d http://ftp.ensembl.org/pub/release-111/fasta/myotis_lucifugus/dna/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz -c 59463 -g http://ftp.ensembl.org/pub/release-111/gtf/myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz

#Morcego offline
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz -c 59463

#Calidris pugnax
#bash  ~/MTD/Customized_host.sh -t 20 -d https://ftp.ensembl.org/pub/release-111/fasta/calidris_pugnax/dna/Calidris_pugnax.ASM143184v1.dna.toplevel.fa.gz -c 198806 -g https://ftp.ensembl.org/pub/release-111/gtf/calidris_pugnax/Calidris_pugnax.ASM143184v1.111.gtf.gz
#Rattus norvegicus
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/-c 10116

#Gallus gallus
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/Gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.gtf.gz -c 9031

unpigz *.fa.gz
mv *.fa genome_${customized}.fa

# Extraia o nome científico da espécie baseado no Taxon_ID
species_name=$(awk -F, -v taxid="$customized" '$1 == taxid {print $3}' "$MTDIR/HostSpecies.csv")

# Verifique se o nome da espécie foi encontrado
if [ -z "$species_name" ]; then
  echo "Error: species name not found for Taxon_ID $customized."
  exit 1
fi

# Extraia o assembly_name do cabeçalho da sequência de entrada
assembly_name=$(grep -m 1 '^>' genome_${customized}.fa | sed -n 's/.*dna:primary_assembly \([^ ]*\).*/\1/p')
echo ''
echo -e "Selected host species:\e[3m $species_name\e[0m"
#echo "Selected host species:$species_name"
echo "Taxon ID: $customized"
echo ''

# Check the header format of the FASTA file
header=$(head -n 1 genome_${customized}.fa)

if [[ "$header" == *"dna:primary_assembly"* || "$header" == *"dna:genescaffold"* ]]; then
    # Ensembl format detected
    echo "Detected Ensembl format. Modifying headers for Ensembl..."
    sed -i "s/^>\(.*\) dna:.* \(.*\):\(.*\):\(.*\):\(.*\) REF$/>kraken:taxid|${customized}|\\3 ${species_name} chromosome \\3, ${assembly_name} Primary Assembly/" genome_${customized}.fa
else
    # NCBI format detected
    echo "Detected NCBI format. Modifying headers for NCBI..."
    sed -i "s/^>\(.*\) \(.*\)/>kraken:taxid|${customized}|\\1 ${species_name} \\2/" genome_${customized}.fa
fi

echo "Header modification complete."
rm -rf $MTDIR/blastdb_$customized
mkdir -p $MTDIR/blastdb_$customized
#cp genome_${customized}.fa $MTDIR/blastdb_$customized 

cd ..
kraken2-build --download-taxonomy --threads $threads --db $DBNAME
kraken2-build --add-to-library $DBNAME/genome_${customized}.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME
# download host GTF
#wget -c $gtf -P ref_${customized} -O ref_${customized}.gtf.gz
echo "Coping GTF file to $ref_${customized}/ref_${customized}.gtf.gz"
rm -rf ref_${customized}
mkdir -p ref_${customized}
cp $gtf ref_${customized}
cd ref_${customized}
mv *.gtf.gz ref_${customized}.gtf.gz

#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Calidris_pugnax.ASM143184v1.111.gtf.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz .
cd ..
echo "Building indexes for hisat2"
rm -rf hisat2_index_${customized}
mkdir -p hisat2_index_${customized}
cd hisat2_index_${customized}
cp ../ref_${customized}/ref_${customized}.gtf.gz .
gzip -d *.gtf.gz
mv *.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
mv ../$DBNAME/genome_${customized}.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo "Creating blast databases for custom reference $customized"
mkdir -p $MTDIR/blastdb_$customized
cd $MTDIR/blastdb_$customized
#cp $DBNAME/genome_${customized}.fa .
cp $download . 
gunzip *.fa.gz 
mv *.fa blastdb_$customized

makeblastdb -in $MTDIR/blastdb_$customized/blastdb_$customized -dbtype nucl -out $MTDIR/blastdb_$customized/blastdb_$customized -parse_seqids


echo "Creating the annotation package for R412"
echo -e "Selected host species:\e[3m $species_name\e[0m"
echo "Taxon ID: $customized"
echo ''
conda activate R412
cd $MTDIR
rm -rf NCBI org.*eg*
################################################################################################

# Run the R script and capture stdout/stderr

RSCRIPT_OUTPUT=$(conda run -n R412 Rscript create_annotation_package.R -t $customized -o $offline_files_folder -c 2>&1)

# Optional: show output
echo "$RSCRIPT_OUTPUT"

# Extract the last line with the .sqlite file name
SQLITE_NAME=$(echo "$RSCRIPT_OUTPUT" | grep -oE 'org\.[A-Za-z]+\.eg\.sqlite' | tail -n1)

# Convert the .sqlite to .db (as expected by install.packages)
PKG_NAME="${SQLITE_NAME/sqlite/db}"

if [ -n "$PKG_NAME" ]; then
    echo "Installing R package: $PKG_NAME"
    R -e "install.packages('$PKG_NAME', repos = NULL, type = 'source')"
    rm -rf NCBI org.*eg*
else
    echo "Error: Could not extract the package name from the R script output."
    exit 1
fi

################################################################################################

exit 1
#the parameter -o 
#/media/me/18TB_BACKUP_LBN/lbn_workspace/RNA-Seq-LBN/viral-rna-seq/MTD/Compressed/MTD/makeOrgPackageFromNCBI/

#Rscript $dir/create_annotation_package.R -t $customized -d $dir
R -e "install.packages('org.Aplatyrhynchos.eg.sqlite', repos = NULL, type = 'source')"

conda deactivate

echo "Customized host reference building is done"
