#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Download specific genomic libraries for use with Kraken 2.
# Supported libraries were chosen based on support from NCBI's FTP site
#   in easily obtaining a good collection of genomic data.  Others may
#   be added upon popular demand.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

dir=~/MTD
LIBRARY_DIR="$dir/$KRAKEN2_DB_NAME/library"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="https://$NCBI_SERVER"
RSYNC_SERVER="rsync://$NCBI_SERVER"
THIS_DIR=$PWD

library_name="$1"
ftp_subdir=$library_name
library_file="library.fna"
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  library_file="library.faa"
fi

function download_file() {
  file="$1"
  if [ -n "$KRAKEN2_USE_FTP" ]
  then
    wget ${FTP_SERVER}${file}
  else
    rsync --no-motd ${RSYNC_SERVER}${file} .
  fi
}

case $library_name in
  "archaea" | "bacteria" | "viral" | "fungi" | "plant" | "human" | "protozoa")
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f assembly_summary.txt
    remote_dir_name=$library_name
    if [ "$library_name" = "human" ]; then
      remote_dir_name="vertebrate_mammalian/Homo_sapiens"
    fi
    if ! download_file "/genomes/refseq/$remote_dir_name/assembly_summary.txt"; then
      1>&2 echo "Error downloading assembly summary file for $library_name, exiting."
      exit 1
    fi
    if [ "$library_name" = "human" ]; then
      grep "Genome Reference Consortium" assembly_summary.txt > x
      mv -f x assembly_summary.txt
    fi
    rm -rf all/ library.f* manifest.txt rsync.err
    rsync_from_ncbi.pl assembly_summary.txt
    scan_fasta_file.pl $library_file >> prelim_map.txt
    ;;
  "plasmid")
    mkdir -p $LIBRARY_DIR/plasmid
    cd $LIBRARY_DIR/plasmid
    rm -f library.f* plasmid.*

    # Use local files instead of downloading via FTP
    echo -n "Using local files for plasmid library..."
    local_download_dir="/media/me/4TB_BACKUP_LBN/Compressed/MTD/Kraken2DB_micro/library/plasmid/"

    if [ -d "$local_download_dir" ]; then
      echo "Local directory found: $local_download_dir"

      # Copy .gz files to the destination directory
      cp "$local_download_dir"/*.gz "$LIBRARY_DIR/plasmid/"

      # Ensure we are in the destination directory
      cd $LIBRARY_DIR/plasmid/

      if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
        ls *.faa.gz > manifest.txt
      else
        ls *.fna.gz > manifest.txt
      fi

      if [ ! -s manifest.txt ]; then
        echo "Error: manifest.txt is empty or was not created correctly."
        exit 1
      fi

      cat manifest.txt | xargs -I{} gunzip -c {} > $library_file
      rm -f plasmid.* .listing
      scan_fasta_file.pl $library_file > prelim_map.txt
      echo " done."
    else
      1>&2 echo "Error: Local directory $local_download_dir not found."
      exit 1
    fi
    ;;
  "nr" | "nt")
    protein_lib=0
    if [ "$library_name" = "nr" ]; then
      protein_lib=1
    fi
    if (( protein_lib == 1 )) && [ -z "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is a protein database, and the Kraken DB specified is nucleotide"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f $library_name.gz
    echo -n "Downloading $library_name database from server... "
    download_file "/blast/db/FASTA/$library_name.gz"
    echo "done."
    echo -n "Uncompressing $library_name database..."
    gunzip $library_name.gz
    mv $library_name $library_file
    echo "done."
    echo -n "Parsing $library_name FASTA file..."
    # The nr/nt files tend to have non-standard sequence IDs, so
    # --lenient is used here.
    scan_fasta_file.pl --lenient $library_file >> prelim_map.txt
    echo "done."
    ;;
  "UniVec" | "UniVec_Core")
    if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is for nucleotide databases only"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    echo -n "Downloading $library_name data from server... "
    download_file "/pub/UniVec/$library_name"
    echo "done."
    # 28384: "other sequences"
    special_taxid=28384
    echo -n "Adding taxonomy ID of $special_taxid to all sequences... "
    sed -e "s/^>/>kraken:taxid|$special_taxid|/" $library_name > library.fna
    scan_fasta_file.pl library.fna > prelim_map.txt
    echo "done."
    ;;
  *)
    1>&2 echo "Unsupported library.  Valid options are: "
    1>&2 echo "  archaea bacteria viral fungi plant protozoa human plasmid"
    1>&2 echo "  nr nt UniVec UniVec_Core"
    exit 1
    ;;
esac

if [ -n "$KRAKEN2_MASK_LC" ]; then
  echo -n "Masking low-complexity regions of downloaded library..."
  mask_low_complexity.sh .
  echo " done."
fi

