#!/bin/bash

# ------------------------------------------------------------
# Colors and text styles
# ------------------------------------------------------------

w=$(tput sgr0 2>/dev/null || true)
r=$(tput setaf 1 2>/dev/null || true)
g=$(tput setaf 2 2>/dev/null || true)
y=$(tput setaf 3 2>/dev/null || true)
p=$(tput setaf 5 2>/dev/null || true)

ital=$(tput sitm 2>/dev/null || printf '\033[3m')
noital=$(tput ritm 2>/dev/null || printf '\033[23m')

echo "${w}"

die() {
    echo "${r}[ERROR] $*${w}" >&2
    exit 1
}

require_file() {
    local f="$1"
    local label="${2:-file}"

    if [[ ! -s "$f" ]]; then
        echo "${r}[MISSING] $label: $f${w}" >&2
        exit 1
    fi
}

run_cmd() {
    echo "${g}[RUN]${w} $*"
    "$@" || die "Command failed: $*"
}

# ------------------------------------------------------------
# Default settings
# ------------------------------------------------------------

pdm="spearman"                         # HAllA method
length=35                              # fastp minimum read length
read_len=75                            # Bracken read length
threads="$(nproc)"                     # CPU threads
blast="hisat"                          # default host alignment method
no_trimm=0                             # default flag
metadata=""                            # optional metadata file

# Optional Kraken2 host-filtering DB override.
# If empty, DB_host is selected automatically from --hostid.
# If provided through --kraken-host-db, only Kraken2 host filtering uses this custom DB.
KRAKEN_HOST_DB=""

# Kraken2/Bracken defaults for microbiome classification
KRAKEN_MICRO_CONF="0.10"
KRAKEN_MICRO_MIN_HIT_GROUPS="3"
BRACKEN_THRESHOLD="10"

# Optional cache folder containing already compressed FASTQ files.
# Used only to speed up --no-trim mode.
# If empty, --no-trim uses FASTQ files from the samplesheet directory.
CUSTOM_PATH=""

show_help() {
cat << EOF
Usage:
  bash $(basename "$0") [options]

Required:
  -i, --input FILE                         Path to samplesheet.csv
  -o, --output DIR                         Output directory
  -h, --hostid TAXID                       Host species taxon ID used for annotation/downstream host analysis

Optional:
  -m, --metadata FILE                      Metadata CSV file
  -p, --pdm METHOD                         HAllA metric: spearman, pearson, mi, nmi, xicor, dcor
                                           Default: ${pdm}
  -l, --trim-length INT                    Minimum read length required by fastp
                                           Default: ${length}
  -r, --bracken-read-len INT               Bracken read length
                                           Default: ${read_len}
      --threads INT                        Number of CPU threads
                                           Default: nproc = ${threads}

Host processing:
  -b, --blast                              Use Magic-BLAST instead of HISAT2
  -t, --no-trim                            Skip fastp trimming
      --custom-raw-path DIR                Folder used when choice="skip"
                                           Default: ${CUSTOM_PATH}

Kraken2 host filtering:
      --kraken-host-db DIR                 Optional custom Kraken2 host-filtering database.
                                           If not provided, DB_host is selected automatically from --hostid.
                                           This does NOT change GTF, BLAST/HISAT2, featureCounts, or host DEG resources.

Kraken2 microbiome classification:
      --kraken-micro-confidence FLOAT      Kraken2 --confidence for microbiome step
                                           Default: ${KRAKEN_MICRO_CONF}
      --kraken-micro-min-hit-groups INT    Kraken2 --minimum-hit-groups for microbiome step
                                           Default: ${KRAKEN_MICRO_MIN_HIT_GROUPS}

Bracken:
      --bracken-threshold INT              Bracken -t minimum read threshold
                                           Default: ${BRACKEN_THRESHOLD}

Other:
      --help                               Show this help message

Examples:
  Default host DB from --hostid:
    bash $(basename "$0") \\
      --input samplesheet.csv \\
      --output MTD_results_Myotis_auto \\
      --hostid 59463 \\
      --blast \\
      --no-trim

  Custom Kraken2 host DB, while keeping --hostid for annotation:
    bash $(basename "$0") \\
      --input samplesheet.csv \\
      --output MTD_results_Carollia_Myotis \\
      --hostid 59463 \\
      --blast \\
      --no-trim \\
      --kraken-host-db /home/me/MTD/kraken2DB_Carollia_Myotis/ \\
      --kraken-micro-confidence 0.10 \\
      --kraken-micro-min-hit-groups 3 \\
      --bracken-threshold 10
EOF
}

# Save original command line before parsing arguments
ORIGINAL_COMMAND="$(printf '%q ' "$0" "$@")"
LAUNCH_DIR="$(pwd)"
# ------------------------------------------------------------
# Parse command-line options
# ------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            inputdr="$2"
            shift 2
            ;;
        --input=*)
            inputdr="${1#*=}"
            shift
            ;;

        -o|--output)
            outputdr="$2"
            shift 2
            ;;
        --output=*)
            outputdr="${1#*=}"
            shift
            ;;

        -h|--hostid)
            hostid="$2"
            shift 2
            ;;
        --hostid=*)
            hostid="${1#*=}"
            shift
            ;;

        -m|--metadata)
            metadata="$2"
            shift 2
            ;;
        --metadata=*)
            metadata="${1#*=}"
            shift
            ;;

        -p|--pdm)
            pdm="$2"
            shift 2
            ;;
        --pdm=*)
            pdm="${1#*=}"
            shift
            ;;

        -l|--trim-length)
            length="$2"
            shift 2
            ;;
        --trim-length=*)
            length="${1#*=}"
            shift
            ;;

        -r|--bracken-read-len)
            read_len="$2"
            shift 2
            ;;
        --bracken-read-len=*)
            read_len="${1#*=}"
            shift
            ;;

        --threads)
            threads="$2"
            shift 2
            ;;
        --threads=*)
            threads="${1#*=}"
            shift
            ;;

        -b|--blast)
            blast="blast"
            shift
            ;;

        -t|--no-trim)
            no_trimm=1
            shift
            ;;

        --custom-raw-path)
            CUSTOM_PATH="$2"
            shift 2
            ;;
        --custom-raw-path=*)
            CUSTOM_PATH="${1#*=}"
            shift
            ;;

        --kraken-host-db|--host-kraken-db)
            KRAKEN_HOST_DB="$2"
            shift 2
            ;;
        --kraken-host-db=*|--host-kraken-db=*)
            KRAKEN_HOST_DB="${1#*=}"
            shift
            ;;

        --kraken-micro-confidence|--kraken-micro-conf)
            KRAKEN_MICRO_CONF="$2"
            shift 2
            ;;
        --kraken-micro-confidence=*|--kraken-micro-conf=*)
            KRAKEN_MICRO_CONF="${1#*=}"
            shift
            ;;

        --kraken-micro-min-hit-groups)
            KRAKEN_MICRO_MIN_HIT_GROUPS="$2"
            shift 2
            ;;
        --kraken-micro-min-hit-groups=*)
            KRAKEN_MICRO_MIN_HIT_GROUPS="${1#*=}"
            shift
            ;;

        --bracken-threshold)
            BRACKEN_THRESHOLD="$2"
            shift 2
            ;;
        --bracken-threshold=*)
            BRACKEN_THRESHOLD="${1#*=}"
            shift
            ;;

        --help)
            show_help
            exit 0
            ;;

        *)
            echo "${r}[ERROR] Unknown option: $1${w}" >&2
            echo
            show_help
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------
# Custom raw path behavior
# ------------------------------------------------------------
# --custom-raw-path is a cache of already compressed FASTQ files.
# It is intended only for --no-trim mode.
# If provided, force no_trimm=1 and copy .gz files directly.

if [[ -n "${CUSTOM_PATH:-}" ]]; then
    if [[ ! -d "$CUSTOM_PATH" ]]; then
        die "--custom-raw-path was provided but directory was not found: $CUSTOM_PATH"
    fi

    if [[ "$no_trimm" != "1" ]]; then
        echo "${y}[WARNING] --custom-raw-path was provided, so --no-trim mode will be enabled automatically.${w}"
        no_trimm=1
    fi
fi

# ------------------------------------------------------------
# Required argument checks
# ------------------------------------------------------------

if [[ -z "${inputdr:-}" ]]; then
    die "Missing required argument: -i or --input samplesheet.csv"
fi

if [[ -z "${outputdr:-}" ]]; then
    die "Missing required argument: -o or --output output_directory"
fi

if [[ -z "${hostid:-}" ]]; then
    die "Missing required argument: -h or --hostid TAXID"
fi

# ------------------------------------------------------------
# Basic value validation
# ------------------------------------------------------------

if ! [[ "$threads" =~ ^[0-9]+$ ]] || [[ "$threads" -lt 1 ]]; then
    die "--threads must be a positive integer. Got: $threads"
fi

if ! [[ "$length" =~ ^[0-9]+$ ]] || [[ "$length" -lt 1 ]]; then
    die "--trim-length must be a positive integer. Got: $length"
fi

if ! [[ "$read_len" =~ ^[0-9]+$ ]] || [[ "$read_len" -lt 1 ]]; then
    die "--bracken-read-len must be a positive integer. Got: $read_len"
fi

if ! [[ "$KRAKEN_MICRO_MIN_HIT_GROUPS" =~ ^[0-9]+$ ]] || [[ "$KRAKEN_MICRO_MIN_HIT_GROUPS" -lt 1 ]]; then
    die "--kraken-micro-min-hit-groups must be a positive integer. Got: $KRAKEN_MICRO_MIN_HIT_GROUPS"
fi

if ! [[ "$BRACKEN_THRESHOLD" =~ ^[0-9]+$ ]] || [[ "$BRACKEN_THRESHOLD" -lt 0 ]]; then
    die "--bracken-threshold must be an integer >= 0. Got: $BRACKEN_THRESHOLD"
fi

if ! awk -v x="$KRAKEN_MICRO_CONF" 'BEGIN { exit !(x >= 0 && x <= 1) }'; then
    die "--kraken-micro-confidence must be between 0 and 1. Got: $KRAKEN_MICRO_CONF"
fi

if [[ "$pdm" != "spearman" && "$pdm" != "pearson" && "$pdm" != "mi" && "$pdm" != "nmi" && "$pdm" != "xicor" && "$pdm" != "dcor" ]]; then
    die "--pdm must be one of: spearman, pearson, mi, nmi, xicor, dcor. Got: $pdm"
fi

# ------------------------------------------------------------
# Output directory behavior
# ------------------------------------------------------------

if [[ -d "$outputdr" ]]; then
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

# ------------------------------------------------------------
# MTD path and Conda environment
# ------------------------------------------------------------

MTDIR=$(dirname "$(readlink -f "$0")")
echo "MTD directory is $MTDIR"

condapath=$(head -n 1 "$MTDIR/condaPath")

source "$condapath/etc/profile.d/conda.sh"
conda deactivate
conda activate MTD

# ------------------------------------------------------------
# Input paths
# ------------------------------------------------------------

samplesheet_file="$inputdr"

if [[ ! -s "$samplesheet_file" ]]; then
    die "Samplesheet not found or empty: $samplesheet_file"
fi

inputdr=$(dirname "$samplesheet_file")

mkdir -p "$outputdr"
mkdir -p "$outputdr/temp"

cd "$outputdr/temp" || die "Could not enter temp directory: $outputdr/temp"

# ------------------------------------------------------------
# Step 0: Host database auto selection from --hostid
# ------------------------------------------------------------

if [[ "$hostid" == 9606 ]]; then
    DB_host="$MTDIR/kraken2DB_human"
    DB_hisat2="$MTDIR/hisat2_index_human/genome_tran"
    DB_blast="$MTDIR/human_blastdb/human_blastdb"
    gtf="$MTDIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz"

elif [[ "$hostid" == 9544 ]]; then
    DB_host="$MTDIR/kraken2DB_rhesus"
    DB_hisat2="$MTDIR/hisat2_index_rhesus/genome_tran"
    DB_blast="$MTDIR/rhesus_blastdb/rhesus_blastdb"
    gtf="$MTDIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz"

elif [[ "$hostid" == 10090 ]]; then
    DB_host="$MTDIR/kraken2DB_mice"
    DB_hisat2="$MTDIR/hisat2_index_mouse/genome_tran"
    DB_blast="$MTDIR/mouse_blastdb/mouse_blastdb"
    gtf="$MTDIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz"

elif [[ -d "$MTDIR/kraken2DB_${hostid}" ]]; then
    DB_host="$MTDIR/kraken2DB_${hostid}"
    DB_hisat2="$MTDIR/hisat2_index_${hostid}/genome_tran"
    DB_blast="$MTDIR/blastdb_${hostid}/blastdb_${hostid}"
    gtf="$MTDIR/ref_${hostid}/ref_${hostid}.gtf.gz"

else
    echo "${r}[ERROR] Host species is not supported for --hostid $hostid.${w}"
    echo "You can use bash Customized_host.sh for building the required resources."
    exit 1
fi

# ------------------------------------------------------------
# Optional Kraken2 host-filtering DB override
# ------------------------------------------------------------
# --hostid still controls downstream host annotation resources:
#   GTF, HISAT2/Magic-BLAST DB, featureCounts and host DEG.
#
# --kraken-host-db, when provided, overrides only DB_host for
# Kraken2 host read filtering.

DB_host_from_hostid="$DB_host"

if [[ -n "${KRAKEN_HOST_DB:-}" ]]; then
    DB_host="$KRAKEN_HOST_DB"
    KRAKEN_HOST_DB_MODE="custom_path_from_--kraken-host-db"
else
    DB_host="$DB_host_from_hostid"
    KRAKEN_HOST_DB_MODE="auto_from_--hostid"
fi

DB_micro="$MTDIR/kraken2DB_micro"

# ------------------------------------------------------------
# Scientific name from HostSpecies.csv
# ------------------------------------------------------------

species_name=$(awk -F, -v taxid="$hostid" '
    NR > 1 && $1 == taxid {
        print $3
        exit
    }
' "$MTDIR/HostSpecies.csv")

if [[ -z "$species_name" ]]; then
    species_name="scientific name not found"
    hostid_display="${hostid} (${y}${species_name}${w})"
else
    hostid_display="${hostid} (${ital}${species_name}${noital})"
fi

# ------------------------------------------------------------
# Opening summary
# ------------------------------------------------------------

echo "${g}Selected pipeline parameters:${w}"
echo "  input samplesheet:              $samplesheet_file"
echo "  input directory:                $inputdr"
echo "  output directory:               $outputdr"
echo "  host annotation taxid:          $hostid_display"
echo "  host Kraken2 DB:                $DB_host"
echo "  host Kraken2 DB mode:           $KRAKEN_HOST_DB_MODE"
echo "  microbiome Kraken2 DB:          $DB_micro"
echo "  metadata:                       ${metadata:-none}"
echo "  HAllA metric:                   $pdm"
echo "  trim length:                    $length"
echo "  Bracken read length:            $read_len"
echo "  threads:                        $threads"
echo "  host alignment mode:            $blast"
echo "  no_trim:                        $no_trimm"
echo "  custom raw cache path:          ${CUSTOM_PATH:-not provided}"
echo "  Kraken micro confidence:        $KRAKEN_MICRO_CONF"
echo "  Kraken micro min hit groups:    $KRAKEN_MICRO_MIN_HIT_GROUPS"
echo "  Bracken threshold:              $BRACKEN_THRESHOLD"

echo "${g}Kraken2 microbiome parameters:${w}"
echo "  --kraken-micro-confidence $KRAKEN_MICRO_CONF"
echo "  --kraken-micro-min-hit-groups $KRAKEN_MICRO_MIN_HIT_GROUPS"

echo "${g}Kraken2 host filtering parameters:${w}"
if [[ "$KRAKEN_HOST_DB_MODE" == "custom_path_from_--kraken-host-db" ]]; then
    echo "  --kraken-host-db $DB_host"
else
    echo "  --kraken-host-db not provided; using DB from --hostid"
fi

echo "${g}============================================"
echo "Selected host annotation species:${w} ${ital}${species_name}${noital}${g}"
echo "Annotation Taxon ID:${w} $hostid ${g}"
echo "Host Kraken2 filtering DB:${w} $DB_host ${g}"
echo "${g}============================================${w}"

# ------------------------------------------------------------
# Export methods and run parameters
# ------------------------------------------------------------

write_methods_log() {
    local methods_dir="$outputdr/methods"
    local methods_csv="$methods_dir/mtd_methods_run_parameters.csv"
    local bracken_dist="$DB_micro/database${read_len}mers.kmer_distrib"

    mkdir -p "$methods_dir"

    csv_escape() {
        local s="${1:-}"
        s="${s//$'\r'/}"
        s="${s//$'\n'/ }"
        s="${s//\"/\"\"}"
        printf '"%s"' "$s"
    }

    csv_row() {
        csv_escape "$1"; printf ","
        csv_escape "$2"; printf ","
        csv_escape "$3"; printf ","
        csv_escape "$4"; printf ","
        csv_escape "$5"; printf "\n"
    }

    get_tool_path() {
        command -v "$1" 2>/dev/null || printf "not_found"
    }

    get_tool_version() {
        local cmd="$1"
        eval "$cmd" 2>&1 | head -n 1 | sed 's/\r//g'
    }

    {
        csv_row "category" "program" "parameter" "value" "description"

        # Run metadata
        csv_row "Run metadata" "MTD_SE.sh" "run_datetime" "$(date -Is)" "Date and time when the run was started"
        csv_row "Run metadata" "MTD_SE.sh" "launch_directory" "$LAUNCH_DIR" "Directory from which the command was launched"
        csv_row "Run metadata" "MTD_SE.sh" "original_command" "$ORIGINAL_COMMAND" "Original command line used to start the run"
        csv_row "Run metadata" "MTD_SE.sh" "MTD_directory" "$MTDIR" "Path to the MTD installation directory"
        csv_row "Run metadata" "MTD_SE.sh" "conda_path" "$condapath" "Path to Conda installation used by the pipeline"

        # Input/output
        csv_row "Input and output" "MTD_SE.sh" "samplesheet_file" "$samplesheet_file" "Input samplesheet CSV"
        csv_row "Input and output" "MTD_SE.sh" "input_directory" "$inputdr" "Directory containing the original FASTQ files or samplesheet"
        csv_row "Input and output" "MTD_SE.sh" "output_directory" "$outputdr" "Main output directory"
        csv_row "Input and output" "MTD_SE.sh" "metadata_file" "${metadata:-none}" "Optional metadata file"

        # Study design
        csv_row "Study design" "MTD_SE.sh" "host_annotation_taxid" "$hostid" "Taxon ID used for host annotation and downstream host analyses"
        csv_row "Study design" "HostSpecies.csv" "host_annotation_species" "$species_name" "Scientific name associated with host annotation taxid"
        csv_row "Study design" "MTD_SE.sh" "HAllA_metric" "$pdm" "Association metric used by HAllA"

        # Read preparation
        csv_row "Read preparation" "fastp / pigz / cp" "no_trim" "$no_trimm" "If 1, fastp trimming is skipped"
        csv_row "Read preparation" "fastp" "trim_length" "$length" "Minimum read length required by fastp"
        csv_row "Read preparation" "MTD_SE.sh" "custom_raw_cache_path" "${CUSTOM_PATH:-not provided}" "Optional cache directory containing already compressed FASTQ files"
        csv_row "Read preparation" "MTD_SE.sh" "prepared_fastq_pattern" "Trimmed_SAMPLE.fq.gz" "FASTQ naming pattern used by downstream Kraken2 host step"

        # Host filtering
        csv_row "Host filtering" "Kraken2" "host_kraken_db" "$DB_host" "Kraken2 database used for host read filtering"
        csv_row "Host filtering" "Kraken2" "host_kraken_db_mode" "$KRAKEN_HOST_DB_MODE" "Whether DB_host came from --hostid or --kraken-host-db"
        csv_row "Host filtering" "Kraken2" "host_kraken_db_from_hostid" "$DB_host_from_hostid" "Default Kraken2 host DB selected from --hostid before optional override"
        csv_row "Host filtering" "Kraken2" "classified_out" "SAMPLE_host.fq" "Reads classified as host"
        csv_row "Host filtering" "Kraken2" "unclassified_out" "SAMPLE_non-host_raw.fq" "Reads not classified as host and used for microbiome classification"

        # Microbiome classification
        csv_row "Microbiome classification" "Kraken2" "microbiome_kraken_db" "$DB_micro" "Kraken2 database used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "confidence" "$KRAKEN_MICRO_CONF" "Kraken2 --confidence used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "minimum_hit_groups" "$KRAKEN_MICRO_MIN_HIT_GROUPS" "Kraken2 --minimum-hit-groups used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "raw_report_pattern" "Report_non-host.raw_SAMPLE.txt" "Raw microbiome Kraken2 report before optional contaminant removal"
        csv_row "Microbiome classification" "Kraken2" "final_report_pattern" "Report_non-host_SAMPLE.txt" "Final microbiome Kraken2 report used by Bracken"

        # Contaminant removal
        csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list" "$MTDIR/conta_ls.txt" "Optional list of contaminant taxids to exclude"
        if [[ -s "$MTDIR/conta_ls.txt" ]]; then
            csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list_status" "present" "conta_ls.txt was found and may be used"
        else
            csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list_status" "not_found_or_empty" "No contaminant list was found or file was empty"
        fi

        # Bracken
        csv_row "Abundance estimation" "Bracken" "bracken_read_length" "$read_len" "Read length used by Bracken -r"
        csv_row "Abundance estimation" "Bracken" "bracken_threshold" "$BRACKEN_THRESHOLD" "Minimum read threshold used by Bracken -t"
        csv_row "Abundance estimation" "Bracken" "bracken_distribution_file" "$bracken_dist" "Expected Bracken read-length-specific kmer distribution file"
        if [[ -s "$bracken_dist" ]]; then
            csv_row "Abundance estimation" "Bracken" "bracken_distribution_file_status" "present" "Bracken distribution file exists"
        else
            csv_row "Abundance estimation" "Bracken" "bracken_distribution_file_status" "missing" "Bracken distribution file not found at run-start logging time"
        fi
        csv_row "Abundance estimation" "Bracken" "taxonomic_levels" "P;G;S" "Bracken is run at phylum, genus and species levels"

        # Host downstream analysis
        csv_row "Host downstream analysis" "$blast" "host_alignment_mode" "$blast" "Host read alignment mode: magicblast if blast, otherwise HISAT2"
        csv_row "Host downstream analysis" "Magic-BLAST" "blast_database" "$DB_blast" "Magic-BLAST database selected from --hostid"
        csv_row "Host downstream analysis" "HISAT2" "hisat2_database" "$DB_hisat2" "HISAT2 database selected from --hostid"
        csv_row "Host downstream analysis" "featureCounts" "gtf_file" "$gtf" "GTF annotation selected from --hostid"
        csv_row "Host downstream analysis" "featureCounts" "output" "$outputdr/host_counts.txt" "Host count matrix output"

        # HUMAnN
        csv_row "Functional profiling" "HUMAnN" "input_pattern" "SAMPLE_non-host.fq" "Final non-host reads used as HUMAnN input"
        csv_row "Functional profiling" "HUMAnN" "threads" "$threads" "Threads used by HUMAnN"
        csv_row "Functional profiling" "HUMAnN" "renormalization_units" "relab" "HUMAnN output is renormalized to relative abundance"

        # HAllA
        csv_row "Association analysis" "HAllA" "metric" "$pdm" "Association metric used for microbiome-host associations"
        csv_row "Association analysis" "HAllA" "RUN_EXTRA_PEARSON" "${RUN_EXTRA_PEARSON:-1}" "Whether extra Pearson HAllA analysis is enabled later in the script"
        csv_row "Association analysis" "HAllA" "RUN_FULL_HALLAGRAM" "${RUN_FULL_HALLAGRAM:-0}" "Whether full hallagram plots are enabled"
        csv_row "Association analysis" "HAllA" "HALLA_DIAGNOSTIC" "${HALLA_DIAGNOSTIC:-0}" "Whether HAllA diagnostic plot is enabled"

        # Outputs
        csv_row "Key outputs" "Kraken2" "raw_global_composition" "$outputdr/kraken/kraken_global_read_composition_raw.tsv" "Global host/microbiome/unclassified composition before optional contaminant removal"
        csv_row "Key outputs" "Kraken2" "final_global_composition" "$outputdr/kraken/kraken_global_read_composition_final.tsv" "Global host/microbiome/unclassified composition after optional contaminant removal"
        csv_row "Key outputs" "Bracken" "species_table" "$outputdr/bracken_species_all" "Combined Bracken species table"
        csv_row "Key outputs" "HAllA" "microbiome_input" "$outputdr/halla/Microbiomes.txt" "Microbiome matrix used by HAllA"

        # Software paths
        csv_row "Software path" "kraken2" "path" "$(get_tool_path kraken2)" "Executable path"
        csv_row "Software path" "bracken" "path" "$(get_tool_path bracken)" "Executable path"
        csv_row "Software path" "fastp" "path" "$(get_tool_path fastp)" "Executable path"
        csv_row "Software path" "pigz" "path" "$(get_tool_path pigz)" "Executable path"
        csv_row "Software path" "python" "path" "$(get_tool_path python)" "Executable path"
        csv_row "Software path" "Rscript" "path" "$(get_tool_path Rscript)" "Executable path"
        csv_row "Software path" "humann" "path" "$(get_tool_path humann)" "Executable path"
        csv_row "Software path" "magicblast" "path" "$(get_tool_path magicblast)" "Executable path"
        csv_row "Software path" "hisat2" "path" "$(get_tool_path hisat2)" "Executable path"
        csv_row "Software path" "featureCounts" "path" "$(get_tool_path featureCounts)" "Executable path"
        csv_row "Software path" "samtools" "path" "$(get_tool_path samtools)" "Executable path"

        # Software versions, best effort
        csv_row "Software version" "kraken2" "version" "$(get_tool_version 'kraken2 --version')" "Best-effort version capture"
        csv_row "Software version" "bracken" "version" "$(get_tool_version 'bracken -v')" "Best-effort version capture"
        csv_row "Software version" "fastp" "version" "$(get_tool_version 'fastp --version')" "Best-effort version capture"
        csv_row "Software version" "python" "version" "$(get_tool_version 'python --version')" "Best-effort version capture"
        csv_row "Software version" "Rscript" "version" "$(get_tool_version 'Rscript --version')" "Best-effort version capture"
        csv_row "Software version" "humann" "version" "$(get_tool_version 'humann --version')" "Best-effort version capture"
        csv_row "Software version" "magicblast" "version" "$(get_tool_version 'magicblast -version')" "Best-effort version capture"
        csv_row "Software version" "hisat2" "version" "$(get_tool_version 'hisat2 --version')" "Best-effort version capture"
        csv_row "Software version" "featureCounts" "version" "$(get_tool_version 'featureCounts -v')" "Best-effort version capture"
        csv_row "Software version" "samtools" "version" "$(get_tool_version 'samtools --version')" "Best-effort version capture"

    } > "$methods_csv"

    echo "${g}[OK] Methods/run parameters exported to:${w}"
    echo "  $methods_csv"
}

write_methods_log

# ------------------------------------------------------------
# for SRR input samples in the samplesheet.csv; download SRR samples
# ------------------------------------------------------------

cd "$inputdr" || die "Could not enter input directory: $inputdr"

if [[ -n "$(cut -f 1 -d ',' "$samplesheet_file" | grep '^SRR' || true)" ]]; then
    for s in $(cut -f 1 -d ',' "$samplesheet_file" | grep '^SRR'); do
        if [[ ! -f "${s}_1.fastq" || ! -f "${s}_2.fastq" ]]; then
            echo "File ${s} fastq files do NOT exist. Start downloading..."
            echo "download SRA files..."
            prefetch -X 999G "$s"
            echo "split SRA files to fastq files..."
            fasterq-dump -p --split-files "$s"
            rm -rf "$s"
        fi
    done
fi

cd "$outputdr/temp" || die "Could not enter temp directory: $outputdr/temp"

# ------------------------------------------------------------
# Extract sample names from input FASTQ files
# Supports .fq.gz, .fastq.gz, .fq, .fastq
# ------------------------------------------------------------

files=$(find "$inputdr" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fastq" \))

lsn=""

for fqfile in $files; do
    fn=$(basename "$fqfile")

    base="$fn"
    base="${base%.gz}"
    base="${base%.fastq}"
    base="${base%.fq}"

    # Remove common single/paired-end suffixes if present
    sn=$(echo "$base" | sed -E 's/_R?1(_001)?$//; s/_1$//')

    lsn="$lsn $sn"
done

# Remove duplicates and normalize spaces
lsn=$(echo "$lsn" | tr " " "\n" | awk 'NF' | sort -u | tr "\n" " " | sed 's/[[:space:]]*$//')

# ------------------------------------------------------------
# Check if input files match samplesheet.csv
# ------------------------------------------------------------

fastq_files=$(echo "$lsn" | tr " " "\n" | sort | tr "\n" " ")
SamplesInSheet=$(cut -f 1 -d ',' "$samplesheet_file" | tail -n +2 | sort | tr "\n" " ")

if [[ "$fastq_files" != "$SamplesInSheet" ]]; then
    echo "${r}[ERROR] The sample FASTQ files in the input folder do not match your samplesheet.csv.${w}"
    echo
    echo "Samples detected from FASTQ files:"
    echo "$fastq_files"
    echo
    echo "Samples listed in samplesheet.csv:"
    echo "$SamplesInSheet"
    echo
    echo "Please ensure no unrelated FASTQ files are under the input folder/subfolders."
    exit 1
fi

# ------------------------------------------------------------
# Study design summary
# ------------------------------------------------------------

echo "${g}============================================"
echo "Main study design:${w}"

awk -F',' '
NR > 1 {
    groups[$2]++
}
END {
    for (g in groups) {
        printf "Group: %s - Number of samples: %d\n", g, groups[g]
    }
}' "$samplesheet_file"

echo "${g}============================================${w}"

# ------------------------------------------------------------
# Optional metadata summary
# ------------------------------------------------------------

if [[ -n "$metadata" ]]; then
    if [[ ! -s "$metadata" ]]; then
        die "Metadata file not found or empty: $metadata"
    fi

    echo "============================================"

    header=$(head -n 1 "$metadata")
    IFS=',' read -ra columns <<< "$header"

    for ((i=3; i<=${#columns[@]}; i++)); do
        col="${columns[$i-1]}"

        echo "Metadata column: $col,"
        echo "Meta-groups:"

        awk -v col_index="$i" -F',' '
        NR > 1 {
            values[$col_index]++
        }
        END {
            for (value in values) {
                printf "  %s: %d\n", value, values[value]
            }
        }' "$metadata"
    done

    echo "============================================"
    echo ""
fi

echo "${g}MTD running  progress:"
echo ">>                  [10%]"

echo "Raw reads preparation${w}"

# ------------------------------------------------------------
# Prepare FASTQ files for the pipeline
# ------------------------------------------------------------
# The host Kraken2 step expects files named:
#   Trimmed_${sample}.fq.gz
#
# Modes:
#   1) Normal mode:
#        no --no-trim and no --custom-raw-path
#        -> run fastp using FASTQ files from samplesheet directory
#
#   2) No-trim mode without custom cache:
#        --no-trim
#        -> copy .gz files or compress uncompressed FASTQ files from samplesheet directory
#
#   3) No-trim mode with custom cache:
#        --custom-raw-path DIR
#        -> copy already compressed FASTQ files from DIR
#        -> no fastp, no pigz compression
# ------------------------------------------------------------

find_fastq_for_sample() {
    local sample="$1"
    local search_dir="$2"
    local compressed_only="${3:-0}"

    if [[ "$compressed_only" == "1" ]]; then
        find "$search_dir" -maxdepth 1 -type f \( \
            -name "Trimmed_${sample}.fq.gz" -o \
            -name "Trimmed_${sample}.fastq.gz" -o \
            -name "${sample}.fq.gz" -o \
            -name "${sample}.fastq.gz" -o \
            -name "${sample}_R1*.fq.gz" -o \
            -name "${sample}_R1*.fastq.gz" -o \
            -name "${sample}_1*.fq.gz" -o \
            -name "${sample}_1*.fastq.gz" -o \
            -name "${sample}*.fq.gz" -o \
            -name "${sample}*.fastq.gz" \
        \) | sort | head -n 1
    else
        find "$search_dir" -maxdepth 1 -type f \( \
            -name "Trimmed_${sample}.fq.gz" -o \
            -name "Trimmed_${sample}.fastq.gz" -o \
            -name "Trimmed_${sample}.fq" -o \
            -name "Trimmed_${sample}.fastq" -o \
            -name "${sample}.fq.gz" -o \
            -name "${sample}.fastq.gz" -o \
            -name "${sample}.fq" -o \
            -name "${sample}.fastq" -o \
            -name "${sample}_R1*.fq.gz" -o \
            -name "${sample}_R1*.fastq.gz" -o \
            -name "${sample}_R1*.fq" -o \
            -name "${sample}_R1*.fastq" -o \
            -name "${sample}_1*.fq.gz" -o \
            -name "${sample}_1*.fastq.gz" -o \
            -name "${sample}_1*.fq" -o \
            -name "${sample}_1*.fastq" -o \
            -name "${sample}*.fq.gz" -o \
            -name "${sample}*.fastq.gz" -o \
            -name "${sample}*.fq" -o \
            -name "${sample}*.fastq" \
        \) | sort | head -n 1
    fi
}

total_cores=$(nproc)

if [[ "$total_cores" -le 4 ]]; then
    threads_per_job=1
elif [[ "$total_cores" -le 8 ]]; then
    threads_per_job=2
elif [[ "$total_cores" -le 16 ]]; then
    threads_per_job=4
else
    threads_per_job=10
fi

max_jobs=$(( total_cores / threads_per_job ))

if [[ "$max_jobs" -lt 1 ]]; then
    max_jobs=1
fi

# ------------------------------------------------------------
# Mode 1: custom cache path provided
# ------------------------------------------------------------

if [[ -n "${CUSTOM_PATH:-}" ]]; then
    echo "${y}[INFO] Using --custom-raw-path cache mode.${w}"
    echo "This mode expects already compressed FASTQ files."
    echo "No fastp and no pigz compression will be performed."
    echo "Custom FASTQ cache:"
    echo "  $CUSTOM_PATH"
    echo

    for i in $lsn; do
        fq=$(find_fastq_for_sample "$i" "$CUSTOM_PATH" 1)

        if [[ -z "$fq" ]]; then
            echo "${r}[ERROR] Could not find compressed FASTQ cache file for sample:${w} $i"
            echo "Expected .fq.gz or .fastq.gz in:"
            echo "  $CUSTOM_PATH"
            exit 1
        fi

        out_fq="$outputdr/temp/Trimmed_${i}.fq.gz"

        echo "============================================================"
        echo "[COPY CACHED FASTQ] Sample: $i"
        echo "Input:  $fq"
        echo "Output: $out_fq"
        echo "============================================================"

        cp "$fq" "$out_fq"

        if [[ ! -s "$out_fq" ]]; then
            echo "${r}[ERROR] Failed to create:${w} $out_fq"
            exit 1
        fi
    done

# ------------------------------------------------------------
# Mode 2: no-trim without custom cache
# ------------------------------------------------------------

elif [[ "$no_trimm" == "1" ]]; then
    echo "${y}[INFO] --no-trim was declared.${w}"
    echo "Using FASTQ files from samplesheet directory:"
    echo "  $inputdr"
    echo "Compressed files will be copied; uncompressed files will be compressed with pigz."
    echo

    for i in $lsn; do
        fq=$(find_fastq_for_sample "$i" "$inputdr" 0)

        if [[ -z "$fq" ]]; then
            echo "${r}[ERROR] Could not find FASTQ file for sample:${w} $i"
            echo "Search directory:"
            echo "  $inputdr"
            exit 1
        fi

        out_fq="$outputdr/temp/Trimmed_${i}.fq.gz"

        echo "============================================================"
        echo "[NO TRIM] Sample: $i"
        echo "Input:  $fq"
        echo "Output: $out_fq"
        echo "============================================================"

        if [[ "$fq" == *.gz ]]; then
            cp "$fq" "$out_fq"
        else
            pigz -p "$threads_per_job" -c "$fq" > "$out_fq"
        fi

        if [[ ! -s "$out_fq" ]]; then
            echo "${r}[ERROR] Failed to create:${w} $out_fq"
            exit 1
        fi
    done

# ------------------------------------------------------------
# Mode 3: normal fastp trimming
# ------------------------------------------------------------

else
    echo "[INFO] Running fastp trimming using FASTQ files from samplesheet directory:"
    echo "  $inputdr"
    echo

    max_fastp_cores=16

    if [[ "$threads" -gt "$max_fastp_cores" ]]; then
        fastp_threads="$max_fastp_cores"
    else
        fastp_threads="$threads"
    fi

    for i in $lsn; do
        fq=$(find_fastq_for_sample "$i" "$inputdr" 0)

        if [[ -z "$fq" ]]; then
            echo "${r}[ERROR] Could not find FASTQ file for sample:${w} $i"
            echo "Search directory:"
            echo "  $inputdr"
            exit 1
        fi

        out_fq="$outputdr/temp/Trimmed_${i}.fq.gz"

        echo "============================================================"
        echo "[FASTP] Sample: $i"
        echo "Input:  $fq"
        echo "Output: $out_fq"
        echo "Threads: $fastp_threads"
        echo "Minimum length: $length"
        echo "============================================================"

        fastp --trim_poly_x \
              --length_required "$length" \
              --thread "$fastp_threads" \
              -i "$fq" \
              -o "$out_fq"

        if [[ ! -s "$out_fq" ]]; then
            echo "${r}[ERROR] fastp failed to create:${w} $out_fq"
            exit 1
        fi
    done
fi
#$MTDIR/MTD_scripts/data_trimming.sh 

echo "${g}MTD running  progress:"
echo '>>>>                [20%]'
echo "Reads classification by kraken2; 1st step for host ${w}"
echo "Host DB: $DB_host"

if [[ ! -d "$DB_host" ]]; then
    echo "[ERROR] Host Kraken2 DB folder not found:"
    echo "$DB_host"
    exit 1
fi

if [[ ! -s "$DB_host/hash.k2d" || ! -s "$DB_host/opts.k2d" || ! -s "$DB_host/taxo.k2d" ]]; then
    echo "[ERROR] Host Kraken2 DB appears incomplete."
    echo "Expected files:"
    echo "  $DB_host/hash.k2d"
    echo "  $DB_host/opts.k2d"
    echo "  $DB_host/taxo.k2d"
    exit 1
fi

summary_file="kraken_host_summary.tsv"
echo -e "sample\thost_classified_reads\thost_classified_pct\thost_unclassified_reads\thost_unclassified_pct" > "$summary_file"

# Threshold para marcar amostras com baixa classificação como host
HOST_LOW_WARN=50

for i in $lsn; do
    echo "============================================================"
    echo "[HOST] Sample: $i"
    echo "Input: Trimmed_${i}.fq.gz"
    echo "============================================================"

    if [[ ! -s "Trimmed_${i}.fq.gz" ]]; then
        echo "[ERROR] Missing input file:"
        echo "Trimmed_${i}.fq.gz"
        exit 1
    fi

    kraken2 --db "$DB_host" --use-names \
        --report "Report_host_${i}.txt" \
        --output "Report_host_${i}.kraken" \
        --threads "$threads" \
        --gzip-compressed \
        --classified-out "${i}_host.fq" \
        --unclassified-out "${i}_non-host_raw.fq" \
        "Trimmed_${i}.fq.gz"

    report="Report_host_${i}.txt"

    host_unclassified_pct=$(awk '$4=="U"{print $1; exit}' "$report")
    host_unclassified_reads=$(awk '$4=="U"{print $2; exit}' "$report")

    host_classified_pct=$(awk '$4=="R" && $5==1{print $1; exit}' "$report")
    host_classified_reads=$(awk '$4=="R" && $5==1{print $2; exit}' "$report")

    # Fallback caso a linha root não apareça por algum motivo
    if [[ -z "$host_classified_pct" ]]; then
        host_classified_pct=$(awk -v u="$host_unclassified_pct" 'BEGIN{printf "%.2f", 100-u}')
    fi

    if [[ -z "$host_classified_reads" ]]; then
        host_classified_reads="NA"
    fi

    echo
    echo "[RESULT] Sample: $i"
    echo "  Classified as host:   ${host_classified_pct}%  (${host_classified_reads} reads)"
    echo "  Unclassified:         ${host_unclassified_pct}%  (${host_unclassified_reads} reads)"

    echo -e "${i}\t${host_classified_reads}\t${host_classified_pct}\t${host_unclassified_reads}\t${host_unclassified_pct}" >> "$summary_file"

    echo
    echo "[HOST] Main taxa with >=1%:"
    awk '
        $4!="U" && !($4=="R" && $5==1) && $1 >= 1 {
            name=$6
            for (j=7; j<=NF; j++) name=name" "$j
            printf "    %7s%%  %12s reads  rank=%-4s taxid=%-10s %s\n", $1, $2, $4, $5, name
        }
    ' "$report" | head -n 20

    # Warning para amostras com pouca classificação como host
    if awk -v p="$host_classified_pct" -v t="$HOST_LOW_WARN" 'BEGIN{exit !(p < t)}'; then
        echo
        echo "  [WARNING] Low host classification for sample $i"
        echo "  Host classified: ${host_classified_pct}%"
        echo "  This sample may contain more microbial reads, contamination, or lower host RNA content."
    fi

    echo
done

echo "============================================================"
echo "[OK] Host Kraken2 summary saved to:"
echo "$summary_file"
echo "============================================================"
column -t "$summary_file"

echo "${g}MTD running  progress:"
echo '>>>>>               [25%]'

# ------------------------------------------------------------
# Reads classification by Kraken2; 2nd step for non-host reads
# Microbiome classification using reads not classified as host
# ------------------------------------------------------------

echo "Reads classification by kraken2; 2nd step for non-host reads ${w}"
echo "Microbiome DB: $DB_micro"

if [[ ! -d "$DB_micro" ]]; then
    echo "[ERROR] Microbiome Kraken2 DB folder not found:"
    echo "$DB_micro"
    exit 1
fi

if [[ ! -s "$DB_micro/hash.k2d" || ! -s "$DB_micro/opts.k2d" || ! -s "$DB_micro/taxo.k2d" ]]; then
    echo "[ERROR] Microbiome Kraken2 DB appears incomplete."
    echo "Expected files:"
    echo "  $DB_micro/hash.k2d"
    echo "  $DB_micro/opts.k2d"
    echo "  $DB_micro/taxo.k2d"
    exit 1
fi

micro_summary="kraken_nonhost_raw_summary.tsv"

echo -e "sample\tmicro_classified_reads\tmicro_classified_pct\tmicro_unclassified_reads\tmicro_unclassified_pct" > "$micro_summary"

# Threshold para alertar amostras com muita classificação microbiana dentro do non-host
MICRO_HIGH_WARN=20

for i in $lsn; do
    echo "============================================================"
    echo "[MICRO RAW] Sample: $i"
    echo "Input: ${i}_non-host_raw.fq"
    echo "============================================================"

    if [[ ! -s "${i}_non-host_raw.fq" ]]; then
        echo "[ERROR] Missing non-host file from host-filtering step:"
        echo "${i}_non-host_raw.fq"
        exit 1
    fi

    kraken2 --db "$DB_micro" --use-names \
    --confidence "$KRAKEN_MICRO_CONF" \
    --minimum-hit-groups "$KRAKEN_MICRO_MIN_HIT_GROUPS" \
    --report "Report_non-host.raw_${i}.txt" \
    --output "Report_non-host_raw_${i}.kraken" \
    --threads "$threads" \
    --classified-out "${i}_raw_cseqs.fq" \
    --unclassified-out "${i}_raw_ucseqs.fq" \
    "${i}_non-host_raw.fq"

    report="Report_non-host.raw_${i}.txt"

    micro_unclassified_pct=$(awk '$4=="U"{print $1; exit}' "$report")
    micro_unclassified_reads=$(awk '$4=="U"{print $2; exit}' "$report")

    micro_classified_pct=$(awk '$4=="R" && $5==1{print $1; exit}' "$report")
    micro_classified_reads=$(awk '$4=="R" && $5==1{print $2; exit}' "$report")

    if [[ -z "$micro_unclassified_pct" ]]; then
        micro_unclassified_pct="NA"
    fi

    if [[ -z "$micro_unclassified_reads" ]]; then
        micro_unclassified_reads="NA"
    fi

    if [[ -z "$micro_classified_pct" ]]; then
        if [[ "$micro_unclassified_pct" != "NA" ]]; then
            micro_classified_pct=$(awk -v u="$micro_unclassified_pct" 'BEGIN{printf "%.2f", 100-u}')
        else
            micro_classified_pct="NA"
        fi
    fi

    if [[ -z "$micro_classified_reads" ]]; then
        micro_classified_reads="NA"
    fi

    echo
    echo "[RESULT] Sample: $i"
    echo "  Classified in DB_micro:   ${micro_classified_pct}%  (${micro_classified_reads} reads)"
    echo "  Unclassified in DB_micro: ${micro_unclassified_pct}%  (${micro_unclassified_reads} reads)"

    echo -e "${i}\t${micro_classified_reads}\t${micro_classified_pct}\t${micro_unclassified_reads}\t${micro_unclassified_pct}" >> "$micro_summary"

    if [[ "$micro_classified_pct" != "NA" ]]; then
        if awk -v p="$micro_classified_pct" -v t="$MICRO_HIGH_WARN" 'BEGIN{exit !(p >= t)}'; then
            echo
            echo "  [WARNING] High DB_micro classification for sample $i"
            echo "  Microbial classification here is ${micro_classified_pct}% of the NON-HOST reads, not of total reads."
            echo
            echo "  Top taxa with >=1% in report:"
            awk '
                $4!="U" && !($4=="R" && $5==1) && $1 >= 1 {
                    name=$6
                    for (j=7; j<=NF; j++) name=name" "$j
                    printf "    %7s%%  %12s reads  rank=%-4s taxid=%-10s %s\n", $1, $2, $4, $5, name
                }
            ' "$report" | head -n 20
        fi
    fi

    echo
done

echo "============================================================"
echo "[OK] Non-host raw Kraken2 summary saved to:"
echo "$micro_summary"
echo "============================================================"
column -s $'\t' -t "$micro_summary"

echo "${g}MTD running  progress:"
echo '>>>>>>              [30%]'

# ------------------------------------------------------------
# Global read composition summary
# Percentages are calculated relative to the original total reads
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Creating global Kraken read composition table"
echo "============================================================"

host_summary="kraken_host_summary.tsv"
micro_summary="kraken_nonhost_raw_summary.tsv"
out_summary="kraken_global_read_composition_raw.tsv"

if [[ ! -s "$host_summary" ]]; then
    echo "[ERROR] Missing file: $host_summary"
    exit 1
fi

if [[ ! -s "$micro_summary" ]]; then
    echo "[ERROR] Missing file: $micro_summary"
    exit 1
fi

echo -e "sample\ttotal_reads\thost\tmicrobiome\tunclassified\tcheck_pct_sum" > "$out_summary"

awk '
BEGIN {
    FS=OFS="\t"
}

NR==FNR {
    if (FNR == 1) next

    sample=$1

    host_reads[sample]=$2
    host_unclassified_reads[sample]=$4

    total_reads[sample]=$2 + $4

    next
}

FNR > 1 {
    sample=$1

    micro_reads=$2
    micro_unclassified_reads=$4

    if (!(sample in total_reads)) {
        print "[WARNING] Sample found in micro summary but not in host summary: " sample > "/dev/stderr"
        next
    }

    total=total_reads[sample]
    host=host_reads[sample]
    micro=micro_reads
    unclassified=micro_unclassified_reads

    host_pct=(host/total)*100
    micro_pct=(micro/total)*100
    unclassified_pct=(unclassified/total)*100

    check_sum=host_pct + micro_pct + unclassified_pct

    host_label=sprintf("%d (%.2f%%)", host, host_pct)
    micro_label=sprintf("%d (%.2f%%)", micro, micro_pct)
    unclassified_label=sprintf("%d (%.2f%%)", unclassified, unclassified_pct)

    printf "%s\t%d\t%s\t%s\t%s\t%.2f\n", \
        sample, total, host_label, micro_label, unclassified_label, check_sum
}
' "$host_summary" "$micro_summary" >> "$out_summary"

echo
echo "[OK] Global read composition saved to:"
echo "$out_summary"
echo

column -s $'\t' -t "$out_summary"

echo
echo "[INFO] Interpretation:"
echo "  total_reads    = original reads after trimming/compression step"
echo "  host           = reads classified as host in Kraken2 host step"
echo "  microbiome     = reads classified by DB_micro after host removal"
echo "  unclassified   = reads not classified as host and not classified by DB_micro"
echo "  check_pct_sum  = should be close to 100.00"
echo "============================================================"

echo "${g}MTD running  progress:"
echo '>>>>>>              [30%]'
mkdir -p "$outputdr/kraken"
mv kraken_global_read_composition_raw.tsv kraken_host_summary.tsv kraken_nonhost_raw_summary.tsv "$outputdr/kraken/"

echo "Decontamination step${w}"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate MTD

export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset PYTHONHOME

# ------------------------------------------------------------
# Optional contaminant removal before final microbiome reports
# ------------------------------------------------------------
# If conta_ls.txt exists and contains taxids, remove those taxa from
# the raw non-host reads and re-run Kraken2 to create final reports.
#
# If conta_ls.txt does not exist, or is empty, use the raw non-host
# Kraken2 outputs as the final non-host outputs.
#
# This guarantees that Bracken always receives:
#   Report_non-host_${sample}.txt
# ------------------------------------------------------------

conta_file="$MTDIR/conta_ls.txt"

# Safety defaults, caso você ainda não tenha definido esses parâmetros antes
KRAKEN_MICRO_CONF="${KRAKEN_MICRO_CONF:-0.10}"
KRAKEN_MICRO_MIN_HIT_GROUPS="${KRAKEN_MICRO_MIN_HIT_GROUPS:-3}"

if [[ -s "$conta_file" ]]; then
    echo "[INFO] Contaminant list found:"
    echo "  $conta_file"

    # Assumes taxids are in the second column, tab-separated.
    conta_ls=$(awk -F '\t' '{gsub(/\r/, "", $2); if ($2 ~ /^[0-9]+$/) print $2}' "$conta_file" | sort -u | paste -sd ' ' -)

    if [[ -n "$conta_ls" ]]; then
        echo "[INFO] TaxIDs to exclude:"
        echo "  $conta_ls"

        for i in $lsn; do
            echo "============================================================"
            echo "[DECONTAMINATION] Sample: $i"
            echo "Input reads: ${i}_non-host_raw.fq"
            echo "Input Kraken output: Report_non-host_raw_${i}.kraken"
            echo "Input Kraken report: Report_non-host.raw_${i}.txt"
            echo "Output reads: ${i}_non-host.fq"
            echo "============================================================"

            if [[ ! -s "Report_non-host_raw_${i}.kraken" ]]; then
                echo "${r}[ERROR] Missing raw Kraken output:${w}"
                echo "  Report_non-host_raw_${i}.kraken"
                exit 1
            fi

            if [[ ! -s "Report_non-host.raw_${i}.txt" ]]; then
                echo "${r}[ERROR] Missing raw Kraken report:${w}"
                echo "  Report_non-host.raw_${i}.txt"
                exit 1
            fi

            if [[ ! -s "${i}_non-host_raw.fq" ]]; then
                echo "${r}[ERROR] Missing raw non-host reads:${w}"
                echo "  ${i}_non-host_raw.fq"
                exit 1
            fi

            python "$MTDIR/Tools/KrakenTools/extract_kraken_reads.py" \
                -k "Report_non-host_raw_${i}.kraken" \
                -s1 "${i}_non-host_raw.fq" \
                -o "${i}_non-host.fq" \
                -r "Report_non-host.raw_${i}.txt" \
                --taxid $conta_ls \
                --exclude \
                --include-children

            if [[ ! -s "${i}_non-host.fq" ]]; then
                echo "${r}[ERROR] Decontaminated non-host file was not created or is empty:${w}"
                echo "  ${i}_non-host.fq"
                exit 1
            fi
        done

        echo "${g}MTD running  progress:"
        echo '>>>>>>>             [35%]'

        echo "Reads classification by kraken2; 3rd step for decontaminated non-host reads to get final reports ${w}"

        for i in $lsn; do
            echo "============================================================"
            echo "[MICRO FINAL] Sample: $i"
            echo "Input: ${i}_non-host.fq"
            echo "Kraken2 confidence: $KRAKEN_MICRO_CONF"
            echo "Kraken2 minimum hit groups: $KRAKEN_MICRO_MIN_HIT_GROUPS"
            echo "============================================================"

            kraken2 --db "$DB_micro" --use-names \
                --confidence "$KRAKEN_MICRO_CONF" \
                --minimum-hit-groups "$KRAKEN_MICRO_MIN_HIT_GROUPS" \
                --report "Report_non-host_${i}.txt" \
                --output "Report_non-host_${i}.kraken" \
                --threads "$threads" \
                --classified-out "${i}_cseqs.fq" \
                --unclassified-out "${i}_ucseqs.fq" \
                "${i}_non-host.fq"

            if [[ ! -s "Report_non-host_${i}.txt" ]]; then
                echo "${r}[ERROR] Final Kraken2 report was not created:${w}"
                echo "  Report_non-host_${i}.txt"
                exit 1
            fi
        done

    else
        echo "${y}[WARNING] conta_ls.txt exists but no valid taxids were found in column 2.${w}"
        echo "[INFO] Using raw non-host Kraken2 outputs as final outputs."

        for i in $lsn; do
            cp "Report_non-host.raw_${i}.txt" "Report_non-host_${i}.txt"
            cp "Report_non-host_raw_${i}.kraken" "Report_non-host_${i}.kraken"
            cp "${i}_non-host_raw.fq" "${i}_non-host.fq"
            if [[ -s "${i}_raw_cseqs.fq" ]]; then
            cp "${i}_raw_cseqs.fq" "${i}_cseqs.fq"
            fi

            if [[ -s "${i}_raw_ucseqs.fq" ]]; then
                cp "${i}_raw_ucseqs.fq" "${i}_ucseqs.fq"
            fi
        done
    fi

else
    echo "[INFO] No contaminant list found:"
    echo "  $conta_file"
    echo "[INFO] Skipping contaminant removal."
    echo "[INFO] Using raw non-host Kraken2 outputs as final outputs."

    for i in $lsn; do
        echo "============================================================"
        echo "[NO DECONTAMINATION] Sample: $i"
        echo "Using raw non-host outputs as final Bracken inputs"
        echo "============================================================"

        if [[ ! -s "Report_non-host.raw_${i}.txt" ]]; then
            echo "${r}[ERROR] Missing raw Kraken2 report:${w}"
            echo "  Report_non-host.raw_${i}.txt"
            exit 1
        fi

        if [[ ! -s "Report_non-host_raw_${i}.kraken" ]]; then
            echo "${r}[ERROR] Missing raw Kraken2 output:${w}"
            echo "  Report_non-host_raw_${i}.kraken"
            exit 1
        fi

        if [[ ! -s "${i}_non-host_raw.fq" ]]; then
            echo "${r}[ERROR] Missing raw non-host reads:${w}"
            echo "  ${i}_non-host_raw.fq"
            exit 1
        fi

        cp "Report_non-host.raw_${i}.txt" "Report_non-host_${i}.txt"
        cp "Report_non-host_raw_${i}.kraken" "Report_non-host_${i}.kraken"
        cp "${i}_non-host_raw.fq" "${i}_non-host.fq"

        if [[ -s "${i}_raw_cseqs.fq" ]]; then
            cp "${i}_raw_cseqs.fq" "${i}_cseqs.fq"
        fi

        if [[ -s "${i}_raw_ucseqs.fq" ]]; then
            cp "${i}_raw_ucseqs.fq" "${i}_ucseqs.fq"
        fi
    done
fi
# ------------------------------------------------------------
# Final Kraken read composition summary
# Uses final non-host reports after optional contaminant removal
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Creating FINAL Kraken read composition table"
echo "============================================================"

final_micro_summary="$outputdr/kraken/kraken_nonhost_final_summary.tsv"
final_out_summary="$outputdr/kraken/kraken_global_read_composition_final.tsv"
host_summary="$outputdr/kraken/kraken_host_summary.tsv"

if [[ ! -s "$host_summary" ]]; then
    echo "${r}[ERROR] Missing host summary:${w}"
    echo "  $host_summary"
    exit 1
fi

echo -e "sample\tmicro_classified_reads\tmicro_classified_pct\tmicro_unclassified_reads\tmicro_unclassified_pct" > "$final_micro_summary"

for i in $lsn; do
    report="Report_non-host_${i}.txt"

    echo "============================================================"
    echo "[FINAL MICRO SUMMARY] Sample: $i"
    echo "Final Kraken2 report: $report"
    echo "============================================================"

    if [[ ! -s "$report" ]]; then
        echo "${r}[ERROR] Missing final Kraken2 report:${w}"
        echo "  $report"
        exit 1
    fi

    micro_unclassified_pct=$(awk '$4=="U"{print $1; exit}' "$report")
    micro_unclassified_reads=$(awk '$4=="U"{print $2; exit}' "$report")

    micro_classified_pct=$(awk '$4=="R" && $5==1{print $1; exit}' "$report")
    micro_classified_reads=$(awk '$4=="R" && $5==1{print $2; exit}' "$report")

    if [[ -z "$micro_unclassified_pct" ]]; then
        micro_unclassified_pct="NA"
    fi

    if [[ -z "$micro_unclassified_reads" ]]; then
        micro_unclassified_reads="NA"
    fi

    if [[ -z "$micro_classified_pct" ]]; then
        if [[ "$micro_unclassified_pct" != "NA" ]]; then
            micro_classified_pct=$(awk -v u="$micro_unclassified_pct" 'BEGIN{printf "%.2f", 100-u}')
        else
            micro_classified_pct="NA"
        fi
    fi

    if [[ -z "$micro_classified_reads" ]]; then
        micro_classified_reads="NA"
    fi

    echo "[RESULT] Final DB_micro classified:   ${micro_classified_pct}%  (${micro_classified_reads} reads)"
    echo "[RESULT] Final DB_micro unclassified: ${micro_unclassified_pct}%  (${micro_unclassified_reads} reads)"

    echo -e "${i}\t${micro_classified_reads}\t${micro_classified_pct}\t${micro_unclassified_reads}\t${micro_unclassified_pct}" >> "$final_micro_summary"
done

echo
echo "[OK] Final non-host Kraken2 summary saved to:"
echo "$final_micro_summary"
echo
column -s $'\t' -t "$final_micro_summary"

echo -e "sample\ttotal_reads\thost\tmicrobiome\tunclassified\tcheck_pct_sum" > "$final_out_summary"

awk '
BEGIN {
    FS=OFS="\t"
}

NR==FNR {
    if (FNR == 1) next

    sample=$1

    host_reads[sample]=$2
    host_unclassified_reads[sample]=$4
    total_reads[sample]=$2 + $4

    next
}

FNR > 1 {
    sample=$1

    micro_reads=$2
    micro_unclassified_reads=$4

    if (!(sample in total_reads)) {
        print "[WARNING] Sample found in final micro summary but not in host summary: " sample > "/dev/stderr"
        next
    }

    total=total_reads[sample]
    host=host_reads[sample]
    micro=micro_reads
    unclassified=micro_unclassified_reads

    host_pct=(host/total)*100
    micro_pct=(micro/total)*100
    unclassified_pct=(unclassified/total)*100

    check_sum=host_pct + micro_pct + unclassified_pct

    host_label=sprintf("%d (%.2f%%)", host, host_pct)
    micro_label=sprintf("%d (%.2f%%)", micro, micro_pct)
    unclassified_label=sprintf("%d (%.2f%%)", unclassified, unclassified_pct)

    printf "%s\t%d\t%s\t%s\t%s\t%.2f\n", \
        sample, total, host_label, micro_label, unclassified_label, check_sum
}
' "$host_summary" "$final_micro_summary" >> "$final_out_summary"

echo
echo "[OK] FINAL global read composition saved to:"
echo "$final_out_summary"
echo
column -s $'\t' -t "$final_out_summary"

echo
echo "[INFO] FINAL interpretation:"
echo "  total_reads    = original reads after trimming/compression step"
echo "  host           = reads classified as host in Kraken2 host step"
echo "  microbiome     = reads classified by DB_micro after optional contaminant removal"
echo "  unclassified   = reads not classified as host and not classified by final DB_micro step"
echo "  check_pct_sum  = should be close to 100.00"
echo "============================================================"

# ------------------------------------------------------------
# Save individual Kraken2 reports and outputs
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Saving individual Kraken2 reports to kraken folder"
echo "============================================================"

mkdir -p "$outputdr/kraken/reports_host"
mkdir -p "$outputdr/kraken/reports_micro_raw"
mkdir -p "$outputdr/kraken/reports_micro_final"

# Host reports
cp Report_host_*.txt Report_host_*.kraken "$outputdr/kraken/reports_host/" 2>/dev/null || true

# Raw microbiome reports
cp Report_non-host.raw_*.txt Report_non-host_raw_*.kraken "$outputdr/kraken/reports_micro_raw/" 2>/dev/null || true

# Final microbiome reports used by Bracken
cp Report_non-host_*.txt Report_non-host_*.kraken "$outputdr/kraken/reports_micro_final/" 2>/dev/null || true

echo "[OK] Kraken2 individual reports copied to:"
echo "  $outputdr/kraken/reports_host"
echo "  $outputdr/kraken/reports_micro_raw"
echo "  $outputdr/kraken/reports_micro_final"
echo "============================================================"

echo "${g}MTD running  progress:"
echo '>>>>>>>>            [40%]'

echo "Bracken analysis ${w}"

# ------------------------------------------------------------
# Check Bracken read-length distribution file
# ------------------------------------------------------------
# Bracken needs a read-length-specific distribution file.
# Example: if read_len=75, the database must contain:
#   database75mers.kmer_distrib
#
# This prevents running Bracken with a read length that was not built
# for the current Kraken2 microbiome database.

BRACKEN_DIST="$DB_micro/database${read_len}mers.kmer_distrib"

echo "Bracken DB: $DB_micro"
echo "Bracken read length: $read_len"
echo "Expected Bracken distribution file: $BRACKEN_DIST"

if [[ ! -s "$BRACKEN_DIST" ]]; then
    echo "${r}[ERROR] Bracken distribution file not found for read length ${read_len}.${w}"
    echo
    echo "Expected file:"
    echo "  $BRACKEN_DIST"
    echo
    echo "Available Bracken distribution files in DB_micro:"
    ls -lh "$DB_micro"/database*mers.kmer_distrib 2>/dev/null || echo "  None found."
    echo
    echo "You probably need to build the Bracken distribution file for read length ${read_len}."
    echo "Example command:"
    echo "  bracken-build -d \"$DB_micro\" -t \"$threads\" -k 35 -l \"$read_len\""
    echo
    exit 1
fi

echo "${g}[OK] Bracken distribution file found:${w} $BRACKEN_DIST"

# Bracken -t is the minimum read threshold, not threads.
# Keep it fixed so results do not change when CPU thread number changes.

echo "Bracken threshold: $BRACKEN_THRESHOLD"

for i in $lsn; do
    echo "============================================================"
    echo "[BRACKEN] Sample: $i"
    echo "Input report: Report_non-host_${i}.txt"
    echo "Read length: $read_len"
    echo "Threshold: $BRACKEN_THRESHOLD"
    echo "============================================================"

    if [[ ! -s "Report_non-host_${i}.txt" ]]; then
        echo "${r}[ERROR] Missing Kraken2 report for Bracken:${w}"
        echo "  Report_non-host_${i}.txt"
        exit 1
    fi

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.phylum.bracken" \
        -r "$read_len" \
        -l P \
        -t "$BRACKEN_THRESHOLD"

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.genus.bracken" \
        -r "$read_len" \
        -l G \
        -t "$BRACKEN_THRESHOLD"

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.species.bracken" \
        -r "$read_len" \
        -l S \
        -t "$BRACKEN_THRESHOLD"
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
sed -i 's/sp\. //g' "$outputdr/bracken_species_all.biom"

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
    --external_annotations 7 \
    --max_clade_size 300

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
        -num_threads $threads #$threads
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
#read -p "PRESS ENTER"
echo "${w}"
echo $MTDIR
echo $outputdr
echo $inputdr
echo $hostid 
echo $metadata
#### BEGIN FUNCTION: prepare_gene_id_cache_from_gtf ####
prepare_gene_id_cache_from_gtf() {
    local outputdr="$1"
    local gtf_file="$2"

    local cache_dir="${outputdr}/Host_DEG"
    local cache_file="${cache_dir}/gene_ID_cache.csv"

    local rscript_bin
    rscript_bin="$(command -v Rscript)"

    echo "------------------------------------------------------------"
    echo "Preparing offline host gene annotation cache from GTF"
    echo "Rscript: $rscript_bin"
    echo "GTF: $gtf_file"
    echo "Cache: $cache_file"
    echo "------------------------------------------------------------"

    mkdir -p "$cache_dir"

    if [[ -s "$cache_file" ]]; then
        echo "[INFO] Existing gene_ID_cache.csv found. Validating..."

        "$rscript_bin" - <<RS
cache_file <- "$cache_file"

gene_ID <- read.csv(cache_file, header = TRUE, check.names = FALSE)

required <- c(
  "gene_name",
  "ensembl_gene_id",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand",
  "gene_biotype",
  "description",
  "gene_length"
)

missing <- setdiff(required, names(gene_ID))

if (length(missing) > 0) {
  stop("Cache is missing columns: ", paste(missing, collapse = ", "))
}

if (nrow(gene_ID) == 0) {
  stop("Cache has zero rows.")
}

cat("[OK] Existing cache is valid. Genes:", nrow(gene_ID), "\\n")
RS

        if [[ "$?" -eq 0 ]]; then
            return 0
        else
            echo "[WARN] Existing cache is invalid. Rebuilding from GTF..."
            rm -f "$cache_file"
        fi
    fi

    if [[ ! -s "$gtf_file" ]]; then
        echo "[ERROR] GTF file not found or empty: $gtf_file"
        return 1
    fi

    if [[ -z "$rscript_bin" || ! -x "$rscript_bin" ]]; then
        echo "[ERROR] Could not find Rscript in PATH."
        return 1
    fi

    "$rscript_bin" - <<RS
gtf_file <- "$gtf_file"
cache_file <- "$cache_file"

cat("[INFO] Reading GTF:", gtf_file, "\\n")

open_gtf <- function(path) {
  if (grepl("\\\\.gz$", path)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- sub(pattern, "\\\\1", regmatches(x, hit)[ok], perl = TRUE)
  out
}

con <- open_gtf(gtf_file)
gtf <- read.delim(
  con,
  header = FALSE,
  sep = "\\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)
close(con)

if (ncol(gtf) < 9) {
  stop("GTF has fewer than 9 columns. Is this a valid GTF?")
}

names(gtf)[1:9] <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

genes <- gtf[gtf\$feature == "gene", , drop = FALSE]

if (nrow(genes) == 0) {
  cat("[WARN] No 'gene' features found. Falling back to exon-derived gene ranges.\\n")

  exons <- gtf[gtf\$feature == "exon", , drop = FALSE]

  if (nrow(exons) == 0) {
    stop("No gene or exon features found in GTF.")
  }

  exons\$gene_id <- extract_attr(exons\$attribute, "gene_id")
  exons\$gene_name <- extract_attr(exons\$attribute, "gene_name")
  exons\$gene_biotype <- extract_attr(exons\$attribute, "gene_biotype")

  if (all(is.na(exons\$gene_biotype))) {
    exons\$gene_biotype <- extract_attr(exons\$attribute, "gene_type")
  }

  exons <- exons[!is.na(exons\$gene_id), , drop = FALSE]

  genes <- aggregate(
    cbind(start, end) ~ gene_id + seqname + strand,
    data = exons,
    FUN = function(z) c(min = min(z), max = max(z))
  )

  genes\$start_position <- genes\$start[, "min"]
  genes\$end_position <- genes\$end[, "max"]
  genes\$start <- NULL
  genes\$end <- NULL

  meta <- exons[!duplicated(exons\$gene_id), c("gene_id", "gene_name", "gene_biotype"), drop = FALSE]
  genes <- merge(genes, meta, by = "gene_id", all.x = TRUE)

} else {
  genes\$gene_id <- extract_attr(genes\$attribute, "gene_id")
  genes\$gene_name <- extract_attr(genes\$attribute, "gene_name")
  genes\$gene_biotype <- extract_attr(genes\$attribute, "gene_biotype")

  if (all(is.na(genes\$gene_biotype))) {
    genes\$gene_biotype <- extract_attr(genes\$attribute, "gene_type")
  }

  genes\$start_position <- genes\$start
  genes\$end_position <- genes\$end
}

genes <- genes[!is.na(genes\$gene_id), , drop = FALSE]
genes <- genes[!duplicated(genes\$gene_id), , drop = FALSE]

genes\$gene_name[is.na(genes\$gene_name) | genes\$gene_name == ""] <- genes\$gene_id[is.na(genes\$gene_name) | genes\$gene_name == ""]
genes\$gene_biotype[is.na(genes\$gene_biotype) | genes\$gene_biotype == ""] <- "unknown"

gene_ID <- data.frame(
  gene_name = genes\$gene_name,
  ensembl_gene_id = genes\$gene_id,
  chromosome_name = genes\$seqname,
  start_position = genes\$start_position,
  end_position = genes\$end_position,
  strand = genes\$strand,
  gene_biotype = genes\$gene_biotype,
  description = genes\$gene_name,
  gene_length = abs(as.numeric(genes\$end_position) - as.numeric(genes\$start_position)) + 1,
  stringsAsFactors = FALSE
)

gene_ID <- gene_ID[!is.na(gene_ID\$ensembl_gene_id), , drop = FALSE]
gene_ID <- gene_ID[!duplicated(gene_ID\$ensembl_gene_id), , drop = FALSE]

write.csv(gene_ID, cache_file, row.names = FALSE, quote = TRUE)

cat("[OK] Offline cache created:", cache_file, "\\n")
cat("[OK] Genes:", nrow(gene_ID), "\\n")
cat("[INFO] First rows:\\n")
print(utils::head(gene_ID, 3))
RS

    local status=$?

    if [[ "$status" -ne 0 || ! -s "$cache_file" ]]; then
        echo "[ERROR] Failed to create gene_ID_cache.csv from GTF."
        return 1
    fi

    echo "[OK] gene_ID_cache.csv is ready."
    return 0
}
#### END FUNCTION: prepare_gene_id_cache_from_gtf ####

#### BEGIN CALL: host gene annotation cache update ####

prepare_gene_id_cache_from_gtf "$outputdr" "$gtf" || {
    echo "${r}[ERROR] Could not prepare gene_ID_cache.csv from local GTF.${w}"
    exit 1
}

# Opcional: só tenta BioMart depois que já existe um cache local seguro.
# Se a internet/DNS falhar, o pipeline continua usando o cache do GTF.
if declare -F update_host_gene_cache_online >/dev/null; then
    update_host_gene_cache_online "$outputdr" "$hostid" "$MTDIR" || true
else
    echo "[INFO] update_host_gene_cache_online function not defined; skipping online cache update."
fi

#### END CALL: host gene annotation cache update ####

Rscript "$MTDIR/DEG_Anno_Plot.R" \
    "$outputdr/host_counts.txt" \
    "$inputdr/samplesheet.csv" \
    "$hostid" \
    "$MTDIR/HostSpecies.csv" \
    $metadata

require_file "$outputdr/Host_DEG/host_counts_TPM.csv" "Host TPM matrix generated by DEG_Anno_Plot.R"

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>     [75%]'

echo "ssGSEA${w}"

require_file "$outputdr/Host_DEG/host_counts_TPM.csv" "Host TPM matrix"

run_cmd Rscript "$MTDIR/gct_making.R" \
    "$outputdr/Host_DEG/host_counts_TPM.csv" \
    "$inputdr/samplesheet.csv"

require_file "$outputdr/ssGSEA/host.gct" "ssGSEA input GCT"

run_cmd Rscript "$MTDIR/Tools/ssGSEA2.0/ssgsea-cli.R" \
    -i "$outputdr/ssGSEA/host.gct" \
    -o "$outputdr/ssGSEA/ssgsea-results" \
    -d "$MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt" \
    -y "$MTDIR/Tools/ssGSEA2.0/config.yaml" \
    -u "$threads"

require_file "$outputdr/ssGSEA/ssgsea-results-scores.gct" "ssGSEA scores"

run_cmd Rscript "$MTDIR/for_halla.R" \
    "$outputdr/ssGSEA/ssgsea-results-scores.gct" \
    "$inputdr/samplesheet.csv" \
    $metadata

require_file "$outputdr/halla/Microbiomes.txt" "HAllA microbiome input"
require_file "$outputdr/halla/Host_gene.txt" "HAllA host gene input"
require_file "$outputdr/halla/Host_score.txt" "HAllA host pathway input"

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>>    [80%]'
echo "MTD DEG analyses are done. Starting microbiome x host association analyses..."

echo "halla: association analysis${w}"

conda deactivate
conda activate halla0820

export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset PYTHONHOME
export MPLBACKEND=Agg
export PYTHONWARNINGS="ignore"

HALLA_THREADS="${threads:-$(nproc)}"
export OMP_NUM_THREADS="$HALLA_THREADS"
export OPENBLAS_NUM_THREADS="$HALLA_THREADS"
export MKL_NUM_THREADS="$HALLA_THREADS"
export NUMEXPR_NUM_THREADS="$HALLA_THREADS"

RUN_EXTRA_PEARSON="${RUN_EXTRA_PEARSON:-1}"
RUN_FULL_HALLAGRAM="${RUN_FULL_HALLAGRAM:-0}"
HALLA_DIAGNOSTIC="${HALLA_DIAGNOSTIC:-0}"

echo "[INFO] HAllA threads: $HALLA_THREADS"
echo "[INFO] RUN_EXTRA_PEARSON: $RUN_EXTRA_PEARSON"
echo "[INFO] RUN_FULL_HALLAGRAM: $RUN_FULL_HALLAGRAM"
echo "[INFO] HALLA_DIAGNOSTIC: $HALLA_DIAGNOSTIC"

run_halla_safe() {
    local xfile="$1"
    local yfile="$2"
    local outdir="$3"
    local metric="$4"
    local xlabel="$5"
    local ylabel="$6"
    local logfile="${outdir}.halla.log"

    echo "============================================================"
    echo "[HALLA] X: $xfile"
    echo "[HALLA] Y: $yfile"
    echo "[HALLA] Output: $outdir"
    echo "[HALLA] Metric: $metric"
    echo "[HALLA] Threads: $HALLA_THREADS"
    echo "============================================================"

    if [[ ! -s "$xfile" ]]; then
        echo "[WARNING] Missing HAllA X input: $xfile"
        return 0
    fi

    if [[ ! -s "$yfile" ]]; then
        echo "[WARNING] Missing HAllA Y input: $yfile"
        return 0
    fi

    mkdir -p "$(dirname "$outdir")"

    if [[ -d "$outdir" ]]; then
        local backup="${outdir}.previous_$(date +%Y%m%d_%H%M%S)"
        echo "[INFO] Existing HAllA output folder found. Moving to: $backup"
        mv "$outdir" "$backup"
    fi

    local halla_diag_args=""
if [[ "${HALLA_DIAGNOSTIC:-0}" == "1" ]]; then
    halla_diag_args="--diagnostic_plot"
fi

halla -x "$xfile" -y "$yfile" -o "$outdir" --x_dataset_label "$xlabel" --y_dataset_label "$ylabel" $halla_diag_args -m "$metric" --num_threads "$HALLA_THREADS" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] HAllA exited with status $status for: $outdir"
        echo "[WARNING] This often happens during report/hallagram plotting after the statistics were computed."
        echo "[WARNING] Log saved at: $logfile"

        mkdir -p "$outdir"

        {
            echo "HAllA finished with warning/error status: $status"
            echo "This may be caused by report/hallagram plotting, especially MatplotlibDeprecationWarning."
            echo "The MTD pipeline continued instead of stopping."
            echo "Log file:"
            echo "$logfile"
            echo
            echo "Relevant log lines:"
            grep -E "Number of significant|significant clusters|Traceback|MatplotlibDeprecationWarning|ERROR|WARNING" "$logfile" | tail -n 100
        } > "$outdir/HAllA_finished_with_warning.txt"

        echo "[INFO] Relevant HAllA log lines:"
        grep -E "Number of significant|significant clusters|Traceback|MatplotlibDeprecationWarning|ERROR|WARNING" "$logfile" | tail -n 40

        return 0
    fi

    echo "[OK] HAllA completed successfully: $outdir"
    echo "[OK] Log saved at: $logfile"
    return 0
}

run_hallagram_safe() {
    local indir="$1"
    local outfile="$2"
    local block_num="$3"
    local xlabel="$4"
    local ylabel="$5"
    local cbar="$6"
    local logfile="${outfile}.log"

    echo "============================================================"
    echo "[HALLAGRAM] Input: $indir"
    echo "[HALLAGRAM] Output: $outfile"
    echo "[HALLAGRAM] block_num: $block_num"
    echo "============================================================"

    if [[ ! -d "$indir" ]]; then
        echo "[WARNING] HAllA folder does not exist; skipping hallagram: $indir"
        return 0
    fi

    mkdir -p "$(dirname "$outfile")"

    hallagram -i "$indir" --cbar_label "$cbar" --x_dataset_label "$xlabel" --y_dataset_label "$ylabel" --output "$outfile" --block_num "$block_num" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] hallagram failed for: $outfile"
        echo "[WARNING] Log saved at: $logfile"
        grep -E "Traceback|Error|WARNING|MatplotlibDeprecationWarning" "$logfile" | tail -n 40
        return 0
    fi

    echo "[OK] hallagram saved: $outfile"
    return 0
}

run_python_plot_safe() {
    local label="$1"
    local cmd="$2"
    local logfile="$3"

    echo "============================================================"
    echo "[PYTHON] $label"
    echo "============================================================"

    eval "$cmd" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] $label failed; continuing."
        echo "[WARNING] Log saved at: $logfile"
        tail -n 40 "$logfile"
        return 0
    fi

    echo "[OK] $label completed."
    return 0
}

if [[ "$pdm" == "spearman" ]]; then
    pdm_name='Pairwise Spearman'
elif [[ "$pdm" == "pearson" ]]; then
    pdm_name='Pairwise Pearson'
elif [[ "$pdm" == "mi" ]]; then
    pdm_name='mi'
elif [[ "$pdm" == "nmi" ]]; then
    pdm_name='nmi'
elif [[ "$pdm" == "xicor" ]]; then
    pdm_name='xicor'
elif [[ "$pdm" == "dcor" ]]; then
    pdm_name='dcor'
else
    pdm_name="$pdm"
fi

echo "${g}Analyzing microbiome x host_genes associations...${w}"

run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_gene.txt" "$outputdr/halla/host_gene" "$pdm" "Microbiomes" "Host_gene"

run_python_plot_safe "PLS-DA microbiome x host_gene" "python $MTDIR/pls_da_analysis.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/pls_da_results.pdf" "$outputdr/halla/pls_da_analysis.log"

run_python_plot_safe "k-means microbiome x host_gene" "python $MTDIR/kmeans_clustering.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/kmeans_results.pdf -k 3" "$outputdr/halla/kmeans_clustering.log"

if [[ "$RUN_EXTRA_PEARSON" == "1" ]]; then
    echo "[INFO] RUN_EXTRA_PEARSON=1, running extra Pearson HAllA analysis."
    run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_gene.txt" "$outputdr/halla/pearson" "pearson" "Microbiomes" "Host_gene"
else
    echo "[INFO] Skipping extra Pearson HAllA run."
    echo "[INFO] To enable it: RUN_EXTRA_PEARSON=1 bash MTD_SE.sh ..."
fi

run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top5.pdf" 5 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top10.pdf" 10 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top25.pdf" 25 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top50.pdf" 50 "Microbiomes" "Host_gene" "$pdm_name"

if [[ "$RUN_FULL_HALLAGRAM" == "1" ]]; then
    run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_all.pdf" -1 "Microbiomes" "Host_gene" "$pdm_name"
else
    echo "[INFO] Skipping full host_gene hallagram by default."
    echo "[INFO] To enable it: RUN_FULL_HALLAGRAM=1 bash MTD_SE.sh ..."
fi

echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'

echo 'Analyzing microbiome x host_pathways associations...'
echo "${w}"

run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_score.txt" "$outputdr/halla/pathway" "$pdm" "Microbiomes" "Host_pathway"

run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top5.pdf" 5 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top10.pdf" 10 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top25.pdf" 25 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top50.pdf" 50 "Microbiomes" "Host_pathway" "$pdm_name"

if [[ "$RUN_FULL_HALLAGRAM" == "1" ]]; then
    run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_all.pdf" -1 "Microbiomes" "Host_pathway" "$pdm_name"
else
    echo "[INFO] Skipping full pathway hallagram by default."
fi

conda deactivate
conda activate MTD

echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "MTD running is finished"
echo -e "${w}"
