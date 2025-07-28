#!/usr/bin/env Rscript

# Definir mirror CRAN
options(repos = c(CRAN = "https://cran.r-project.org"))

#Apagar $offline
if (dir.exists("NCBI")) {
  unlink("NCBI", recursive = TRUE)
}

# Carregar pacotes necessários (instala se precisar)
if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
    BiocManager::install('AnnotationForge', force = TRUE, update = FALSE)
}
if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse")
}
if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt", update = FALSE)
}
if (!requireNamespace("httr", quietly = TRUE)) {
    install.packages("httr")
}

library(AnnotationForge)
library(optparse)
library(readr)
library(biomaRt)
library(httr)


# Definir opções de linha de comando
option_list <- list(
    make_option(c("-t", "--taxid"), type = "character", default = NULL,
                help = "NCBI Taxonomy ID of the species (required)", metavar = "character"),
    make_option(c("-d", "--dest_dir"), type = "character", default = getwd(),
                help = "Output directory for annotation package [default: current directory]", metavar = "character"),
    make_option(c("-o", "--offline"), type = "character", default = NULL,
                help = "Path to directory with pre-downloaded NCBI files (optional)", metavar = "character"),
    make_option(c("-c", "--copy"), action = "store_true", default = FALSE,
                help = "Copy offline files to './NCBI' folder before creating package")
)

# Parsear argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validar taxid obrigatório
if (is.null(opt$taxid)) {
    stop("The taxid must be provided using -t or --taxid")
}

# Informar destino
message("Output directory for annotation package: ", opt$dest_dir)

# Caminho do arquivo CSV (HostSpecies.csv) esperado na pasta atual (script rodando lá)
csv_file_path <- file.path(getwd(), "HostSpecies.csv")

if (!file.exists(csv_file_path)) {
    stop(paste("The file HostSpecies.csv was not found in current directory:", getwd()))
}

# Ler CSV
species_data <- read_csv(csv_file_path, show_col_types = FALSE)

if (!all(c("Taxon_ID", "Scientific_name") %in% colnames(species_data))) {
    stop("The CSV file must contain 'Taxon_ID' and 'Scientific_name' columns.")
}

# Função para extrair gênero e espécie
get_species_info_from_csv <- function(taxid, species_data) {
    matched <- species_data[species_data$Taxon_ID == as.numeric(taxid), ]
    if (nrow(matched) == 0) {
        stop(paste("No species found for TaxID:", taxid))
    }
    sp_name <- matched$Scientific_name[1]
    sp_split <- strsplit(sp_name, " ")[[1]]
    if (length(sp_split) != 2) stop("Scientific_name must be 'Genus species'")
    genus <- sp_split[1]
    species <- sp_split[2]
    message("Taxid: ", taxid)
    message("Genus: ", genus)
    message("Species: ", species)
    list(genus = genus, species = species)
}

species_info <- get_species_info_from_csv(opt$taxid, species_data)

# Copiar arquivos offline se solicitado
if (!is.null(opt$offline)) {
    message("Offline mode enabled.")
    if (!dir.exists(opt$offline)) {
        stop("Offline directory does not exist: ", opt$offline)
    }
    if (opt$copy) {
        message("Copying files from ", opt$offline, " to './NCBI' ...")
        dir.create("NCBI", showWarnings = FALSE)
        success <- file.copy(
            from = list.files(opt$offline, full.names = TRUE),
            to = "NCBI",
            overwrite = TRUE,
            recursive = TRUE
        )
        if (!all(success)) {
            warning("Some files could not be copied!")
        } else {
            message("Files copied successfully.")
        }
    } else {
        message("Note: You did not set --copy; offline files will NOT be copied to './NCBI'.")
        message("Ensure the files exist where makeOrgPackageFromNCBI expects them.")
    }
}


#Verificar se os arquivos de cache offline estao iguais ao NCBI 
# Function to compare file sizes and suggest updating if mismatch is found
check_ncbi_file_consistency <- function(local_path, remote_url) {
    if (!file.exists(local_path)) {
        warning("Local file not found: ", local_path)
        return(FALSE)
    }

    # Local file size
    local_size <- file.info(local_path)$size

    # Remote file size via HEAD request
    r <- httr::HEAD(remote_url)
    remote_size <- as.numeric(httr::headers(r)[["content-length"]])

    if (is.na(remote_size)) {
        warning("Could not retrieve remote file size for: ", remote_url)
        return(FALSE)
    }

    if (local_size != remote_size) {
        warning(paste0(
            "File mismatch detected: ", basename(local_path), "\n",
            "  - Local size: ", local_size, " bytes\n",
            "  - Remote size: ", remote_size, " bytes\n",
            "Suggestion: Your local file may be outdated or corrupted.\n",
            "Consider re-downloading it from: ", remote_url
        ))
        return(FALSE)
    }

    message("File is up to date: ", basename(local_path), " (OK)")
    return(TRUE)
}

# List of files and base path
ncbi_files <- c("gene2pubmed.gz", "gene2accession.gz", "gene2refseq.gz",
                "gene_info.gz", "gene2go.gz")

ncbi_base_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/"
local_dir <- file.path(getwd(), "NCBI")

# Check each file
for (f in ncbi_files) {
    local_file <- file.path(local_dir, f)
    remote_file <- paste0(ncbi_base_url, f)
    check_ncbi_file_consistency(local_file, remote_file)
}
# Rodar a criação do pacote
message("Creating annotation package...")

httr::set_config(httr::config(ssl_verifypeer = FALSE))

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Generic Author",
                       maintainer = "Generic Maintainer <maintainer@example.com>",
                       outputDir = opt$dest_dir,
                       tax_id = opt$taxid,
                       genus = species_info$genus,
                       species = species_info$species)

message("Annotation package created successfully.")
