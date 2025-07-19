#!/usr/bin/env Rscript

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
biomartCacheClear()

# Definir mirror CRAN
options(repos = c(CRAN = "https://cran.r-project.org"))

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
