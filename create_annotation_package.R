
# Load necessary packages
if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
    BiocManager::install('AnnotationForge', force = TRUE, update = FALSE)
}
if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse")
}
if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
}
library(AnnotationForge)
library(optparse)
library(readr)
library(biomaRt)
biomartCacheClear()
# Set a CRAN mirror to ensure package installation
options(repos = c(CRAN = "https://cran.r-project.org"))

# Define command-line options
option_list <- list(
    make_option(c("-t", "--taxid"), type = "character", default = NULL,
                help = "NCBI Taxonomy ID of the species", metavar = "character"),
    make_option(c("-d", "--dir"), type = "character", default = NULL,
                help = "Directory where the HostSpecies.csv is located", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if taxid and directory were provided
if (is.null(opt$taxid)) {
    stop("The taxid must be provided using the -t or --taxid option")
}
if (is.null(opt$dir)) {
    stop("The directory must be provided using the -d or --dir option")
}

# Define the path to the CSV file
csv_file_path <- file.path(opt$dir, "HostSpecies.csv")

# Check if the file exists
if (!file.exists(csv_file_path)) {
    stop(paste("The file HostSpecies.csv was not found in directory:", opt$dir))
}

# Read the CSV file
species_data <- read_csv(csv_file_path, show_col_types = FALSE)

# Check if required columns exist
if (!all(c("Taxon_ID", "Scientific_name") %in% colnames(species_data))) {
    stop("The CSV file must contain 'Taxon_ID' and 'Scientific_name' columns.")
}

# Function to retrieve genus and species names from the CSV file based on taxid
get_species_info_from_csv <- function(taxid, species_data) {
    # Filter the data to match the provided taxid (column 1)
    matched_species <- species_data[species_data$Taxon_ID == as.numeric(taxid), ]
    
    if (nrow(matched_species) == 0) {
        stop(paste("No species found for the provided taxid:", taxid))
    }
    
    # The species name is in the 'Scientific_name' column (formatted as "Genus species")
    species_name <- matched_species$Scientific_name[1]
    
    # Split the species name into genus and species
    species_split <- strsplit(species_name, " ")
    
    if (length(species_split[[1]]) != 2) {
        stop("The species name format should be 'Genus species'.")
    }
    
    genus <- species_split[[1]][1]
    species <- species_split[[1]][2]
    
    # Display informational messages
    message("Taxid: ", taxid)
    message("Genus: ", genus)
    message("Species: ", species)
    
    return(list(genus = genus, species = species))
}

# Retrieve species information using the provided taxid
species_info <- get_species_info_from_csv(opt$taxid, species_data)

# Ensure that species and genus are single strings
if (length(species_info$species) != 1 || length(species_info$genus) != 1) {
    stop("'species' and 'genus' must be single strings.")
}

# Create the annotation package using AnnotationForge
message("
Creating annotation package...
")
httr::set_config(httr::config(ssl_verifypeer = FALSE))
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Generic Author",
                       maintainer = "Generic Maintainer <maintainer@example.com>",
                       outputDir = ".",
                       tax_id = opt$taxid,
                       genus = species_info$genus,
                       species = species_info$species)

message("Annotation package created successfully.")

