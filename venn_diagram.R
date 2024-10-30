# Set a default CRAN mirror (you can adjust to a closer mirror if necessary)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to check and install packages if necessary
install_if_missing <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = TRUE)
    }
}

# Check and install required packages
install_if_missing("VennDiagram")
install_if_missing("argparser")

# Function to extract species list from .krona files
extract_species_list <- function(krona_files) {
    combined_species <- c()
    for (krona_file in krona_files) {
        temp_file <- paste0(krona_file, "_species_list.txt")
        grep_command <- sprintf("grep 's__' %s | awk -F 's__' '{print $2}' | sort | uniq > %s", krona_file, temp_file)
        system(grep_command)
        species <- scan(temp_file, what = "", sep = "\n", quiet = TRUE)
        combined_species <- unique(c(combined_species, species))
        file.remove(temp_file)
    }
    return(combined_species)
}

# Create argument parser
p <- arg_parser("Generate a Venn diagram based on multiple input krona files.")
p <- add_argument(p, "--krona_files1", help="Paths to first group krona file(s)", nargs = Inf)
p <- add_argument(p, "--krona_files2", help="Paths to second group krona file(s)", nargs = Inf)
p <- add_argument(p, "--group_label1", help="Label for first group")
p <- add_argument(p, "--group_label2", help="Label for second group")

# Parse arguments
args <- parse_args(p)

# Extract species from the two groups of krona files
species_group1 <- extract_species_list(args$krona_files1)
species_group2 <- extract_species_list(args$krona_files2)

# Calculate unique and shared species
group1_unique <- length(setdiff(species_group1, species_group2))          # Unique species for group 1
group2_unique <- length(setdiff(species_group2, species_group1))          # Unique species for group 2
intersection <- length(intersect(species_group1, species_group2))         # Shared species

# Define a scale factor for the circles
scale_factor <- 0.75  # Reduce the size of the circles proportionally

# Define circle radii based on adjusted factors and apply the scale factor
radius1 <- scale_factor * (0.15 + (group1_unique + intersection) / 2000)
radius2 <- scale_factor * (0.15 + (group2_unique + intersection) / 2000)

# Calculate automatic coordinates for positioning
# Circle 1 (left)
center_x1 <- 0.35
center_y1 <- 0.55

# Circle 2 (right)
center_x2 <- 0.55
center_y2 <- 0.55

# Intersection (center)
center_intersection_x <- (center_x1 + center_x2) / 2
center_intersection_y <- center_y1  # They share the same Y axis

# Define additional horizontal offset
offset <- 0.05  # Small offset to separate further to the sides

# Draw circles manually with adjusted sizes
draw_venn_diagram <- function() {
    # Draw circle for group 1 (left)
    grid.circle(x = center_x1, y = center_y1, r = radius1, gp = gpar(fill = "#fc9272", alpha = 0.5, col = "black"))
    
    # Draw circle for group 2 (right)
    grid.circle(x = center_x2, y = center_y2, r = radius2, gp = gpar(fill = "#a1d99b", alpha = 0.5, col = "black"))
    
    # Center the numbers in the corresponding areas with horizontal offset
    grid.text(label = as.character(group1_unique), x = center_x1 - radius1 / 2 - offset, y = center_y1, gp = gpar(fontsize = 18, fontface = "plain"))  # Position in red
    grid.text(label = as.character(group2_unique), x = center_x2 + radius2 / 2 + offset, y = center_y2, gp = gpar(fontsize = 18, fontface = "plain"))  # Position in green
    grid.text(label = as.character(intersection), x = center_intersection_x, y = center_intersection_y, gp = gpar(fontsize = 18, fontface = "plain"))   # Position in the center (yellow)
    
    # Title at the top (larger than the numbers and labels)
    grid.text("Total number of species detected in the microbiome", y = 0.98, gp = gpar(fontsize = 22, fontface = "bold"))
    
    # Group labels below the circles (centered and italicized)
    grid.text(args$group_label1, x = 0.35, y = 0.1, gp = gpar(fontsize = 20, fontface = "italic"))  # Centered, moved down, and italicized
    grid.text(args$group_label2, x = 0.55, y = 0.1, gp = gpar(fontsize = 20, fontface = "italic"))  # Centered, moved down, and italicized
}

# Save the diagram in various formats, with an increased canvas
output_base <- sprintf("venn_diagram_%s_vs_%s", args$group_label1, args$group_label2)

# Adjustment for PNG, with a larger canvas
png(paste0(output_base, ".png"), width = 3500, height = 3000, res = 300)  # Increase the canvas
draw_venn_diagram()
dev.off()

pdf(paste0(output_base, ".pdf"), width = 16, height = 14)  # Large canvas for PDF
draw_venn_diagram()
dev.off()

svg(paste0(output_base, ".svg"), width = 16, height = 14)  # Large canvas for SVG
draw_venn_diagram()
dev.off()

# Remove the Rplots.pdf file if it exists
if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
}