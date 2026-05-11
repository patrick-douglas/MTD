#!/usr/bin/env Rscript

# ============================================================
# Volcano plot generator for microbiome differential abundance
# Labels are taken from the first column of the input table
# ============================================================

# -----------------------------
# Help message
# -----------------------------
print_help <- function() {
    cat("
Volcano plot generator for microbiome differential abundance results

Usage:
  Rscript EV.volcano_microbiome.R --de_results FILE.csv [options]

Required:
  --de_results FILE
      Input table containing differential abundance/expression results.

      The first column will be used as the label source.
      Required columns:
        log2FoldChange
        padj

Main options:
  -h, --help
      Show this help message and exit.

  --labels
      Show labels for all significant features.

      Significance is defined by:
        padj < --padj
        abs(log2FoldChange) >= --logfc

  --label_top N
      Show labels only for the top N features with the lowest padj.
      If this option is used, it overrides --labels.

  --padj VALUE
      Adjusted p-value / FDR cutoff.
      Default: 0.05

  --logfc VALUE
      Absolute log2 fold-change cutoff.
      Default: 2

  --group_label1 TEXT
      Name of the first group shown in the plot title.
      Default: Group1

  --group_label2 TEXT
      Name of the second group shown in the plot title.
      Default: Group2

  --output FILE.pdf
      Output PDF file name.
      Default: input file name with .pdf extension.

Input format options:
  --sep auto|comma|tab|semicolon
      Table separator.
      Default: auto

Plot size options:
  --width VALUE
      PDF width in inches.
      Default: 10

  --height VALUE
      PDF height in inches.
      Default: 8

Examples:

  1) Basic volcano plot without labels:
     Rscript volcano_microbiome.R \\
       --de_results microbiome_DE_results.csv

  2) Show labels for all significant features:
     Rscript volcano_microbiome.R \\
       --de_results microbiome_DE_results.csv \\
       --labels \\
       --padj 0.05 \\
       --logfc 2

  3) Show only the top 20 most significant features:
     Rscript volcano_microbiome.R \\
       --de_results microbiome_DE_results.csv \\
       --label_top 20

  4) Add custom group names:
     Rscript volcano_microbiome.R \\
       --de_results microbiome_DE_results.csv \\
       --label_top 20 \\
       --group_label1 Liver \\
       --group_label2 Telencephalon

  5) Save with a custom output name:
     Rscript volcano_microbiome.R \\
       --de_results microbiome_DE_results.csv \\
       --label_top 30 \\
       --output volcano_liver_vs_telencephalon.pdf

Notes:
  - Labels are taken from the first column of the input table.
  - The first column may contain species, genus, ASV, OTU, or feature names.
  - The x-axis is log2FoldChange.
  - The y-axis is -log10(FDR), using the padj column.

")
}

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    print_help()

    if (length(args) == 0) {
        quit(status = 1)
    } else {
        quit(status = 0)
    }
}

de_results <- NULL
output_file <- NULL

show_labels <- FALSE
label_top <- NULL

padj_cutoff <- 0.05
logfc_cutoff <- 2

group_label1 <- "Group1"
group_label2 <- "Group2"

sep_option <- "auto"

pdf_width <- 10
pdf_height <- 8

i <- 1

while (i <= length(args)) {

    arg <- args[i]

    if (arg == "--de_results") {
        if (i + 1 <= length(args)) {
            de_results <- args[i + 1]
            i <- i + 2
            next
        } else {
            stop("Error: Missing file name after --de_results.", call. = FALSE)
        }
    }

    if (arg == "--output") {
        if (i + 1 <= length(args)) {
            output_file <- args[i + 1]
            i <- i + 2
            next
        } else {
            stop("Error: Missing file name after --output.", call. = FALSE)
        }
    }

    if (arg == "--labels") {
        show_labels <- TRUE
        i <- i + 1
        next
    }

    if (arg == "--label_top") {
        if (i + 1 <= length(args)) {
            label_top <- as.numeric(args[i + 1])
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --label_top.", call. = FALSE)
        }
    }

    if (arg == "--padj") {
        if (i + 1 <= length(args)) {
            padj_cutoff <- as.numeric(args[i + 1])
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --padj.", call. = FALSE)
        }
    }

    if (arg == "--logfc") {
        if (i + 1 <= length(args)) {
            logfc_cutoff <- as.numeric(args[i + 1])
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --logfc.", call. = FALSE)
        }
    }

    if (arg == "--group_label1") {
        if (i + 1 <= length(args)) {
            group_label1 <- args[i + 1]
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --group_label1.", call. = FALSE)
        }
    }

    if (arg == "--group_label2") {
        if (i + 1 <= length(args)) {
            group_label2 <- args[i + 1]
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --group_label2.", call. = FALSE)
        }
    }

    if (arg == "--sep") {
        if (i + 1 <= length(args)) {
            sep_option <- args[i + 1]
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --sep.", call. = FALSE)
        }
    }

    if (arg == "--width") {
        if (i + 1 <= length(args)) {
            pdf_width <- as.numeric(args[i + 1])
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --width.", call. = FALSE)
        }
    }

    if (arg == "--height") {
        if (i + 1 <= length(args)) {
            pdf_height <- as.numeric(args[i + 1])
            i <- i + 2
            next
        } else {
            stop("Error: Missing value after --height.", call. = FALSE)
        }
    }

    stop(
        paste0(
            "Error: Unknown option: ", arg, "\n\n",
            "Use:\n",
            "  Rscript volcano_microbiome.R --help\n\n",
            "to see all available options."
        ),
        call. = FALSE
    )
}

# -----------------------------
# Validate arguments
# -----------------------------
if (is.null(de_results)) {
    cat("\nError: Missing required argument --de_results\n\n")
    print_help()
    quit(status = 1)
}

if (!file.exists(de_results)) {
    stop(
        paste0("Error: Input file not found: ", de_results),
        call. = FALSE
    )
}

if (!is.null(label_top)) {
    if (is.na(label_top) || label_top <= 0) {
        stop("Error: --label_top must be a positive number.", call. = FALSE)
    }
}

if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
    stop("Error: --padj must be a number between 0 and 1.", call. = FALSE)
}

if (is.na(logfc_cutoff) || logfc_cutoff < 0) {
    stop("Error: --logfc must be a positive number or zero.", call. = FALSE)
}

if (!(sep_option %in% c("auto", "comma", "tab", "semicolon"))) {
    stop(
        "Error: --sep must be one of: auto, comma, tab, semicolon.",
        call. = FALSE
    )
}

if (is.na(pdf_width) || pdf_width <= 0) {
    stop("Error: --width must be a positive number.", call. = FALSE)
}

if (is.na(pdf_height) || pdf_height <= 0) {
    stop("Error: --height must be a positive number.", call. = FALSE)
}

# -----------------------------
# Check and install packages
# -----------------------------
necessary_packages <- c(
    "dplyr",
    "ggplot2",
    "ggrepel",
    "EnhancedVolcano"
)

install_if_needed <- function(packages) {

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages(
            "BiocManager",
            repos = "https://cloud.r-project.org/"
        )
    }

    for (pkg in packages) {

        if (!requireNamespace(pkg, quietly = TRUE)) {

            message("Installing missing package: ", pkg)

            if (pkg == "EnhancedVolcano") {
                BiocManager::install(
                    pkg,
                    ask = FALSE,
                    update = FALSE
                )
            } else {
                install.packages(
                    pkg,
                    dependencies = TRUE,
                    repos = "https://cloud.r-project.org/"
                )
            }
        }

        suppressPackageStartupMessages(
            library(pkg, character.only = TRUE)
        )
    }
}

install_if_needed(necessary_packages)

# -----------------------------
# Detect separator
# -----------------------------
detect_separator <- function(file) {

    first_lines <- readLines(file, n = 5, warn = FALSE)
    first_lines <- first_lines[nchar(first_lines) > 0]

    if (length(first_lines) == 0) {
        stop("Error: Input file appears to be empty.", call. = FALSE)
    }

    first_line <- first_lines[1]

    n_tabs <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
    n_commas <- lengths(regmatches(first_line, gregexpr(",", first_line)))
    n_semicolons <- lengths(regmatches(first_line, gregexpr(";", first_line)))

    if (n_tabs >= n_commas && n_tabs >= n_semicolons && n_tabs > 0) {
        return("\t")
    }

    if (n_semicolons >= n_commas && n_semicolons > 0) {
        return(";")
    }

    return(",")
}

if (sep_option == "auto") {
    sep <- detect_separator(de_results)
} else if (sep_option == "comma") {
    sep <- ","
} else if (sep_option == "tab") {
    sep <- "\t"
} else if (sep_option == "semicolon") {
    sep <- ";"
}

message("Input file: ", de_results)
message("Detected/selected separator: ", ifelse(sep == "\t", "tab", sep))

# -----------------------------
# Load DE results
# -----------------------------
res2 <- read.table(
    de_results,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    fill = TRUE
)

if (ncol(res2) < 3) {
    stop(
        "Error: Input file has too few columns. Please check the separator using --sep.",
        call. = FALSE
    )
}

# The first column contains microbiome/taxon/feature names
colnames(res2)[1] <- "Feature_ID"

# Clean labels
res2$Feature_ID <- trimws(as.character(res2$Feature_ID))

# Check required columns
required_cols <- c("log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(res2))

if (length(missing_cols) > 0) {
    stop(
        paste(
            "Error: Missing required column(s):",
            paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
    )
}

# Convert important columns to numeric
res2$log2FoldChange <- suppressWarnings(
    as.numeric(res2$log2FoldChange)
)

res2$padj <- suppressWarnings(
    as.numeric(res2$padj)
)

# Remove rows without usable values
n_before <- nrow(res2)

res2 <- res2 |>
    dplyr::filter(
        !is.na(log2FoldChange),
        !is.na(padj),
        !is.na(Feature_ID),
        Feature_ID != ""
    )

n_after <- nrow(res2)

if (n_after == 0) {
    stop(
        "Error: No valid rows left after filtering NA values in log2FoldChange, padj, or first column.",
        call. = FALSE
    )
}

if (n_after < n_before) {
    message("Removed ", n_before - n_after, " row(s) with missing/invalid values.")
}

# Avoid problems with padj = 0
res2$padj_plot <- ifelse(
    res2$padj <= 0,
    .Machine$double.xmin,
    res2$padj
)

# -----------------------------
# Output file
# -----------------------------
if (is.null(output_file)) {

    output_file <- sub("\\.[^.]+$", ".pdf", basename(de_results))

    if (output_file == basename(de_results)) {
        output_file <- paste0(basename(de_results), ".pdf")
    }
}

if (!grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    output_file <- paste0(output_file, ".pdf")
}

message("Output PDF: ", output_file)

# -----------------------------
# Dynamic axis limits
# -----------------------------
x_min <- min(res2$log2FoldChange, na.rm = TRUE)
x_max <- max(res2$log2FoldChange, na.rm = TRUE)

x_abs <- max(abs(c(x_min, x_max)), na.rm = TRUE)
x_padding <- x_abs * 0.08

if (x_padding == 0 || is.na(x_padding)) {
    x_padding <- 1
}

x_range <- c(
    x_min - x_padding,
    x_max + x_padding
)

y_values <- -log10(res2$padj_plot)

y_max <- max(y_values, na.rm = TRUE)

if (is.infinite(y_max) || is.na(y_max)) {
    y_max <- 1
}

y_range <- c(
    0,
    y_max * 1.12
)

# -----------------------------
# Significant features
# -----------------------------
significant_features <- sum(
    res2$padj < padj_cutoff &
        abs(res2$log2FoldChange) >= logfc_cutoff,
    na.rm = TRUE
)

message("Total variables: ", nrow(res2))
message("Significant variables: ", significant_features)

# -----------------------------
# Select labels from first column
# -----------------------------
if (!is.null(label_top)) {

    labels_to_show <- res2 |>
        dplyr::arrange(
            padj,
            dplyr::desc(abs(log2FoldChange))
        ) |>
        dplyr::slice_head(n = label_top) |>
        dplyr::pull(Feature_ID)

    message("Label mode: top ", label_top, " by lowest padj.")

} else if (show_labels) {

    labels_to_show <- res2 |>
        dplyr::filter(
            padj < padj_cutoff,
            abs(log2FoldChange) >= logfc_cutoff
        ) |>
        dplyr::pull(Feature_ID)

    message("Label mode: all significant features.")

} else {

    labels_to_show <- character(0)
    message("Label mode: no labels.")
}

labels_to_show <- labels_to_show[
    !is.na(labels_to_show) &
        labels_to_show != ""
]

# -----------------------------
# Make volcano plot
# -----------------------------
pdf(
    file = output_file,
    width = pdf_width,
    height = pdf_height
)

volcano_plot <- EnhancedVolcano(
    res2,

    lab = res2$Feature_ID,
    selectLab = labels_to_show,

    x = "log2FoldChange",
    y = "padj_plot",

    title = bquote(.(group_label1)~italic(versus)~.(group_label2)),
    subtitle = "Differential abundance",

    xlab = bquote(~Log[2]~"Fold Change"),
    ylab = bquote(-1*Log[10]~"FDR"),

    pCutoff = padj_cutoff,
    FCcutoff = logfc_cutoff,

    pointSize = 2.0,
    labSize = 3.0,
    labFace = "bold",
    boxedLabels = TRUE,

    colAlpha = 1,

    cutoffLineType = "blank",

    hline = c(padj_cutoff),
    hlineCol = c("grey75"),
    hlineType = "longdash",
    hlineWidth = 0.4,

    vline = c(-logfc_cutoff, logfc_cutoff),
    vlineCol = c("grey75"),
    vlineType = "longdash",
    vlineWidth = 0.4,

    drawConnectors = TRUE,
    colConnectors = "grey50",
    maxoverlapsConnectors = Inf,
    widthConnectors = 0.5,

    border = "full",
    borderColour = "black",
    borderWidth = 0.5,

    gridlines.major = TRUE,
    gridlines.minor = FALSE,

    legendPosition = "top",
    legendLabSize = 12,
    legendIconSize = 3.0,
    legendLabels = c(
        "NS",
        expression(Log[2]~FC),
        "FDR",
        expression(FDR~and~Log[2]~FC)
    )
) +
    ggplot2::coord_cartesian(
        xlim = x_range,
        ylim = y_range,
        clip = "off"
    ) +
    ggplot2::scale_x_continuous(
        breaks = seq(
            floor(x_range[1]),
            ceiling(x_range[2]),
            length.out = 5
        )
    ) +
    ggplot2::scale_y_continuous(
        breaks = seq(
            0,
            ceiling(y_range[2]),
            length.out = 5
        )
    ) +
    ggplot2::annotate(
        "text",
        x = logfc_cutoff + 1,
        y = max(y_range) * 0.95,
        label = paste0("+", logfc_cutoff),
        color = "grey75",
        hjust = 0,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::annotate(
        "text",
        x = -logfc_cutoff - 1,
        y = max(y_range) * 0.95,
        label = paste0("-", logfc_cutoff),
        color = "grey75",
        hjust = 1,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::annotate(
        "text",
        x = min(x_range),
        y = -log10(padj_cutoff) + 0.5,
        label = as.character(padj_cutoff),
        color = "grey75",
        hjust = 0,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::labs(
        caption = paste(
            "Total =",
            nrow(res2),
            "variables | Significant =",
            significant_features
        )
    ) +
    ggplot2::theme(
        plot.margin = ggplot2::margin(10, 35, 10, 10)
    )

print(volcano_plot)

dev.off()

message("Done.")
message("Volcano plot saved as: ", output_file)
