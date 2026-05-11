#!/usr/bin/env Rscript

# ============================================================
# Krona species Venn + Euler diagrams - paper style
# ------------------------------------------------------------
# - Does not depend on Conda
# - Checks required R packages
# - Tries to install missing packages automatically
# - Generates:
#     1) Venn diagram using ggVennDiagram
#     2) Proportional Euler diagram using eulerr
# - Outputs only PDF and SVG
# ============================================================

# ============================================================
# Help
# ============================================================

print_help <- function(exit_status = 0) {
  cat("\nGenerate publication-style Venn and Euler diagrams from Krona/MPA-like files.\n\nThe script extracts species names after 's__' from one or more files per group,\ncalculates unique/shared species, and generates both:\n\n  1) Venn diagram\n  2) Proportional Euler diagram\n\nThe Euler diagram is useful when you want the circle areas to approximately\nrepresent the number of species in each group.\n\nThis script does NOT depend on any Conda environment.\nIt uses regular R packages from CRAN.\n\nRequired R packages:\n  ggplot2\n  ggVennDiagram\n  eulerr\n\nIf packages are missing, the script will try to install them automatically.\nIf installation fails, it will print the command you need to run manually.\n\nUsage:\n  Rscript krona_venn_paper.R \\\n    --krona_files1 file1.krona file2.krona \\\n    --krona_files2 file3.krona file4.krona \\\n    --group_label1 Liver \\\n    --group_label2 Telencephalon\n\nRequired arguments:\n  --krona_files1       One or more Krona files for group 1\n  --krona_files2       One or more Krona files for group 2\n  --group_label1       Label for group 1\n  --group_label2       Label for group 2\n\nOptional arguments:\n  --title              Plot title\n                       Default: 'Total species detected in the microbiomes'\n\n  --subtitle           Plot subtitle\n                       Default: automatic summary\n\n  --output_prefix      Output file prefix\n                       Default: species_overlap_<group1>_vs_<group2>\n\n  --plot_type          Which plot to generate:\n                         both\n                         venn\n                         euler\n                       Default: both\n\n  --width              Figure width in inches\n                       Default: 10\n\n  --height             Figure height in inches\n                       Default: 8\n\n  --label_type         Label type for plot regions:\n                         count\n                         percent\n                         both\n                       Default: count\n\n  --write_lists        Write TXT files with species lists:\n                         <prefix>_group1_species.txt\n                         <prefix>_group2_species.txt\n                         <prefix>_unique_group1.txt\n                         <prefix>_unique_group2.txt\n                         <prefix>_shared_species.txt\n\n  --replace_underscores\n                       Replace underscores with spaces in species names\n\n  --no_install         Do not try to install missing packages automatically\n\n  -h, --help           Show this help message\n\nExamples:\n  Rscript krona_venn_paper.R \\\n    --krona_files1 ../krona/*LIVER* \\\n    --group_label1 Liver \\\n    --krona_files2 ../krona/*TEL* \\\n    --group_label2 Telencephalon\n\n  Rscript krona_venn_paper.R \\\n    --krona_files1 ../krona/*LIVER* \\\n    --group_label1 Liver \\\n    --krona_files2 ../krona/*TEL* \\\n    --group_label2 Telencephalon \\\n    --label_type both \\\n    --write_lists\n\n  Rscript krona_venn_paper.R \\\n    --krona_files1 ../krona/*LIVER* \\\n    --group_label1 Liver \\\n    --krona_files2 ../krona/*TEL* \\\n    --group_label2 Telencephalon \\\n    --plot_type euler\n\nOutputs:\n  <prefix>_venn.pdf\n  <prefix>_venn.svg\n\n  <prefix>_euler.pdf\n  <prefix>_euler.svg\n\n  <prefix>_summary.tsv\n\nIf --write_lists is used:\n  <prefix>_group1_species.txt\n  <prefix>_group2_species.txt\n  <prefix>_unique_group1.txt\n  <prefix>_unique_group2.txt\n  <prefix>_shared_species.txt\n\n")
  quit(status = exit_status)
}

# ============================================================
# Argument parser
# ============================================================

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    print_help(0)
  }

  opts <- list(
    krona_files1 = NULL,
    krona_files2 = NULL,
    group_label1 = NULL,
    group_label2 = NULL,
    title = "Total species detected",
    subtitle = NULL,
    output_prefix = NULL,
    plot_type = "both",
    width = 10,
    height = 8,
    label_type = "count",
    write_lists = FALSE,
    replace_underscores = FALSE,
    auto_install = TRUE
  )

  known_options <- c(
    "krona_files1",
    "krona_files2",
    "group_label1",
    "group_label2",
    "title",
    "subtitle",
    "output_prefix",
    "plot_type",
    "width",
    "height",
    "label_type",
    "write_lists",
    "replace_underscores",
    "no_install"
  )

  i <- 1

  while (i <= length(args)) {
    key <- args[i]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    name <- sub("^--", "", key)

    if (!name %in% known_options) {
      stop("Unknown option: --", name, call. = FALSE)
    }

    if (name == "write_lists") {
      opts$write_lists <- TRUE
      i <- i + 1
      next
    }

    if (name == "replace_underscores") {
      opts$replace_underscores <- TRUE
      i <- i + 1
      next
    }

    if (name == "no_install") {
      opts$auto_install <- FALSE
      i <- i + 1
      next
    }

    i <- i + 1
    values <- character()

    while (i <= length(args) && !startsWith(args[i], "--")) {
      values <- c(values, args[i])
      i <- i + 1
    }

    if (length(values) == 0) {
      stop("Missing value for option --", name, call. = FALSE)
    }

    if (name %in% c("krona_files1", "krona_files2")) {
      opts[[name]] <- values
    } else if (name %in% c("width", "height")) {
      opts[[name]] <- as.numeric(values[1])
    } else {
      opts[[name]] <- paste(values, collapse = " ")
    }
  }

  required <- c("krona_files1", "krona_files2", "group_label1", "group_label2")

  for (r in required) {
    value <- opts[[r]]

    missing_value <- is.null(value) ||
      length(value) == 0 ||
      all(is.na(value)) ||
      all(trimws(as.character(value)) == "")

    if (missing_value) {
      stop("Missing required argument --", r, call. = FALSE)
    }
  }

  if (!opts$plot_type %in% c("both", "venn", "euler")) {
    stop("--plot_type must be one of: both, venn, euler", call. = FALSE)
  }

  if (!opts$label_type %in% c("count", "percent", "both")) {
    stop("--label_type must be one of: count, percent, both", call. = FALSE)
  }

  if (is.na(opts$width) || opts$width <= 0) {
    stop("--width must be a positive number.", call. = FALSE)
  }

  if (is.na(opts$height) || opts$height <= 0) {
    stop("--height must be a positive number.", call. = FALSE)
  }

  opts
}

# ============================================================
# Package check / install
# ============================================================

setup_user_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER")

  if (user_lib == "") {
    user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
    Sys.setenv(R_LIBS_USER = user_lib)
  }

  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }

  .libPaths(unique(c(user_lib, .libPaths())))

  invisible(user_lib)
}

check_and_install_packages <- function(pkgs, auto_install = TRUE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) == 0) {
    message("[OK] Required R packages are already installed.")
    return(invisible(TRUE))
  }

  message("[INFO] Missing R package(s): ", paste(missing_pkgs, collapse = ", "))

  install_cmd <- paste0(
    "Rscript -e 'install.packages(c(",
    paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "),
    "), repos=\"https://cloud.r-project.org\")'"
  )

  if (!auto_install) {
    stop(
      "Missing required R package(s): ",
      paste(missing_pkgs, collapse = ", "),
      "\n\nInstall manually with:\n  ",
      install_cmd,
      "\n",
      call. = FALSE
    )
  }

  message("[INSTALL] Trying to install missing packages from CRAN...")

  setup_user_library()

  install_ok <- TRUE

  tryCatch(
    {
      install.packages(
        missing_pkgs,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
      )
    },
    error = function(e) {
      install_ok <<- FALSE
      message("[ERROR] Automatic package installation failed.")
      message("[ERROR] ", conditionMessage(e))
    }
  )

  still_missing <- missing_pkgs[
    !vapply(missing_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (!install_ok || length(still_missing) > 0) {
    stop(
      "Could not install/load the following package(s): ",
      paste(still_missing, collapse = ", "),
      "\n\nPlease install manually with:\n  ",
      install_cmd,
      "\n\nIf you do not have write permission to the system R library, try:\n",
      "  mkdir -p ~/R/library\n",
      "  echo 'R_LIBS_USER=~/R/library' >> ~/.Renviron\n",
      "  ",
      install_cmd,
      "\n",
      call. = FALSE
    )
  }

  message("[OK] Missing packages installed successfully.")

  invisible(TRUE)
}

# ============================================================
# Utility functions
# ============================================================

sanitize_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x
}

expand_files <- function(files) {
  expanded <- unlist(lapply(files, function(f) {
    hit <- Sys.glob(f)

    if (length(hit) == 0) {
      return(f)
    }

    hit
  }))

  expanded <- unique(expanded)

  missing_files <- expanded[!file.exists(expanded)]

  if (length(missing_files) > 0) {
    stop(
      "The following file(s) do not exist:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  expanded
}

clean_species_name <- function(x, replace_underscores = FALSE) {
  x <- trimws(x)

  # Remove content after common separators used in Krona/MPA-like files.
  x <- sub("\t.*$", "", x)
  x <- sub(";.*$", "", x)
  x <- sub("\\|.*$", "", x)

  # Remove repeated spaces.
  x <- gsub("[[:space:]]+", " ", x)

  if (replace_underscores) {
    x <- gsub("_", " ", x)
    x <- gsub("[[:space:]]+", " ", x)
  }

  trimws(x)
}

extract_species_list <- function(krona_files, replace_underscores = FALSE) {
  combined_species <- character()

  for (krona_file in krona_files) {
    message("[READING] ", krona_file)

    lines <- readLines(krona_file, warn = FALSE)

    species_lines <- grep("s__", lines, value = TRUE, fixed = TRUE)

    if (length(species_lines) == 0) {
      warning("No lines containing 's__' found in: ", krona_file)
      next
    }

    # Extract text after the LAST occurrence of s__.
    species <- sub("^.*s__", "", species_lines)

    species <- clean_species_name(
      species,
      replace_underscores = replace_underscores
    )

    # Remove empty and uninformative entries.
    species <- species[species != ""]
    species <- species[
      !grepl(
        "^(NA|NaN|null|unknown|unclassified|uncultured)$",
        species,
        ignore.case = TRUE
      )
    ]

    combined_species <- c(combined_species, species)
  }

  sort(unique(combined_species))
}

write_species_file <- function(x, file) {
  writeLines(sort(unique(x)), con = file)
}

safe_metric <- function(x, field) {
  tryCatch(
    {
      if (!is.null(x[[field]])) {
        as.numeric(x[[field]])[1]
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
}

# ============================================================
# Save functions
# ============================================================

save_ggplot_pdf_svg <- function(plot, prefix, width, height) {
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  grDevices::svg(
    filename = svg_file,
    width = width,
    height = height,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  message("[OK] Files generated:")
  message("  ", pdf_file)
  message("  ", svg_file)
}

make_euler_grob <- function(fit, quantities_arg) {
  grid::grid.grabExpr({
    euler_plot <- plot(
      fit,
      fills = list(
  fill = c("#E9C9CF", "#BFE5B8"),
  alpha = 0.65
),
      edges = list(
        col = "black",
        lwd = 2
      ),
      labels = list(
        cex = 1.25,
        font = 2
      ),
      quantities = quantities_arg,
      main = NULL
    )

    # Some versions of eulerr draw directly; others return a grob.
    # This makes the function robust across package versions.
    if (!is.null(euler_plot)) {
      try(grid::grid.draw(euler_plot), silent = TRUE)
    }
  })
}

save_euler_pdf_svg <- function(fit,
                               prefix,
                               width,
                               height,
                               title,
                               subtitle,
                               caption,
                               label_type) {
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  if (label_type == "percent") {
    quantities_arg <- list(type = "percent", cex = 1.1, font = 2)
  } else if (label_type == "both") {
    quantities_arg <- list(type = c("counts", "percent"), cex = 1.0, font = 2)
  } else {
    quantities_arg <- list(type = "counts", cex = 1.2, font = 2)
  }

  euler_grob <- make_euler_grob(
    fit = fit,
    quantities_arg = quantities_arg
  )

  draw_one <- function() {
    grid::grid.newpage()

    grid::pushViewport(
      grid::viewport(
        x = 0.5,
        y = 0.50,
        width = 0.88,
        height = 0.72
      )
    )
    grid::grid.draw(euler_grob)
    grid::popViewport()

    grid::grid.text(
      title,
      x = 0.5,
      y = 0.97,
      gp = grid::gpar(fontsize = 18, fontface = "bold")
    )

    grid::grid.text(
      subtitle,
      x = 0.5,
      y = 0.925,
      gp = grid::gpar(fontsize = 10)
    )

    grid::grid.text(
      caption,
      x = 0.5,
      y = 0.035,
      gp = grid::gpar(fontsize = 8, col = "grey35")
    )
  }

  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    bg = "white"
  )
  draw_one()
  grDevices::dev.off()

  grDevices::svg(
    filename = svg_file,
    width = width,
    height = height,
    bg = "white"
  )
  draw_one()
  grDevices::dev.off()

  message("[OK] Files generated:")
  message("  ", pdf_file)
  message("  ", svg_file)
}

# ============================================================
# Main
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(args)

required_packages <- c("ggplot2", "ggVennDiagram", "eulerr")

check_and_install_packages(
  pkgs = required_packages,
  auto_install = opt$auto_install
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(eulerr)
  library(grid)
})

krona_files1 <- expand_files(opt$krona_files1)
krona_files2 <- expand_files(opt$krona_files2)

label1 <- opt$group_label1
label2 <- opt$group_label2

if (is.null(opt$output_prefix)) {
  output_prefix <- paste0(
    "species_overlap_",
    sanitize_name(label1),
    "_vs_",
    sanitize_name(label2)
  )
} else {
  output_prefix <- opt$output_prefix
}

message("============================================================")
message("Krona species overlap diagrams")
message("No Conda environment required")
message("Plot type: ", opt$plot_type)
message("Group 1: ", label1)
message("Group 2: ", label2)
message("Files group 1: ", length(krona_files1))
message("Files group 2: ", length(krona_files2))
message("Output prefix: ", output_prefix)
message("============================================================")

species1 <- extract_species_list(
  krona_files1,
  replace_underscores = opt$replace_underscores
)

species2 <- extract_species_list(
  krona_files2,
  replace_underscores = opt$replace_underscores
)

unique1 <- setdiff(species1, species2)
unique2 <- setdiff(species2, species1)
shared <- intersect(species1, species2)

n1 <- length(species1)
n2 <- length(species2)
n_unique1 <- length(unique1)
n_unique2 <- length(unique2)
n_shared <- length(shared)
n_union <- length(union(species1, species2))

message("")
message("Summary:")
message("  ", label1, " total species: ", n1)
message("  ", label2, " total species: ", n2)
message("  Unique to ", label1, ": ", n_unique1)
message("  Unique to ", label2, ": ", n_unique2)
message("  Shared species: ", n_shared)
message("  Total union: ", n_union)
message("")

if (n1 == 0 && n2 == 0) {
  stop(
    "No species were detected in either group. ",
    "Check your Krona files and whether they contain 's__'.",
    call. = FALSE
  )
}

if (is.null(opt$subtitle)) {
  subtitle_text <- paste0(
    label1, ": ", n1,
    " species | ",
    label2, ": ", n2,
    " species | Shared: ", n_shared,
    " | Union: ", n_union
  )
} else {
  subtitle_text <- opt$subtitle
}

caption_text <- paste0(
  "Species were extracted from Krona/MPA-like entries containing 's__'. ",
  "Unique ", label1, ": ", n_unique1,
  "; unique ", label2, ": ", n_unique2,
  "; shared: ", n_shared, "."
)

# ============================================================
# Euler fit and diagnostics
# ============================================================

euler_fit <- NULL
euler_stress <- NA_real_
euler_diag_error <- NA_real_

if (n1 > 0 && n2 > 0) {
  euler_counts <- numeric(3)

  names(euler_counts) <- c(
    label1,
    label2,
    paste(label1, label2, sep = "&")
  )

  # eulerr expects disjoint region counts:
  # A only, B only, and A&B shared.
  euler_counts[label1] <- n_unique1
  euler_counts[label2] <- n_unique2
  euler_counts[paste(label1, label2, sep = "&")] <- n_shared

  euler_fit <- eulerr::euler(euler_counts)

  euler_stress <- safe_metric(euler_fit, "stress")
  euler_diag_error <- safe_metric(euler_fit, "diagError")
}

# ============================================================
# Summary table
# ============================================================

summary_df <- data.frame(
  comparison = paste(label1, "vs", label2),
  group1 = label1,
  group2 = label2,
  group1_total_species = n1,
  group2_total_species = n2,
  unique_to_group1 = n_unique1,
  unique_to_group2 = n_unique2,
  shared_species = n_shared,
  union_species = n_union,
  euler_stress = euler_stress,
  euler_diag_error = euler_diag_error,
  stringsAsFactors = FALSE
)

summary_file <- paste0(output_prefix, "_summary.tsv")

write.table(
  summary_df,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[OK] Summary table written:")
message("  ", summary_file)

if (opt$write_lists) {
  write_species_file(species1, paste0(output_prefix, "_group1_species.txt"))
  write_species_file(species2, paste0(output_prefix, "_group2_species.txt"))
  write_species_file(unique1, paste0(output_prefix, "_unique_group1.txt"))
  write_species_file(unique2, paste0(output_prefix, "_unique_group2.txt"))
  write_species_file(shared, paste0(output_prefix, "_shared_species.txt"))

  message("[OK] Species TXT files written:")
  message("  ", output_prefix, "_group1_species.txt")
  message("  ", output_prefix, "_group2_species.txt")
  message("  ", output_prefix, "_unique_group1.txt")
  message("  ", output_prefix, "_unique_group2.txt")
  message("  ", output_prefix, "_shared_species.txt")
}

# ============================================================
# Venn diagram
# ============================================================

if (opt$plot_type %in% c("both", "venn")) {
  message("")
  message("[PLOT] Generating Venn diagram...")

  venn_input <- list()
  venn_input[[label1]] <- species1
  venn_input[[label2]] <- species2

  p_venn <- suppressMessages({
    ggVennDiagram::ggVennDiagram(
      venn_input,
      label = opt$label_type,
      label_alpha = 0,
      label_geom = "label",
      label_size = 5,
      edge_size = 1.1,
      set_size = 5.2
    ) +
      ggplot2::scale_fill_gradientn(
  colours = c("#FFF8F9", "#F1D5DB", "#D9ECD2", "#BFE5B8"),
  name = "Species\ncount"
) +
      ggplot2::labs(
        title = opt$title,
        subtitle = subtitle_text,
        caption = caption_text
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.10, 0.06))
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::coord_fixed(clip = "off") +
      ggplot2::theme_void(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 22,
          margin = ggplot2::margin(b = 8)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 13,
          margin = ggplot2::margin(b = 20)
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0.5,
          size = 10,
          color = "grey35",
          margin = ggplot2::margin(t = 18)
        ),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        plot.margin = ggplot2::margin(20, 25, 20, 55)
      )
  })

  save_ggplot_pdf_svg(
    plot = p_venn,
    prefix = paste0(output_prefix, "_venn"),
    width = opt$width,
    height = opt$height
  )
}

# ============================================================
# Euler diagram
# ============================================================

if (opt$plot_type %in% c("both", "euler")) {
  message("")
  message("[PLOT] Generating proportional Euler diagram...")

  if (is.null(euler_fit)) {
    warning(
      "Euler diagram was skipped because one of the groups has zero detected species."
    )
  } else {
    euler_caption <- paste0(
      "Euler diagram: circle areas are fitted to approximate species counts. ",
      "Stress: ", signif(euler_stress, 4),
      "; diagnostic error: ", signif(euler_diag_error, 4), "."
    )

    save_euler_pdf_svg(
      fit = euler_fit,
      prefix = paste0(output_prefix, "_euler"),
      width = opt$width,
      height = opt$height,
      title = opt$title,
      subtitle = subtitle_text,
      caption = euler_caption,
      label_type = opt$label_type
    )
  }
}

message("")
message("[DONE] Finished successfully.")
message("")
message("Suggested interpretation:")
message("  - Use the Venn diagram when you want a clean conceptual overlap figure.")
message("  - Use the Euler diagram when you want the circle areas to reflect species counts.")
message("  - Check ", output_prefix, "_summary.tsv for counts and Euler fit diagnostics.")

