#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  ix <- which(args == flag)
  if (length(ix) == 0) return(default)
  if (ix[1] == length(args)) stop("Missing value for argument: ", flag)
  args[ix[1] + 1]
}

scores_file <- get_arg("--scores")
samplesheet_file <- get_arg("--samplesheet")
outdir <- get_arg("--outdir")
top_var <- as.integer(get_arg("--top-var", "50"))
top_diff <- as.integer(get_arg("--top-diff", "20"))
pca_top_var <- as.integer(get_arg("--pca-top-var", "500"))

if (is.null(scores_file) || is.null(samplesheet_file) || is.null(outdir)) {
  stop(
    "Usage:\n",
    "Rscript plot_ssgsea_results.R ",
    "--scores ssgsea-results-scores.gct ",
    "--samplesheet samplesheet.csv ",
    "--outdir ssGSEA/plots ",
    "[--top-var 50] [--top-diff 20] [--pca-top-var 500]\n"
  )
}

if (!file.exists(scores_file)) {
  stop("Scores file not found: ", scores_file)
}

if (!file.exists(samplesheet_file)) {
  stop("Samplesheet not found: ", samplesheet_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] ssGSEA scores: ", scores_file)
message("[INFO] samplesheet: ", samplesheet_file)
message("[INFO] output directory: ", outdir)

read_gct_scores <- function(file) {
  header_lines <- readLines(file, n = 2)

  if (!grepl("^#1\\.", header_lines[1])) {
    stop("Input file does not look like a GCT file: ", file)
  }

  dims <- strsplit(header_lines[2], "[\t ]+")[[1]]
  dims <- dims[dims != ""]

  if (length(dims) < 2) {
    stop("Could not parse GCT dimensions from second line: ", header_lines[2])
  }

  n_samples <- as.integer(dims[2])

  tab <- read.delim(
    file,
    skip = 2,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(tab) < n_samples + 1) {
    stop("GCT file has fewer columns than expected.")
  }

  term_ids <- as.character(tab[[1]])
  sample_cols <- tail(colnames(tab), n_samples)

  mat <- as.matrix(tab[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- term_ids

  desc_col <- NULL
  if (ncol(tab) >= 2) {
    desc_col <- as.character(tab[[2]])
  } else {
    desc_col <- term_ids
  }

  desc <- data.frame(
    term = term_ids,
    description = desc_col,
    stringsAsFactors = FALSE
  )

  list(mat = mat, desc = desc, sample_cols = sample_cols)
}

read_samplesheet <- function(file, score_samples) {
  ss <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(ss) < 2) {
    stop("Samplesheet must have at least two columns: sample and group.")
  }

  cn <- colnames(ss)
  cn_lower <- tolower(cn)

  sample_col <- cn[1]

  group_candidates <- c("group", "condition", "treatment", "grupo")
  group_col <- NULL

  for (candidate in group_candidates) {
    hit <- which(cn_lower == candidate)
    if (length(hit) > 0) {
      group_col <- cn[hit[1]]
      break
    }
  }

  if (is.null(group_col)) {
    group_col <- cn[2]
  }

  meta <- data.frame(
    sample = as.character(ss[[sample_col]]),
    group = as.character(ss[[group_col]]),
    stringsAsFactors = FALSE
  )

  meta$sample <- trimws(meta$sample)
  meta$group <- trimws(meta$group)

  meta <- meta[meta$sample != "" & meta$group != "", , drop = FALSE]
  meta <- meta[!duplicated(meta$sample), , drop = FALSE]

  missing_in_scores <- setdiff(meta$sample, score_samples)
  missing_in_sheet <- setdiff(score_samples, meta$sample)

  if (length(missing_in_scores) > 0) {
    warning(
      "Samples present in samplesheet but absent from ssGSEA scores: ",
      paste(missing_in_scores, collapse = ", ")
    )
  }

  if (length(missing_in_sheet) > 0) {
    warning(
      "Samples present in ssGSEA scores but absent from samplesheet: ",
      paste(missing_in_sheet, collapse = ", ")
    )
  }

  meta <- meta[meta$sample %in% score_samples, , drop = FALSE]

  if (nrow(meta) < 2) {
    stop("Fewer than two samples matched between samplesheet and ssGSEA scores.")
  }

  meta$group <- factor(meta$group, levels = unique(meta$group))

  meta
}

read_contrasts <- function(file, meta) {
  ss <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )

  groups <- levels(meta$group)

  empty_contrasts <- data.frame(
    group1 = character(0),
    group2 = character(0),
    source = character(0),
    stringsAsFactors = FALSE
  )

  if (length(groups) < 2L) {
    message(
      "[INFO] Differential ssGSEA: fewer than two groups; ",
      "no pairwise contrasts are possible."
    )
    return(empty_contrasts)
  }

  cn_lower <- tolower(colnames(ss))
  group1_index <- match("group1", cn_lower, nomatch = 0L)
  group2_index <- match("group2", cn_lower, nomatch = 0L)

  if (xor(group1_index > 0L, group2_index > 0L)) {
    stop(
      "Samplesheet must contain both group1 and group2 columns, ",
      "or neither of them."
    )
  }

  if (group1_index > 0L && group2_index > 0L) {
    group1 <- trimws(as.character(ss[[group1_index]]))
    group2 <- trimws(as.character(ss[[group2_index]]))

    group1[is.na(group1)] <- ""
    group2[is.na(group2)] <- ""

    incomplete <- xor(group1 != "", group2 != "")

    if (any(incomplete)) {
      bad_rows <- which(incomplete) + 1L
      stop(
        "Incomplete group1/group2 comparison in samplesheet row(s): ",
        paste(bad_rows, collapse = ", ")
      )
    }

    explicit <- data.frame(
      group1 = group1[group1 != "" & group2 != ""],
      group2 = group2[group1 != "" & group2 != ""],
      stringsAsFactors = FALSE
    )

    explicit <- unique(explicit)

    if (nrow(explicit) > 0L) {
      invalid_group <- (
        !(explicit$group1 %in% groups) |
          !(explicit$group2 %in% groups)
      )

      if (any(invalid_group)) {
        invalid_pairs <- paste0(
          explicit$group1[invalid_group],
          " vs ",
          explicit$group2[invalid_group]
        )

        stop(
          "Samplesheet contains differential ssGSEA contrasts with groups ",
          "that are not present among matched samples: ",
          paste(invalid_pairs, collapse = "; ")
        )
      }

      same_group <- explicit$group1 == explicit$group2

      if (any(same_group)) {
        invalid_pairs <- paste0(
          explicit$group1[same_group],
          " vs ",
          explicit$group2[same_group]
        )

        stop(
          "Differential ssGSEA contrasts must compare different groups: ",
          paste(invalid_pairs, collapse = "; ")
        )
      }

      explicit$source <- "samplesheet_group1_group2"

      message(
        "[INFO] Differential ssGSEA contrasts from samplesheet: ",
        nrow(explicit)
      )

      return(explicit)
    }
  }

  pair_matrix <- t(combn(groups, 2L))

  automatic <- data.frame(
    group1 = pair_matrix[, 1],
    group2 = pair_matrix[, 2],
    source = "automatic_all_pairwise",
    stringsAsFactors = FALSE
  )

  message(
    "[INFO] No explicit group1/group2 contrasts were supplied; ",
    "using all pairwise group comparisons: ",
    nrow(automatic)
  )

  automatic
}

safe_path_component <- function(x) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  value <- gsub("^_+|_+$", "", value)

  if (!nzchar(value)) {
    value <- "group"
  }

  value
}

clean_matrix <- function(mat) {
  keep <- rowSums(is.finite(mat)) >= 2
  mat <- mat[keep, , drop = FALSE]

  for (i in seq_len(nrow(mat))) {
    bad <- !is.finite(mat[i, ])
    if (any(bad)) {
      mat[i, bad] <- mean(mat[i, !bad], na.rm = TRUE)
    }
  }

  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[is.finite(vars) & vars > 0, , drop = FALSE]

  mat
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

save_plot <- function(plot, filename_base, width = 10, height = 8) {
  png_file <- paste0(filename_base, ".png")
  pdf_file <- paste0(filename_base, ".pdf")

  ggsave(png_file, plot = plot, width = width, height = height, dpi = 180)
  ggsave(pdf_file, plot = plot, width = width, height = height)

  message("[OK] Saved: ", png_file)
  message("[OK] Saved: ", pdf_file)
}

plot_heatmap <- function(mat, meta, outfile_base, title) {
  z <- row_zscore(mat)

  sample_order <- meta$sample[order(meta$group, meta$sample)]
  z <- z[, sample_order, drop = FALSE]

  term_order <- rownames(z)

  df <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
  colnames(df) <- c("Term", "Sample", "Z_score")

  df$Sample <- factor(df$Sample, levels = sample_order)
  df$Term <- factor(df$Term, levels = rev(term_order))

  p <- ggplot(df, aes(x = Sample, y = Term, fill = Z_score)) +
    geom_tile() +
    theme_bw() +
    labs(
      title = title,
      x = "Sample",
      y = "GO term",
      fill = "Row z-score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(hjust = 0.5)
    )

  h <- max(7, 0.18 * nrow(z) + 3)
  save_plot(p, outfile_base, width = 11, height = h)
}

plot_sample_correlation <- function(mat, meta, outfile_base) {
  sample_order <- meta$sample[order(meta$group, meta$sample)]
  mat <- mat[, sample_order, drop = FALSE]

  cmat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
  df <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
  colnames(df) <- c("Sample1", "Sample2", "Correlation")

  df$Sample1 <- factor(df$Sample1, levels = sample_order)
  df$Sample2 <- factor(df$Sample2, levels = rev(sample_order))

  p <- ggplot(df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile() +
    theme_bw() +
    labs(
      title = "Sample correlation based on ssGSEA scores",
      x = "",
      y = "",
      fill = "Spearman r"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  save_plot(p, outfile_base, width = 9, height = 8)
}

plot_pca <- function(mat, meta, outfile_base, pca_top_var = 500) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(pca_top_var, length(vars)))]

  pca_mat <- mat[selected, meta$sample, drop = FALSE]
  pca <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)

  var_exp <- summary(pca)$importance[2, ] * 100

  df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  df <- merge(df, meta, by = "sample", all.x = TRUE)

  p <- ggplot(df, aes(x = PC1, y = PC2, shape = group)) +
    geom_point(size = 4) +
    geom_text(aes(label = sample), vjust = -0.8, size = 3) +
    theme_bw() +
    labs(
      title = "PCA based on ssGSEA scores",
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      shape = "Group"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  save_plot(p, outfile_base, width = 9, height = 7)
}

run_differential <- function(
  mat,
  meta,
  desc,
  contrasts,
  outdir,
  top_diff = 20
) {
  combined_file <- file.path(
    outdir,
    "ssGSEA_differential_scores.tsv"
  )

  skipped_file <- file.path(
    outdir,
    "ssGSEA_differential_skipped.txt"
  )

  comparison_file <- file.path(
    outdir,
    "ssGSEA_differential_comparisons.tsv"
  )

  unlink(
    c(
      combined_file,
      skipped_file,
      comparison_file
    ),
    force = TRUE
  )

  if (nrow(contrasts) == 0L) {
    msg <- paste0(
      "Differential ssGSEA skipped because fewer than two groups ",
      "were available after matching the samplesheet to the score matrix."
    )

    writeLines(msg, skipped_file)
    message("[INFO] ", msg)
    return(invisible(NULL))
  }

  write.table(
    contrasts,
    comparison_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("[OK] Saved: ", comparison_file)

  differential_root <- file.path(
    outdir,
    "differential"
  )

  dir.create(
    differential_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  all_results <- vector(
    "list",
    nrow(contrasts)
  )

  for (comparison_index in seq_len(nrow(contrasts))) {
    g1 <- contrasts$group1[comparison_index]
    g2 <- contrasts$group2[comparison_index]
    comparison_source <- contrasts$source[comparison_index]
    comparison_label <- paste0(g1, "_vs_", g2)

    x1_idx <- meta$sample[as.character(meta$group) == g1]
    x2_idx <- meta$sample[as.character(meta$group) == g2]

    if (length(x1_idx) == 0L || length(x2_idx) == 0L) {
      stop(
        "Differential ssGSEA contrast has no matched samples: ",
        g1,
        " vs ",
        g2
      )
    }

    comparison_dir <- file.path(
      differential_root,
      paste0(
        sprintf("%02d", comparison_index),
        "_",
        safe_path_component(g1),
        "_vs_",
        safe_path_component(g2)
      )
    )

    dir.create(
      comparison_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    message(
      "[DIFFERENTIAL ssGSEA] [",
      comparison_index,
      "/",
      nrow(contrasts),
      "] ",
      g1,
      " vs ",
      g2,
      " (",
      length(x1_idx),
      " vs ",
      length(x2_idx),
      " samples)"
    )

    res_list <- lapply(rownames(mat), function(term) {
      x1 <- as.numeric(mat[term, x1_idx])
      x2 <- as.numeric(mat[term, x2_idx])

      p_t <- tryCatch(
        t.test(x1, x2)$p.value,
        error = function(e) NA_real_
      )

      p_w <- tryCatch(
        wilcox.test(x1, x2, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      data.frame(
        comparison_index = comparison_index,
        comparison = comparison_label,
        comparison_source = comparison_source,
        group1 = g1,
        group2 = g2,
        n_group1 = length(x1),
        n_group2 = length(x2),
        term = term,
        mean_group1 = mean(x1, na.rm = TRUE),
        mean_group2 = mean(x2, na.rm = TRUE),
        diff_group1_minus_group2 = (
          mean(x1, na.rm = TRUE) -
            mean(x2, na.rm = TRUE)
        ),
        pvalue_ttest = p_t,
        pvalue_wilcox = p_w,
        stringsAsFactors = FALSE
      )
    })

    res <- do.call(rbind, res_list)

    # Preserve the historical MTD behavior: BH correction is performed
    # independently within each biological contrast across all tested terms.
    res$padj_ttest <- p.adjust(
      res$pvalue_ttest,
      method = "BH"
    )

    res$padj_wilcox <- p.adjust(
      res$pvalue_wilcox,
      method = "BH"
    )

    res <- merge(
      res,
      desc,
      by = "term",
      all.x = TRUE,
      sort = FALSE
    )

    res <- res[
      order(
        res$padj_ttest,
        -abs(res$diff_group1_minus_group2),
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    comparison_tsv <- file.path(
      comparison_dir,
      "ssGSEA_differential_scores.tsv"
    )

    write.table(
      res,
      comparison_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    message("[OK] Saved: ", comparison_tsv)

    all_results[[comparison_index]] <- res

    top_terms <- res$term[
      is.finite(res$padj_ttest)
    ]

    top_terms <- head(
      top_terms,
      min(top_diff, length(top_terms))
    )

    if (length(top_terms) == 0L) {
      writeLines(
        paste0(
          "No valid differential ssGSEA terms found for ",
          g1,
          " vs ",
          g2,
          "."
        ),
        file.path(
          comparison_dir,
          "ssGSEA_differential_no_valid_terms.txt"
        )
      )

      next
    }

    pair_samples <- meta$sample[
      as.character(meta$group) %in% c(g1, g2)
    ]

    pair_meta <- meta[
      match(pair_samples, meta$sample),
      ,
      drop = FALSE
    ]

    pair_meta$group <- factor(
      as.character(pair_meta$group),
      levels = c(g1, g2)
    )

    long_list <- lapply(top_terms, function(term) {
      data.frame(
        term = term,
        sample = pair_meta$sample,
        group = pair_meta$group,
        score = as.numeric(
          mat[term, pair_meta$sample]
        ),
        stringsAsFactors = FALSE
      )
    })

    df <- do.call(rbind, long_list)

    df$term <- factor(
      df$term,
      levels = rev(top_terms)
    )

    df$group <- factor(
      df$group,
      levels = c(g1, g2)
    )

    p <- ggplot(df, aes(x = group, y = score)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
      facet_wrap(~ term, scales = "free_y", ncol = 4) +
      theme_bw() +
      labs(
        title = paste0(
          "Top differential ssGSEA GO scores: ",
          g1,
          " vs ",
          g2
        ),
        x = "Group",
        y = "ssGSEA score"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        ),
        strip.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)
      )

    h <- max(
      7,
      ceiling(length(top_terms) / 4) * 2.5
    )

    save_plot(
      p,
      file.path(
        comparison_dir,
        "ssGSEA_top_differential_boxplots"
      ),
      width = 13,
      height = h
    )

    # Backward-compatible top-level boxplot names for a run that
    # contains exactly one differential contrast.
    if (nrow(contrasts) == 1L) {
      save_plot(
        p,
        file.path(
          outdir,
          "ssGSEA_top_differential_boxplots"
        ),
        width = 13,
        height = h
      )
    }
  }

  combined <- do.call(
    rbind,
    all_results
  )

  write.table(
    combined,
    combined_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message(
    "[OK] Saved combined differential ssGSEA table: ",
    combined_file
  )

  message(
    "[OK] Differential ssGSEA comparisons completed: ",
    nrow(contrasts)
  )

  invisible(combined)
}

gct <- read_gct_scores(scores_file)
mat <- clean_matrix(gct$mat)

meta <- read_samplesheet(samplesheet_file, colnames(mat))
mat <- mat[, meta$sample, drop = FALSE]

contrasts <- read_contrasts(
  samplesheet_file,
  meta
)

desc <- gct$desc[gct$desc$term %in% rownames(mat), , drop = FALSE]

summary_file <- file.path(outdir, "ssGSEA_plot_summary.txt")
summary_lines <- c(
  paste0("scores_file\t", scores_file),
  paste0("samplesheet_file\t", samplesheet_file),
  paste0("n_terms_after_filter\t", nrow(mat)),
  paste0("n_samples\t", ncol(mat)),
  paste0("samples\t", paste(colnames(mat), collapse = ",")),
  paste0("groups\t", paste(levels(meta$group), collapse = ",")),
  paste0("n_differential_contrasts\t", nrow(contrasts)),
  paste0(
    "differential_contrasts\t",
    if (nrow(contrasts) == 0L) {
      "none"
    } else {
      paste(
        paste0(contrasts$group1, "_vs_", contrasts$group2),
        collapse = ","
      )
    }
  ),
  paste0("top_var\t", top_var),
  paste0("top_diff\t", top_diff),
  paste0("pca_top_var\t", pca_top_var)
)
writeLines(summary_lines, summary_file)
message("[OK] Saved: ", summary_file)

vars <- apply(mat, 1, var, na.rm = TRUE)
top_var_terms <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_var, length(vars)))]

plot_heatmap(
  mat[top_var_terms, , drop = FALSE],
  meta,
  file.path(outdir, "ssGSEA_top_variable_heatmap"),
  paste0("Top ", length(top_var_terms), " most variable ssGSEA GO scores")
)

plot_sample_correlation(
  mat,
  meta,
  file.path(outdir, "ssGSEA_sample_correlation_heatmap")
)

plot_pca(
  mat,
  meta,
  file.path(outdir, "ssGSEA_PCA_samples"),
  pca_top_var = pca_top_var
)

run_differential(
  mat,
  meta,
  desc,
  contrasts,
  outdir,
  top_diff = top_diff
)

message("[OK] ssGSEA plotting step finished.")
