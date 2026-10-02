# symmetry_tests.R
#
# Ancient vs modern substitution symmetry across taxa (Lemmon-Kishi et al. 2026,
# Sections 4.2.1 and 4.2.3; Tables 2-3, Fig. 3).
#
# Input: a directory of per-read CSVs from analyze_substitutions_patterns_aggregate.py -r
# (one per library), and a list of ancient tax IDs. Reads whose tax ID is not ancient are
# modern. Optionally, contaminant taxa found in control samples are removed first.
#
# For two sets of substitution pairs:
#   Complementary (symmetric) pairs, e.g. C>T vs T>C: ~1:1 under time-reversible evolution
#   Watson-Crick complement pairs,   e.g. C>T vs G>A: ~1:1 in double-stranded libraries
# it computes:
#   1. Ratio plots of each pair in ancient and modern reads (Fig. 3)
#   2. Chi-squared tests (ancient/modern x pair) with Cramér's V (Table 2)
#   3. Permutation tests of the ancient - modern ratio difference, permuting
#      ancient/modern labels across reads; p = (b + 1) / (m + 1)
#
# Command line (run with --help for all options):
#   Rscript symmetry_tests.R READS_DIR --ancient-taxa ancient_taxa_id.csv \
#     --control-taxa control_taxa.csv --modern-taxa nonControl_modern_taxa.csv -o results/
#
# Interactive (R / RStudio):
#   source("symmetry_tests.R")
#   res <- run_symmetry_tests("reads/", "ancient_taxa_id.csv",
#                             control_taxa = "control_taxa.csv", outdir = "results/")
#   res$ratio_plots$combined

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(parallel)
})

# =============================================================================
# SUBSTITUTION PAIRS
# =============================================================================

sub_cols <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
              "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")

ratio_sets <- list(
  complementary = list(
    label = "Complementary",
    # Ratio definitions: numerator / denominator
    defs = list(
      "C>A / A>C" = c("C>A", "A>C"),
      "G>A / A>G" = c("G>A", "A>G"),
      "A>T / T>A" = c("A>T", "T>A"),
      "C>G / G>C" = c("C>G", "G>C"),
      "C>T / T>C" = c("C>T", "T>C"),
      "G>T / T>G" = c("G>T", "T>G")
    ),
    # Chi-squared pair definitions: col1, col2 in contingency table
    pair_defs = list(
      AC_CA = c("A>C", "C>A"),
      AG_GA = c("A>G", "G>A"),
      AT_TA = c("A>T", "T>A"),
      CG_GC = c("C>G", "G>C"),
      CT_TC = c("C>T", "T>C"),
      GT_TG = c("G>T", "T>G")
    ),
    colors = c("C>A / A>C" = "#D55E00", "G>A / A>G" = "#0072B2",
               "A>T / T>A" = "#009E73", "C>G / G>C" = "#CC79A7",
               "C>T / T>C" = "#E69F00", "G>T / T>G" = "#56B4E9")
  ),
  watson_crick = list(
    label = "Watson-Crick",
    defs = list(
      "T>G / A>C" = c("T>G", "A>C"),
      "A>G / T>C" = c("A>G", "T>C"),
      "A>T / T>A" = c("A>T", "T>A"),
      "C>G / G>C" = c("C>G", "G>C"),
      "C>T / G>A" = c("C>T", "G>A"),
      "G>T / C>A" = c("G>T", "C>A")
    ),
    pair_defs = list(
      AC_TG = c("A>C", "T>G"),
      AG_TC = c("A>G", "T>C"),
      AT_TA = c("A>T", "T>A"),
      CG_GC = c("C>G", "G>C"),
      CT_GA = c("C>T", "G>A"),
      GT_CA = c("G>T", "C>A")
    ),
    colors = c("T>G / A>C" = "#56B4E9", "A>G / T>C" = "#0072B2",
               "A>T / T>A" = "#009E73", "C>G / G>C" = "#CC79A7",
               "C>T / G>A" = "#E69F00", "G>T / C>A" = "#D55E00")
  )
)

# =============================================================================
# DATA LOADING
# =============================================================================

read_taxa <- function(path) {
  if (!file.exists(path)) stop("Taxa file not found: ", path, call. = FALSE)
  taxa <- read.csv(path, check.names = FALSE)
  if (!"taxa_id" %in% names(taxa)) {
    stop("Taxa file needs a 'taxa_id' column: ", path, call. = FALSE)
  }
  as.character(taxa$taxa_id)
}

load_reads <- function(reads_dir, pattern) {
  files <- list.files(reads_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files matching '", pattern, "' in ", reads_dir, call. = FALSE)
  }
  names(files) <- basename(files)
  cols <- c("read_name", "tax_id", sub_cols)
  reads <- lapply(files, function(f) {
    missing <- setdiff(cols, names(fread(f, nrows = 0)))
    if (length(missing) > 0) {
      stop(f, " is missing column(s): ", paste(missing, collapse = ", "),
           ". Expected per-read CSVs from analyze_substitutions_patterns_aggregate.py -r.",
           call. = FALSE)
    }
    fread(f, select = cols, colClasses = list(character = "tax_id"))
  })
  as.data.frame(bind_rows(reads, .id = "filename"))
}

# Compare the modern taxa in the data with the expected list. With controls removed,
# modern should be exactly the listed taxa, and any other taxon is warned about. With
# controls kept, the extra taxa are expected (they are the control taxa) and only noted.
check_modern_taxa <- function(data, modern_taxa, remove_controls) {
  observed <- unique(data$tax_id[data$taxa_category == "modern"])
  unexpected <- setdiff(observed, modern_taxa)
  absent <- setdiff(modern_taxa, observed)
  if (length(unexpected) > 0) {
    if (remove_controls) {
      warning("Modern reads include ", length(unexpected), " taxa not in the modern taxa list: ",
              paste(head(sort(unexpected), 10), collapse = ", "),
              if (length(unexpected) > 10) ", ..." else "", call. = FALSE)
    } else {
      cat("  Modern pool has", length(unexpected),
          "taxa beyond the modern taxa list (controls kept)\n")
    }
  }
  if (length(absent) > 0) {
    warning("Listed modern taxa with no reads: ", paste(sort(absent), collapse = ", "),
            call. = FALSE)
  }
  if (length(absent) == 0 && (length(unexpected) == 0 || !remove_controls)) {
    cat("  All listed modern taxa are present\n")
  }
}

# =============================================================================
# PART 1: RATIO PLOTS
# =============================================================================

theme_ratio <- theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )

# Pooled ratio of each pair in the set, one row per category and ratio
compute_ratios_long <- function(summary_df, rset) {
  ratios <- summary_df
  for (rname in names(rset$defs)) {
    cols_pair <- rset$defs[[rname]]
    ratios[[rname]] <- ifelse(ratios[[cols_pair[2]]] > 0,
                              ratios[[cols_pair[1]]] / ratios[[cols_pair[2]]], NA)
  }
  ratios %>%
    select(taxa_category, all_of(names(rset$defs))) %>%
    pivot_longer(cols = all_of(names(rset$defs)),
                 names_to = "ratio_type", values_to = "ratio_value") %>%
    filter(!is.na(ratio_value) & is.finite(ratio_value) & ratio_value > 0)
}

# y_range zooms the view (coord_cartesian) rather than dropping points outside it
make_ratio_plot <- function(ratios_long, rset, y_range) {
  ggplot(ratios_long, aes(x = taxa_category, y = ratio_value)) +
    geom_point(aes(color = ratio_type, shape = taxa_category),
               alpha = 1, size = 4, stroke = 1.5,
               # Fixed seed so the jitter (and the saved figure) is the same every run
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8,
                                               seed = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_color_manual(values = rset$colors, name = "Substitution Ratio") +
    coord_cartesian(ylim = y_range) +
    labs(x = "Taxa Category", y = "Substitution Ratio", shape = "Taxa Category") +
    theme_ratio
}

# =============================================================================
# PART 2: CHI-SQUARED TESTS
# =============================================================================

# 2x2 contingency tables (ancient/modern x pair) from summarized counts
create_pair_tables <- function(summary_df, pair_defs) {
  lapply(pair_defs, function(cols) {
    matrix(c(summary_df[summary_df$taxa_category == "ancient", cols[1]],
             summary_df[summary_df$taxa_category == "ancient", cols[2]],
             summary_df[summary_df$taxa_category == "modern", cols[1]],
             summary_df[summary_df$taxa_category == "modern", cols[2]]),
           nrow = 2, byrow = TRUE,
           dimnames = list(c("ancient", "modern"), cols))
  })
}

analyze_pairs_with_effect <- function(pair_tables, analysis_type) {
  bind_rows(lapply(names(pair_tables), function(pair_name) {
    tbl <- pair_tables[[pair_name]]
    test <- chisq.test(tbl, correct = TRUE)
    ancient_ratio <- tbl[1, 1] / tbl[1, 2]
    modern_ratio  <- tbl[2, 1] / tbl[2, 2]
    tibble(
      analysis_type = analysis_type,
      pair = pair_name,
      chi_squared = test$statistic,
      p_value = test$p.value,
      ancient_ratio = ancient_ratio,
      modern_ratio = modern_ratio,
      ratio_difference = ancient_ratio - modern_ratio,
      percent_difference = ((ancient_ratio - modern_ratio) / modern_ratio) * 100,
      ancient_deviation_from_symmetry = abs(ancient_ratio - 1),
      modern_deviation_from_symmetry = abs(modern_ratio - 1),
      cramer_v = sqrt(test$statistic / sum(tbl))
    )
  }))
}

write_chisq_details <- function(pair_tables, label) {
  cat("\n\n####################", label, "####################\n")
  for (pair_name in names(pair_tables)) {
    cat("\n==========", pair_name, "==========\n")
    test <- chisq.test(pair_tables[[pair_name]])
    print(test)
    cat("\nObserved:\n")
    print(test$observed)
    cat("\nExpected:\n")
    print(test$expected)
    cat("\nStandardized Residuals:\n")
    print(test$stdres)
    n <- sum(test$observed)
    k <- min(nrow(test$observed), ncol(test$observed))
    v <- sqrt(test$statistic / (n * (k - 1)))
    cat("\nCramér's V:", round(v, 4), "\n")
  }
}

# =============================================================================
# PART 3: PERMUTATION TESTS
# =============================================================================

calc_ratio_diffs_fast <- function(is_ancient, sub_matrix, total_sums, rdefs_idx) {
  anc_sums <- colSums(sub_matrix[is_ancient, , drop = FALSE])
  mod_sums <- total_sums - anc_sums
  sapply(rdefs_idx, function(idx) {
    (anc_sums[idx[1]] / anc_sums[idx[2]]) - (mod_sums[idx[1]] / mod_sums[idx[2]])
  })
}

# Permute ancient/modern labels across reads; one pass covers all ratios.
# Results depend on seed and ncores (each worker gets seed + worker index).
run_permutation_test <- function(df, rdefs, n_perm, seed, ncores, min_sub_count) {
  original_categories <- df$taxa_category
  sub_matrix <- as.matrix(df[, sub_cols])
  n_reads <- nrow(sub_matrix)

  for (cat_name in c("ancient", "modern")) {
    cat_sums <- colSums(sub_matrix[original_categories == cat_name, , drop = FALSE])
    if (any(cat_sums < min_sub_count)) {
      warning("Category '", cat_name, "' has substitution counts below ", min_sub_count,
              call. = FALSE)
    }
  }

  total_sums <- colSums(sub_matrix)
  rdefs_idx <- lapply(rdefs, function(cols) match(cols, colnames(sub_matrix)))

  is_ancient <- original_categories == "ancient"
  n_ancient <- sum(is_ancient)

  obs_diffs <- calc_ratio_diffs_fast(is_ancient, sub_matrix, total_sums, rdefs_idx)
  names(obs_diffs) <- names(rdefs)

  chunk_size <- ceiling(n_perm / ncores)
  worker_seeds <- seed + seq_len(ncores)

  # mclapply forks the process; workers share sub_matrix copy-on-write.
  # Forking is not available on Windows, where only ncores = 1 works.
  perm_results <- mclapply(worker_seeds, function(s) {
    set.seed(s)
    local_diffs <- matrix(NA_real_, nrow = chunk_size, ncol = length(rdefs))
    for (i in seq_len(chunk_size)) {
      perm_flag <- logical(n_reads)
      perm_flag[sample.int(n_reads, n_ancient)] <- TRUE
      local_diffs[i, ] <- calc_ratio_diffs_fast(perm_flag, sub_matrix, total_sums, rdefs_idx)
    }
    local_diffs
  }, mc.cores = ncores)

  perm_diffs <- do.call(rbind, perm_results)
  perm_diffs <- perm_diffs[seq_len(n_perm), , drop = FALSE]
  colnames(perm_diffs) <- names(rdefs)

  # Two-sided p = (b + 1) / (m + 1), counting the observed labelling as one permutation
  # (Phipson & Smyth 2010); the smallest possible value is 1 / (n_perm + 1).
  p_values <- sapply(seq_along(rdefs), function(j) {
    (sum(abs(perm_diffs[, j]) >= abs(obs_diffs[j])) + 1) / (n_perm + 1)
  })
  names(p_values) <- names(rdefs)

  list(
    observed_diffs = obs_diffs,
    perm_diffs     = perm_diffs,
    p_values       = p_values,
    n_ancient      = n_ancient,
    n_modern       = n_reads - n_ancient
  )
}

# Load cached permutations if they were made with the same settings; otherwise run them.
get_permutations <- function(data, rdefs, settings, perm_cache) {
  if (!is.null(perm_cache) && file.exists(perm_cache)) {
    cached <- readRDS(perm_cache)
    if (identical(cached$settings, settings)) {
      cat("  Loaded cached permutations from", perm_cache, "\n")
      return(cached)
    }
    warning("Ignoring ", perm_cache, ": it was made with different settings or data",
            call. = FALSE)
  }
  cat("  Running", settings$n_perm, "permutations on", settings$ncores, "core(s)...\n")
  result <- run_permutation_test(data, rdefs, settings$n_perm, settings$seed,
                                 settings$ncores, settings$min_sub_count)
  result$settings <- settings
  if (!is.null(perm_cache)) {
    saveRDS(result, perm_cache)
    cat("  Saved permutations to", perm_cache, "\n")
  }
  result
}

create_perm_plot <- function(perm_result, ratio_col) {
  plot_data <- data.frame(perm_diff = perm_result$perm_diffs[, ratio_col])
  obs_diff  <- perm_result$observed_diffs[ratio_col]
  p_val     <- perm_result$p_values[ratio_col]

  ggplot(plot_data, aes(x = perm_diff)) +
    geom_histogram(bins = 50, fill = "#56B4E9", alpha = 0.7, color = "white") +
    geom_vline(xintercept = obs_diff, color = "#D55E00",
               linetype = "dashed", linewidth = 1.2) +
    labs(
      title = ratio_col,
      subtitle = sprintf("Observed diff = %.4f, p = %.4f", obs_diff, p_val),
      x = "Permuted Ratio Difference (Ancient - Modern)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40")
    )
}

# =============================================================================
# MAIN
# =============================================================================

# Run all three analyses and write results to outdir.
#   reads_dir      directory of per-read CSVs from analyze_substitutions_patterns_aggregate.py -r
#   pattern        regex selecting the per-read CSVs in reads_dir
#   ancient_taxa   CSV with a taxa_id column; all other taxa are modern
#   control_taxa   CSV of contaminant taxa to remove from the modern pool; NULL keeps them
#   modern_taxa    CSV of expected modern taxa, used only as a consistency check
#   perm_cache     .rds file to save permutations to; reused when the data and settings
#                  match, overwritten otherwise
#   ylim           c(min, max) y-axis range for the ratio plots; NULL fits the data
# Output files are prefixed controls_removed_ or controls_included_.
# Returns the summarized counts, test results and plots invisibly.
run_symmetry_tests <- function(reads_dir, ancient_taxa, control_taxa = NULL,
                               modern_taxa = NULL, outdir = ".", pattern = "_reads\\.csv$",
                               min_sub_count = 100, n_perm = 10000, ncores = 4,
                               seed = 42, perm_cache = NULL, ylim = NULL) {
  if (!is.null(ylim) && !(is.numeric(ylim) && length(ylim) == 2 && ylim[1] < ylim[2])) {
    stop("ylim must be two increasing numbers, e.g. c(0.5, 8)", call. = FALSE)
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  remove_controls <- !is.null(control_taxa)
  prefix <- if (remove_controls) "controls_removed" else "controls_included"
  out <- function(name) file.path(outdir, paste0(prefix, "_", name))
  saved <- function(path) cat("  Saved:", path, "\n")

  # ---------------------------------------------------------------------------
  # Data
  # ---------------------------------------------------------------------------
  cat("=== Loading data ===\n")
  ancient_ids <- read_taxa(ancient_taxa)
  control_ids <- if (remove_controls) read_taxa(control_taxa) else NULL
  modern_ids  <- if (!is.null(modern_taxa)) read_taxa(modern_taxa) else NULL

  data <- load_reads(reads_dir, pattern) %>%
    mutate(taxa_category = if_else(tax_id %in% ancient_ids, "ancient", "modern"))

  if (remove_controls) {
    n_before <- nrow(data)
    data <- data %>% filter(!tax_id %in% control_ids)
    cat("  Removed", n_before - nrow(data), "reads from control taxa\n")
  } else {
    cat("  Control taxa kept in the modern pool\n")
  }

  cat("  Reads:", nrow(data), "\n")
  cat("    Ancient:", sum(data$taxa_category == "ancient"), "\n")
  cat("    Modern:", sum(data$taxa_category == "modern"), "\n")
  if (!all(c("ancient", "modern") %in% data$taxa_category)) {
    stop("Need both ancient and modern reads; check the taxa lists.", call. = FALSE)
  }
  if (!is.null(modern_ids)) {
    check_modern_taxa(data, modern_ids, remove_controls)
  }

  data_summarized <- data %>%
    group_by(taxa_category) %>%
    summarise(across(all_of(sub_cols), sum), .groups = "drop") %>%
    as.data.frame()

  # ---------------------------------------------------------------------------
  # Part 1: ratio plots
  # ---------------------------------------------------------------------------
  cat("\n=== PART 1: Ratio plots ===\n")
  # A category is plotted only if every substitution type has at least min_sub_count
  data_summarized_filtered <- data_summarized %>%
    filter(if_all(all_of(sub_cols), ~ . >= min_sub_count))
  dropped <- setdiff(data_summarized$taxa_category, data_summarized_filtered$taxa_category)
  if (length(dropped) > 0) {
    warning("Not plotting ", paste(dropped, collapse = ", "),
            ": some substitution counts are below ", min_sub_count, " (see --min-sub-count)",
            call. = FALSE)
  }

  ratios_long <- lapply(ratio_sets, compute_ratios_long, summary_df = data_summarized_filtered)

  # One y-range for both sets so the panels are comparable: by default all plotted ratios
  # plus 1 (the symmetry line), with ggplot's usual padding
  y_range <- if (is.null(ylim)) {
    range(c(1, unlist(lapply(ratios_long, `[[`, "ratio_value"))))
  } else {
    ylim
  }
  cat("  y-axis range:", paste(signif(y_range, 3), collapse = " to "), "\n")

  ratio_plots <- list()
  for (set_name in names(ratio_sets)) {
    p <- make_ratio_plot(ratios_long[[set_name]], ratio_sets[[set_name]], y_range)
    path <- out(sprintf("ratios_%s.png", set_name))
    ggsave(path, p, width = 16, height = 7, dpi = 300)
    saved(path)
    ratio_plots[[set_name]] <- p
  }

  ratio_plots$combined <- ratio_plots$complementary / ratio_plots$watson_crick +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14))
  path <- out("ratios_combined.png")
  ggsave(path, ratio_plots$combined, width = 12, height = 10, dpi = 300)
  saved(path)

  # ---------------------------------------------------------------------------
  # Part 2: chi-squared tests
  # ---------------------------------------------------------------------------
  cat("\n=== PART 2: Chi-squared tests ===\n")
  details_path <- out("chisq_detailed.txt")
  chisq_all <- list()
  sink(details_path)
  tryCatch({
    for (set_name in names(ratio_sets)) {
      rset <- ratio_sets[[set_name]]
      pair_tables <- create_pair_tables(data_summarized, rset$pair_defs)
      chisq_all[[set_name]] <- analyze_pairs_with_effect(pair_tables, rset$label)
      write_chisq_details(pair_tables, rset$label)
    }
  }, finally = sink())
  chisq_results <- bind_rows(chisq_all)
  path <- out("chisq_results.csv")
  write.csv(chisq_results, path, row.names = FALSE)
  saved(path)
  saved(details_path)

  # ---------------------------------------------------------------------------
  # Part 3: permutation tests
  # ---------------------------------------------------------------------------
  cat("\n=== PART 3: Permutation tests ===\n")
  # All ratio definitions, deduplicated (A>T / T>A and C>G / G>C are in both sets)
  all_ratio_defs <- c(ratio_sets$complementary$defs, ratio_sets$watson_crick$defs)
  all_ratio_defs <- all_ratio_defs[!duplicated(names(all_ratio_defs))]

  settings <- list(prefix = prefix, n_perm = n_perm, seed = seed, ncores = ncores,
                   min_sub_count = min_sub_count,
                   n_reads = nrow(data), n_ancient = sum(data$taxa_category == "ancient"),
                   sub_totals = colSums(as.matrix(data[, sub_cols])))
  perm_result <- get_permutations(data, all_ratio_defs, settings, perm_cache)

  sub_matrix_all <- as.matrix(data[, sub_cols])
  is_ancient_all <- data$taxa_category == "ancient"
  is_modern_all  <- data$taxa_category == "modern"

  perm_tables <- list()
  perm_plots <- list()
  for (set_name in names(ratio_sets)) {
    rset <- ratio_sets[[set_name]]
    cat("\n---", rset$label, "permutation results ---\n")

    results <- bind_rows(lapply(names(rset$defs), function(ratio_name) {
      cols_pair <- rset$defs[[ratio_name]]
      cat("  ", ratio_name, ": p =", perm_result$p_values[ratio_name], "\n")
      data.frame(
        ratio_type    = ratio_name,
        n_ancient     = perm_result$n_ancient,
        n_modern      = perm_result$n_modern,
        ancient_ratio = sum(sub_matrix_all[is_ancient_all, cols_pair[1]]) /
                        sum(sub_matrix_all[is_ancient_all, cols_pair[2]]),
        modern_ratio  = sum(sub_matrix_all[is_modern_all,  cols_pair[1]]) /
                        sum(sub_matrix_all[is_modern_all,  cols_pair[2]]),
        observed_diff = perm_result$observed_diffs[[ratio_name]],
        p_value       = perm_result$p_values[[ratio_name]],
        stringsAsFactors = FALSE
      )
    })) %>%
      mutate(significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ))

    path <- out(sprintf("permutation_results_%s.csv", set_name))
    write.csv(results, path, row.names = FALSE)
    saved(path)

    plots <- lapply(names(rset$defs), function(r) create_perm_plot(perm_result, r))
    perm_plots[[set_name]] <- wrap_plots(plots, ncol = 2)
    path <- out(sprintf("permutation_distributions_%s.png", set_name))
    ggsave(path, perm_plots[[set_name]], width = 12, height = 9, dpi = 300)
    saved(path)

    perm_tables[[set_name]] <- results
  }

  cat("\n=== All analyses complete ===\n")
  invisible(list(
    data_summarized   = data_summarized,
    chisq_results     = chisq_results,
    permutation       = perm_tables,
    perm_result       = perm_result,
    ratio_plots       = ratio_plots,
    permutation_plots = perm_plots
  ))
}

# =============================================================================
# COMMAND LINE
# =============================================================================

usage <- "Usage: Rscript symmetry_tests.R READS_DIR --ancient-taxa FILE [options]

Ancient vs modern substitution symmetry: ratio plots, chi-squared tests with
Cramer's V, and label-permutation tests, for complementary and Watson-Crick pairs.

Arguments:
  READS_DIR               Directory of per-read CSVs from
                          analyze_substitutions_patterns_aggregate.py -r

Options:
  --ancient-taxa FILE     CSV with a taxa_id column of ancient taxa (required);
                          all other taxa are modern
  --control-taxa FILE     CSV of contaminant taxa to remove from the modern pool
                          (default: keep all taxa)
  --modern-taxa FILE      CSV of expected modern taxa; warns on any mismatch
  -o, --outdir DIR        Output directory (default: current directory)
  --pattern REGEX         Which files in READS_DIR to read (default: _reads\\.csv$)
  --min-sub-count NUM     Minimum count of every substitution type for a category
                          to be plotted (default: 100)
  --permutations NUM      Number of permutations (default: 10000)
  --cores NUM             Cores for permutations; 1 on Windows (default: 4)
  --seed NUM              Random seed (default: 42)
  --perm-cache FILE       .rds file to save permutations to; reused when the data
                          and settings match, overwritten otherwise
  --ylim MIN,MAX          y-axis range for the ratio plots (default: fit all
                          ratios and 1; points outside a set range are not dropped,
                          only out of view)
  -h, --help              Show this message
"

parse_cli <- function(args) {
  opts <- list(reads_dir = NULL, ancient_taxa = NULL, control_taxa = NULL,
               modern_taxa = NULL, outdir = ".", pattern = "_reads\\.csv$",
               min_sub_count = 100, n_perm = 10000, ncores = 4, seed = 42,
               perm_cache = NULL, ylim = NULL)

  # Accept both "--opt value" and "--opt=value"
  args <- unlist(lapply(args, function(a) {
    if (startsWith(a, "--") && grepl("=", a)) strsplit(sub("=", "\001", a), "\001")[[1]] else a
  }))

  value_of <- function(i, flag) {
    if (i + 1 > length(args)) stop("Missing value for ", flag, call. = FALSE)
    args[i + 1]
  }
  integer_of <- function(i, flag, min = 1) {
    x <- suppressWarnings(as.numeric(value_of(i, flag)))
    if (is.na(x) || x < min || x != round(x)) {
      stop(flag, " must be an integer >= ", min, call. = FALSE)
    }
    as.integer(x)
  }

  string_opts <- c("--ancient-taxa" = "ancient_taxa", "--control-taxa" = "control_taxa",
                   "--modern-taxa" = "modern_taxa", "-o" = "outdir", "--outdir" = "outdir",
                   "--pattern" = "pattern", "--perm-cache" = "perm_cache")
  integer_opts <- c("--min-sub-count" = "min_sub_count", "--permutations" = "n_perm",
                    "--cores" = "ncores", "--seed" = "seed")

  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      cat(usage)
      quit(status = 0)
    } else if (a %in% names(string_opts)) {
      opts[[string_opts[[a]]]] <- value_of(i, a); i <- i + 2
    } else if (a %in% names(integer_opts)) {
      opts[[integer_opts[[a]]]] <- integer_of(i, a, min = if (a == "--min-sub-count") 0 else 1)
      i <- i + 2
    } else if (a == "--ylim") {
      y <- suppressWarnings(as.numeric(strsplit(value_of(i, a), ",")[[1]]))
      if (length(y) != 2 || anyNA(y) || y[1] >= y[2]) {
        stop("--ylim must be two increasing numbers, e.g. --ylim 0.5,8", call. = FALSE)
      }
      opts$ylim <- y; i <- i + 2
    } else if (startsWith(a, "-")) {
      stop("Unknown option: ", a, call. = FALSE)
    } else if (is.null(opts$reads_dir)) {
      opts$reads_dir <- a; i <- i + 1
    } else {
      stop("Unexpected argument: ", a, call. = FALSE)
    }
  }

  if (is.null(opts$reads_dir) || is.null(opts$ancient_taxa)) {
    cat(usage)
    quit(status = 1)
  }
  if (!dir.exists(opts$reads_dir)) {
    stop("Reads directory not found: ", opts$reads_dir, call. = FALSE)
  }
  opts
}

# Run only via Rscript, not when source()d
if (!interactive() && sys.nframe() == 0) {
  tryCatch({
    opts <- parse_cli(commandArgs(trailingOnly = TRUE))
    do.call(run_symmetry_tests, opts)
  }, error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
  })
}
