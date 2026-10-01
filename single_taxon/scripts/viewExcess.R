# viewExcess.R
#
# Positional substitution patterns for a single taxon (Lemmon-Kishi et al. 2026, Fig. 2).
#
# Input:  the --pos-output CSV from calc_subs.py.
# Output: a four-panel figure of the first/last N read positions (5' and 3' facets):
#   A  observed rate (πQ), all 12 substitution types
#   B  observed rate (πQ), excluded pairs (default: deamination, C<>T and A<>G) removed
#   C  observed-to-expected ratio, all 12 substitution types
#   D  observed-to-expected ratio, excluded pairs removed
# Panels C and D need expected rates, i.e. calc_subs.py run with --gtr-params or
# --unrest-params. For counts-only input (no model), only panels A and B are drawn.
#
# Each symmetric pair shares a color; the solid line is one direction and the dashed line
# its reverse (e.g. solid C>T vs dashed T>C). Under symmetry, matched solid and dashed
# lines should overlap in A/B, and all lines should overlap in C/D.
#
# Command line (saves the figure; run with --help for all options):
#   Rscript viewExcess.R ex_gtr_pos.csv                  # -> ex_gtr_figure.png
#   Rscript viewExcess.R ex_gtr_pos.csv -o fig2.pdf      # format from extension
#
# Interactive (R / RStudio; shows the figure, saves only if output is given):
#   source("viewExcess.R")
#   view_excess("ex_gtr_pos.csv")
#   view_excess("ex_gtr_pos.csv", output = "fig2.pdf")
#   p <- view_excess("ex_gtr_pos.csv", show = FALSE)     # ggplot object to modify

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# =============================================================================
# SUBSTITUTION STYLES
# =============================================================================

# Okabe-Ito (colorblind-safe) color per symmetric pair
pair_colors <- c(
  "A<>C" = "#D55E00",  # vermillion
  "A<>G" = "#0072B2",  # blue
  "A<>T" = "#009E73",  # bluish green
  "C<>G" = "#CC79A7",  # reddish purple
  "C<>T" = "#E69F00",  # orange
  "G<>T" = "#56B4E9"   # sky blue
)

# Both directions of each pair, in legend order: forward (solid) then reverse (dashed)
sub_order <- unlist(lapply(names(pair_colors), function(p) {
  b <- strsplit(p, "<>")[[1]]
  c(paste0(b[1], ">", b[2]), paste0(b[2], ">", b[1]))
}))

sub_pair      <- setNames(rep(names(pair_colors), each = 2), sub_order)
sub_colors    <- setNames(rep(pair_colors, each = 2), sub_order)
sub_linetypes <- setNames(rep(c("solid", "dashed"), length(pair_colors)), sub_order)
sub_shapes    <- setNames(rep(c(16, 17), length(pair_colors)), sub_order)

# Pairs removed from panels B and D by default (deamination: C>T and its complement G>A)
default_exclude <- c("C<>T", "A<>G")

required_columns <- c("position", "position_type", "substitution", "pi_q",
                      "expected_rate", "obs_exp_ratio")

# calc_subs.py writes expected_rate = 0 for every row when run without a model
is_counts_only <- function(data) {
  all(data$expected_rate == 0)
}

# =============================================================================
# PLOTTING
# =============================================================================

theme_clean <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks       = element_line(color = "grey40", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    legend.position  = "bottom",
    legend.key       = element_blank(),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )

# One panel: same mapping and scales for every panel so the legends merge
make_panel <- function(df, yvar, title, ylab) {
  ggplot(df, aes(x = position, y = .data[[yvar]],
                 color = substitution,
                 linetype = substitution,
                 shape = substitution,
                 group = substitution)) +
    geom_line(linewidth = 0.75, alpha = 0.9) +
    geom_point(size = 2, stroke = 0.4) +
    facet_wrap(~position_type, scales = "free_x") +
    scale_color_manual(values = sub_colors) +
    scale_linetype_manual(values = sub_linetypes) +
    scale_shape_manual(values = sub_shapes) +
    theme_clean +
    labs(title = title, x = "Position in Read", y = ylab)
}

# Normalize pair names so "C<>T", "T<>C", "CT" and "ct" all become "C<>T".
# "none" (or an empty vector) excludes nothing.
normalize_pairs <- function(pairs) {
  pairs <- toupper(trimws(pairs))
  pairs <- pairs[pairs != "" & pairs != "NONE"]
  normalized <- vapply(pairs, function(p) {
    bases <- sort(strsplit(gsub("[^ACGT]", "", p), "")[[1]])
    paste(bases, collapse = "<>")
  }, character(1), USE.NAMES = FALSE)
  invalid <- pairs[!normalized %in% names(pair_colors)]
  if (length(invalid) > 0) {
    stop("Invalid substitution pair(s): ", paste(invalid, collapse = ", "),
         ". Use two different bases, e.g. CT or C<>T.", call. = FALSE)
  }
  unique(normalized)
}

# Build the figure from a calc_subs.py --pos-output CSV.
#   input   path to the CSV
#   output  file to save to (format from extension: .png, .pdf, .svg, ...); NULL = don't save
#   exclude pairs removed from panels B and D; "none" or character(0) removes nothing
#   height  NULL = 10 inches for four panels, 5 for counts-only input (panels A and B)
#   show    draw the figure on the current graphics device
# Returns the patchwork/ggplot object invisibly.
view_excess <- function(input, output = NULL, exclude = default_exclude,
                        width = 14, height = NULL, dpi = 300, show = interactive()) {
  if (!file.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }
  data <- read.csv(input)
  missing <- setdiff(required_columns, names(data))
  if (length(missing) > 0) {
    stop("Input is missing column(s): ", paste(missing, collapse = ", "),
         ". Expected the --pos-output CSV from calc_subs.py.", call. = FALSE)
  }
  exclude <- normalize_pairs(exclude)

  counts_only <- is_counts_only(data)
  if (counts_only) {
    warning("No expected rates in ", input, " (calc_subs.py was run without ",
            "--gtr-params or --unrest-params); plotting observed rates (panels A, B) only.",
            call. = FALSE)
  }
  if (is.null(height)) height <- if (counts_only) 5 else 10

  data <- data %>%
    mutate(
      position_type = factor(position_type, levels = c("5prime", "3prime")),
      pair          = sub_pair[substitution],
      substitution  = factor(substitution, levels = sub_order)
    )

  data_filtered <- data %>%
    filter(!pair %in% exclude)

  pi_q_label  <- expression(pi*" Q")
  ratio_label <- "Observed/Expected Ratio"

  # Panel A carries the single legend for the whole figure
  plot_a <- make_panel(data, "pi_q", "A", pi_q_label) +
    labs(color = "Substitution", linetype = "Substitution", shape = "Substitution") +
    guides(
      color    = guide_legend(nrow = 2, override.aes = list(size = 2.4, linewidth = 0.9)),
      linetype = guide_legend(nrow = 2),
      shape    = guide_legend(nrow = 2)
    )

  no_legend <- guides(color = "none", linetype = "none", shape = "none")

  plot_b <- make_panel(data_filtered, "pi_q", "B", pi_q_label) + no_legend

  if (counts_only) {
    layout <- plot_a | plot_b
  } else {
    plot_c <- make_panel(data,          "obs_exp_ratio", "C", ratio_label) + no_legend
    plot_d <- make_panel(data_filtered, "obs_exp_ratio", "D", ratio_label) + no_legend
    layout <- (plot_a | plot_b) / (plot_c | plot_d)
  }

  combined_plot <- layout +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (show) print(combined_plot)

  if (!is.null(output)) {
    ggsave(output, plot = combined_plot, width = width, height = height, dpi = dpi)
    message("Wrote figure to: ", output)
  }

  invisible(combined_plot)
}

# =============================================================================
# COMMAND LINE
# =============================================================================

usage <- "Usage: Rscript viewExcess.R INPUT [options]

Plot positional substitution patterns (paper Fig. 2) from the --pos-output CSV
of calc_subs.py, and save the figure. Observed-to-expected panels (C, D) need
calc_subs.py run with a model; counts-only input gives panels A and B only.

Arguments:
  INPUT                 --pos-output CSV from calc_subs.py

Options:
  -o, --output FILE     Output file; format from extension (.png, .pdf, .svg, ...)
                        (default: INPUT with _pos.csv replaced by _figure.png)
  -e, --exclude PAIRS   Comma-separated pairs removed from panels B and D, written
                        as CT or C<>T; 'none' removes nothing (default: CT,AG)
  --width NUM           Figure width in inches (default: 14)
  --height NUM          Figure height in inches (default: 10, or 5 for counts-only)
  --dpi NUM             Resolution for raster formats (default: 300)
  -h, --help            Show this message

To view the figure without saving, use it interactively:
  source(\"viewExcess.R\"); view_excess(\"INPUT\")
"

# ex_gtr_pos.csv -> ex_gtr_figure.png, next to the input
default_output <- function(input) {
  stem <- sub("_pos$", "", tools::file_path_sans_ext(basename(input)))
  file.path(dirname(input), paste0(stem, "_figure.png"))
}

parse_cli <- function(args) {
  opts <- list(input = NULL, output = NULL, exclude = default_exclude,
               width = 14, height = NULL, dpi = 300)

  # Accept both "--opt value" and "--opt=value"
  args <- unlist(lapply(args, function(a) {
    if (startsWith(a, "--") && grepl("=", a)) strsplit(sub("=", "\001", a), "\001")[[1]] else a
  }))

  value_of <- function(i, flag) {
    if (i + 1 > length(args)) stop("Missing value for ", flag, call. = FALSE)
    args[i + 1]
  }
  number_of <- function(i, flag) {
    x <- suppressWarnings(as.numeric(value_of(i, flag)))
    if (is.na(x) || x <= 0) stop(flag, " must be a positive number", call. = FALSE)
    x
  }

  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      cat(usage)
      quit(status = 0)
    } else if (a %in% c("-o", "--output")) {
      opts$output <- value_of(i, a); i <- i + 2
    } else if (a %in% c("-e", "--exclude")) {
      opts$exclude <- strsplit(value_of(i, a), ",")[[1]]; i <- i + 2
    } else if (a == "--width") {
      opts$width <- number_of(i, a); i <- i + 2
    } else if (a == "--height") {
      opts$height <- number_of(i, a); i <- i + 2
    } else if (a == "--dpi") {
      opts$dpi <- number_of(i, a); i <- i + 2
    } else if (startsWith(a, "-")) {
      stop("Unknown option: ", a, call. = FALSE)
    } else if (is.null(opts$input)) {
      opts$input <- a; i <- i + 1
    } else {
      stop("Unexpected argument: ", a, call. = FALSE)
    }
  }

  if (is.null(opts$input)) {
    cat(usage)
    quit(status = 1)
  }
  if (is.null(opts$output)) opts$output <- default_output(opts$input)
  opts
}

# Run only via Rscript, not when source()d
if (!interactive() && sys.nframe() == 0) {
  tryCatch({
    opts <- parse_cli(commandArgs(trailingOnly = TRUE))
    view_excess(opts$input, output = opts$output, exclude = opts$exclude,
                width = opts$width, height = opts$height, dpi = opts$dpi,
                show = FALSE)
  }, error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1)
  })
}
