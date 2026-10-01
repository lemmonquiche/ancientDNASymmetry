library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)

custom_colors <- c(
  "A>C" = "#e41a1c", "A>G" = "#377eb8", "A>T" = "#4daf4a",  # A subs: red, blue, green
  "C>A" = "#ff7f00", "C>G" = "#ffff33", "C>T" = "#a65628",  # C subs: orange, yellow, brown
  "G>A" = "#f781bf", "G>C" = "#999999", "G>T" = "#984ea3",  # G subs: pink, gray, purple
  "T>A" = "#00ffff", "T>C" = "#8dd3c7", "T>G" = "#00ff00"   # T subs: cyan, teal, lime
)

custom_colors_pair <- c(
  "A<>C" = "#e41a1c", "A<>G" = "#377eb8", "A<>T" = "#4daf4a",
  "C<>G" = "#999999", "C<>T" = "#a65628", "G<>T" = "#984ea3"
)

# custom_colors_pair_2 <- c(
# "A<>C" = "#FF6B35",  # Bright Orange
#   "A<>G" = "#004E89",  # Deep Blue
#   "A<>T" = "#00A878",  # Emerald Green
#   "C<>G" = "#FFC400",  # Golden Yellow
#   "C<>T" = "#DC143C",  # Crimson Red
#   "G<>T" = "#9D02D7"   # Vivid Purple
#   )

custom_colors_pair_2 <- c(
  "A<>C" = "#FF6B35",  # Bright Orange (kept vibrant)
  "A<>G" = "#004E89",  # Deep Blue (kept vibrant)
  "A<>T" = "#7FA99B",  # Desaturated Green (still recognizably green)
  "C<>G" = "#D4C48A",  # Desaturated Yellow (still recognizably yellow)
  "C<>T" = "#DC143C",  # Crimson Red (kept vibrant)
  "G<>T" = "#9D02D7"   # Vivid Purple (kept vibrant)
)

custom_colors_sym <- c(
  "A>C" = "#e41a1c", "A>G" = "#377eb8", "A>T" = "#4daf4a",  # A subs: red, blue, green
  "C>A" = "#e41a1c", "C>G" = "#999999", "C>T" = "#a65628",  # C subs: orange, yellow, brown
  "G>A" = "#377eb8", "G>C" = "#999999", "G>T" = "#984ea3",  # G subs: pink, gray, purple
  "T>A" = "#4daf4a", "T>C" = "#a65628", "T>G" = "#984ea3"   # T subs: cyan, teal, lime
)



# Function to standardize ratios by reference substitutions
# Can specify different references for 5' and 3' ends
standardize_by_substitution <- function(df, ref_sub_5prime, ref_sub_3prime = NULL) {
  # If only one reference provided, use it for both ends
  if (is.null(ref_sub_3prime)) {
    ref_sub_3prime <- ref_sub_5prime
  }
  
  # Get reference substitution data for 5' end
  ref_data_5prime <- df %>%
    filter(substitution == ref_sub_5prime, position_type == "5prime") %>%
    select(position, position_type, obs_exp_ratio) %>%
    rename(ref_ratio = obs_exp_ratio)
  
  # Get reference substitution data for 3' end
  ref_data_3prime <- df %>%
    filter(substitution == ref_sub_3prime, position_type == "3prime") %>%
    select(position, position_type, obs_exp_ratio) %>%
    rename(ref_ratio = obs_exp_ratio)
  
  # Combine reference data
  ref_data <- bind_rows(ref_data_5prime, ref_data_3prime)
  
  # Merge with main dataframe
  df_merged <- df %>%
    left_join(ref_data, by = c("position", "position_type")) %>%
    mutate(standardized_ratio = obs_exp_ratio / ref_ratio,
           ref_substitution = ifelse(position_type == "5prime", ref_sub_5prime, ref_sub_3prime))
  
  return(df_merged)
}

# Function to trim first and last x positions from the data
trim_positions <- function(df, trim_amount) {
  # Filter data:
  # For 5prime: keep positions > trim_amount
  # For 3prime: keep positions < -trim_amount (since they're negative)
  df_trimmed <- df %>%
    filter(
      (position_type == "5prime" & position > trim_amount) |
      (position_type == "3prime" & position < -trim_amount)
    )
  
  return(df_trimmed)
}

get_symmetric_pair <- function(sub) {
  bases <- strsplit(sub, ">")[[1]]
  paste(bases[2], bases[1], sep = ">")
}

# Read data
data <- read.csv("./combined_positional.csv")

file_pref <- "betula_combined_positional_"

# Reorder position_type factor so 5prime comes first
data$position_type <- factor(data$position_type, levels = c("5prime",
                                                            "3prime"))

# Create symmetric pair labels with <>
data <- data %>%
  mutate(
    symmetric_pair = pmin(substitution, sapply(substitution, get_symmetric_pair)),
    is_primary = substitution == symmetric_pair,
    # Create the <> label for legend
    pair_label = paste(substr(symmetric_pair, 1, 1), 
                      substr(symmetric_pair, 3, 3), 
                      sep = "<>")
  )


df_ratios <- data %>%
  # Define symmetric pairs
  mutate(
    pair = case_when(
      substitution %in% c("A>C", "C>A") ~ "A<>C",
      substitution %in% c("A>G", "G>A") ~ "A<>G",
      substitution %in% c("A>T", "T>A") ~ "A<>T",
      substitution %in% c("C>G", "G>C") ~ "C<>G",
      substitution %in% c("C>T", "T>C") ~ "C<>T",
      substitution %in% c("G>T", "T>G") ~ "G<>T"
    )
  ) %>%
  # Group by position, position_type, and pair
  group_by(position, position_type, pair) %>%
  # Calculate ratio (max/min to ensure >= 1)
  summarise(
    pi_q_max = max(pi_q),
    pi_q_min = min(pi_q),
    ratio = max(pi_q) / min(pi_q),
    .groups = "drop"
  )


# Choose your reference substitution types for standardization
# You can use different references for 5' and 3' ends
#REFERENCE_SUBSTITUTION_5PRIME <- 'C>G'
#REFERENCE_SUBSTITUTION_3PRIME <- 'C>G'  # Change this if you want a different reference for 3' end

# Create standardized data
#data_standardized <- standardize_by_substitution(data, 
#                                                  REFERENCE_SUBSTITUTION_5PRIME, 
#                                                  REFERENCE_SUBSTITUTION_3PRIME)

#data_trimmed <- trim_positions(data_standardized, trim_amount = 5)

data_trim <- trim_positions(data, trim_amount = 4)

df_ratios_trim <- data_trim %>%
  # Define symmetric pairs
  mutate(
    pair = case_when(
      substitution %in% c("A>C", "C>A") ~ "A<>C",
      substitution %in% c("A>G", "G>A") ~ "A<>G",
      substitution %in% c("A>T", "T>A") ~ "A<>T",
      substitution %in% c("C>G", "G>C") ~ "C<>G",
      substitution %in% c("C>T", "T>C") ~ "C<>T",
      substitution %in% c("G>T", "T>G") ~ "G<>T"
    )
  ) %>%
  # Group by position, position_type, and pair
  group_by(position, position_type, pair) %>%
  # Calculate ratio (max/min to ensure >= 1)
  summarise(
    pi_q_max = max(pi_q),
    pi_q_min = min(pi_q),
    ratio = max(pi_q) / min(pi_q),
    .groups = "drop"
  )

# # =============================================================================
# # ORIGINAL PLOTS (obs_exp_ratio)
# # =============================================================================
# 
# # Plot all substitution types
# ggplot(data, aes(x = position, y = obs_exp_ratio, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns (Obs/Exp Ratios)",
#        x = "Position in Read",
#        y = "Observed/Expected Ratio",
#        color = "Substitution") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_ratio.png"))
# 
# # Plot all substitution types - observed rate
# ggplot(data, aes(x = position, y = pi_q, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = expression(pi[i] * Q[ij]),
#        color = "Substitution") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_rate.png"))
# 
# # Plot all substitution types - observed rate
# ggplot(df_ratios, aes(x = position, y = ratio, color = pair, group = pair)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair) +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = expression(pi[i] * Q[ij] * " Ratio"),
#        color = "Substitution") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_rate.png"))
# 
# # Highlight specific substitutions (e.g., transitions vs transversions)
# data_annotated <- data %>%
#   mutate(sub_class = case_when(
#     substitution %in% c("C>T", "T>C", "A>G", "G>A") ~ "Transition",
#     TRUE ~ "Transversion"
#   ))
# 
# # Make sure position_type is ordered
# data_annotated$position_type <- factor(data_annotated$position_type, 
#                                        levels = c("5prime", "3prime"))
# 
# ggplot(data_annotated, aes(x = position, y = obs_exp_ratio, 
#                            color = substitution, linetype = sub_class)) +
#   geom_line(linewidth = 0.8) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_linetype_manual(values = c("Transition" = "solid", "Transversion" = "dotted")) +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Substitution Patterns: Transitions vs Transversions",
#        x = "Position in Read",
#        y = "Observed/Expected Ratio",
#        color = "Substitution",
#        linetype = "Type") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_ratio_type.png"))
# 
# # Create an absolute position for plotting (5' positions stay as is, 3' positions become positive)
# data <- data %>%
#   mutate(
#     abs_position = ifelse(position < 0, abs(position), position),
#     position_type = factor(position_type, levels = c("5prime", "3prime"))
#   )
# 
# # Alternative: Separate panels for each substitution type
# ggplot(data, aes(x = abs_position, y = obs_exp_ratio, color = position_type)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~substitution, scales = "free_y", ncol = 4) +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns by Type",
#        x = "Position from End",
#        y = "Observed/Expected Ratio",
#        color = "End") +
#   theme(legend.position = "bottom")
# ggsave(paste0(file_pref, "_all_ratio_panel.png"))
# 
# 
# 
# # =============================================================================
# # ORIGINAL TRIMMED PLOTS 
# # =============================================================================
# 
# # Plot all substitution types - STANDARDIZED
# ggplot(data_trim, aes(x = position, y = obs_exp_ratio, color = substitution,
#                       group = substitution)) + geom_line(linewidth = 0.8) +
#                    geom_point(size = 1.5) + facet_wrap(~position_type, scales =
#                                                        "free_x") +
#                    scale_color_manual(values = custom_colors) + theme_bw() +
#                    labs(title = "Trimmed Substitution Patterns",
#                              x = "Position in
#                                Read", y = "Observed/Expected Ratio", color =
#                                "Substitution") + theme(legend.position =
#                            "right", plot.subtitle = element_text(size = 9,
#                                                                  color =
#                                                                      "gray40"))
#                         ggsave(paste0(file_pref, "_all_ratio_trim.png"))
# 
# # Plot all substitution types - observed rate
# ggplot(data_trim, aes(x = position, y = pi_q, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Trimmed Positional Substitution Patterns",
#        x = "Position in Read",
#        y = expression(pi[i] * Q[ij]),
#        color = "Substitution") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_rate_trim.png"))
# 
# # Plot all substitution types - observed rate
# ggplot(df_ratios_trim, aes(x = position, y = ratio, color = pair, group = pair)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair) +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = expression(pi[i] * Q[ij] * " Ratio"),
#        color = "Substitution") +
#   theme(legend.position = "right")
# ggsave(paste0(file_pref, "_all_rate.png"))
# 
# 
# # Standardized with transitions vs transversions annotation
# data_annotated_trimmed <- data_trim %>%
#   mutate(sub_class = case_when(
#     substitution %in% c("C>T", "T>C", "A>G", "G>A") ~ "Transition",
#     TRUE ~ "Transversion"
#   ))
# 
# data_annotated_trimmed$position_type <- factor(data_annotated_trimmed$position_type, 
#                                                     levels = c("5prime", "3prime"))
# 
# ggplot(data_annotated_trimmed, aes(x = position, y = obs_exp_ratio, 
#                                         color = substitution, linetype = sub_class)) +
#   geom_line(linewidth = 0.8) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_linetype_manual(values = c("Transition" = "solid", "Transversion" = "dotted")) +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Trimmed Patterns: Transitions vs Transversions",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position in Read",
#        y = "Observed/Expected Ratio",
#        color = "Substitution",
#        linetype = "Type") +
#   theme(legend.position = "right",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave(paste0(file_pref,"_all_type_trim.png"))
# 
# # Add absolute position to standardized data
# data_trim_abs <- data_trim %>%
#   mutate(
#     abs_position = ifelse(position < 0, abs(position), position),
#     position_type = factor(position_type, levels = c("5prime", "3prime"))
#   )
# 
# # Alternative: Separate panels for each substitution type - STANDARDIZED
# ggplot(data_trim_abs, aes(x = abs_position, y = obs_exp_ratio, color = position_type)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~substitution, scales = "free_y", ncol = 4) +
#   theme_bw() +
#   labs(title = "Trimmed Substitution Patterns by Type",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position from End",
#        y = "Observed/Expected Ratio",
#        color = "End") +
#   theme(legend.position = "bottom",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave(paste0(file_pref, "_all_panel_trim.png"))


# # =============================================================================
# # STANDARDIZED PLOTS (standardized_ratio relative to reference substitution)
# # =============================================================================
# 
# # Plot all substitution types - STANDARDIZED
# ggplot(data_standardized, aes(x = position, y = standardized_ratio, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Standardized Substitution Patterns",
#        subtitle = paste0("5' relative to ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' relative to ", REFERENCE_SUBSTITUTION_3PRIME, 
#                         " | Values > 1 = excess, < 1 = depletion"),
#        x = "Position in Read",
#        y = "Standardized Ratio",
#        color = "Substitution") +
#   theme(legend.position = "right",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized.png")
# 
# # Standardized with transitions vs transversions annotation
# data_standardized_annotated <- data_standardized %>%
#   mutate(sub_class = case_when(
#     substitution %in% c("C>T", "T>C", "A>G", "G>A") ~ "Transition",
#     TRUE ~ "Transversion"
#   ))
# 
# data_standardized_annotated$position_type <- factor(data_standardized_annotated$position_type, 
#                                                     levels = c("5prime", "3prime"))
# 
# ggplot(data_standardized_annotated, aes(x = position, y = standardized_ratio, 
#                                         color = substitution, linetype = sub_class)) +
#   geom_line(linewidth = 0.8) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_linetype_manual(values = c("Transition" = "solid", "Transversion" = "dotted")) +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Standardized Patterns: Transitions vs Transversions",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position in Read",
#        y = "Standardized Ratio",
#        color = "Substitution",
#        linetype = "Type") +
#   theme(legend.position = "right",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized_type.png")
# 
# # Add absolute position to standardized data
# data_standardized <- data_standardized %>%
#   mutate(
#     abs_position = ifelse(position < 0, abs(position), position),
#     position_type = factor(position_type, levels = c("5prime", "3prime"))
#   )
# 
# # Alternative: Separate panels for each substitution type - STANDARDIZED
# ggplot(data_standardized, aes(x = abs_position, y = standardized_ratio, color = position_type)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
#   facet_wrap(~substitution, scales = "free_y", ncol = 4) +
#   theme_bw() +
#   labs(title = "Standardized Substitution Patterns by Type",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position from End",
#        y = "Standardized Ratio",
#        color = "End") +
#   theme(legend.position = "bottom",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized_panel.png")
# 
# # =============================================================================
# # STANDARDIZED TRIMMED PLOTS (standardized_ratio relative to reference substitution for trimmed data)
# # =============================================================================
# 
# # Plot all substitution types - STANDARDIZED
# ggplot(data_trimmed, aes(x = position, y = standardized_ratio, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Standardized Substitution Patterns",
#        subtitle = paste0("5' relative to ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' relative to ", REFERENCE_SUBSTITUTION_3PRIME, 
#                         " | Values > 1 = excess, < 1 = depletion"),
#        x = "Position in Read",
#        y = "Standardized Ratio",
#        color = "Substitution") +
#   theme(legend.position = "right",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized_trimmed.png")
# 
# # Standardized with transitions vs transversions annotation
# data_standardized_annotated_trimmed <- data_trimmed %>%
#   mutate(sub_class = case_when(
#     substitution %in% c("C>T", "T>C", "A>G", "G>A") ~ "Transition",
#     TRUE ~ "Transversion"
#   ))
# 
# data_standardized_annotated_trimmed$position_type <- factor(data_standardized_annotated_trimmed$position_type, 
#                                                     levels = c("5prime", "3prime"))
# 
# ggplot(data_standardized_annotated_trimmed, aes(x = position, y = standardized_ratio, 
#                                         color = substitution, linetype = sub_class)) +
#   geom_line(linewidth = 0.8) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_linetype_manual(values = c("Transition" = "solid", "Transversion" = "dotted")) +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "Standardized Patterns: Transitions vs Transversions",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position in Read",
#        y = "Standardized Ratio",
#        color = "Substitution",
#        linetype = "Type") +
#   theme(legend.position = "right",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized_type_trimmed.png")
# 
# # Add absolute position to standardized data
# data_standardized_trimmed <- data_trimmed %>%
#   mutate(
#     abs_position = ifelse(position < 0, abs(position), position),
#     position_type = factor(position_type, levels = c("5prime", "3prime"))
#   )
# 
# # Alternative: Separate panels for each substitution type - STANDARDIZED
# ggplot(data_standardized_trimmed, aes(x = abs_position, y = standardized_ratio, color = position_type)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
#   facet_wrap(~substitution, scales = "free_y", ncol = 4) +
#   theme_bw() +
#   labs(title = "Standardized Substitution Patterns by Type",
#        subtitle = paste0("5' vs ", REFERENCE_SUBSTITUTION_5PRIME, 
#                         ", 3' vs ", REFERENCE_SUBSTITUTION_3PRIME),
#        x = "Position from End",
#        y = "Standardized Ratio",
#        color = "End") +
#   theme(legend.position = "bottom",
#         plot.subtitle = element_text(size = 9, color = "gray40"))
# ggsave("./0_betula_10_all_standardized_panel_trimmed.png")
# 
# 
# # =============================================================================
# # Export standardized data
# # =============================================================================
# write.csv(data_standardized, "./betula_mine_standardized_data.csv", row.names = FALSE)
#
# # Print summary
# cat("\n=== Standardization Complete ===\n")
# cat("Reference Substitution (5' end):", REFERENCE_SUBSTITUTION_5PRIME, "\n")
# cat("Reference Substitution (3' end):", REFERENCE_SUBSTITUTION_3PRIME, "\n")
# cat("\nOriginal plots saved:\n")
# cat("  - betula_mine_all_ratio.png\n")
# cat("  - betula_mine_all_rate.png\n")
# cat("  - betula_mine_all_ratio_type.png\n")
# cat("  - betula_mine_all_ratio_panel.png\n")
# cat("\nStandardized plots saved:\n")
# cat("  - betula_mine_all_standardized.png\n")
# cat("  - betula_mine_all_standardized_type.png\n")
# cat("  - betula_mine_all_standardized_panel.png\n")
# cat("\nStandardized data saved:\n")
# cat("  - betula_mine_standardized_data.csv\n")

# =============================================================================
# COMBINED FIGURE
# =============================================================================

  # Remove C<>T and G<>A substitutions
  #data_trim
data_filtered_trim <- data %>%
  filter(!pair_label %in% c("C<>T", "A<>G"))



# Create individual plots for combination
# Plot a: Observed rate - full length
# plot_a <- ggplot(data, aes(x = position, y = observed_rate, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "a",
#        x = "Position in Read",
#        y = "Observed Rate",
#        color = "Substitution") +
#   theme(legend.position = "none")

# Plot b: Observed rate - trimmed
#plot_b <- ggplot(data_trim, aes(x = position, y = observed_rate, color = substitution, group = substitution)) +
#  geom_line(linewidth = 0.8) +
#  geom_point(size = 1.5) +
#  facet_wrap(~position_type, scales = "free_x") +
#  scale_color_manual(values = custom_colors) +
#  theme_bw() +
#  labs(title = "b",
#       x = "Position in Read",
#       y = "Observed Rate",
#       color = "Substitution") +
#  theme(legend.position = "none")

# Combine plots with patchwork
# combined_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")

# Create color and linetype mappings for all 12 substitution types
sub_order <- c("A>C","C>A", "A>G","G>A", "A>T","T>A",
               "C>G","G>C", "C>T","T>C", "G>T","T>G")
data$substitution <- factor(data$substitution, levels = sub_order)
data_filtered_trim$substitution <- factor(data_filtered_trim$substitution, levels = sub_order)

sub_linetypes <- c(
  "A>C" = "solid", "C>A" = "dashed",
  "A>G" = "solid", "G>A" = "dashed",
  "A>T" = "solid", "T>A" = "dashed",
  "C>G" = "solid", "G>C" = "dashed",
  "C>T" = "solid", "T>C" = "dashed",
  "G>T" = "solid", "T>G" = "dashed"
)

sub_shapes <- c(
  "A>C" = 16, "C>A" = 17,
  "A>G" = 16, "G>A" = 17,
  "A>T" = 16, "T>A" = 17,
  "C>G" = 16, "G>C" = 17,
  "C>T" = 16, "T>C" = 17,
  "G>T" = 16, "T>G" = 17
)

custom_colors_pair_2 <- c(
  "A<>C" = "#D55E00",  # vermillion
  "A<>G" = "#0072B2",  # blue
  "A<>T" = "#009E73",  # bluish green
  "C<>G" = "#CC79A7",  # reddish purple
  "C<>T" = "#E69F00",  # orange
  "G<>T" = "#56B4E9"   # sky blue
)

sub_colors <- c(
  "A>C" = custom_colors_pair_2[["A<>C"]],
  "C>A" = custom_colors_pair_2[["A<>C"]],
  "A>G" = custom_colors_pair_2[["A<>G"]],
  "G>A" = custom_colors_pair_2[["A<>G"]],
  "A>T" = custom_colors_pair_2[["A<>T"]],
  "T>A" = custom_colors_pair_2[["A<>T"]],
  "C>G" = custom_colors_pair_2[["C<>G"]],
  "G>C" = custom_colors_pair_2[["C<>G"]],
  "C>T" = custom_colors_pair_2[["C<>T"]],
  "T>C" = custom_colors_pair_2[["C<>T"]],
  "G>T" = custom_colors_pair_2[["G<>T"]],
  "T>G" = custom_colors_pair_2[["G<>T"]]
)

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
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),       
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 14)  
  )

# Shared aesthetic layer: same mapping + scales for every panel
sub_aes <- function(yvar) {
  list(
    aes(x = position, y = .data[[yvar]],
        color = substitution,
        linetype = substitution,
        shape = substitution,
        group = substitution),
    geom_line(linewidth = 0.75, alpha = 0.9),
    geom_point(size = 2, stroke = 0.4),
    facet_wrap(~position_type, scales = "free_x"),
    scale_color_manual(values = sub_colors),
    scale_linetype_manual(values = sub_linetypes),
    scale_shape_manual(values = sub_shapes),
    theme_clean
  )
}

# Plot a: carries the single legend for the whole figure
plot_a <- ggplot(data) +
  sub_aes("pi_q") +
  labs(title = "A",
       x = "Position in Read",
       y = expression(pi*" Q"),
       color = "Substitution",
       linetype = "Substitution",
       shape = "Substitution") +
  guides(
    color    = guide_legend(nrow = 2, override.aes = list(size = 2.4, linewidth = 0.9)),
    linetype = guide_legend(nrow = 2),
    shape    = guide_legend(nrow = 2)
  )

# Plot b
plot_b <- ggplot(data_filtered_trim) +
  sub_aes("pi_q") +
  labs(title = "B",
       x = "Position in Read",
       y = expression(pi*" Q")) +
  guides(color = "none", linetype = "none", shape = "none")

# Plot c
plot_c <- ggplot(data) +
  sub_aes("obs_exp_ratio") +
  labs(title = "D",
       x = "Position in Read",
       y = "Observed/Expected Ratio") +
  guides(color = "none", linetype = "none", shape = "none")

# Plot d
plot_d <- ggplot(data_filtered_trim) +
  sub_aes("obs_exp_ratio") +
  labs(title = "D",
       x = "Position in Read",
       y = "Observed/Expected Ratio") +
  guides(color = "none", linetype = "none", shape = "none")

# Combine
combined_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

ggsave(paste0(file_pref, "_combined_figure.png"),
       plot = combined_plot,
       width = 14,
       height = 10,
       dpi = 300)

# # Plot a: Keep legend
# plot_a <- ggplot(data, aes(x = position, y = pi_q,
#                  color = substitution,
#                  linetype = substitution,
#                  shape = substitution,
#                  group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = sub_colors) +
#   scale_linetype_manual(values = sub_linetypes) +
#   scale_shape_manual(values = sub_shapes) +
#   theme_bw() +
#   labs(title = "a",
#        x = "Position in Read",
#        y = expression(pi*" Q"),
#        color = "Substitution (a, b)",
#        linetype = "Substitution (a, b)",
#        shape = "Substitution (a, b)")
# 
# # Plot b: Remove legend (same as before)
# plot_b <- ggplot(data_filtered_trim, aes(x = position, y = pi_q,
#                  color = pair_label,
#                  group = substitution,
#                  linetype = is_primary,
#                  shape = is_primary)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair_2) +
#   scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
#                         guide = "none") +
#   scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17),
#                      guide = "none") +
#   theme_bw() +
#   labs(title = "b",
#        x = "Position in Read",
#        y = expression(pi*" Q")) +
#   guides(color = "none")
# 
# # Plot c: Keep legend
# plot_c <- ggplot(data, aes(x = position, y = obs_exp_ratio, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "c",
#        x = "Position in Read",
#        y = "Observed/Expected Ratio",
#        color = "Substitution (c,d)")  # Different name
# 
# # Plot d: Remove legend
# plot_d <- ggplot(data_filtered_trim, aes(x = position, y = obs_exp_ratio, color = substitution, group = substitution)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors) +
#   theme_bw() +
#   labs(title = "d",
#        x = "Position in Read",
#        y = "Observed/Expected Ratio") +
#   guides(color = "none")  # Explicitly remove guide
# 
# # Combine plots
# combined_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# combined_plot
# 
# # Save combined figure
# ggsave(paste0(file_pref, "_combined_figure.png"),
#        plot = combined_plot,
#        width = 14,
#        height = 10,
#        dpi = 300)
# 
# 
# # Create symmetric pair labels with <>
# data <- data %>%
#   mutate(
#     symmetric_pair = pmin(substitution, sapply(substitution, get_symmetric_pair)),
#     is_primary = substitution == symmetric_pair,
#     # Create the <> label for legend
#     pair_label = paste(substr(symmetric_pair, 1, 1), 
#                       substr(symmetric_pair, 3, 3), 
#                       sep = "<>")
#   )
# 
# # Plot with <> labels and no linetype legend
# ggplot(data, aes(x = position, y = pi_q, 
#                  color = pair_label, 
#                  group = substitution,
#                  linetype = is_primary)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair) +
#   scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
#                         guide = "none") +  # This removes the linetype legend
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = "Observed Rate (pi_q)",
#        color = "Substitution") +
#   theme(legend.position = "right")
# 
# data_trim <- data_trim %>%
#   mutate(
#     symmetric_pair = pmin(substitution, sapply(substitution, get_symmetric_pair)),
#     is_primary = substitution == symmetric_pair,
#     # Create the <> label for legend
#     pair_label = paste(substr(symmetric_pair, 1, 1), 
#                       substr(symmetric_pair, 3, 3), 
#                       sep = "<>")
#   )
# 
# # Plot with <> labels and no linetype legend
# ggplot(data_trim, aes(x = position, y = pi_q, 
#                  color = pair_label, 
#                  group = substitution,
#                  linetype = is_primary)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair_2) +
#   scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
#                         guide = "none") +  # This removes the linetype legend
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = "Observed Rate (pi_q)",
#        color = "Substitution") +
#   theme(legend.position = "right")
# 
#   # Remove C<>T and G<>A substitutions
# data_filtered_trim <- data_trim %>%
#   filter(!pair_label %in% c("C<>T", "A<>G"))
# 
# # Then plot with the filtered data
# ggplot(data_filtered_trim, aes(x = position, y = pi_q, 
#                  color = pair_label, 
#                  group = substitution,
#                  linetype = is_primary)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 1.5) +
#   facet_wrap(~position_type, scales = "free_x") +
#   scale_color_manual(values = custom_colors_pair_2) +
#   scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
#                         guide = "none") +
#   theme_bw() +
#   labs(title = "Positional Substitution Patterns",
#        x = "Position in Read",
#        y = "Observed Rate (pi_q)",
#        color = "Substitution") +
#   theme(legend.position = "right")
