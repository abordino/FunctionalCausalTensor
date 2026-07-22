setwd("~/Desktop/code")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# Load saved results ------------------------------------------------------------------

results_dir = "Results/masked_nonmasked_rank3_bootstrap"
plots_dir = "Plots"

target_name = "robbery"
target_label = "Robbery rate, log"

result_file = file.path(
  results_dir,
  "castle_masked_nonmasked_bootstrap_results_robbery.rds"
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(result_file)) {
  stop("Result file not found: ", result_file)
}

results = readRDS(result_file)

functional_order = results$metadata$functional_order
rank_value = results$metadata$rank_value
B = results$metadata$B

method_colors = c(
  Tensor = "#1f77b4",
  Matrix = "#ff7f0e"
)

design_levels = c(
  "Masked: 5 full rows",
  "Masked: 10 full rows",
  "Masked: 15 full rows",
  "Non-masked"
)

# Helpers -----------------------------------------------------------------

design_short = function(design_id, design_label) {
  case_when(
    str_detect(design_id, "masked_full5_") ~ "Masked: 5 full rows",
    str_detect(design_id, "masked_full10_") ~ "Masked: 10 full rows",
    str_detect(design_id, "masked_full15_") ~ "Masked: 15 full rows",
    design_id == "nonmasked_original" ~ "Non-masked",
    TRUE ~ design_label
  )
}

symmetric_limits = function(data, padding = 1.08) {
  values = c(data$ci_low, data$ci_high, data$point_estimate)
  values = values[is.finite(values)]
  
  limit = if (length(values)) max(abs(values)) * padding else 1
  
  if (!is.finite(limit) || limit == 0) {
    limit = 1
  }
  
  c(-limit, limit)
}

# Psi0 plot data ----------------------------------------------------------

plot_data = results$results_with_ci %>%
  filter(
    crime == target_name,
    quantity %in% c("Psi0", "Psi0_matrix")
  ) %>%
  mutate(
    method = recode(
      quantity,
      Psi0 = "Tensor",
      Psi0_matrix = "Matrix"
    ),
    method = factor(
      method,
      levels = c("Tensor", "Matrix")
    ),
    design_short = factor(
      design_short(design_id, design_label),
      levels = design_levels
    ),
    functional = factor(
      functional,
      levels = functional_order
    ),
    functional_num = as.numeric(
      factor(
        functional,
        levels = rev(functional_order)
      )
    ),
    y_pos = functional_num +
      if_else(method == "Tensor", -0.13, 0.13)
  )

# Display and save Psi0 plot ----------------------------------------------

functional_breaks = seq_along(rev(functional_order))
functional_labels = rev(functional_order)

psi0_figure = ggplot(plot_data) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = ci_low,
      xend = ci_high,
      y = y_pos,
      yend = y_pos,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = ci_low,
      xend = ci_low,
      y = y_pos - 0.055,
      yend = y_pos + 0.055,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = ci_high,
      xend = ci_high,
      y = y_pos - 0.055,
      yend = y_pos + 0.055,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      x = point_estimate,
      y = y_pos,
      color = method,
      shape = method
    ),
    size = 2.6,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ design_short,
    ncol = 2,
    scales = "fixed"
  ) +
  coord_cartesian(
    xlim = symmetric_limits(plot_data)
  ) +
  scale_x_continuous(
    breaks = breaks_extended(n = 5),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    breaks = functional_breaks,
    labels = functional_labels,
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(
    values = c(
      Tensor = 16,
      Matrix = 17
    )
  ) +
  labs(
    title = paste0(
      "Castle Doctrine | Psi0: ",
      target_label
    ),
    subtitle = paste0(
      "Rank r = ",
      rank_value,
      "; B = ",
      B,
      "; intervals are point estimate +/- 1.96 x bootstrap SE"
    ),
    x = expression(
      "Estimated quantity with bootstrap-SE 95% CI"
    ),
    y = "Functional type",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

plot_file = file.path(
  plots_dir,
  "castle_psi0_robbery_rate_log_robbery.png"
)

print(psi0_figure)

ggsave(
  filename = plot_file,
  plot = psi0_figure,
  width = 11,
  height = 7.5,
  dpi = 320
)

cat("\nSaved plot:\n", plot_file, "\n", sep = "")

# Print non-masked Delta table ---------------------------------------

delta_table = results$results_with_ci %>%
  filter(
    crime == target_name,
    design_id == "nonmasked_original",
    quantity %in% c(
      "Delta_h",
      "Delta_h_matrix"
    )
  ) %>%
  mutate(
    method = recode(
      quantity,
      Delta_h = "\u0394(h)",
      Delta_h_matrix = "\u0394(h)^mat"
    ),
    functional = factor(
      functional,
      levels = functional_order
    ),
    result = if_else(
      is.na(ci_low) | is.na(ci_high),
      sprintf(
        "%.4f\n(NA, NA)",
        point_estimate
      ),
      sprintf(
        "%.4f\n(%.4f, %.4f)",
        point_estimate,
        ci_low,
        ci_high
      )
    )
  ) %>%
  select(
    method,
    functional,
    result
  ) %>%
  pivot_wider(
    names_from = functional,
    values_from = result
  ) %>%
  arrange(
    factor(
      method,
      levels = c(
        "\u0394(h)",
        "\u0394(h)^mat"
      )
    )
  )

cat("\nNon-masked robbery results:\n")
print(
  delta_table,
  width = Inf
)
