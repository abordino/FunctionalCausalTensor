setwd("~/Desktop/code")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# Configuration -----------------------------------------------------------

results_dir = "Results/oxford_deaths_target_layer_bootstrap/"
plots_dir = "Plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

results_files = list.files(
  results_dir,
  pattern = "^oxford_deaths_target_layer_bootstrap_.*\\.rds$",
  full.names = TRUE
)

if (!length(results_files)) {
  stop("No saved analysis RDS found in: ", results_dir)
}

results_file = results_files[which.max(file.info(results_files)$mtime)]
analysis = readRDS(results_file)

Y = analysis$Y
Omega_policy = analysis$Omega_policy
results_with_ci = analysis$results_with_ci
main_r = analysis$metadata$r_grid[1]
B = analysis$metadata$B

tensor_data = readRDS(analysis$metadata$tensor_file)

# Plot 1: two outcome panels ----------------------------------------------

layers = dimnames(Y)[[3]]

if (length(layers) != 2) {
  stop("The saved Y tensor must contain exactly two layers.")
}

layer_map = analysis$layer_map %>%
  transmute(
    layer = as.character(outcome),
    policy = as.character(policy),
    panel_label = as.character(panel_label)
  )

plot_layer_map = tibble(layer = layers) %>%
  left_join(layer_map, by = "layer") %>%
  mutate(
    outcome_cap = unname(tensor_data$outcome_caps[layer])
  )

if (anyNA(plot_layer_map$outcome_cap)) {
  stop("An outcome cap is missing for at least one plotted Y layer.")
}

outcome_plot_data = as.data.frame.table(
  Y,
  responseName = "value",
  stringsAsFactors = FALSE
)

names(outcome_plot_data)[1:3] = c(
  "country",
  "time",
  "layer"
)

policy_plot_data = as.data.frame.table(
  Omega_policy,
  responseName = "policy_on",
  stringsAsFactors = FALSE
)

names(policy_plot_data)[1:3] = c(
  "country",
  "time",
  "layer"
)

blue_pal = grDevices::colorRampPalette(
  c("lightblue", "blue")
)(101)

red_pal = grDevices::colorRampPalette(
  c("pink", "red")
)(101)

country_levels = dimnames(Y)[[1]]

country_labels = country_levels
country_labels[seq_along(country_labels) %% 1 == 1] = ""
names(country_labels) = country_levels

y_plot_df = as_tibble(outcome_plot_data) %>%
  left_join(
    as_tibble(policy_plot_data),
    by = c("country", "time", "layer")
  ) %>%
  left_join(
    plot_layer_map,
    by = "layer"
  ) %>%
  drop_na(value, policy_on) %>%
  mutate(
    date = as.Date(time),
    policy_on = as.logical(policy_on),
    
    value_capped = pmin(value, outcome_cap),
    value_scaled = value_capped / outcome_cap,
    
    intensity_index =
      as.integer(round(value_scaled * 100)) + 1L,
    
    intensity_index =
      pmin(101L, pmax(1L, intensity_index)),
    
    fill_colour = if_else(
      policy_on,
      red_pal[intensity_index],
      blue_pal[intensity_index]
    ),
    
    country = factor(
      country,
      levels = country_levels
    ),
    
    panel_label = factor(
      panel_label,
      levels = plot_layer_map$panel_label
    )
  )

p_y = ggplot(
  y_plot_df,
  aes(
    x = date,
    y = country,
    fill = fill_colour
  )
) +
  geom_tile() +
  facet_wrap(
    ~ panel_label,
    nrow = 1
  ) +
  scale_y_discrete(
    breaks = country_levels,
    labels = country_labels,
    drop = FALSE
  ) +
  scale_fill_identity() +
  labs(
    title = paste0(
      "Outcome tensor values between ",
      tensor_data$start.date,
      " to ",
      tensor_data$time_horizon,
      " \nusing outcomes at t + ",
      tensor_data$delay_days,
      " days"
    ),
    x = "Policy date",
    y = "Country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.position = "none"
  )

panel_file = file.path(
  plots_dir,
  "data_panels.png"
)

ggsave(
  filename = panel_file,
  plot = p_y,
  width = 18,
  height = 8,
  units = "in",
  dpi = 600
)

print(p_y)

# Plot 2: final ATE and Trend results -------------------------------------

method_colors = c(
  "Tensor" = "#1f77b4",
  "Matrix" = "#ff7f0e"
)

get_symmetric_ylim = function(
    dat,
    pad_mult = 1.08
) {
  y_vals = c(
    dat$ci_low,
    dat$ci_high,
    dat$point_estimate
  )
  
  y_vals = y_vals[is.finite(y_vals)]
  
  if (length(y_vals) == 0) {
    return(c(-1, 1))
  }
  
  y_abs_max = max(
    abs(y_vals),
    na.rm = TRUE
  )
  
  if (!is.finite(y_abs_max) || y_abs_max == 0) {
    return(c(-1, 1))
  }
  
  y_abs_max = y_abs_max * pad_mult
  
  c(-y_abs_max, y_abs_max)
}

psi_delta_plot_data = results_with_ci %>%
  filter(
    r == main_r,
    functional %in% c("ATE", "Trend"),
    quantity %in% c(
      "Psi0",
      "Psi0_matrix",
      "Delta_h",
      "Delta_h_matrix"
    )
  ) %>%
  mutate(
    estimand = case_when(
      quantity %in%
        c("Psi0", "Psi0_matrix") ~ "Psi0",
      
      quantity %in%
        c("Delta_h", "Delta_h_matrix") ~ "Delta"
    ),
    
    method = case_when(
      quantity %in%
        c("Psi0", "Delta_h") ~ "Tensor",
      
      quantity %in%
        c("Psi0_matrix", "Delta_h_matrix") ~ "Matrix"
    ),
    
    method = factor(
      method,
      levels = c("Tensor", "Matrix")
    ),
    
    x_group = paste(
      functional,
      estimand,
      sep = "\n"
    ),
    
    x_group = factor(
      x_group,
      levels = c(
        "ATE\nPsi0",
        "ATE\nDelta",
        "Trend\nPsi0",
        "Trend\nDelta"
      )
    ),
    
    x_num = as.numeric(x_group),
    
    method_offset = if_else(
      method == "Tensor",
      -0.13,
      0.13
    ),
    
    x_pos = x_num + method_offset
  )

y_limits = get_symmetric_ylim(
  psi_delta_plot_data
)

final_plot = ggplot(
  psi_delta_plot_data
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = x_pos,
      xend = x_pos,
      y = ci_low,
      yend = ci_high,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = x_pos - 0.055,
      xend = x_pos + 0.055,
      y = ci_low,
      yend = ci_low,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = x_pos - 0.055,
      xend = x_pos + 0.055,
      y = ci_high,
      yend = ci_high,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      x = x_pos,
      y = point_estimate,
      color = method,
      shape = method
    ),
    size = 2.6,
    na.rm = TRUE
  ) +
  coord_cartesian(
    ylim = y_limits
  ) +
  scale_x_continuous(
    breaks = seq_along(
      levels(psi_delta_plot_data$x_group)
    ),
    labels = levels(
      psi_delta_plot_data$x_group
    ),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 5),
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  scale_color_manual(
    values = method_colors
  ) +
  scale_shape_manual(
    values = c(
      "Tensor" = 16,
      "Matrix" = 17
    )
  ) +
  labs(
    title =
      "CovidOx | Deaths target layer: Psi0 and Delta",
    
    subtitle = paste0(
      "Rank r = ",
      main_r,
      "; B = ",
      B,
      "; intervals are point estimate +/- 1.96 x bootstrap SE"
    ),
    
    x = "",
    
    y =
      "Estimated quantity with bootstrap-SE 95% CI",
    
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9)
  )

final_file = file.path(
  plots_dir,
  "final_ATE_Trend.png"
)

ggsave(
  filename = final_file,
  plot = final_plot,
  width = 11,
  height = 7.5,
  dpi = 320
)

print(final_plot)

cat("Loaded:", results_file, "\n")
cat("Saved:", panel_file, "\n")
cat("Saved:", final_file, "\n")