setwd("~/Documents/phd/projects/causalMatrix/code/synthetic/4Block")

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(grid)

res_dir = "result"
fig_dir = "figure"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

read_result = function(file) {
  path = file.path(res_dir, file)
  
  if (!file.exists(path)) {
    stop(
      "File not found: ", path,
      "\n\nAvailable summary files are:\n",
      paste(list.files(res_dir, pattern = "summary.*\\.csv$", recursive = TRUE), collapse = "\n")
    )
  }
  
  read.csv(path)
}

clean_methods = function(df) {
  df %>%
    mutate(
      method_label = case_when(
        method == "oracle_pool" ~ "Oracle pooled",
        method == "oracle_nopool" ~ "Oracle no-pool",
        method == "oracle_local" ~ "Oracle local",
        method == "real_pool" ~ "Estimated pooled",
        method == "real_nopool" ~ "Estimated no-pool",
        method == "real_pool_misspec" ~ "Estimated pooled",
        method == "real_nopool_misspec" ~ "Estimated no-pool",
        TRUE ~ method
      ),
      method_label = factor(
        method_label,
        levels = c(
          "Oracle pooled",
          "Oracle no-pool",
          "Oracle local",
          "Estimated pooled",
          "Estimated no-pool"
        )
      )
    )
}

# ------------------------------------------------------------
# Load all four summaries
# ------------------------------------------------------------

A = read_result(
  "summary_compare_pooling_oracle_varyK_randomXY_SNR1_N1win70-70_T1win60-60.csv"
)

B = read_result(
  "Zoomsummary_compare_pooling_oracle_varyK_randomXY_SNR1_N1win70-70_T1win60-60.csv"
)

C = read_result(
  "summary_compare_pooling_oracle_varyK_SNR1_N1win30-70_T1win30-60_leadingEigXY_trueRank6_estRank11.csv"
)

D = read_result(
  "Zoomsummary_compare_pooling_oracleLocal_varyK_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.csv"
)

A = clean_methods(A)
B = clean_methods(B)
C = clean_methods(C)
D = clean_methods(D)

# ------------------------------------------------------------
# Keep methods for each panel
# ------------------------------------------------------------

A = A %>%
  filter(method %in% c(
    "oracle_pool",
    "oracle_nopool",
    "oracle_local",
    "real_pool",
    "real_nopool"
  ))

B = B %>%
  filter(method %in% c(
    "oracle_pool",
    "oracle_local",
    "real_pool"
  ))

C = C %>%
  filter(method %in% c(
    "oracle_pool",
    "oracle_nopool",
    "oracle_local",
    "real_pool_misspec",
    "real_nopool_misspec"
  ))

D = D %>%
  filter(method %in% c(
    "oracle_local",
    "real_pool"
  )) %>%
  mutate(
    SNR_label = factor(
      SNR_target,
      levels = c(1, 100, 10000),
      labels = c("SNR = 1", "SNR = 100", "SNR = 10000")
    ),
    logMSE = log10(MSE),
    log_low = log10(pmax(MSE - CI95, .Machine$double.eps)),
    log_high = log10(MSE + CI95)
  )

method_cols = c(
  "Oracle pooled" = "#7570b3",
  "Oracle no-pool" = "#e7298a",
  "Oracle local" = "#f95919",
  "Estimated pooled" = "#1b9e77",
  "Estimated no-pool" = "#d95f02"
)

method_shapes = c(
  "Oracle pooled" = 15,
  "Oracle no-pool" = 18,
  "Oracle local" = 20,
  "Estimated pooled" = 16,
  "Estimated no-pool" = 17
)

snr_linetypes = c(
  "SNR = 1" = "solid",
  "SNR = 100" = "dashed",
  "SNR = 10000" = "dotted"
)

base_theme = theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "gray88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.35, "lines"),
    legend.spacing.x = unit(0.8, "lines"),
    legend.spacing.y = unit(0.10, "lines"),
    plot.margin = margin(6, 6, 6, 6)
  )

# ------------------------------------------------------------
# Generic MSE panel
# ------------------------------------------------------------

mse_plot = function(df, title, subtitle = NULL, tag = NULL, show_legend = FALSE) {
  ggplot(
    df,
    aes(
      x = K,
      y = MSE,
      color = method_label,
      shape = method_label,
      group = method_label
    )
  ) +
    geom_errorbar(
      aes(ymin = pmax(0, MSE - CI95), ymax = MSE + CI95),
      width = 0,
      linewidth = 0.35,
      alpha = 0.55,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.75, show.legend = show_legend) +
    geom_point(size = 2.1, show.legend = show_legend) +
    scale_color_manual(
      name = NULL,
      values = method_cols,
      drop = TRUE
    ) +
    scale_shape_manual(
      name = NULL,
      values = method_shapes,
      drop = TRUE
    ) +
    scale_x_continuous(breaks = sort(unique(df$K))) +
    scale_y_continuous(labels = scientific) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(K),
      y = "MSE",
      tag = tag
    ) +
    guides(
      color = guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0.9, size = 2.4)
      ),
      shape = guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    base_theme +
    theme(
      plot.tag = element_text(face = "bold", size = 13),
      plot.tag.position = c(0.02, 0.98),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# ------------------------------------------------------------
# Panels
# ------------------------------------------------------------

pA = mse_plot(
  A,
  title = "Random-XY design",
  subtitle = expression(N[1] == 70 ~ "," ~ T[1] == 60 ~ "," ~ SNR == 1),
  tag = "(a)",
  show_legend = TRUE
)

pB = mse_plot(
  B,
  title = "Random-XY design, large K",
  subtitle = expression(N[1] == 70 ~ "," ~ T[1] == 60 ~ "," ~ SNR == 1),
  tag = "(b)",
  show_legend = FALSE
)

pC = mse_plot(
  C,
  title = "Leading-eigenvector XY, misspecified rank",
  subtitle = expression(true~rank == 6 ~ "," ~ estimated~rank == 11 ~ "," ~ N[1] %in% "[30,70]" ~ "," ~ T[1] %in% "[30,60]"),
  tag = "(c)",
  show_legend = FALSE
)

pD = ggplot(
  D,
  aes(
    x = K,
    y = logMSE,
    color = method_label,
    shape = method_label,
    linetype = SNR_label,
    group = interaction(method_label, SNR_label)
  )
) +
  geom_errorbar(
    aes(ymin = log_low, ymax = log_high),
    width = 0,
    linewidth = 0.35,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.75, show.legend = TRUE) +
  geom_point(size = 2.1, show.legend = FALSE) +
  scale_color_manual(
    name = NULL,
    values = method_cols,
    drop = TRUE,
    guide = "none"
  ) +
  scale_shape_manual(
    name = NULL,
    values = method_shapes,
    drop = TRUE,
    guide = "none"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = snr_linetypes
  ) +
  scale_x_continuous(breaks = sort(unique(D$K))) +
  labs(
    title = expression(log[10](MSE) ~ "vs." ~ K),
    subtitle = expression(N[1] == 70 ~ "," ~ T[1] == 60),
    x = expression(K),
    y = expression(log[10](MSE)),
    tag = "(d)"
  ) +
  guides(
    linetype = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(color = "black", shape = NA, linewidth = 0.9)
    )
  ) +
  base_theme +
  theme(
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "bottom"
  )


combined = (pA + pB) / (pC + pD) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank(),
    legend.spacing.y = unit(0.1, "lines"),
    legend.spacing.x = unit(0.8, "lines")
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  filename = file.path(fig_dir, "combined_2x2_ggplot.png"),
  plot = combined,
  width = 12,
  height = 8,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "combined_2x2_ggplot.pdf"),
  plot = combined,
  width = 12,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

print(combined)

message("Saved:")
message(file.path(fig_dir, "combined_2x2_ggplot.png"))
message(file.path(fig_dir, "combined_2x2_ggplot.pdf"))
