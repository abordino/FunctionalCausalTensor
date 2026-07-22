setwd("~/Desktop/code")

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(grid)

res_dir = "Results"
plot_dir = "Plots"
runtime_adjustment_power = 1

runtime_summary = read.csv(
  file.path(
    res_dir,
    "ATE_runtime_summary_quadraticPsi_vs_linearReducedAnchor.csv"
  )
)

accuracy_summary = read.csv(
  file.path(
    res_dir,
    "ATE_accuracy_summary_quadraticPsi_vs_linearReducedAnchor.csv"
  )
)

runtime_summary = runtime_summary %>%
  mutate(
    mean_runtime_adjusted = mean_runtime_sec / o_k^runtime_adjustment_power,
    se_runtime_adjusted = se_runtime_sec / o_k^runtime_adjustment_power,
    method = factor(
      method,
      levels = c("Quadratic Psi", "Linear reduced-anchor")
    )
  )

accuracy_summary = accuracy_summary %>%
  mutate(
    method = factor(
      method,
      levels = c("Quadratic Psi", "Linear reduced-anchor")
    )
  )

runtime_adjustment_label = if (runtime_adjustment_power == 1) {
  expression("Runtime / " * o[k])
} else if (runtime_adjustment_power == 2) {
  expression("Runtime / " * o[k]^2)
} else {
  bquote("Runtime / " * o[k]^.(runtime_adjustment_power))
}

runtime_adjustment_subtitle = if (runtime_adjustment_power == 1) {
  expression("Average runtime divided by " * o[k])
} else if (runtime_adjustment_power == 2) {
  expression("Average runtime divided by " * o[k]^2)
} else {
  bquote("Average runtime divided by " * o[k]^.(runtime_adjustment_power))
}

method_cols = c(
  "Quadratic Psi" = "steelblue",
  "Linear reduced-anchor" = "orange"
)

method_shapes = c(
  "Quadratic Psi" = 17,
  "Linear reduced-anchor" = 16
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
    legend.spacing.y = unit(0.1, "lines"),
    plot.margin = margin(6, 6, 6, 6),
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.02, 0.98)
  )

p_runtime = ggplot(
  runtime_summary,
  aes(
    x = o_k,
    y = mean_runtime_adjusted,
    color = method,
    shape = method,
    group = method
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_runtime_adjusted - se_runtime_adjusted),
      ymax = mean_runtime_adjusted + se_runtime_adjusted
    ),
    width = 0,
    linewidth = 0.35,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.1) +
  scale_color_manual(values = method_cols) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = sort(unique(runtime_summary$o_k))) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Adjusted runtime comparison",
    subtitle = runtime_adjustment_subtitle,
    x = expression(o[k]),
    y = runtime_adjustment_label,
    tag = "(a)"
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
  base_theme

p_accuracy = ggplot(
  accuracy_summary,
  aes(
    x = o_k,
    y = mean_abs_error,
    color = method,
    shape = method,
    group = method
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_abs_error - se_abs_error),
      ymax = mean_abs_error + se_abs_error
    ),
    width = 0,
    linewidth = 0.35,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.75, show.legend = FALSE) +
  geom_point(size = 2.1, show.legend = FALSE) +
  scale_color_manual(values = method_cols) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = sort(unique(accuracy_summary$o_k))) +
  scale_y_continuous(labels = scientific) +
  labs(
    title = "Statistical accuracy",
    subtitle = expression("Mean absolute error " %+-% " standard error"),
    x = expression(o[k]),
    y = "Mean absolute error",
    tag = "(b)"
  ) +
  base_theme +
  theme(legend.position = "none")

combined_fig = p_runtime + p_accuracy +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.spacing.y = unit(0.1, "lines"),
    legend.spacing.x = unit(0.8, "lines")
  )

figure_path = file.path(
  plot_dir,
  "ATE_runtime_accuracy_adjusted_1x2.png"
)

ggsave(
  filename = figure_path,
  plot = combined_fig,
  width = 10.5,
  height = 4.4,
  dpi = 400,
  bg = "white"
)

print(combined_fig)
message("Saved: ", figure_path)
