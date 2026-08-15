# Setup -----------------------------------------------------------------------

setwd("~/Desktop/code")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

save_plot_png = function(plot, filename, width, height, dpi = 320) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(filename)
}

# Settings --------------------------------------------------------------------

JUMP = 80L
STEP_SIZE = 3L # 3L
REP = 30L

RUN_CY = TRUE

N_UNITS = 50L
N_PERIODS = 11L
TUCKER_RANK_UNIT = 2L
TUCKER_RANK_TIME = 2L

target_name = "layer2"

PLOT_GRID = crossing(
  N_LAYERS = 4L,
  rank_value = c(1L, 2L, 3L)
)

cy_file_tag = if (RUN_CY) "" else "_noCY"
rep_file_tag = if (REP == 1L) "" else paste0("_rep", REP)

method_levels = c(
  # "Tensor: all layers masked",
  "Tensor: target layer only",
  "Matrix"
)

if (RUN_CY) {
  method_levels = c(
    method_levels,
    "CY biased",
    "CY debiased"
  )
}

method_colors = c(
  `Tensor: all layers masked` = "#1f77b4",
  `Tensor: target layer only` = "#d62728",
  Matrix = "#ff7f0e",
  `CY biased` = "#2ca02c",
  `CY debiased` = "#9467bd"
)

method_shapes = c(
  `Tensor: all layers masked` = 16,
  `Tensor: target layer only` = 17,
  Matrix = 15,
  `CY biased` = 18,
  `CY debiased` = 8
)

# Summary plot ----------------------------------------------------------------

plot_one_result = function(N_LAYERS, rank_value) {
  N_LAYERS = as.integer(N_LAYERS)
  rank_value = as.integer(rank_value)
  
  DATASET_TAG = sprintf(
    "synthetic_tucker2_%dx%dx%d",
    N_UNITS,
    N_PERIODS,
    N_LAYERS
  )
  
  RESULTS_STUDY_TAG = paste0(
    DATASET_TAG,
    "_staggered_masking_signalrank",
    TUCKER_RANK_UNIT,
    "x",
    TUCKER_RANK_TIME
  )
  
  results_dir = file.path(
    "Results",
    RESULTS_STUDY_TAG
  )
  
  plots_dir = file.path(
    "Plots",
    RESULTS_STUDY_TAG
  )
  
  dir.create(
    plots_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  result_file = file.path(
    results_dir,
    paste0(
      DATASET_TAG,
      "_staggered_masking_results_",
      target_name,
      "_rank",
      rank_value,
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      ".rds"
    )
  )
  
  if (!file.exists(result_file)) {
    stop("Result file not found: ", result_file)
  }
  
  results = readRDS(result_file)
  
  if (results$metadata$REP != REP) {
    stop("REP does not match the saved results.")
  }
  
  target_label = results$metadata$target_label
  tensor_dimensions =
    as.integer(results$metadata$tensor_dimensions)
  
  # plot_data = results$average_squared_error_by_step %>%
  #   filter(n_artificial_masked >= JUMP) %>%
  #   mutate(
  plot_data = results$average_squared_error_by_step %>%
        filter(
          n_artificial_masked >= JUMP,
          method != "Tensor: all layers masked"
        ) %>%
        mutate(
      method = factor(method, levels = method_levels),
      ribbon_lower = pmax(
        0,
        average_squared_error -
          coalesce(se_squared_error, 0)
      ),
      ribbon_upper =
        average_squared_error +
        coalesce(se_squared_error, 0)
    ) %>%
    arrange(method, n_artificial_masked)
  
  if (nrow(plot_data) == 0L) {
    stop("No summary data are available.")
  }
  
  first_mask_fraction = min(
    plot_data$masking_fraction,
    na.rm = TRUE
  )
  
  last_mask_fraction = max(
    plot_data$masking_fraction,
    na.rm = TRUE
  )
  
  summary_figure = ggplot(
    plot_data,
    aes(
      x = masking_fraction,
      y = average_squared_error,
      color = method,
      fill = method,
      shape = method,
      group = method
    )
  ) +
    geom_ribbon(
      aes(
        ymin = ribbon_lower,
        ymax = ribbon_upper
      ),
      alpha = 0.12,
      color = NA,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2.1,
      na.rm = TRUE
    ) +
    scale_x_continuous(
      limits = c(
        first_mask_fraction,
        last_mask_fraction
      ),
      breaks = breaks_extended(n = 6),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 0.1),
      breaks = breaks_extended(n = 5),
      labels = label_number(accuracy = 0.0001)
    ) +
    scale_color_manual(
      values = method_colors[method_levels],
      breaks = method_levels,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = method_colors[method_levels],
      breaks = method_levels,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = method_shapes[method_levels],
      breaks = method_levels,
      drop = FALSE
    ) +
    guides(
      fill = "none",
      shape = "none"
    ) +
    labs(
      title = paste0(
        "Staggered masking error: ",
        target_label
      ),
      subtitle = paste0(
        REP,
        " replicates; shaded bands show mean +/- 1 SE; estimation rank = ",
        rank_value
      ),
      x = "Share of terminal artificial mask",
      y = "Mean squared error",
      color = "Method"
    ) +
    theme_minimal(base_size = 22) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  plot_file = file.path(
    plots_dir,
    paste0(
      DATASET_TAG,
      "_staggered_masking_summary_",
      target_name,
      "_rank",
      rank_value,
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      ".png"
    )
  )
  
  save_plot_png(
    plot = summary_figure,
    filename = plot_file,
    width = 11,
    height = 7,
    dpi = 320
  )
  
  tibble(
    N_LAYERS = N_LAYERS,
    rank_value = rank_value,
    REP = REP,
    result_file = result_file,
    plot_file = plot_file,
    first_mask_fraction = first_mask_fraction,
    last_mask_fraction = last_mask_fraction,
    tensor_dimensions = paste(
      tensor_dimensions,
      collapse = " x "
    )
  )
}

# Generate summary plots ------------------------------------------------------

plot_summary = pmap_dfr(
  PLOT_GRID,
  ~ plot_one_result(
    N_LAYERS = ..1,
    rank_value = ..2
  )
)

dir.create(
  "Plots",
  recursive = TRUE,
  showWarnings = FALSE
)