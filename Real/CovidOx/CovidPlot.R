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

JUMP = 30L
STEP_SIZE = 1L # 3L
REP = 30L 

RUN_CY = FALSE

N_LAYERS = 2L

target_name = "deaths"

DATASET_TAG = "oxford_deaths"
RESULTS_STUDY_TAG = "oxford_deaths_target_layer_bootstrap"

tensor_file =
  "Real/CovidOx/data/Omega_Y_until_2020-04-05_delay_28.rds"

`%||%` = function(x, y) {
  if (is.null(x)) y else x
}

tensor_data = readRDS(tensor_file)

slice_label = paste0(
  "until_",
  tensor_data$time_horizon %||% "unknown",
  "_delay_",
  tensor_data$delay_days %||% "unknown"
)

slice_label = stringr::str_replace_all(
  slice_label,
  "[^A-Za-z0-9._-]+",
  "_"
)

slice_file_tag = paste0("_", slice_label)

PLOT_GRID = crossing(
  N_LAYERS = N_LAYERS,
  rank_value = c(1L, 2L)
)

cy_file_tag = if (RUN_CY) "" else "_noCY"
rep_file_tag = if (REP == 1L) "" else paste0("_rep", REP)

method_levels = c(
  "Tensor: all layers masked",
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

plot_method_levels = setdiff(
  method_levels,
  "Tensor: all layers masked"
)

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

rank_line_types = c(
  `Rank 3` = "solid",
  `Rank 4` = "longdash"
)

rank_shapes = c(
  `Rank 3` = 17,
  `Rank 4` = 17
)

rank_colors = c(
  `Rank 3` = "#d62728",
  `Rank 4` = "#e377c2"
)

# Summary plot: ranks 1 and 2 -------------------------------------------------

plot_one_result = function(N_LAYERS, rank_value) {
  N_LAYERS = as.integer(N_LAYERS)
  rank_value = as.integer(rank_value)
  
  results_dir = file.path(
    "Results",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
  )
  
  plots_dir = file.path(
    "Plots",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
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
      slice_file_tag,
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
  
  if (results$metadata$JUMP != JUMP) {
    stop("JUMP does not match the saved results.")
  }
  
  if (results$metadata$STEP_SIZE != STEP_SIZE) {
    stop("STEP_SIZE does not match the saved results.")
  }
  
  if (results$metadata$N_LAYERS != N_LAYERS) {
    stop("N_LAYERS does not match the saved results.")
  }
  
  if (results$metadata$RUN_CY != RUN_CY) {
    stop("RUN_CY does not match the saved results.")
  }
  
  if (results$metadata$slice_label != slice_label) {
    stop("The Oxford data slice does not match the saved results.")
  }
  
  if (
    normalizePath(
      results$metadata$tensor_file,
      winslash = "/",
      mustWork = FALSE
    ) !=
    normalizePath(
      tensor_file,
      winslash = "/",
      mustWork = FALSE
    )
  ) {
    stop("tensor_file does not match the saved results.")
  }
  
  target_label = results$metadata$target_label
  tensor_dimensions =
    as.integer(results$metadata$tensor_dimensions)
  
  plot_data = results$average_squared_error_by_step %>%
    filter(
      n_artificial_masked >= JUMP,
      method != "Tensor: all layers masked"
    ) %>%
    mutate(
      method = factor(
        method,
        levels = plot_method_levels
      ),
      # ribbon_lower = pmax(
      #   0,
      #   average_squared_error -
      #     coalesce(se_squared_error, 0)
      # ),
      ribbon_lower = pmax(
        1e-10,
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
    # scale_y_continuous(
    #   limits = c(0, 20),
    #   breaks = breaks_extended(n = 5),
    #   labels = label_number(accuracy = 0.0001)
    # ) +
    scale_y_log10(
      breaks = breaks_log(n = 6),
      labels = label_number()
    ) +
    coord_cartesian(
      ylim = c(1e-2, 20)
    ) +
    scale_color_manual(
      values = method_colors[plot_method_levels],
      breaks = plot_method_levels,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = method_colors[plot_method_levels],
      breaks = plot_method_levels,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = method_shapes[plot_method_levels],
      breaks = plot_method_levels,
      drop = FALSE
    ) +
    guides(
      fill = "none",
      shape = "none"
    ) +
    labs(
      x = "Share of terminal artificial mask",
      y = "Mean squared error (log scale)",
      color = "Method"
    ) +
    theme_minimal(base_size = 20) +
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
      slice_file_tag,
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

# Summary plot: ranks 3 and 4 -------------------------------------------------

plot_rank_3_4_tensor = function(N_LAYERS) {
  N_LAYERS = as.integer(N_LAYERS)
  
  results_dir = file.path(
    "Results",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
  )
  
  plots_dir = file.path(
    "Plots",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
  )
  
  dir.create(
    plots_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  rank_values = c(3L, 4L)
  
  rank_results = map(
    rank_values,
    function(rank_value) {
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
          slice_file_tag,
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
      
      if (results$metadata$JUMP != JUMP) {
        stop("JUMP does not match the saved results.")
      }
      
      if (results$metadata$STEP_SIZE != STEP_SIZE) {
        stop("STEP_SIZE does not match the saved results.")
      }
      
      if (results$metadata$N_LAYERS != N_LAYERS) {
        stop("N_LAYERS does not match the saved results.")
      }
      
      if (results$metadata$RUN_CY != RUN_CY) {
        stop("RUN_CY does not match the saved results.")
      }
      
      if (results$metadata$slice_label != slice_label) {
        stop("The Oxford data slice does not match the saved results.")
      }
      
      if (
        normalizePath(
          results$metadata$tensor_file,
          winslash = "/",
          mustWork = FALSE
        ) !=
        normalizePath(
          tensor_file,
          winslash = "/",
          mustWork = FALSE
        )
      ) {
        stop("tensor_file does not match the saved results.")
      }
      
      plot_data = results$average_squared_error_by_step %>%
        filter(
          n_artificial_masked >= JUMP,
          method == "Tensor: target layer only"
        ) %>%
        mutate(
          rank_value = rank_value,
          rank_label = paste0("Rank ", rank_value),
          # ribbon_lower = pmax(
          #   0,
          #   average_squared_error -
          #     coalesce(se_squared_error, 0)
          # ),
          ribbon_lower = pmax(
            1e-10,
            average_squared_error -
              coalesce(se_squared_error, 0)
          ),
          ribbon_upper =
            average_squared_error +
            coalesce(se_squared_error, 0)
        ) %>%
        arrange(n_artificial_masked)
      
      list(
        results = results,
        result_file = result_file,
        plot_data = plot_data
      )
    }
  )
  
  target_labels = map_chr(
    rank_results,
    ~ .x$results$metadata$target_label
  )
  
  if (length(unique(target_labels)) != 1L) {
    stop("Target labels differ across ranks 3 and 4.")
  }
  
  target_label = target_labels[[1]]
  
  tensor_dimensions = as.integer(
    rank_results[[1]]$results$metadata$tensor_dimensions
  )
  
  plot_data = map_dfr(
    rank_results,
    "plot_data"
  ) %>%
    mutate(
      rank_label = factor(
        rank_label,
        levels = c("Rank 3", "Rank 4")
      )
    )
  
  if (nrow(plot_data) == 0L) {
    stop("No rank 3/4 target-layer Tensor data are available.")
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
      color = rank_label,
      fill = rank_label,
      linetype = rank_label,
      shape = rank_label,
      group = rank_label
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
    scale_color_manual(
      values = rank_colors,
      breaks = c("Rank 3", "Rank 4"),
      drop = FALSE
    ) +
    scale_fill_manual(
      values = rank_colors,
      breaks = c("Rank 3", "Rank 4"),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = rank_line_types,
      breaks = c("Rank 3", "Rank 4"),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = rank_shapes,
      breaks = c("Rank 3", "Rank 4"),
      drop = FALSE
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
    # scale_y_continuous(
    #   limits = c(0, 20),
    #   breaks = breaks_extended(n = 5),
    #   labels = label_number(accuracy = 0.0001)
    # ) +
    scale_y_log10(
      breaks = breaks_log(n = 6),
      labels = label_number()
    ) +
    coord_cartesian(
      ylim = c(1e-2, 20)
    ) +
    guides(
      fill = "none",
      shape = "none"
    ) +
    labs(
      x = "Share of terminal artificial mask",
      y = "Mean squared error",
      color = "Estimation rank",
      linetype = "Estimation rank"
    ) +
    theme_minimal(base_size = 20) +
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
      "_rank3_rank4_target_tensor",
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      slice_file_tag,
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
    rank_value = "3 & 4",
    REP = REP,
    result_file = paste(
      map_chr(rank_results, "result_file"),
      collapse = "; "
    ),
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

plot_summary_rank_1_2 = pmap_dfr(
  PLOT_GRID,
  ~ plot_one_result(
    N_LAYERS = ..1,
    rank_value = ..2
  )
)

plot_summary_rank_3_4 = plot_rank_3_4_tensor(
  N_LAYERS = N_LAYERS
)

plot_summary = bind_rows(
  plot_summary_rank_1_2 %>%
    mutate(rank_value = as.character(rank_value)),
  plot_summary_rank_3_4
)

plots_dir = file.path(
  "Plots",
  RESULTS_STUDY_TAG,
  paste0("step", STEP_SIZE)
)

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_summary_file = file.path(
  plots_dir,
  paste0(
    DATASET_TAG,
    "_staggered_masking_plot_summary_",
    target_name,
    "_jump",
    JUMP,
    "_step",
    STEP_SIZE,
    rep_file_tag,
    cy_file_tag,
    slice_file_tag,
    ".csv"
  )
)

print(
  plot_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nSaved plot summary:\n",
  plot_summary_file,
  "\n",
  sep = ""
)
