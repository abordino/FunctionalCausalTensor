setwd("~/Desktop/code")

source("bilinearTensorAllFunction.R")

library(ggplot2)
library(dplyr)
library(scales)
library(grid)

results_dir = "Results"
plots_dir = "Plots"

required_files = c(
  "tau_sensitivity_results.rds",
  "tau_sensitivity_config.csv"
)

missing_files = required_files[
  !file.exists(file.path(results_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing result files: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

results = readRDS(file.path(results_dir, "tau_sensitivity_results.rds"))
config = read.csv(file.path(results_dir, "tau_sensitivity_config.csv"))

diagnostics = results$diagnostics
summary_data = results$summary
targets = results$targets

targets = targets %>%
  mutate(
    target_id = row_number(),
    target_label = paste0(
      "B", target_id,
      " (", a - 1, ",", b, ")"
    )
  )

summary_data = summary_data %>%
  left_join(
    targets %>% select(target_id, target_label),
    by = "target_id"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Tensor-pooled", "Matrix-only")
    ),
    target_label = factor(
      target_label,
      levels = targets$target_label
    )
  )

simulate_missingness = function(
    N,
    Tt,
    K,
    r,
    sigma,
    n_adopt_times,
    p_never,
    p_initial,
    seed
) {
  set.seed(seed)
  
  rand_orth = function(n, rank) {
    qr.Q(qr(matrix(rnorm(n * rank), nrow = n)))[
      , seq_len(rank), drop = FALSE
    ]
  }
  
  U = rand_orth(N, r)
  V = rand_orth(Tt, r)
  C_core = array(rnorm(r * r * K), dim = c(r, r, K))
  
  M = array(0, dim = c(N, Tt, K))
  for (k in seq_len(K)) {
    M[, , k] = U %*% C_core[, , k] %*% t(V)
  }
  
  invisible(array(rnorm(N * Tt * K, sd = sigma), dim = c(N, Tt, K)))
  
  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(FALSE, dim = c(N, Tt, K))
  first_adopt_time = max(2, floor(p_initial * Tt) + 1)
  
  for (k in seq_len(K)) {
    adopt_grid = unique(round(seq(
      first_adopt_time,
      Tt,
      length.out = n_adopt_times[k]
    )))
    
    for (i in seq_len(N)) {
      A[i, k] = if (runif(1) < p_never) Inf else sample(adopt_grid, 1)
      
      if (is.infinite(A[i, k])) {
        Omega[i, , k] = TRUE
      } else {
        Omega[i, seq_len(A[i, k] - 1), k] = TRUE
      }
    }
  }
  
  list(Omega = Omega, A = A)
}

get_target_blocks = function(sim, k) {
  Tt = dim(sim$Omega)[2]
  obs_len = ifelse(is.infinite(sim$A[, k]), Tt, sim$A[, k] - 1)
  row_perm = order(-obs_len, seq_along(obs_len))
  obs_len_perm = obs_len[row_perm]
  m_desc = unique(obs_len_perm)
  
  row_parts = lapply(m_desc, function(m) which(obs_len_perm == m))
  T_part = diff(c(0, rev(m_desc)))
  ends = cumsum(T_part)
  starts = c(1, head(ends, -1) + 1)
  
  list(
    row_perm = row_perm,
    row_parts = row_parts,
    col_parts = Map(seq, starts, ends),
    N_part = vapply(row_parts, length, integer(1)),
    T_part = T_part
  )
}

save_ggplot = function(plot, name, width, height) {
  ggsave(
    file.path(plots_dir, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white"
  )
  
  ggsave(
    file.path(plots_dir, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

save_base_plot = function(draw, name, width, height) {
  png(
    file.path(plots_dir, paste0(name, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 400
  )
  draw()
  dev.off()
  
  pdf(
    file.path(plots_dir, paste0(name, ".pdf")),
    width = width,
    height = height
  )
  draw()
  dev.off()
}

make_tau_plot = function(summary_data, add_standard_error = TRUE) {
  method_cols = c(
    "Tensor-pooled" = "steelblue",
    "Matrix-only" = "orange"
  )
  
  block_shapes = setNames(
    c(16, 17, 15, 18, 8, 4),
    levels(summary_data$target_label)
  )
  
  tau_breaks = 10^seq(
    floor(log10(min(summary_data$tau))),
    ceiling(log10(max(summary_data$tau))),
    by = 1
  )
  
  p = ggplot(
    summary_data,
    aes(
      x = tau,
      y = mean_abs_error,
      color = method,
      shape = target_label,
      group = interaction(method, target_label)
    )
  )
  
  if (add_standard_error) {
    p = p +
      geom_errorbar(
        aes(
          ymin = pmax(0, mean_abs_error - se_abs_error),
          ymax = mean_abs_error + se_abs_error
        ),
        width = 0,
        linewidth = 0.30,
        alpha = 0.35,
        show.legend = FALSE
      )
  }
  
  p +
    geom_line(linewidth = 0.78) +
    geom_point(size = 2.15) +
    scale_color_manual(
      name = "Method",
      values = method_cols,
      drop = FALSE
    ) +
    scale_shape_manual(
      name = "Missing block",
      values = block_shapes,
      drop = FALSE
    ) +
    scale_x_log10(
      breaks = tau_breaks,
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(labels = scientific) +
    labs(
      title = "Sensitivity to the regularisation parameter",
      subtitle = if (add_standard_error) {
        "Mean absolute error ± standard error across queries"
      } else {
        "Mean absolute error across queries"
      },
      x = expression(tau),
      y = "Mean absolute error"
    ) +
    guides(
      color = guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(shape = 16, linewidth = 0.9, size = 2.4)
      ),
      shape = guide_legend(
        order = 2,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(color = "gray25", linewidth = 0, size = 2.8)
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "gray88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.key.width = unit(1.25, "lines"),
      legend.spacing.x = unit(0.7, "lines"),
      legend.spacing.y = unit(0.15, "lines"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

make_missingness_plots = function(sim, targets, target_k = 1) {
  K = dim(sim$Omega)[3]
  blocks = get_target_blocks(sim, target_k)
  row_perm = blocks$row_perm
  layers = unique(c(target_k, if (K >= 2) 2 else integer(0), K))
  
  target_rows = targets[targets$k == target_k, , drop = FALSE]
  
  target_rects = bind_rows(lapply(seq_len(nrow(target_rows)), function(i) {
    a = target_rows$a[i]
    b = target_rows$b[i]
    rows = blocks$row_parts[[a]]
    cols = blocks$col_parts[[b]]
    y_rows = dim(sim$Omega)[1] - rows + 1
    
    data.frame(
      xmin = min(cols) - 0.5,
      xmax = max(cols) + 0.5,
      ymin = min(y_rows) - 0.5,
      ymax = max(y_rows) + 0.5,
      x = mean(range(cols)),
      y = mean(range(y_rows)),
      label = paste0("B", i, "\n(", a - 1, ",", b, ")")
    )
  }))
  
  row_cuts = cumsum(blocks$N_part)
  col_cuts = cumsum(blocks$T_part)
  
  lapply(layers, function(layer) {
    Omega = sim$Omega[row_perm, , layer]
    N = nrow(Omega)
    Tt = ncol(Omega)
    
    plot_data = expand.grid(
      row_order = seq_len(N),
      time = seq_len(Tt),
      KEEP.OUT.ATTRS = FALSE
    )
    plot_data$status = factor(
      ifelse(as.vector(Omega), "Observed", "Missing"),
      levels = c("Missing", "Observed")
    )
    plot_data$unit_plot = N - plot_data$row_order + 1
    
    unit_labels = pretty(seq_len(N))
    unit_labels = unit_labels[unit_labels >= 1 & unit_labels <= N]
    
    p = ggplot(plot_data, aes(time, unit_plot, fill = status)) +
      geom_raster() +
      scale_fill_manual(
        values = c("Missing" = "red", "Observed" = "blue"),
        guide = "none"
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(
        breaks = N - unit_labels + 1,
        labels = unit_labels,
        expand = c(0, 0)
      ) +
      coord_cartesian(
        xlim = c(0.5, Tt + 0.5),
        ylim = c(0.5, N + 0.5),
        expand = FALSE
      ) +
      labs(
        title = if (layer == target_k) {
          paste0("Target layer ", layer)
        } else {
          paste0("Layer ", layer)
        },
        subtitle = if (layer == target_k) {
          "Blue = observed; red = missing; yellow = evaluated blocks"
        } else {
          NULL
        },
        x = "Time",
        y = "Unit"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.5, hjust = 0.5),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    if (layer == target_k) {
      p = p +
        geom_vline(xintercept = col_cuts + 0.5, linewidth = 0.32) +
        geom_hline(yintercept = N - row_cuts + 0.5, linewidth = 0.32) +
        geom_rect(
          data = target_rects,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = NA,
          color = "yellow",
          linewidth = 0.85
        ) +
        geom_text(
          data = target_rects,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          color = "yellow",
          size = 2.7,
          fontface = "bold",
          lineheight = 0.88
        )
    }
    
    p
  })
}

draw_plot_row = function(plots) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, length(plots))))
  
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      newpage = FALSE,
      vp = viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }
  
  popViewport()
}

draw_combined_plot = function(missingness_plots, tau_plot) {
  n_top = length(missingness_plots)
  
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = n_top,
        heights = unit(c(1, 1.65), "null")
      )
    )
  )
  
  for (i in seq_len(n_top)) {
    print(
      missingness_plots[[i]],
      newpage = FALSE,
      vp = viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }
  
  print(
    tau_plot,
    newpage = FALSE,
    vp = viewport(layout.pos.row = 2, layout.pos.col = seq_len(n_top))
  )
  
  popViewport()
}

plot_missingness_layer = function(sim, layer, row_perm = NULL, main = NULL) {
  Omega = sim$Omega[, , layer]
  if (!is.null(row_perm)) {
    Omega = Omega[row_perm, , drop = FALSE]
  }
  
  N = nrow(Omega)
  Tt = ncol(Omega)
  
  image(
    x = seq_len(Tt),
    y = seq_len(N),
    z = t(matrix(as.numeric(Omega), N, Tt)[N:1, ]),
    col = c("red", "blue"),
    axes = FALSE,
    xlab = "Time index t",
    ylab = "Unit index i",
    main = main
  )
  axis(1)
  axis(2, at = pretty(seq_len(N)), labels = rev(pretty(seq_len(N))))
  box()
}

draw_before_after = function(sim, target_k = 1) {
  K = dim(sim$Omega)[3]
  row_perm = get_target_blocks(sim, target_k)$row_perm
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, K), mar = c(4, 4, 3, 1))
  
  for (layer in seq_len(K)) {
    plot_missingness_layer(
      sim,
      layer,
      main = paste0("Before, layer ", layer)
    )
  }
  
  for (layer in seq_len(K)) {
    plot_missingness_layer(
      sim,
      layer,
      row_perm = row_perm,
      main = paste0("After target k=", target_k, ", layer ", layer)
    )
  }
}

make_accuracy_diagnostics = function(diagnostics, summary_data, target_id = 1) {
  selected_tau = summary_data %>%
    filter(.data$target_id == target_id) %>%
    group_by(tau) %>%
    summarise(mean_abs_error = mean(mean_abs_error), .groups = "drop") %>%
    slice_min(mean_abs_error, n = 1, with_ties = FALSE) %>%
    pull(tau)
  
  selected = diagnostics %>%
    filter(
      .data$target_id == target_id,
      near(tau, selected_tau)
    )
  
  tensor = selected %>%
    filter(method == "Tensor-pooled") %>%
    select(query, mu_true, mu_hat_tensor = mu_hat, abs_error_tensor = abs_error)
  
  matrix = selected %>%
    filter(method == "Matrix-only") %>%
    select(query, mu_hat_matrix = mu_hat, abs_error_matrix = abs_error)
  
  list(
    tau = selected_tau,
    data = inner_join(tensor, matrix, by = "query")
  )
}

draw_accuracy_diagnostics = function(accuracy) {
  data = accuracy$data
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  lims = range(c(data$mu_true, data$mu_hat_tensor, data$mu_hat_matrix))
  
  plot(
    data$mu_true,
    data$mu_hat_tensor,
    xlab = expression(mu[true]),
    ylab = expression(mu[hat]),
    main = paste0("Estimates at tau = ", signif(accuracy$tau, 4)),
    pch = 19,
    col = "steelblue",
    xlim = lims,
    ylim = lims
  )
  points(data$mu_true, data$mu_hat_matrix, pch = 17, col = "orange")
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  legend(
    "topleft",
    legend = c("tensor-pooled", "matrix-only", "45-degree line"),
    col = c("steelblue", "orange", "red"),
    pch = c(19, 17, NA),
    lty = c(NA, NA, 2),
    lwd = c(NA, NA, 2),
    bty = "n"
  )
  
  ymax = max(data$abs_error_tensor, data$abs_error_matrix)
  plot(
    data$query,
    data$abs_error_tensor,
    type = "b",
    pch = 19,
    col = "steelblue",
    ylim = c(0, ymax),
    xlab = "Query index",
    ylab = "Absolute error",
    main = "Absolute error by query"
  )
  lines(
    data$query,
    data$abs_error_matrix,
    type = "b",
    pch = 17,
    col = "orange"
  )
  abline(h = mean(data$abs_error_tensor), col = "steelblue", lwd = 2, lty = 2)
  abline(h = mean(data$abs_error_matrix), col = "orange", lwd = 2, lty = 2)
  
  boxplot(
    data$abs_error_tensor,
    data$abs_error_matrix,
    names = c("tensor", "matrix"),
    col = c("steelblue", "orange"),
    ylab = "Absolute error",
    main = "Absolute error comparison"
  )
  
  relative_tensor = data$abs_error_tensor / pmax(abs(data$mu_true), 1e-12)
  relative_matrix = data$abs_error_matrix / pmax(abs(data$mu_true), 1e-12)
  
  plot(
    abs(data$mu_true),
    relative_tensor,
    pch = 19,
    col = "steelblue",
    ylim = c(0, max(relative_tensor, relative_matrix)),
    xlab = expression(abs(mu[true])),
    ylab = "Relative error",
    main = "Relative error vs signal size"
  )
  points(abs(data$mu_true), relative_matrix, pch = 17, col = "orange")
  legend(
    "topright",
    legend = c("tensor-pooled", "matrix-only"),
    col = c("steelblue", "orange"),
    pch = c(19, 17),
    bty = "n"
  )
}

n_adopt_times = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5)

sim = simulate_missingness(
  N = config$N[1],
  Tt = config$Tt[1],
  K = config$K[1],
  r = config$r[1],
  sigma = config$sigma[1],
  n_adopt_times = n_adopt_times,
  p_never = config$p_never[1],
  p_initial = config$p_initial[1],
  seed = config$experiment_seed[1]
)

tau_plot = make_tau_plot(summary_data)
missingness_plots = make_missingness_plots(sim, targets, target_k = 1)
accuracy = make_accuracy_diagnostics(diagnostics, summary_data, target_id = 1)

save_ggplot(
  tau_plot,
  "tau_sensitivity_six_missing_blocks_12_lines",
  width = 10.5,
  height = 6.2
)

save_base_plot(
  function() draw_plot_row(missingness_plots),
  "missingness_panels",
  width = 12,
  height = 4.2
)

# save_base_plot(
#   function() draw_combined_plot(missingness_plots, tau_plot),
#   "missingness_and_tau_sensitivity",
#   width = 12,
#   height = 9
# )
# 
# save_base_plot(
#   function() draw_before_after(sim, target_k = 1),
#   "missingness_before_after_rearrangement",
#   width = 24,
#   height = 8
# )
# 
# save_base_plot(
#   function() draw_accuracy_diagnostics(accuracy),
#   "accuracy_diagnostics_target_1_best_tau",
#   width = 10,
#   height = 8
# )

message("Plots saved to: ", normalizePath(plots_dir))



# Show the final plot
# print(missingness_plots)
print(tau_plot)
