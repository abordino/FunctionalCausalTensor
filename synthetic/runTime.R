setwd("~/Documents/phd/projects/causalMatrix/code")

source("bilinearTensorStaggered.R")
source("bilinearTensorStaggeredPsi.R")
source("LinearAggregate.R")

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(grid)

res_dir = "result"
fig_dir = "figure"

if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
}

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}


# ===============================================================
# 1. Simulator with fixed N, fixed Tt, varying o_k
# ===============================================================

simulate_staggered_tucker2_fixed_ok = function(
    N = 150,
    Tt = 200,
    K = 5,
    r = 5,
    sigma = 0,
    o_k_target = 6,
    target_k = 1,
    p_initial = 0.10,
    seed = 123
) {
  set.seed(seed)
  
  stopifnot(r <= min(N, Tt))
  stopifnot(o_k_target >= 2)
  stopifnot(o_k_target <= min(N, Tt))
  stopifnot(target_k >= 1, target_k <= K)
  
  rand_orth = function(n, r) {
    qr.Q(qr(matrix(rnorm(n * r), nrow = n, ncol = r)))[, seq_len(r), drop = FALSE]
  }
  
  U = rand_orth(N, r)
  V = rand_orth(Tt, r)
  
  C_core = array(
    rnorm(r * r * K),
    dim = c(r, r, K)
  )
  
  M = array(0, dim = c(N, Tt, K))
  
  for (kk in seq_len(K)) {
    M[, , kk] = U %*% C_core[, , kk] %*% t(V)
  }
  
  E = array(
    rnorm(N * Tt * K, mean = 0, sd = sigma),
    dim = c(N, Tt, K)
  )
  
  Y_full = M + E
  
  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(TRUE, dim = c(N, Tt, K))
  Y_obs = Y_full
  
  # -------------------------------------------------------------
  # Target layer: exactly o_k_target row/time blocks.
  # -------------------------------------------------------------
  
  min_obs_len = max(1, floor(p_initial * Tt))
  
  finite_obs_lens = unique(round(seq(
    min_obs_len,
    Tt - 1,
    length.out = o_k_target - 1
  )))
  
  group_labels = rep(seq_len(o_k_target), length.out = N)
  group_labels = sample(group_labels)
  
  A[, target_k] = Inf
  
  for (i in seq_len(N)) {
    g = group_labels[i]
    
    if (g == 1) {
      A[i, target_k] = Inf
    } else {
      obs_len_i = finite_obs_lens[g - 1]
      A[i, target_k] = obs_len_i + 1
    }
  }
  
  Omega[, , target_k] = FALSE
  
  for (i in seq_len(N)) {
    if (is.infinite(A[i, target_k])) {
      Omega[i, , target_k] = TRUE
    } else if (A[i, target_k] > 1) {
      Omega[i, seq_len(A[i, target_k] - 1), target_k] = TRUE
    }
  }
  
  Y_obs[, , target_k][!Omega[, , target_k]] = NA_real_
  
  for (kk in setdiff(seq_len(K), target_k)) {
    A[, kk] = Inf
    Omega[, , kk] = TRUE
    Y_obs[, , kk] = Y_full[, , kk]
  }
  
  out = list(
    Y_obs = Y_obs,
    Y_full = Y_full,
    M = M,
    E = E,
    Omega = Omega,
    A = A,
    U = U,
    V = V,
    C_core = C_core,
    r = r,
    sigma = sigma,
    o_k_target = o_k_target,
    target_k = target_k
  )
  
  class(out) = "fixed_ok_staggered_tucker2_sim"
  
  out
}


# ===============================================================
# 2. Truth and estimator wrappers
# ===============================================================

get_target_blocks_sim = function(sim, k) {
  A_k = sim$A[, k]
  Tt = dim(sim$Y_obs)[2]
  
  obs_len = ifelse(is.infinite(A_k), Tt, pmin(Tt, A_k - 1))
  row_perm = order(-obs_len, seq_along(obs_len))
  
  obs_len_perm = obs_len[row_perm]
  m_desc = unique(obs_len_perm)
  
  o_k = length(m_desc)
  
  row_parts = lapply(m_desc, function(m) which(obs_len_perm == m))
  
  m_asc = rev(m_desc)
  T_part = diff(c(0, m_asc))
  
  ends = cumsum(T_part)
  starts = c(1, head(ends, -1) + 1)
  col_parts = Map(seq, starts, ends)
  
  list(
    row_perm = row_perm,
    obs_len_perm = obs_len_perm,
    o_k = o_k,
    row_parts = row_parts,
    col_parts = col_parts,
    N_part = vapply(row_parts, length, integer(1)),
    T_part = T_part
  )
}


true_missing_region_ATE = function(sim, k = 1) {
  missing = !sim$Omega[, , k]
  
  if (!any(missing)) {
    stop("No missing entries in target layer.")
  }
  
  mean(sim$M[, , k][missing])
}


true_block_formula_ATE = function(sim, k = 1) {
  blocks = get_target_blocks_sim(sim, k)
  
  o_k = blocks$o_k
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      if (a + b <= o_k + 1) next
      
      rows_sorted = blocks$row_parts[[a]]
      rows_original = blocks$row_perm[rows_sorted]
      cols = blocks$col_parts[[b]]
      
      N_a = length(rows_original)
      T_b = length(cols)
      
      x_a = rep(1 / sqrt(N_a), N_a)
      y_b = rep(1 / sqrt(T_b), T_b)
      
      M_ab = sim$M[rows_original, cols, k, drop = FALSE][, , 1]
      
      mu_ab = as.numeric(t(x_a) %*% M_ab %*% y_b)
      
      weighted_sum = weighted_sum + sqrt(N_a * T_b) * mu_ab
      normalizer = normalizer + N_a * T_b
    }
  }
  
  weighted_sum / normalizer
}


estimate_quadratic_Psi_ATE = function(sim, k = 1, tau = 1e-2) {
  bilinearTensorStaggeredATE(
    Y = sim$Y_obs,
    k = k,
    r = sim$r,
    tau = tau,
    A = sim$A,
    Omega = sim$Omega
  )
}


estimate_linear_reduced_anchor_ATE = function(sim, k = 1, tau = 1e-2) {
  bilinearTensorStaggeredATELinearReducedAnchor(
    Y = sim$Y_obs,
    k = k,
    r = sim$r,
    tau = tau,
    A = sim$A,
    Omega = sim$Omega
  )
}


# ===============================================================
# 3. Study 1: noiseless comparison
# ===============================================================

run_noiseless_comparison = function(
    N = 150,
    Tt = 200,
    K = 5,
    r = 5,
    o_k_target = 8,
    k = 1,
    tau = 1e-2,
    seed = 123
) {
  sim = simulate_staggered_tucker2_fixed_ok(
    N = N,
    Tt = Tt,
    K = K,
    r = r,
    sigma = 0,
    o_k_target = o_k_target,
    target_k = k,
    seed = seed
  )
  
  blocks = get_target_blocks_sim(sim, k)
  
  true_ate_direct = true_missing_region_ATE(sim, k)
  true_ate_formula = true_block_formula_ATE(sim, k)
  
  psi_quad = estimate_quadratic_Psi_ATE(sim, k, tau)
  linear_red = estimate_linear_reduced_anchor_ATE(sim, k, tau)
  
  out = data.frame(
    method = c(
      "Truth: direct missing-region mean",
      "Truth: block formula",
      "Quadratic Psi ATE",
      "Linear reduced-anchor ATE"
    ),
    estimate = c(
      true_ate_direct,
      true_ate_formula,
      psi_quad,
      linear_red
    )
  )
  
  out$abs_error_vs_truth = abs(out$estimate - true_ate_direct)
  
  cat("\nNoiseless comparison\n")
  cat("--------------------------------------------------\n")
  cat("Requested o_k:", o_k_target, "\n")
  cat("Realized o_k :", blocks$o_k, "\n")
  cat("sigma        :", sim$sigma, "\n\n")
  
  print(out)
  
  cat("\nKey differences\n")
  cat("--------------------------------------------------\n")
  cat("Block formula - direct truth:",
      true_ate_formula - true_ate_direct, "\n")
  cat("Quadratic Psi - truth       :",
      psi_quad - true_ate_direct, "\n")
  cat("Linear reduced - truth      :",
      linear_red - true_ate_direct, "\n")
  cat("Quadratic Psi - linear red  :",
      psi_quad - linear_red, "\n")
  
  invisible(list(
    sim = sim,
    blocks = blocks,
    results = out
  ))
}


# ===============================================================
# 4. Study 2: runtime comparison
# ===============================================================

run_runtime_comparison = function(
    N = 150,
    Tt = 200,
    K = 5,
    r = 5,
    o_grid = c(4, 6, 8, 10, 12, 15, 20),
    n_runs = 100,
    k = 1,
    tau = 1e-2,
    sigma = 0,
    base_seed = 1000
) {
  results = vector("list", length(o_grid) * n_runs)
  counter = 1
  
  for (oo in o_grid) {
    cat("\nRuntime study: o_k =", oo, "\n")
    
    for (rrun in seq_len(n_runs)) {
      seed_run = base_seed + 10000 * oo + rrun
      
      sim = simulate_staggered_tucker2_fixed_ok(
        N = N,
        Tt = Tt,
        K = K,
        r = r,
        sigma = sigma,
        o_k_target = oo,
        target_k = k,
        seed = seed_run
      )
      
      blocks = get_target_blocks_sim(sim, k)
      
      gc()
      
      time_quad = system.time({
        psi_quad = estimate_quadratic_Psi_ATE(sim, k, tau)
      })[["elapsed"]]
      
      gc()
      
      time_linear = system.time({
        linear_red = estimate_linear_reduced_anchor_ATE(sim, k, tau)
      })[["elapsed"]]
      
      results[[counter]] = data.frame(
        o_k = oo,
        run = rrun,
        method = c("Quadratic Psi", "Linear reduced-anchor"),
        runtime_sec = c(time_quad, time_linear),
        estimate = c(psi_quad, linear_red),
        difference_quad_minus_linear = c(
          psi_quad - linear_red,
          psi_quad - linear_red
        )
      )
      
      counter = counter + 1
      
      if (rrun %% 10 == 0) {
        cat("  completed", rrun, "of", n_runs, "\n")
      }
    }
  }
  
  timings = do.call(rbind, results)
  
  summary = timings %>%
    group_by(o_k, method) %>%
    summarise(
      mean_runtime_sec = mean(runtime_sec),
      sd_runtime_sec = sd(runtime_sec),
      n_runs = n(),
      se_runtime_sec = sd_runtime_sec / sqrt(n_runs),
      .groups = "drop"
    )
  
  list(
    timings = timings,
    summary = summary
  )
}


# ===============================================================
# 5. Study 3: noisy accuracy comparison
# ===============================================================

run_noisy_accuracy_comparison = function(
    N = 150,
    Tt = 200,
    K = 5,
    r = 5,
    o_grid = c(4, 6, 8, 10, 12, 15, 20),
    n_runs = 100,
    k = 1,
    tau = 1e-2,
    sigma = 0.05,
    base_seed = 5000
) {
  results = vector("list", length(o_grid) * n_runs)
  counter = 1
  
  for (oo in o_grid) {
    cat("\nNoisy accuracy study: o_k =", oo, "\n")
    
    for (rrun in seq_len(n_runs)) {
      seed_run = base_seed + 10000 * oo + rrun
      
      sim = simulate_staggered_tucker2_fixed_ok(
        N = N,
        Tt = Tt,
        K = K,
        r = r,
        sigma = sigma,
        o_k_target = oo,
        target_k = k,
        seed = seed_run
      )
      
      blocks = get_target_blocks_sim(sim, k)
      
      true_ate = true_missing_region_ATE(sim, k)
      
      psi_quad = estimate_quadratic_Psi_ATE(sim, k, tau)
      linear_red = estimate_linear_reduced_anchor_ATE(sim, k, tau)
      
      err_quad = psi_quad - true_ate
      err_linear = linear_red - true_ate
      
      results[[counter]] = data.frame(
        o_k = oo,
        run = rrun,
        sigma = sigma,
        method = c("Quadratic Psi", "Linear reduced-anchor"),
        true_ate = true_ate,
        estimate = c(psi_quad, linear_red),
        error = c(err_quad, err_linear),
        abs_error = abs(c(err_quad, err_linear)),
        sq_error = c(err_quad, err_linear)^2
      )
      
      counter = counter + 1
      
      if (rrun %% 10 == 0) {
        cat("  completed", rrun, "of", n_runs, "\n")
      }
    }
  }
  
  diagnostics = do.call(rbind, results)
  
  summary = diagnostics %>%
    group_by(o_k, method) %>%
    summarise(
      mean_abs_error = mean(abs_error),
      sd_abs_error = sd(abs_error),
      n_runs = n(),
      se_abs_error = sd_abs_error / sqrt(n_runs),
      median_abs_error = median(abs_error),
      rmse = sqrt(mean(sq_error)),
      bias = mean(error),
      sd_error = sd(error),
      se_error = sd_error / sqrt(n_runs),
      .groups = "drop"
    )
  
  list(
    diagnostics = diagnostics,
    summary = summary
  )
}


# ===============================================================
# 6. 1x2 plot: runtime and accuracy
# ===============================================================

plot_runtime_accuracy_1x2 = function(runtime_out,
                                     accuracy_out,
                                     filename_prefix = "ATE_runtime_accuracy_1x2",
                                     runtime_adjustment_power = 1) {
  
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
  
  runtime_summary = runtime_out$summary %>%
    mutate(
      mean_runtime_adjusted = mean_runtime_sec / (o_k^runtime_adjustment_power),
      se_runtime_adjusted = se_runtime_sec / (o_k^runtime_adjustment_power),
      method_label = case_when(
        method == "Quadratic Psi" ~ "Quadratic Psi",
        method == "Linear reduced-anchor" ~ "Linear reduced-anchor",
        TRUE ~ method
      ),
      method_label = factor(
        method_label,
        levels = c("Quadratic Psi", "Linear reduced-anchor")
      )
    )
  
  accuracy_summary = accuracy_out$summary %>%
    mutate(
      method_label = case_when(
        method == "Quadratic Psi" ~ "Quadratic Psi",
        method == "Linear reduced-anchor" ~ "Linear reduced-anchor",
        TRUE ~ method
      ),
      method_label = factor(
        method_label,
        levels = c("Quadratic Psi", "Linear reduced-anchor")
      )
    )
  
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
      legend.spacing.y = unit(0.10, "lines"),
      plot.margin = margin(6, 6, 6, 6),
      plot.tag = element_text(face = "bold", size = 13),
      plot.tag.position = c(0.02, 0.98)
    )
  
  # -------------------------------------------------------------
  # Panel (a): slope-adjusted runtime
  # -------------------------------------------------------------
  
  pA = ggplot(
    runtime_summary,
    aes(
      x = o_k,
      y = mean_runtime_adjusted,
      color = method_label,
      shape = method_label,
      group = method_label
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
    geom_line(linewidth = 0.75, show.legend = TRUE) +
    geom_point(size = 2.1, show.legend = TRUE) +
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
    scale_x_continuous(
      breaks = sort(unique(runtime_summary$o_k))
    ) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.001)
    ) +
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
  
  # -------------------------------------------------------------
  # Panel (b): accuracy, mean absolute error +/- standard error
  # -------------------------------------------------------------
  
  pB = ggplot(
    accuracy_summary,
    aes(
      x = o_k,
      y = mean_abs_error,
      color = method_label,
      shape = method_label,
      group = method_label
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
    scale_x_continuous(
      breaks = sort(unique(accuracy_summary$o_k))
    ) +
    scale_y_continuous(
      labels = scientific
    ) +
    labs(
      title = "Statistical accuracy",
      subtitle = expression("Mean absolute error " %+-% " standard error"),
      x = expression(o[k]),
      y = "Mean absolute error",
      tag = "(b)"
    ) +
    base_theme +
    theme(
      legend.position = "none"
    )
  
  combined = pA + pB +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.spacing.y = unit(0.1, "lines"),
      legend.spacing.x = unit(0.8, "lines")
    )
  
  png_path = file.path(fig_dir, paste0(filename_prefix, ".png"))
  
  ggsave(
    filename = png_path,
    plot = combined,
    width = 10.5,
    height = 4.4,
    dpi = 400,
    bg = "white"
  )
  
  print(combined)
  
  message("Saved:")
  message(png_path)
  
  invisible(combined)
}

# ===============================================================
# 7. Run all studies
# ===============================================================

N_main = 150
Tt_main = 200
K_main = 5
r_main = 5
k_target = 1
tau_main = 1e-3

o_grid_main = c(4, 6, 8, 10, 12, 15, 20)

n_runs_runtime = 1000
n_runs_accuracy = 1000

sigma_runtime = 0
sigma_accuracy = 0.01


# ---------------------------------------------------------------
# Study 1: noiseless comparison
# ---------------------------------------------------------------

noiseless_out = run_noiseless_comparison(
  N = N_main,
  Tt = Tt_main,
  K = K_main,
  r = r_main,
  o_k_target = 4,
  k = k_target,
  tau = tau_main,
  seed = 123
)

write.csv(
  noiseless_out$results,
  file = file.path(res_dir, "ATE_noiseless_comparison.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------------
# Study 2: runtime comparison
# ---------------------------------------------------------------

runtime_out = run_runtime_comparison(
  N = N_main,
  Tt = Tt_main,
  K = K_main,
  r = r_main,
  o_grid = o_grid_main,
  n_runs = n_runs_runtime,
  k = k_target,
  tau = tau_main,
  sigma = sigma_runtime,
  base_seed = 1000
)

print(runtime_out$summary)

write.csv(
  runtime_out$timings,
  file = file.path(res_dir, "ATE_runtime_timings_quadraticPsi_vs_linearReducedAnchor.csv"),
  row.names = FALSE
)

write.csv(
  runtime_out$summary,
  file = file.path(res_dir, "ATE_runtime_summary_quadraticPsi_vs_linearReducedAnchor.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------------
# Study 3: noisy accuracy comparison
# ---------------------------------------------------------------

accuracy_out = run_noisy_accuracy_comparison(
  N = N_main,
  Tt = Tt_main,
  K = K_main,
  r = r_main,
  o_grid = o_grid_main,
  n_runs = n_runs_accuracy,
  k = k_target,
  tau = tau_main,
  sigma = sigma_accuracy,
  base_seed = 5000
)

print(accuracy_out$summary)

write.csv(
  accuracy_out$diagnostics,
  file = file.path(res_dir, "ATE_accuracy_diagnostics_quadraticPsi_vs_linearReducedAnchor.csv"),
  row.names = FALSE
)

write.csv(
  accuracy_out$summary,
  file = file.path(res_dir, "ATE_accuracy_summary_quadraticPsi_vs_linearReducedAnchor.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------------
# Final figure: runtime and accuracy side by side
# ---------------------------------------------------------------

combined_fig = plot_runtime_accuracy_1x2(
  runtime_out = runtime_out,
  accuracy_out = accuracy_out,
  filename_prefix = "ATE_runtime_accuracy_adjusted_1x2",
  runtime_adjustment_power = 1
)