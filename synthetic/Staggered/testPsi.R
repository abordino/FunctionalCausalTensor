setwd("~/Documents/phd/projects/causalMatrix/code")

source("bilinearTensorStaggered.R")
source("bilinearMatrixStaggered.R")
source("bilinearTensorStaggeredPsi.R")
source("bilinearMatrixStaggeredPsi.R")
source("pluginPsi_c1.R")


# ===============================================================
# 0. Global parameters
# ===============================================================

seed = 221198
N_REP = 30

set.seed(seed)

K_values = c(1, 2, 5)
sigma_values = c(0, 0.01)

N = 80
Tt = 100
r = 3
target_k = 1
tau = 1e-3

PLOT_REPS = c(1)

MAKE_DETAILED_PLOTS = TRUE


# ===============================================================
# 1. Direct oracle calculation of Psi from known signal M with staircase miissingness
# ===============================================================

directPsi = function(M, k = 1, N_parts, T_parts,
                     functional = c("ATE", "RowHet", "Local", "Trend"),
                     eta = NULL,
                     row_index = NULL) {
  functional = match.arg(functional)
  
  make_partition_indices = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  row_parts = make_partition_indices(N_parts[[k]])
  col_parts = make_partition_indices(T_parts[[k]])
  
  o_k = length(N_parts[[k]])
  
  if (functional == "RowHet") {
    stopifnot(!is.null(eta))
  }
  
  if (functional == "Local") {
    stopifnot(!is.null(row_index))
    
    a0 = which(vapply(
      row_parts,
      function(idx) row_index %in% idx,
      logical(1)
    ))
    
    stopifnot(length(a0) == 1)
    
    local_pos = match(row_index, row_parts[[a0]])
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      
      if (a + b <= o_k + 1) {
        next
      }
      
      if (functional == "Local" && a != a0) {
        next
      }
      
      rows = row_parts[[a]]
      cols = col_parts[[b]]
      
      Nik = length(rows)
      Tbk = length(cols)
      
      M_ab = M[rows, cols, drop = FALSE]
      
      if (functional == "ATE") {
        
        x = rep(1 / sqrt(Nik), Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "RowHet") {
        
        x = eta[rows] / sqrt(Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "Local") {
        
        x = rep(0, Nik)
        x[local_pos] = 1
        
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Tbk)
        normalizer_increment = Tbk
        
      } else if (functional == "Trend") {
        
        if (Tbk <= 1) {
          next
        }
        
        z = seq_len(Tbk)
        z_centered = z - mean(z)
        
        x = rep(1 / sqrt(Nik), Nik)
        y = z_centered / sqrt(sum(z_centered^2))
        
        weight = 1 / sqrt(Nik * Tbk * (Tbk^2 - 1) / 12)
        normalizer_increment = 1
      }
      
      mu_ab = as.numeric(t(x) %*% M_ab %*% y)
      
      weighted_sum = weighted_sum + weight * mu_ab
      normalizer = normalizer + normalizer_increment
    }
  }
  
  if (normalizer == 0) {
    stop(
      paste0(
        "No target blocks were included in directPsi for functional = ",
        functional)
    )
  }
  
  weighted_sum / normalizer
}


# ===============================================================
# 2. Random adoption-time generator
# ===============================================================

draw_middle_adoption_times = function(
    n,
    Tt,
    center_frac = 0.50,
    sd_frac = 0.12,
    min_frac = 0.20,
    max_frac = 0.80
) {
  center = center_frac * Tt
  sd = sd_frac * Tt
  
  lower = max(2, round(min_frac * Tt))
  upper = min(Tt, round(max_frac * Tt))
  
  A = round(rnorm(n, mean = center, sd = sd))
  A = pmax(lower, pmin(upper, A))
  
  as.integer(A)
}


# ===============================================================
# 3. General staggered Tucker2 simulation
# ===============================================================

simulate_general_staggered_tucker2 = function(
    N = 80,
    Tt = 100,
    K = 5,
    r = 3,
    sigma = 0,
    p_never = 0.20,
    adoption_center_frac = 0.50,
    adoption_sd_frac = 0.12,
    adoption_min_frac = 0.20,
    adoption_max_frac = 0.80,
    seed = 123
) {
  set.seed(seed)
  
  stopifnot(r <= min(N, Tt))
  stopifnot(p_never >= 0, p_never < 1)
  
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
  
  for (j in seq_len(K)) {
    M[, , j] = U %*% C_core[, , j] %*% t(V)
  }
  
  E = array(
    rnorm(N * Tt * K, mean = 0, sd = sigma),
    dim = c(N, Tt, K)
  )
  
  Y_full = M + E
  
  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(FALSE, dim = c(N, Tt, K))
  Y_obs = array(NA_real_, dim = c(N, Tt, K))
  
  for (j in seq_len(K)) {
    
    n_never = max(1, round(p_never * N))
    never_rows = sample(seq_len(N), size = n_never)
    
    adopter_rows = setdiff(seq_len(N), never_rows)
    
    A[adopter_rows, j] = draw_middle_adoption_times(
      n = length(adopter_rows),
      Tt = Tt,
      center_frac = adoption_center_frac,
      sd_frac = adoption_sd_frac,
      min_frac = adoption_min_frac,
      max_frac = adoption_max_frac
    )
    
    A[never_rows, j] = Inf
    
    for (i in seq_len(N)) {
      if (is.infinite(A[i, j])) {
        Omega[i, , j] = TRUE
      } else if (A[i, j] > 1) {
        Omega[i, seq_len(A[i, j] - 1), j] = TRUE
      }
    }
    
    Y_obs[, , j][Omega[, , j]] = Y_full[, , j][Omega[, , j]]
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
    p_never = p_never,
    adoption_center_frac = adoption_center_frac,
    adoption_sd_frac = adoption_sd_frac,
    adoption_min_frac = adoption_min_frac,
    adoption_max_frac = adoption_max_frac,
    seed = seed
  )
  
  class(out) = "general_staggered_tucker2_sim"
  
  out
}


# ===============================================================
# 4. Staircase diagnostics
# ===============================================================

is_staircase_layer = function(Omega_k) {
  N = nrow(Omega_k)
  Tt = ncol(Omega_k)
  
  obs_len = integer(N)
  
  for (i in seq_len(N)) {
    obs = Omega_k[i, ]
    
    if (any(diff(as.integer(obs)) > 0)) {
      return(FALSE)
    }
    
    first_unobs = which(!obs)[1]
    
    if (is.na(first_unobs)) {
      obs_len[i] = Tt
    } else {
      obs_len[i] = first_unobs - 1
    }
  }
  
  all(diff(obs_len) <= 0)
}


get_target_plugin_inputs = function(A, Tt, k) {
  A_k = A[, k]
  
  obs_len = ifelse(
    is.infinite(A_k),
    Tt,
    pmin(Tt, A_k - 1)
  )

  row_perm = order(-obs_len, seq_along(obs_len))
  
  obs_len_sorted = obs_len[row_perm]
  
  if (any(diff(obs_len_sorted) > 0)) {
    stop("Target-layer row sorting failed.")
  }
  
  m_desc = unique(obs_len_sorted)
  o_k = length(m_desc)
  
  if (m_desc[1] != Tt) {
    stop("Target layer must contain at least one fully observed row block.")
  }
  
  if (tail(m_desc, 1) <= 0) {
    stop("Target layer must contain a non-empty initial observed time period.")
  }
  
  row_parts = lapply(m_desc, function(m) which(obs_len_sorted == m))
  
  m_asc = rev(m_desc)
  T_part = diff(c(0, m_asc))
  
  col_ends = cumsum(T_part)
  col_starts = c(1, head(col_ends, -1) + 1)
  col_parts = Map(seq, col_starts, col_ends)
  
  N_part = vapply(row_parts, length, integer(1))
  
  list(
    row_perm = row_perm,
    obs_len_sorted = obs_len_sorted,
    row_parts = row_parts,
    col_parts = col_parts,
    N_part = N_part,
    T_part = T_part,
    N_parts = list(N_part),
    T_parts = list(T_part),
    o_k = o_k
  )
}


choose_valid_local_row = function(sim, k = 1) {
  A_k = sim$A[, k]
  Tt = dim(sim$Y_obs)[2]
  
  obs_len = ifelse(
    is.infinite(A_k),
    Tt,
    pmin(Tt, A_k - 1)
  )
  
  candidate_rows = which(is.finite(A_k) & obs_len < Tt)
  
  if (length(candidate_rows) == 0) {
    stop("No valid adopter rows for Local functional.")
  }
  
  candidate_rows[which.min(obs_len[candidate_rows])]
}


print_staggered_diagnostics = function(sim, k) {
  K = dim(sim$Y_obs)[3]
  Tt = dim(sim$Y_obs)[2]
  
  cat("\nMissingness diagnostics\n")
  cat("--------------------------------------------------\n")
  
  blocks = get_target_plugin_inputs(
    A = sim$A,
    Tt = Tt,
    k = k
  )
  
  before = logical(K)
  after_target_perm = logical(K)
  
  for (j in seq_len(K)) {
    before[j] = is_staircase_layer(sim$Omega[, , j])
    after_target_perm[j] = is_staircase_layer(sim$Omega[blocks$row_perm, , j])
  }
  
  out = data.frame(
    layer = seq_len(K),
    staircase_before_sorting = before,
    staircase_after_target_permutation = after_target_perm
  )
  
  print(out)
  
  cat("\nTarget layer:", k, "\n")
  cat("Target staircase after sorting:", after_target_perm[k], "\n")
  cat("N_parts for plug-in:", paste(blocks$N_part, collapse = " "), "\n")
  cat("T_parts for plug-in:", paste(blocks$T_part, collapse = " "), "\n")
  
  finite_A = sim$A[is.finite(sim$A[, k]), k]
  
  cat("\nTarget-layer adoption-time summary\n")
  cat("--------------------------------------------------\n")
  print(summary(finite_A))
  cat("Mean finite adoption time:", mean(finite_A), "\n")
  cat("Tt / 2:", Tt / 2, "\n")
  
  invisible(out)
}


# ===============================================================
# 5. Evaluate three methods for one simulation
# ===============================================================

evaluate_three_methods = function(
    sim,
    k = 1,
    tau = 1e-3,
    eta = NULL,
    row_index = NULL
) {
  N = dim(sim$Y_obs)[1]
  Tt = dim(sim$Y_obs)[2]
  
  if (is.null(eta)) {
    eta = rep(c(1, -1), length.out = N)
  }
  
  if (is.null(row_index)) {
    row_index = choose_valid_local_row(sim, k = k)
  }
  
  stopifnot(length(eta) == N)
  stopifnot(row_index >= 1, row_index <= N)
  
  blocks = get_target_plugin_inputs(
    A = sim$A,
    Tt = Tt,
    k = k
  )
  
  # -------------------------------------------------------------
  # Reorder target layer for plug-in and oracle calculations.
  # -------------------------------------------------------------
  
  row_perm = blocks$row_perm
  
  M_target_sorted = sim$M[row_perm, , k]
  Y_full_target_sorted = sim$Y_full[row_perm, , k]
  
  eta_sorted = eta[row_perm]
  row_index_sorted = match(row_index, row_perm)

  
  N_parts_plugin = blocks$N_parts
  T_parts_plugin = blocks$T_parts
  k_plugin = 1
  
  functionals = c("ATE", "RowHet", "Local", "Trend")
  
  out = vector("list", length(functionals))
  
  for (h_idx in seq_along(functionals)) {
    h = functionals[h_idx]
    
    truth_signal = directPsi(
      M = M_target_sorted,
      k = k_plugin,
      N_parts = N_parts_plugin,
      T_parts = T_parts_plugin,
      functional = h,
      eta = if (h == "RowHet") eta_sorted else NULL,
      row_index = if (h == "Local") row_index_sorted else NULL
    )
    
    plugin_signal = pluginPsi_c1(
      Y = M_target_sorted,
      k = k_plugin,
      N_parts = N_parts_plugin,
      T_parts = T_parts_plugin,
      functional = h,
      eta = if (h == "RowHet") eta_sorted else NULL,
      row_index = if (h == "Local") row_index_sorted else NULL
    )
    
    plugin_fullY = pluginPsi_c1(
      Y = Y_full_target_sorted,
      k = k_plugin,
      N_parts = N_parts_plugin,
      T_parts = T_parts_plugin,
      functional = h,
      eta = if (h == "RowHet") eta_sorted else NULL,
      row_index = if (h == "Local") row_index_sorted else NULL
    )
    
    tensor_est = switch(
      h,
      ATE = bilinearTensorStaggeredATE(
        Y = sim$Y_obs,
        k = k,
        r = sim$r,
        tau = tau,
        A = sim$A,
        Omega = sim$Omega
      ),
      RowHet = bilinearTensorStaggeredRowHet(
        Y = sim$Y_obs,
        k = k,
        r = sim$r,
        eta = eta,
        tau = tau,
        A = sim$A,
        Omega = sim$Omega
      ),
      Local = bilinearTensorStaggeredLocal(
        Y = sim$Y_obs,
        k = k,
        r = sim$r,
        row_index = row_index,
        tau = tau,
        A = sim$A,
        Omega = sim$Omega
      ),
      Trend = bilinearTensorStaggeredTrend(
        Y = sim$Y_obs,
        k = k,
        r = sim$r,
        tau = tau,
        A = sim$A,
        Omega = sim$Omega
      )
    )
    
    matrix_est = switch(
      h,
      ATE = bilinearMatrixStaggeredATE(
        Y_mat = sim$Y_obs[, , k],
        r = sim$r,
        tau = tau,
        A = sim$A[, k],
        Omega = sim$Omega[, , k]
      ),
      RowHet = bilinearMatrixStaggeredRowHet(
        Y_mat = sim$Y_obs[, , k],
        r = sim$r,
        eta = eta,
        tau = tau,
        A = sim$A[, k],
        Omega = sim$Omega[, , k]
      ),
      Local = bilinearMatrixStaggeredLocal(
        Y_mat = sim$Y_obs[, , k],
        r = sim$r,
        row_index = row_index,
        tau = tau,
        A = sim$A[, k],
        Omega = sim$Omega[, , k]
      ),
      Trend = bilinearMatrixStaggeredTrend(
        Y_mat = sim$Y_obs[, , k],
        r = sim$r,
        tau = tau,
        A = sim$A[, k],
        Omega = sim$Omega[, , k]
      )
    )
    
    out[[h_idx]] = data.frame(
      sigma = sim$sigma,
      K = dim(sim$Y_obs)[3],
      k = k,
      functional = h,
      
      truth_signal = truth_signal,
      
      plugin_signal = plugin_signal,
      plugin_fullY = plugin_fullY,
      tensor_est = tensor_est,
      matrix_est = matrix_est,
      
      error_plugin_signal = plugin_signal - truth_signal,
      error_plugin_fullY = plugin_fullY - truth_signal,
      error_tensor = tensor_est - truth_signal,
      error_matrix = matrix_est - truth_signal,
      
      abs_error_plugin_signal = abs(plugin_signal - truth_signal),
      abs_error_plugin_fullY = abs(plugin_fullY - truth_signal),
      abs_error_tensor = abs(tensor_est - truth_signal),
      abs_error_matrix = abs(matrix_est - truth_signal)
    )
  }
  
  results = do.call(rbind, out)
  
  attr(results, "N_parts") = blocks$N_part
  attr(results, "T_parts") = blocks$T_part
  attr(results, "row_perm") = row_perm
  attr(results, "row_index_original") = row_index
  attr(results, "row_index_sorted") = row_index_sorted
  
  results
}


# ===============================================================
# 6. Visual checks for one simulation
# ===============================================================

plot_missingness_layer = function(
    Omega_k,
    row_perm = NULL,
    main = NULL
) {
  if (!is.null(row_perm)) {
    Omega_k = Omega_k[row_perm, , drop = FALSE]
  }
  
  N = nrow(Omega_k)
  Tt = ncol(Omega_k)
  
  Z = matrix(as.numeric(Omega_k), nrow = N, ncol = Tt)
  
  image(
    x = seq_len(Tt),
    y = seq_len(N),
    z = t(Z[N:1, ]),
    col = c("red", "blue"),
    axes = FALSE,
    xlab = "Time index t",
    ylab = "Unit index i",
    main = main
  )
  
  axis(1)
  axis(2, at = pretty(seq_len(N)), labels = rev(pretty(seq_len(N))))
  box()
  
  invisible(NULL)
}


plot_before_after_target_sorting = function(sim, k = 1) {
  blocks = get_target_plugin_inputs(
    A = sim$A,
    Tt = dim(sim$Y_obs)[2],
    k = k
  )
  
  row_perm = blocks$row_perm
  K = dim(sim$Y_obs)[3]
  
  layers_to_plot = unique(c(k, if (K >= 2) 2 else integer(0), K))
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, length(layers_to_plot)), mar = c(4, 4, 3, 1))
  
  for (j in layers_to_plot) {
    plot_missingness_layer(
      Omega_k = sim$Omega[, , j],
      main = paste0("Before sorting, layer ", j)
    )
  }
  
  for (j in layers_to_plot) {
    plot_missingness_layer(
      Omega_k = sim$Omega[, , j],
      row_perm = row_perm,
      main = paste0("After target k=", k, " sorting, layer ", j)
    )
  }
  
  invisible(NULL)
}


plot_target_blocks = function(sim, k = 1) {
  blocks = get_target_plugin_inputs(
    A = sim$A,
    Tt = dim(sim$Y_obs)[2],
    k = k
  )
  
  row_perm = blocks$row_perm
  Omega_k = sim$Omega[row_perm, , k]
  
  N = nrow(Omega_k)
  Tt = ncol(Omega_k)
  
  Z = matrix(as.numeric(Omega_k), nrow = N, ncol = Tt)
  
  image(
    x = seq_len(Tt),
    y = seq_len(N),
    z = t(Z[N:1, ]),
    col = c("red", "blue"),
    axes = FALSE,
    xlab = "Time index t",
    ylab = "Sorted unit index",
    main = paste0("Target layer ", k, " staircase blocks")
  )
  
  axis(1)
  axis(2, at = pretty(seq_len(N)), labels = rev(pretty(seq_len(N))))
  box()
  
  row_cuts = cumsum(blocks$N_part)
  col_cuts = cumsum(blocks$T_part)
  
  abline(v = col_cuts + 0.5, col = "black", lwd = 1)
  abline(h = N - row_cuts + 0.5, col = "black", lwd = 1)
  
  o_k = blocks$o_k
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      if (a + b <= o_k + 1) {
        next
      }
      
      rows = blocks$row_parts[[a]]
      cols = blocks$col_parts[[b]]
      
      xleft = min(cols) - 0.5
      xright = max(cols) + 0.5
      ybottom = N - max(rows) + 0.5
      ytop = N - min(rows) + 0.5
      
      rect(
        xleft = xleft,
        ybottom = ybottom,
        xright = xright,
        ytop = ytop,
        border = "yellow",
        lwd = 2
      )
      
      text(
        x = mean(c(xleft, xright)),
        y = mean(c(ybottom, ytop)),
        labels = paste0("(", a, ",", b, ")"),
        col = "yellow",
        font = 2,
        cex = 0.7
      )
    }
  }
  
  invisible(NULL)
}


plot_estimates_vs_truth = function(results, main = NULL) {
  if (is.null(main)) {
    main = paste0(
      "Estimates vs truth, K = ",
      unique(results$K),
      ", sigma = ",
      unique(results$sigma),
      ", rep = ",
      unique(results$rep)
    )
  }
  
  truth = results$truth_signal
  
  vals = c(
    truth,
    results$plugin_signal,
    results$plugin_fullY,
    results$tensor_est,
    results$matrix_est
  )
  
  lims = range(vals, finite = TRUE)
  
  plot(
    truth,
    results$tensor_est,
    xlim = lims,
    ylim = lims,
    pch = 19,
    col = "steelblue",
    xlab = "Truth from signal M",
    ylab = "Estimate",
    main = main
  )
  
  points(truth, results$matrix_est, pch = 17, col = "orange")
  points(truth, results$plugin_fullY, pch = 15, col = "darkgreen")
  points(truth, results$plugin_signal, pch = 1, col = "black")
  
  abline(0, 1, lty = 2, lwd = 2, col = "red")
  
  text(
    truth,
    results$tensor_est,
    labels = results$functional,
    pos = 3,
    cex = 0.75,
    col = "steelblue"
  )
  
  legend(
    "topleft",
    legend = c(
      "tensor-pooled",
      "matrix-only",
      "plugin full Y",
      "plugin signal",
      "45-degree line"
    ),
    col = c("steelblue", "orange", "darkgreen", "black", "red"),
    pch = c(19, 17, 15, 1, NA),
    lty = c(NA, NA, NA, NA, 2),
    lwd = c(NA, NA, NA, NA, 2),
    bty = "n"
  )
  
  invisible(NULL)
}


plot_abs_error_barplot = function(results, main = NULL) {
  if (is.null(main)) {
    main = paste0(
      "Absolute errors, K = ",
      unique(results$K),
      ", sigma = ",
      unique(results$sigma),
      ", rep = ",
      unique(results$rep)
    )
  }
  
  err_mat = rbind(
    plugin_signal = results$abs_error_plugin_signal,
    plugin_fullY = results$abs_error_plugin_fullY,
    tensor = results$abs_error_tensor,
    matrix = results$abs_error_matrix
  )
  
  colnames(err_mat) = results$functional
  
  barplot(
    err_mat,
    beside = TRUE,
    names.arg = results$functional,
    col = c("gray80", "darkgreen", "steelblue", "orange"),
    ylab = "Absolute error",
    main = main,
    legend.text = rownames(err_mat),
    args.legend = list(x = "topright", bty = "n")
  )
  
  invisible(NULL)
}


plot_adoption_times = function(sim, k = 1) {
  A_k = sim$A[, k]
  finite_A = A_k[is.finite(A_k)]
  Tt = dim(sim$Y_obs)[2]
  
  hist(
    finite_A,
    breaks = 20,
    col = "gray80",
    border = "white",
    main = paste0("Finite adoption times, layer ", k),
    xlab = "Adoption time A[i,k]",
    xlim = c(1, Tt)
  )
  
  abline(v = Tt / 2, col = "red", lwd = 2, lty = 2)
  abline(v = mean(finite_A), col = "blue", lwd = 2)
  
  legend(
    "topright",
    legend = c("Tt / 2", "mean finite adoption time"),
    col = c("red", "blue"),
    lty = c(2, 1),
    lwd = 2,
    bty = "n"
  )
  
  invisible(NULL)
}


plot_single_rep_diagnostics = function(sim, results, k = 1) {
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  plot_adoption_times(sim, k = k)
  plot_before_after_target_sorting(sim, k = k)
  plot_target_blocks(sim, k = k)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot_estimates_vs_truth(results)
  plot_abs_error_barplot(results)
  
  invisible(NULL)
}


# ===============================================================
# 7. One repetition for one setting
# ===============================================================

run_single_rep = function(
    K,
    sigma,
    rep_id,
    seed,
    N = 80,
    Tt = 100,
    r = 3,
    k = 1,
    tau = 1e-3,
    make_plots = FALSE
) {
  cat("\n--------------------------------------------------\n")
  cat("K =", K, "| sigma =", sigma, "| rep =", rep_id, "| seed =", seed, "\n")
  cat("--------------------------------------------------\n")
  
  sim = simulate_general_staggered_tucker2(
    N = N,
    Tt = Tt,
    K = K,
    r = r,
    sigma = sigma,
    p_never = 0.20,
    
    adoption_center_frac = 0.50,
    adoption_sd_frac = 0.12,
    adoption_min_frac = 0.20,
    adoption_max_frac = 0.80,
    
    seed = seed
  )
  
  if (make_plots) {
    print_staggered_diagnostics(sim, k = k)
  }
  
  eta = rep(c(1, -1), length.out = N)
  row_index = choose_valid_local_row(sim, k = k)
  
  results = evaluate_three_methods(
    sim = sim,
    k = k,
    tau = tau,
    eta = eta,
    row_index = row_index
  )
  
  results$rep = rep_id
  results$seed = seed
  
  if (make_plots) {
    cat("\nLocal row choice\n")
    cat("--------------------------------------------------\n")
    cat("Original row_index:", row_index, "\n")
    cat("A[row_index, k]:", sim$A[row_index, k], "\n")
    cat(
      "Observed prefix length:",
      ifelse(
        is.infinite(sim$A[row_index, k]),
        Tt,
        sim$A[row_index, k] - 1
      ),
      "\n"
    )
    
    cat("\nResults for plotted repetition\n")
    cat("--------------------------------------------------\n")
    print(results)
    
    plot_single_rep_diagnostics(sim, results, k = k)
  }
  
  invisible(list(
    sim = sim,
    results = results
  ))
}


# ===============================================================
# 8. Run n_rep repetitions for one setting
# ===============================================================

run_setting_reps = function(
    K,
    sigma,
    sigma_id,
    n_rep = 30,
    master_seed = 221198,
    N = 80,
    Tt = 100,
    r = 3,
    k = 1,
    tau = 1e-3,
    plot_reps = c(1),
    make_detailed_plots = TRUE
) {
  cat("\n==================================================\n")
  cat("Running setting: K =", K, ", sigma =", sigma, "\n")
  cat("Number of repetitions:", n_rep, "\n")
  cat("==================================================\n")
  
  rep_results = vector("list", n_rep)
  
  for (rep_id in seq_len(n_rep)) {
    
    seed = master_seed +
      1000000 * K +
      100000 * sigma_id +
      rep_id
    
    make_plots = make_detailed_plots && rep_id %in% plot_reps
    
    run = run_single_rep(
      K = K,
      sigma = sigma,
      rep_id = rep_id,
      seed = seed,
      N = N,
      Tt = Tt,
      r = r,
      k = k,
      tau = tau,
      make_plots = make_plots
    )
    
    rep_results[[rep_id]] = run$results
  }
  
  setting_results = do.call(rbind, rep_results)
  
  cat("\nSetting summary: K =", K, ", sigma =", sigma, "\n")
  cat("--------------------------------------------------\n")
  
  setting_summary = aggregate(
    cbind(
      abs_error_plugin_signal,
      abs_error_plugin_fullY,
      abs_error_tensor,
      abs_error_matrix
    ) ~ functional,
    data = setting_results,
    FUN = mean
  )
  
  print(setting_summary)
  
  invisible(setting_results)
}


# ===============================================================
# 9. Aggregate Monte Carlo plots
# ===============================================================

plot_method_summary_by_setting = function(all_results) {
  summary_df = aggregate(
    cbind(
      abs_error_plugin_signal,
      abs_error_plugin_fullY,
      abs_error_tensor,
      abs_error_matrix
    ) ~ K + sigma,
    data = all_results,
    FUN = mean
  )
  
  labels = paste0(
    "K=", summary_df$K,
    "\n",
    "sigma=", summary_df$sigma
  )
  
  err_mat = t(as.matrix(summary_df[, c(
    "abs_error_plugin_signal",
    "abs_error_plugin_fullY",
    "abs_error_tensor",
    "abs_error_matrix"
  )]))
  
  rownames(err_mat) = c(
    "plugin signal",
    "plugin full Y",
    "tensor",
    "matrix"
  )
  
  barplot(
    err_mat,
    beside = TRUE,
    names.arg = labels,
    col = c("gray80", "darkgreen", "steelblue", "orange"),
    ylab = "Mean absolute error",
    main = paste0("Mean absolute error by setting, n_rep = ", N_REP),
    legend.text = rownames(err_mat),
    args.legend = list(x = "topright", bty = "n")
  )
  
  invisible(summary_df)
}


plot_method_summary_by_functional = function(all_results) {
  summary_df = aggregate(
    cbind(
      abs_error_plugin_signal,
      abs_error_plugin_fullY,
      abs_error_tensor,
      abs_error_matrix
    ) ~ functional,
    data = all_results,
    FUN = mean
  )
  
  err_mat = t(as.matrix(summary_df[, c(
    "abs_error_plugin_signal",
    "abs_error_plugin_fullY",
    "abs_error_tensor",
    "abs_error_matrix"
  )]))
  
  rownames(err_mat) = c(
    "plugin signal",
    "plugin full Y",
    "tensor",
    "matrix"
  )
  
  barplot(
    err_mat,
    beside = TRUE,
    names.arg = summary_df$functional,
    col = c("gray80", "darkgreen", "steelblue", "orange"),
    ylab = "Mean absolute error",
    main = paste0("Mean absolute error by functional, n_rep = ", N_REP),
    legend.text = rownames(err_mat),
    args.legend = list(x = "topright", bty = "n")
  )
  
  invisible(summary_df)
}


plot_error_boxplots = function(all_results, method = c("tensor", "matrix", "plugin_fullY")) {
  method = match.arg(method)
  
  err_col = switch(
    method,
    tensor = "abs_error_tensor",
    matrix = "abs_error_matrix",
    plugin_fullY = "abs_error_plugin_fullY"
  )
  
  boxplot(
    all_results[[err_col]] ~ all_results$K + all_results$sigma,
    xlab = "K.sigma setting",
    ylab = "Absolute error",
    main = paste0("Distribution of absolute errors: ", method),
    col = "gray80"
  )
  
  invisible(NULL)
}


plot_tensor_vs_matrix_errors = function(all_results) {
  lims = range(
    c(all_results$abs_error_tensor, all_results$abs_error_matrix),
    finite = TRUE
  )
  
  plot(
    all_results$abs_error_tensor,
    all_results$abs_error_matrix,
    xlim = lims,
    ylim = lims,
    pch = 19,
    col = "gray40",
    xlab = "Tensor-pooled absolute error",
    ylab = "Matrix-only absolute error",
    main = "Matrix-only error vs tensor-pooled error"
  )
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  invisible(NULL)
}


# ===============================================================
# 10. Main Monte Carlo experiment
# ===============================================================

all_setting_results = list()
counter = 1

sigma_ids = seq_along(sigma_values)

for (K in K_values) {
  for (s_idx in sigma_ids) {
    sigma = sigma_values[s_idx]
    
    setting_results = run_setting_reps(
      K = K,
      sigma = sigma,
      sigma_id = s_idx,
      n_rep = N_REP,
      master_seed = MASTER_SEED,
      N = N,
      Tt = Tt,
      r = r,
      k = target_k,
      tau = tau,
      plot_reps = PLOT_REPS,
      make_detailed_plots = MAKE_DETAILED_PLOTS
    )
    
    all_setting_results[[counter]] = setting_results
    counter = counter + 1
  }
}

all_results = do.call(rbind, all_setting_results)


cat("\n==================================================\n")
cat("All Monte Carlo results\n")
cat("==================================================\n")

print(all_results)


cat("\n==================================================\n")
cat("Monte Carlo summary by K, sigma, functional\n")
cat("==================================================\n")

summary_by_functional = aggregate(
  cbind(
    abs_error_plugin_signal,
    abs_error_plugin_fullY,
    abs_error_tensor,
    abs_error_matrix
  ) ~ K + sigma + functional,
  data = all_results,
  FUN = mean
)

print(summary_by_functional)


cat("\n==================================================\n")
cat("Monte Carlo summary by K and sigma\n")
cat("==================================================\n")

summary_by_setting = aggregate(
  cbind(
    abs_error_plugin_signal,
    abs_error_plugin_fullY,
    abs_error_tensor,
    abs_error_matrix
  ) ~ K + sigma,
  data = all_results,
  FUN = mean
)

print(summary_by_setting)


cat("\n==================================================\n")
cat("Monte Carlo standard deviations by K and sigma\n")
cat("==================================================\n")

sd_by_setting = aggregate(
  cbind(
    abs_error_plugin_signal,
    abs_error_plugin_fullY,
    abs_error_tensor,
    abs_error_matrix
  ) ~ K + sigma,
  data = all_results,
  FUN = sd
)

print(sd_by_setting)


cat("\n==================================================\n")
cat("Overall method comparison\n")
cat("==================================================\n")

overall_summary = data.frame(
  method = c("plugin_signal", "plugin_fullY", "tensor", "matrix"),
  mean_abs_error = c(
    mean(all_results$abs_error_plugin_signal),
    mean(all_results$abs_error_plugin_fullY),
    mean(all_results$abs_error_tensor),
    mean(all_results$abs_error_matrix)
  ),
  median_abs_error = c(
    median(all_results$abs_error_plugin_signal),
    median(all_results$abs_error_plugin_fullY),
    median(all_results$abs_error_tensor),
    median(all_results$abs_error_matrix)
  ),
  sd_abs_error = c(
    sd(all_results$abs_error_plugin_signal),
    sd(all_results$abs_error_plugin_fullY),
    sd(all_results$abs_error_tensor),
    sd(all_results$abs_error_matrix)
  )
)

print(overall_summary)


# ===============================================================
# 11. Visual summaries
# ===============================================================

old_par = par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

plot_method_summary_by_setting(all_results)

plot_method_summary_by_functional(all_results)

par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))
plot_error_boxplots(all_results, method = "tensor")
plot_error_boxplots(all_results, method = "matrix")
plot_error_boxplots(all_results, method = "plugin_fullY")

par(mfrow = c(1, 1), mar = c(5, 4, 3, 1))
plot_tensor_vs_matrix_errors(all_results)
