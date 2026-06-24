setwd("~/Documents/phd/projects/causalMatrix/code")

source("bilinearTensorStaggered.R")
source("bilinearMatrixStaggered.R")

library(ggplot2)
library(dplyr)
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
# 1. Simulate Tucker2 tensor with staggered missingness
# ===============================================================

simulate_general_staggered_tucker2 = function(
    N = 150,
    Tt = 200,
    K = 10,
    r = 5,
    sigma = 1,
    n_adopt_times = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5),
    p_never = 0.20,
    p_initial = 0.10,
    seed = 123
) {
  set.seed(seed)
  
  stopifnot(length(n_adopt_times) == K)
  stopifnot(r <= min(N, Tt))
  
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
  
  for (k in seq_len(K)) {
    M[, , k] = U %*% C_core[, , k] %*% t(V)
  }
  
  E = array(
    rnorm(N * Tt * K, mean = 0, sd = sigma),
    dim = c(N, Tt, K)
  )
  
  Y_full = M + E
  
  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(FALSE, dim = c(N, Tt, K))
  Y_obs = array(NA_real_, dim = c(N, Tt, K))
  
  first_adopt_time = max(2, floor(p_initial * Tt) + 1)
  
  for (k in seq_len(K)) {
    adopt_grid = unique(round(seq(
      first_adopt_time,
      Tt,
      length.out = n_adopt_times[k]
    )))
    
    for (i in seq_len(N)) {
      if (runif(1) < p_never) {
        A[i, k] = Inf
      } else {
        A[i, k] = sample(adopt_grid, size = 1)
      }
      
      if (is.infinite(A[i, k])) {
        Omega[i, , k] = TRUE
      } else if (A[i, k] > 1) {
        Omega[i, seq_len(A[i, k] - 1), k] = TRUE
      }
    }
    
    Y_obs[, , k][Omega[, , k]] = Y_full[, , k][Omega[, , k]]
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
    n_adopt_times = n_adopt_times,
    p_never = p_never,
    p_initial = p_initial
  )
  
  class(out) = "general_staggered_tucker2_sim"
  
  out
}


#------------------------------------------------------------------------------
# 2. Target-layer staircase blocks after row rearrangement
# ===============================================================

get_target_blocks = function(sim, k) {
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


# ===============================================================
# 3. Choose valid missing target blocks
# ===============================================================

choose_missing_target = function(sim, k = 1, prefer = "middle") {
  blocks = get_target_blocks(sim, k)
  o_k = blocks$o_k
  
  targets = expand.grid(
    a = seq_len(o_k),
    b = seq_len(o_k)
  )
  
  targets = targets[targets$a + targets$b > o_k + 1, , drop = FALSE]
  
  if (nrow(targets) == 0) {
    stop("No missing target blocks found.")
  }
  
  if (prefer == "deep") {
    idx = which.max(targets$a + targets$b)
  } else {
    center = (o_k + 1) / 2
    idx = which.min((targets$a - center)^2 + (targets$b - center)^2)
  }
  
  list(
    k = k,
    a = targets$a[idx],
    b = targets$b[idx],
    o_k = o_k,
    blocks = blocks
  )
}


make_missing_targets = function(sim, k = 1, n_targets = 3) {
  blocks = get_target_blocks(sim, k)
  o_k = blocks$o_k
  
  targets = expand.grid(
    k = k,
    a = seq_len(o_k),
    b = seq_len(o_k)
  )
  
  targets = targets[targets$a + targets$b > o_k + 1, , drop = FALSE]
  
  if (nrow(targets) < n_targets) {
    stop("Not enough missing target blocks.")
  }
  
  targets[seq_len(n_targets), , drop = FALSE]
}


# ===============================================================
# 4. Check staircase property
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


print_general_staggered_missingness = function(sim) {
  K = dim(sim$Y_obs)[3]
  
  for (k in seq_len(K)) {
    A_k = sim$A[, k]
    
    cat("\n")
    cat("============================================\n")
    cat("Layer k =", k, "\n")
    cat("Distinct finite adoption times:",
        paste(sort(unique(A_k[is.finite(A_k)])), collapse = " "),
        "\n")
    cat("Never adopters:", sum(is.infinite(A_k)), "\n")
    cat("Observed entries:",
        sum(sim$Omega[, , k]),
        "out of",
        length(sim$Omega[, , k]),
        "\n")
    cat("Staircase before target sorting:",
        is_staircase_layer(sim$Omega[, , k]),
        "\n")
  }
  
  invisible(NULL)
}


test_target_rearrangement = function(sim, k) {
  blocks = get_target_blocks(sim, k)
  row_perm = blocks$row_perm
  
  K = dim(sim$Y_obs)[3]
  
  before = logical(K)
  after = logical(K)
  
  for (j in seq_len(K)) {
    before[j] = is_staircase_layer(sim$Omega[, , j])
    after[j] = is_staircase_layer(sim$Omega[row_perm, , j])
  }
  
  out = data.frame(
    layer = seq_len(K),
    staircase_before = before,
    staircase_after_common_target_permutation = after
  )
  
  cat("\nTarget layer:", k, "\n")
  cat("Target layer staircase after permutation:", after[k], "\n\n")
  
  print(out)
  
  invisible(out)
}


# ===============================================================
# 5. Plot missingness
# ===============================================================

plot_missingness_layer = function(
    sim,
    layer = 1,
    row_perm = NULL,
    main = NULL
) {
  Omega_k = sim$Omega[, , layer]
  
  if (!is.null(row_perm)) {
    Omega_k = Omega_k[row_perm, , drop = FALSE]
  }
  
  N = nrow(Omega_k)
  Tt = ncol(Omega_k)
  
  Z = matrix(as.numeric(Omega_k), nrow = N, ncol = Tt)
  
  if (is.null(main)) {
    main = paste0(
      "Missingness, layer ",
      layer,
      " | blue = observed, red = missing"
    )
  }
  
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


plot_before_after_rearrangement = function(sim, k) {
  K = dim(sim$Y_obs)[3]
  blocks = get_target_blocks(sim, k)
  row_perm = blocks$row_perm
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, K), mar = c(4, 4, 3, 1))
  
  for (j in seq_len(K)) {
    plot_missingness_layer(
      sim = sim,
      layer = j,
      main = paste0("Before, layer ", j)
    )
  }
  
  for (j in seq_len(K)) {
    plot_missingness_layer(
      sim = sim,
      layer = j,
      row_perm = row_perm,
      main = paste0("After target k=", k, ", layer ", j)
    )
  }
  
  invisible(NULL)
}


# ===============================================================
# 6. Single comparison
# ===============================================================

test_bilinearStaggered_compare = function(
    N = 150,
    Tt = 200,
    K = 10,
    r = 5,
    sigma = 0,
    n_adopt_times = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5),
    p_never = 0.20,
    p_initial = 0.10,
    k = 1,
    a = NULL,
    b = NULL,
    tau = 1e-2,
    seed = 123
) {
  sim = simulate_general_staggered_tucker2(
    N = N,
    Tt = Tt,
    K = K,
    r = r,
    sigma = sigma,
    n_adopt_times = n_adopt_times,
    p_never = p_never,
    p_initial = p_initial,
    seed = seed
  )
  
  if (is.null(a) || is.null(b)) {
    target = choose_missing_target(sim, k = k)
    a = target$a
    b = target$b
  }
  
  blocks = get_target_blocks(sim, k)
  
  if (a + b <= blocks$o_k + 1) {
    stop("Chosen target block is observed. Pick a + b > o_k + 1.")
  }
  
  set.seed(seed + 1)
  
  x = rnorm(blocks$N_part[a])
  x = x / sqrt(sum(x^2))
  
  y = rnorm(blocks$T_part[b])
  y = y / sqrt(sum(y^2))
  
  mu_hat_tensor = bilinearTensorStaggered(
    Y = sim$Y_obs,
    k = k,
    a = a,
    b = b,
    r = r,
    x = x,
    y = y,
    A = sim$A,
    Omega = sim$Omega,
    tau = tau
  )
  
  mu_hat_matrix = bilinearMatrixStaggered(
    Y_mat = sim$Y_obs[, , k],
    a = a,
    b = b,
    r = r,
    x = x,
    y = y,
    A = sim$A[, k],
    Omega = sim$Omega[, , k],
    tau = tau
  )
  
  row_idx_perm = blocks$row_parts[[a]]
  col_idx = blocks$col_parts[[b]]
  
  row_idx_original = blocks$row_perm[row_idx_perm]
  
  mu_true = as.numeric(
    t(x) %*% sim$M[row_idx_original, col_idx, k] %*% y
  )
  
  out = list(
    mu_true = mu_true,
    mu_hat_tensor = mu_hat_tensor,
    mu_hat_matrix = mu_hat_matrix,
    abs_error_tensor = abs(mu_hat_tensor - mu_true),
    abs_error_matrix = abs(mu_hat_matrix - mu_true),
    rel_error_tensor = abs(mu_hat_tensor - mu_true) / max(abs(mu_true), 1e-12),
    rel_error_matrix = abs(mu_hat_matrix - mu_true) / max(abs(mu_true), 1e-12),
    target = list(k = k, a = a, b = b, o_k = blocks$o_k),
    blocks = blocks,
    sim = sim
  )
  
  print(out[c(
    "mu_true",
    "mu_hat_tensor",
    "mu_hat_matrix",
    "abs_error_tensor",
    "abs_error_matrix",
    "rel_error_tensor",
    "rel_error_matrix",
    "target"
  )])
  
  invisible(out)
}


# ===============================================================
# 7. Multiple-query accuracy check
# ===============================================================

run_accuracy_check_compare = function(
    sim,
    k = 1,
    a = NULL,
    b = NULL,
    n_queries = 100,
    tau = 1e-2,
    seed = 999
) {
  blocks = get_target_blocks(sim, k)
  
  if (is.null(a) || is.null(b)) {
    target = choose_missing_target(sim, k = k)
    a = target$a
    b = target$b
  }
  
  if (a + b <= blocks$o_k + 1) {
    stop("Chosen target block is observed. Need a + b > o_k + 1.")
  }
  
  set.seed(seed)
  
  Nx = blocks$N_part[a]
  Ty = blocks$T_part[b]
  
  row_idx_perm = blocks$row_parts[[a]]
  col_idx = blocks$col_parts[[b]]
  
  row_idx_original = blocks$row_perm[row_idx_perm]
  
  M_target = sim$M[row_idx_original, col_idx, k, drop = FALSE][, , 1]
  
  mu_true_vec = numeric(n_queries)
  
  mu_hat_tensor_vec = numeric(n_queries)
  mu_hat_matrix_vec = numeric(n_queries)
  
  abs_err_tensor_vec = numeric(n_queries)
  abs_err_matrix_vec = numeric(n_queries)
  
  rel_err_tensor_vec = numeric(n_queries)
  rel_err_matrix_vec = numeric(n_queries)
  
  for (q in seq_len(n_queries)) {
    x = rnorm(Nx)
    x = x / sqrt(sum(x^2))
    
    y = rnorm(Ty)
    y = y / sqrt(sum(y^2))
    
    mu_hat_tensor = bilinearTensorStaggered(
      Y = sim$Y_obs,
      k = k,
      a = a,
      b = b,
      r = sim$r,
      x = x,
      y = y,
      A = sim$A,
      Omega = sim$Omega,
      tau = tau
    )
    
    mu_hat_matrix = bilinearMatrixStaggered(
      Y_mat = sim$Y_obs[, , k],
      a = a,
      b = b,
      r = sim$r,
      x = x,
      y = y,
      A = sim$A[, k],
      Omega = sim$Omega[, , k],
      tau = tau
    )
    
    mu_true = as.numeric(t(x) %*% M_target %*% y)
    
    mu_true_vec[q] = mu_true
    
    mu_hat_tensor_vec[q] = mu_hat_tensor
    mu_hat_matrix_vec[q] = mu_hat_matrix
    
    abs_err_tensor_vec[q] = abs(mu_hat_tensor - mu_true)
    abs_err_matrix_vec[q] = abs(mu_hat_matrix - mu_true)
    
    rel_err_tensor_vec[q] = abs(mu_hat_tensor - mu_true) / max(abs(mu_true), 1e-12)
    rel_err_matrix_vec[q] = abs(mu_hat_matrix - mu_true) / max(abs(mu_true), 1e-12)
  }
  
  diagnostics = data.frame(
    query = seq_len(n_queries),
    mu_true = mu_true_vec,
    
    mu_hat_tensor = mu_hat_tensor_vec,
    mu_hat_matrix = mu_hat_matrix_vec,
    
    abs_error_tensor = abs_err_tensor_vec,
    abs_error_matrix = abs_err_matrix_vec,
    
    rel_error_tensor = rel_err_tensor_vec,
    rel_error_matrix = rel_err_matrix_vec,
    
    target_k = k,
    target_a = a,
    target_b = b,
    target_o_k = blocks$o_k
  )
  
  cat("\nAccuracy summary: tensor-pooled vs matrix-only\n")
  cat("--------------------------------------------------\n")
  cat("Target: k =", k, ", a =", a, ", b =", b, ", o_k =", blocks$o_k, "\n")
  cat("Noise sigma:", sim$sigma, "\n\n")
  
  cat("Tensor-pooled method\n")
  cat("Mean absolute error  :", mean(abs_err_tensor_vec), "\n")
  
  cat("Matrix-only method\n")
  cat("Mean absolute error  :", mean(abs_err_matrix_vec), "\n")
  
  cat("Error ratio: matrix / tensor\n")
  cat("Mean absolute error ratio  :",
      mean(abs_err_matrix_vec) / max(mean(abs_err_tensor_vec), 1e-12),
      "\n")
  
  invisible(diagnostics)
}


# ===============================================================
# 8. Plot comparison diagnostics
# ===============================================================

plot_accuracy_diagnostics_compare = function(diagnostics) {
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  lims = range(
    c(
      diagnostics$mu_true,
      diagnostics$mu_hat_tensor,
      diagnostics$mu_hat_matrix
    )
  )
  
  plot(
    diagnostics$mu_true,
    diagnostics$mu_hat_tensor,
    xlab = expression(mu[true]),
    ylab = expression(mu[hat]),
    main = expression(paste(mu[hat], " vs ", mu[true])),
    pch = 19,
    col = "steelblue",
    xlim = lims,
    ylim = lims
  )
  
  points(
    diagnostics$mu_true,
    diagnostics$mu_hat_matrix,
    pch = 17,
    col = "orange"
  )
  
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
  
  ymax = max(
    diagnostics$abs_error_tensor,
    diagnostics$abs_error_matrix
  )
  
  plot(
    diagnostics$query,
    diagnostics$abs_error_tensor,
    type = "b",
    pch = 19,
    col = "steelblue",
    ylim = c(0, ymax),
    xlab = "Query index",
    ylab = "Absolute error",
    main = "Absolute error by query"
  )
  
  lines(
    diagnostics$query,
    diagnostics$abs_error_matrix,
    type = "b",
    pch = 17,
    col = "orange"
  )
  
  abline(
    h = mean(diagnostics$abs_error_tensor),
    col = "steelblue",
    lwd = 2,
    lty = 2
  )
  
  abline(
    h = mean(diagnostics$abs_error_matrix),
    col = "orange",
    lwd = 2,
    lty = 2
  )
  
  legend(
    "topright",
    legend = c("tensor-pooled", "matrix-only"),
    col = c("steelblue", "orange"),
    pch = c(19, 17),
    lty = 1,
    bty = "n"
  )
  
  boxplot(
    diagnostics$abs_error_tensor,
    diagnostics$abs_error_matrix,
    names = c("tensor", "matrix"),
    col = c("steelblue", "orange"),
    ylab = "Absolute error",
    main = "Absolute error comparison"
  )
  
  ymax_rel = max(
    diagnostics$rel_error_tensor,
    diagnostics$rel_error_matrix
  )
  
  plot(
    abs(diagnostics$mu_true),
    diagnostics$rel_error_tensor,
    pch = 19,
    col = "steelblue",
    ylim = c(0, ymax_rel),
    xlab = expression(abs(mu[true])),
    ylab = "Relative error",
    main = "Relative error vs signal size"
  )
  
  points(
    abs(diagnostics$mu_true),
    diagnostics$rel_error_matrix,
    pch = 17,
    col = "orange"
  )
  
  legend(
    "topright",
    legend = c("tensor-pooled", "matrix-only"),
    col = c("steelblue", "orange"),
    pch = c(19, 17),
    bty = "n"
  )
  
  invisible(NULL)
}



# ===============================================================
# 9. Tau-sensitivity study for six missing target blocks
# ===============================================================

prepare_target_queries = function(
    sim,
    k,
    a,
    b,
    n_queries = 100,
    seed = 999
) {
  blocks = get_target_blocks(sim, k)
  
  if (a + b <= blocks$o_k + 1) {
    stop("Chosen target block is observed. Need a + b > o_k + 1.")
  }
  
  Nx = blocks$N_part[a]
  Ty = blocks$T_part[b]
  
  row_idx_perm = blocks$row_parts[[a]]
  row_idx_original = blocks$row_perm[row_idx_perm]
  col_idx = blocks$col_parts[[b]]
  
  M_target = sim$M[row_idx_original, col_idx, k, drop = FALSE][, , 1]
  
  set.seed(seed)
  
  x_list = vector("list", n_queries)
  y_list = vector("list", n_queries)
  mu_true = numeric(n_queries)
  
  for (q in seq_len(n_queries)) {
    x = rnorm(Nx)
    x = x / sqrt(sum(x^2))
    
    y = rnorm(Ty)
    y = y / sqrt(sum(y^2))
    
    x_list[[q]] = x
    y_list[[q]] = y
    mu_true[q] = as.numeric(t(x) %*% M_target %*% y)
  }
  
  list(
    x = x_list,
    y = y_list,
    mu_true = mu_true,
    blocks = blocks,
    row_idx_original = row_idx_original,
    col_idx = col_idx
  )
}


run_tau_sensitivity_compare = function(
    sim,
    targets = NULL,
    k = 1,
    n_targets = 6,
    tau_seq = 10^seq(-6, 0, by = 0.5),
    n_queries = 100,
    query_seed = 999,
    verbose = TRUE
) {
  tau_seq = sort(unique(as.numeric(tau_seq)))
  
  if (length(tau_seq) < 2) {
    stop("tau_seq must contain at least two distinct values.")
  }
  
  if (any(!is.finite(tau_seq)) || any(tau_seq <= 0)) {
    stop("All values in tau_seq must be finite and strictly positive.")
  }
  
  if (is.null(targets)) {
    targets = make_missing_targets(
      sim = sim,
      k = k,
      n_targets = n_targets
    )
  }
  
  if (nrow(targets) != n_targets) {
    stop("targets must contain exactly n_targets rows.")
  }
  
  diagnostics = vector(
    "list",
    length = nrow(targets) * length(tau_seq)
  )
  
  counter = 1
  
  for (m in seq_len(nrow(targets))) {
    k_m = targets$k[m]
    a_m = targets$a[m]
    b_m = targets$b[m]
    
    target_label = paste0(
      "B", m, " (", a_m - 1, ",", b_m, ")"
    )
    
    query_data = prepare_target_queries(
      sim = sim,
      k = k_m,
      a = a_m,
      b = b_m,
      n_queries = n_queries,
      seed = query_seed + m
    )
    
    if (verbose) {
      cat("\nMissing block", m, "of", nrow(targets), ":", target_label, "\n")
    }
    
    for (tt in seq_along(tau_seq)) {
      tau = tau_seq[tt]
      
      if (verbose) {
        cat("  tau =", format(tau, scientific = TRUE), "\n")
      }
      
      mu_hat_tensor = numeric(n_queries)
      mu_hat_matrix = numeric(n_queries)
      
      for (q in seq_len(n_queries)) {
        x = query_data$x[[q]]
        y = query_data$y[[q]]
        
        mu_hat_tensor[q] = tryCatch(
          bilinearTensorStaggered(
            Y = sim$Y_obs,
            k = k_m,
            a = a_m,
            b = b_m,
            r = sim$r,
            x = x,
            y = y,
            A = sim$A,
            Omega = sim$Omega,
            tau = tau
          ),
          error = function(e) {
            stop(
              "Tensor estimator failed for block ", m,
              ", tau = ", tau,
              ", query = ", q,
              ": ", conditionMessage(e),
              call. = FALSE
            )
          }
        )
        
        mu_hat_matrix[q] = tryCatch(
          bilinearMatrixStaggered(
            Y_mat = sim$Y_obs[, , k_m],
            a = a_m,
            b = b_m,
            r = sim$r,
            x = x,
            y = y,
            A = sim$A[, k_m],
            Omega = sim$Omega[, , k_m],
            tau = tau
          ),
          error = function(e) {
            stop(
              "Matrix estimator failed for block ", m,
              ", tau = ", tau,
              ", query = ", q,
              ": ", conditionMessage(e),
              call. = FALSE
            )
          }
        )
      }
      
      truth = query_data$mu_true
      err_tensor = mu_hat_tensor - truth
      err_matrix = mu_hat_matrix - truth
      
      diagnostics[[counter]] = rbind(
        data.frame(
          target_id = m,
          target_label = target_label,
          target_k = k_m,
          target_a = a_m,
          target_b = b_m,
          tau = tau,
          query = seq_len(n_queries),
          method = "Tensor-pooled",
          mu_true = truth,
          mu_hat = mu_hat_tensor,
          error = err_tensor,
          abs_error = abs(err_tensor),
          sq_error = err_tensor^2
        ),
        data.frame(
          target_id = m,
          target_label = target_label,
          target_k = k_m,
          target_a = a_m,
          target_b = b_m,
          tau = tau,
          query = seq_len(n_queries),
          method = "Matrix-only",
          mu_true = truth,
          mu_hat = mu_hat_matrix,
          error = err_matrix,
          abs_error = abs(err_matrix),
          sq_error = err_matrix^2
        )
      )
      
      counter = counter + 1
    }
  }
  
  diagnostics = do.call(rbind, diagnostics)
  
  diagnostics$method = factor(
    diagnostics$method,
    levels = c("Tensor-pooled", "Matrix-only")
  )
  
  diagnostics$target_label = factor(
    diagnostics$target_label,
    levels = unique(diagnostics$target_label)
  )
  
  summary = diagnostics %>%
    group_by(
      target_id,
      target_label,
      target_k,
      target_a,
      target_b,
      tau,
      method
    ) %>%
    summarise(
      mean_abs_error = mean(abs_error),
      sd_abs_error = sd(abs_error),
      n_queries = n(),
      se_abs_error = sd_abs_error / sqrt(n_queries),
      median_abs_error = median(abs_error),
      rmse = sqrt(mean(sq_error)),
      bias = mean(error),
      .groups = "drop"
    )
  
  list(
    diagnostics = diagnostics,
    summary = summary,
    targets = targets,
    tau_seq = tau_seq
  )
}


make_missingness_panel_plots = function(
    sim,
    targets,
    target_k = 1
) {
  K = dim(sim$Omega)[3]
  blocks = get_target_blocks(sim, target_k)
  row_perm = blocks$row_perm
  
  layers_to_plot = unique(c(
    target_k,
    if (K >= 2) 2 else integer(0),
    K
  ))
  
  target_rows = targets[targets$k == target_k, , drop = FALSE]
  
  if (nrow(target_rows) == 0) {
    stop("No target blocks belong to target_k.")
  }
  
  target_rects = do.call(
    rbind,
    lapply(seq_len(nrow(target_rows)), function(m) {
      aa = target_rows$a[m]
      bb = target_rows$b[m]
      
      rows = blocks$row_parts[[aa]]
      cols = blocks$col_parts[[bb]]
      y_rows = nrow(sim$Omega[, , target_k]) - rows + 1
      
      data.frame(
        xmin = min(cols) - 0.5,
        xmax = max(cols) + 0.5,
        ymin = min(y_rows) - 0.5,
        ymax = max(y_rows) + 0.5,
        x = mean(range(cols)),
        y = mean(range(y_rows)),
        label = paste0(
          "B", m, "\n(", aa - 1, ",", bb, ")"
        )
      )
    })
  )
  
  row_cuts = cumsum(blocks$N_part)
  col_cuts = cumsum(blocks$T_part)
  
  lapply(seq_along(layers_to_plot), function(ell) {
    j = layers_to_plot[ell]
    Omega_j = sim$Omega[row_perm, , j, drop = FALSE][, , 1]
    
    N = nrow(Omega_j)
    Tt = ncol(Omega_j)
    
    missing_data = expand.grid(
      row_order = seq_len(N),
      time = seq_len(Tt),
      KEEP.OUT.ATTRS = FALSE
    )
    
    missing_data$status = factor(
      ifelse(as.vector(Omega_j), "Observed", "Missing"),
      levels = c("Missing", "Observed")
    )
    
    # Put the first row after rearrangement at the top of the image.
    missing_data$unit_plot = N - missing_data$row_order + 1
    
    time_labels = pretty(seq_len(Tt))
    time_labels = time_labels[time_labels >= 1 & time_labels <= Tt]
    
    unit_labels = pretty(seq_len(N))
    unit_labels = unit_labels[unit_labels >= 1 & unit_labels <= N]
    unit_breaks = N - unit_labels + 1
    
    panel_title = if (j == target_k) {
      paste0("Target layer ", j)
    } else {
      paste0("Layer ", j)
    }
    
    panel_subtitle = if (j == target_k) {
      "Blue = observed; red = missing; yellow = evaluated blocks"
    } else {
      NULL
    }
    
    p_missing = ggplot(
      missing_data,
      aes(x = time, y = unit_plot, fill = status)
    ) +
      geom_raster(interpolate = FALSE) +
      scale_fill_manual(
        values = c("Missing" = "red", "Observed" = "blue"),
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = time_labels,
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = unit_breaks,
        labels = unit_labels,
        expand = c(0, 0)
      ) +
      coord_cartesian(
        xlim = c(0.5, Tt + 0.5),
        ylim = c(0.5, N + 0.5),
        expand = FALSE
      ) +
      labs(
        title = panel_title,
        subtitle = panel_subtitle,
        x = "Time",
        y = "Unit"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(
          size = 10.5,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          size = 7.5,
          hjust = 0.5
        ),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    if (j == target_k) {
      p_missing = p_missing +
        geom_vline(
          xintercept = col_cuts + 0.5,
          color = "black",
          linewidth = 0.32
        ) +
        geom_hline(
          yintercept = N - row_cuts + 0.5,
          color = "black",
          linewidth = 0.32
        ) +
        geom_rect(
          data = target_rects,
          aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax
          ),
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
    
    p_missing
  })
}


draw_missingness_and_tau_figure = function(
    missingness_plots,
    tau_plot
) {
  n_top = length(missingness_plots)
  
  grid.newpage()
  
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = n_top,
        heights = unit(c(1.0, 1.65), "null"),
        widths = unit(rep(1, n_top), "null")
      )
    )
  )
  
  for (ell in seq_len(n_top)) {
    print(
      missingness_plots[[ell]],
      newpage = FALSE,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = ell
      )
    )
  }
  
  print(
    tau_plot,
    newpage = FALSE,
    vp = viewport(
      layout.pos.row = 2,
      layout.pos.col = seq_len(n_top)
    )
  )
  
  popViewport()
  
  invisible(NULL)
}


plot_tau_sensitivity_12_lines = function(
    tau_out,
    filename_prefix = "tau_sensitivity_six_missing_blocks_12_lines",
    add_standard_error = TRUE
) {
  plot_data = tau_out$summary %>%
    mutate(
      method = factor(
        method,
        levels = c("Tensor-pooled", "Matrix-only")
      ),
      target_label = factor(
        target_label,
        levels = unique(target_label)
      )
    )
  
  if (nlevels(plot_data$target_label) != 6) {
    stop("The plot expects exactly six missing target blocks.")
  }
  
  method_cols = c(
    "Tensor-pooled" = "steelblue",
    "Matrix-only" = "orange"
  )
  
  block_shapes = setNames(
    c(16, 17, 15, 18, 8, 4),
    levels(plot_data$target_label)
  )
  
  base_theme = theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(
        color = "gray88",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5
      ),
      
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(
        size = 9,
        face = "bold"
      ),
      legend.text = element_text(size = 8.5),
      legend.key.width = unit(1.25, "lines"),
      legend.spacing.x = unit(0.7, "lines"),
      legend.spacing.y = unit(0.15, "lines"),
      
      plot.margin = margin(6, 6, 6, 6)
    )
  
  p_tau = ggplot(
    plot_data,
    aes(
      x = tau,
      y = mean_abs_error,
      color = method,
      shape = target_label,
      group = interaction(method, target_label)
    )
  )
  
  if (add_standard_error) {
    p_tau = p_tau +
      geom_errorbar(
        aes(
          ymin = pmax(
            0,
            mean_abs_error - se_abs_error
          ),
          ymax = mean_abs_error + se_abs_error
        ),
        width = 0,
        linewidth = 0.30,
        alpha = 0.35,
        show.legend = FALSE
      )
  }
  
  tau_breaks = 10^seq(
    floor(log10(min(plot_data$tau))),
    ceiling(log10(max(plot_data$tau))),
    by = 1
  )
  
  p_tau = p_tau +
    geom_line(
      linewidth = 0.78,
      linetype = "solid"
    ) +
    
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
      labels = trans_format(
        "log10",
        math_format(10^.x)
      )
    ) +
    
    scale_y_continuous(
      labels = scientific
    ) +
    
    labs(
      title = "Sensitivity to the regularisation parameter",
      subtitle = if (add_standard_error) {
        expression(
          "Mean absolute error " %+-% 
            " standard error across queries"
        )
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
        override.aes = list(
          shape = 16,
          linetype = 1,
          linewidth = 0.9,
          size = 2.4
        )
      ),
      
      shape = guide_legend(
        order = 2,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          color = "gray25",
          linetype = 0,
          linewidth = 0,
          size = 2.8
        )
      )
    ) +
    
    base_theme
  
  png_path = file.path(
    fig_dir,
    paste0(filename_prefix, ".png")
  )
  
  pdf_path = file.path(
    fig_dir,
    paste0(filename_prefix, ".pdf")
  )
  
  ggsave(
    filename = png_path,
    plot = p_tau,
    width = 10.5,
    height = 6.2,
    units = "in",
    dpi = 400,
    bg = "white"
  )
  
  ggsave(
    filename = pdf_path,
    plot = p_tau,
    width = 10.5,
    height = 6.2,
    units = "in",
    bg = "white",
    device = cairo_pdf
  )
  
  print(p_tau)
  
  message("Saved:")
  message(png_path)
  message(pdf_path)
  
  invisible(list(
    tau_plot = p_tau,
    png_path = png_path,
    pdf_path = pdf_path
  ))
}

# ===============================================================
# 10. Run the six-block tau experiment
# ===============================================================

experiment_seed = 221198
query_seed = 222198

K_main = 10
n_adopt_times_main = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5)

tau_grid_main = 10^seq(-0.5, 0, by = 0.025)

n_queries_main = 100
sigma_main = 0.02

sim_tau = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = K_main,
  r = 5,
  sigma = sigma_main,
  n_adopt_times = n_adopt_times_main,
  p_never = 0.20,
  p_initial = 0.10,
  seed = experiment_seed
)

targets_tau = make_missing_targets(
  sim = sim_tau,
  k = 1,
  n_targets = 6
)

print(targets_tau)

tau_out = run_tau_sensitivity_compare(
  sim = sim_tau,
  targets = targets_tau,
  k = 1,
  n_targets = 6,
  tau_seq = tau_grid_main,
  n_queries = n_queries_main,
  query_seed = query_seed,
  verbose = TRUE
)

write.csv(
  tau_out$diagnostics,
  file = file.path(res_dir, "tau_sensitivity_diagnostics_six_missing_blocks.csv"),
  row.names = FALSE
)

write.csv(
  tau_out$summary,
  file = file.path(res_dir, "tau_sensitivity_summary_six_missing_blocks.csv"),
  row.names = FALSE
)

combined_tau_fig = plot_tau_sensitivity_12_lines(
  tau_out = tau_out,
  filename_prefix = "tau_sensitivity_six_missing_blocks_12_lines",
  add_standard_error = TRUE
)

cat("\nOverall mean absolute error across tau values and blocks\n")
cat("--------------------------------------------------\n")
print(
  tau_out$diagnostics %>%
    group_by(method) %>%
    summarise(
      mean_abs_error = mean(abs_error),
      median_abs_error = median(abs_error),
      rmse = sqrt(mean(sq_error)),
      .groups = "drop"
    )
)
