setwd("~/Documents/phd/projects/causalMatrix/code")

source("bilinearTensorStaggered.R")
source("bilinearMatrixStaggered.R")


# ===============================================================
# 1. Simulate Tucker2 tensor with general staggered adoption
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
  cat("Median absolute error:", median(abs_err_tensor_vec), "\n")
  cat("Max absolute error   :", max(abs_err_tensor_vec), "\n")
  cat("Mean relative error  :", mean(rel_err_tensor_vec), "\n")
  cat("Median relative error:", median(rel_err_tensor_vec), "\n")
  cat("Correlation          :", cor(mu_hat_tensor_vec, mu_true_vec), "\n\n")
  
  cat("Matrix-only method\n")
  cat("Mean absolute error  :", mean(abs_err_matrix_vec), "\n")
  cat("Median absolute error:", median(abs_err_matrix_vec), "\n")
  cat("Max absolute error   :", max(abs_err_matrix_vec), "\n")
  cat("Mean relative error  :", mean(rel_err_matrix_vec), "\n")
  cat("Median relative error:", median(rel_err_matrix_vec), "\n")
  cat("Correlation          :", cor(mu_hat_matrix_vec, mu_true_vec), "\n\n")
  
  cat("Error ratio: matrix / tensor\n")
  cat("Mean absolute error ratio  :",
      mean(abs_err_matrix_vec) / max(mean(abs_err_tensor_vec), 1e-12),
      "\n")
  cat("Median absolute error ratio:",
      median(abs_err_matrix_vec) / max(median(abs_err_tensor_vec), 1e-12),
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
# 9. Run examples with K = 10
# ===============================================================

K = 10
n_adopt_times = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5)
# n_adopt_times = c(2, 5, 6, 4, 5, 6, 4, 5, 6, 5)

# ---------------------------------------------------------------
# Noiseless simulation
# ---------------------------------------------------------------

sim0 = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = K,
  r = 5,
  sigma = 0,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  seed = 123
)

print_general_staggered_missingness(sim0)

test_target_rearrangement(sim0, k = 1)

plot_before_after_rearrangement(sim0, k = 1)

target = choose_missing_target(sim0, k = 1)

k = target$k
a = target$a
b = target$b

cat("\nChosen target\n")
cat("--------------------------------------------------\n")
cat("k =", k, ", a =", a, ", b =", b, ", o_k =", target$o_k, "\n")

diagnostics0_compare = run_accuracy_check_compare(
  sim = sim0,
  k = k,
  a = a,
  b = b,
  n_queries = 100,
  tau = 1e-2,
  seed = 999
)

plot_accuracy_diagnostics_compare(diagnostics0_compare)


# ---------------------------------------------------------------
# Noisy simulation
# ---------------------------------------------------------------

sim = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = K,
  r = 5,
  sigma = 0.01,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  seed = 231198
)

print_general_staggered_missingness(sim)

test_target_rearrangement(sim, k = 1)

diagnostics_compare = run_accuracy_check_compare(
  sim = sim,
  k = k,
  a = a,
  b = b,
  n_queries = 100,
  tau = 1e-2,
  seed = 091132
)

plot_accuracy_diagnostics_compare(diagnostics_compare)


# ---------------------------------------------------------------
# Single scalar checks
# ---------------------------------------------------------------

res0 = test_bilinearStaggered_compare(
  N = 150,
  Tt = 200,
  K = K,
  r = 5,
  sigma = 0,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  k = k,
  a = a,
  b = b,
  tau = 1e-2,
  seed = 231198
)

res = test_bilinearStaggered_compare(
  N = 150,
  Tt = 200,
  K = K,
  r = 5,
  sigma = 0.01,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  k = k,
  a = a,
  b = b,
  tau = 1e-2,
  seed = 231198
)


# ===============================================================
# 10. Compare across multiple target blocks
# ===============================================================

experiment_seed = 221198

sim = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = K,
  r = 5,
  sigma = 0.01,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  seed = experiment_seed
)

# n_targets should be number of missing blocks
# targets = make_missing_targets(sim, k = 1, n_targets = 1)
targets = make_missing_targets(sim, k = 1, n_targets = 6)

diagnostics_targets = vector("list", nrow(targets))

for (m in seq_len(nrow(targets))) {
  cat("\n==================================================\n")
  cat("Running target", m, "\n")
  cat("k =", targets$k[m], "\n")
  cat("a =", targets$a[m], "\n")
  cat("b =", targets$b[m], "\n")
  cat("==================================================\n")
  
  diagnostics_targets[[m]] = run_accuracy_check_compare(
    sim = sim,
    k = targets$k[m],
    a = targets$a[m],
    b = targets$b[m],
    n_queries = 100,
    tau = 1e-2,
    seed = experiment_seed + 1000 + m
  )
  
  diagnostics_targets[[m]]$target_id = m
}

diagnostics_all = do.call(rbind, diagnostics_targets)


# ===============================================================
# Plot top and bottom figures
# ===============================================================

plot_missingness_and_abs_error_barplot = function(
    sim,
    diagnostics_all,
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
  
  n_top = length(layers_to_plot)
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  layout_mat = rbind(
    seq_len(n_top),
    rep(n_top + 1, n_top)
  )
  
  layout(layout_mat, heights = c(1, 1.15))
  
  par(mar = c(4, 4, 3, 1))
  
  for (ell in seq_along(layers_to_plot)) {
    j = layers_to_plot[ell]
    
    Omega_j = sim$Omega[row_perm, , j]
    
    N = nrow(Omega_j)
    Tt = ncol(Omega_j)
    
    Z = matrix(as.numeric(Omega_j), nrow = N, ncol = Tt)
    
    main_j = if (j == target_k) {
      paste0("Target layer ", j)
    } else if (j == 2 && K > 3) {
      paste0("Layer ", j)
    } else {
      paste0("Layer ", j)
    }
    
    image(
      x = seq_len(Tt),
      y = seq_len(N),
      z = t(Z[N:1, ]),
      col = c("red", "blue"),
      axes = FALSE,
      xlab = "Time",
      ylab = "Unit",
      main = main_j
    )
    
    axis(1)
    axis(2, at = pretty(seq_len(N)), labels = rev(pretty(seq_len(N))))
    box()
    
    if (j == target_k) {
      row_cuts = cumsum(blocks$N_part)
      col_cuts = cumsum(blocks$T_part)
      
      abline(v = col_cuts + 0.5, col = "black", lwd = 1)
      abline(h = N - row_cuts + 0.5, col = "black", lwd = 1)
      
      for (m in seq_len(nrow(targets))) {
        aa = targets$a[m] 
        bb = targets$b[m]
        
        rows = blocks$row_parts[[aa]]
        cols = blocks$col_parts[[bb]]
        
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
          labels = paste0("(", aa-1, ",", bb, ")"),
          col = "yellow",
          cex = 0.9,
          font = 2
        )
      }
    }
  }

  
  mean_err_by_target = aggregate(
    cbind(abs_error_tensor, abs_error_matrix) ~ target_id,
    data = diagnostics_all,
    FUN = mean
  )
  
  mean_err_by_target = mean_err_by_target[order(mean_err_by_target$target_id), ]
  
  target_labels = paste0(
    "(",
    targets$a[mean_err_by_target$target_id] -1,
    ",",
    targets$b[mean_err_by_target$target_id],
    ")"
  )
  
  par(mar = c(5, 4, 3, 1))
  
  barplot(
    t(as.matrix(mean_err_by_target[, c("abs_error_tensor", "abs_error_matrix")])),
    beside = TRUE,
    names.arg = target_labels,
    col = c("steelblue", "orange"),
    ylab = "Mean absolute error",
    xlab = "Target block (a,b)",
    main = "Mean absolute error by target block",
    ylim = c(
      0,
      1.15 * max(
        mean_err_by_target$abs_error_tensor,
        mean_err_by_target$abs_error_matrix
      )
    )
  )
  
  legend(
    "top",
    legend = c("tensor-pooled", "matrix-only"),
    fill = c("steelblue", "orange"),
    bty = "n"
  )
  
  invisible(NULL)
}


plot_missingness_and_abs_error_barplot(
  sim = sim,
  diagnostics_all = diagnostics_all,
  targets = targets,
  target_k = 1
)