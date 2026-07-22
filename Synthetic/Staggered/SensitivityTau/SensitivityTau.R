setwd("~/Desktop/code")

source("bilinearTensorAllFunction.R")

library(dplyr)

results_dir = "Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


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
  
  rand_orth = function(n, rank) {
    qr.Q(qr(matrix(rnorm(n * rank), nrow = n)))[, seq_len(rank), drop = FALSE]
  }
  
  U = rand_orth(N, r)
  V = rand_orth(Tt, r)
  C_core = array(rnorm(r * r * K), dim = c(r, r, K))
  
  M = array(0, dim = c(N, Tt, K))
  for (k in seq_len(K)) {
    M[, , k] = U %*% C_core[, , k] %*% t(V)
  }
  
  E = array(rnorm(N * Tt * K, sd = sigma), dim = c(N, Tt, K))
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
      A[i, k] = if (runif(1) < p_never) Inf else sample(adopt_grid, 1)
      
      if (is.infinite(A[i, k])) {
        Omega[i, , k] = TRUE
      } else {
        Omega[i, seq_len(A[i, k] - 1), k] = TRUE
      }
    }
    
    Y_obs[, , k][Omega[, , k]] = Y_full[, , k][Omega[, , k]]
  }
  
  list(
    Y_obs = Y_obs,
    M = M,
    Omega = Omega,
    A = A,
    r = r,
    sigma = sigma
  )
}


get_target_blocks = function(sim, k) {
  Tt = dim(sim$Y_obs)[2]
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
    o_k = length(m_desc),
    row_parts = row_parts,
    col_parts = Map(seq, starts, ends),
    N_part = vapply(row_parts, length, integer(1)),
    T_part = T_part
  )
}


make_missing_targets = function(sim, k = 1, n_targets = 6) {
  o_k = get_target_blocks(sim, k)$o_k
  targets = expand.grid(k = k, a = seq_len(o_k), b = seq_len(o_k))
  targets = targets[targets$a + targets$b > o_k + 1, , drop = FALSE]
  
  if (nrow(targets) < n_targets) {
    stop("Not enough missing target blocks.")
  }
  
  targets[seq_len(n_targets), , drop = FALSE]
}


prepare_target_queries = function(sim, k, a, b, n_queries, seed) {
  blocks = get_target_blocks(sim, k)
  row_idx = blocks$row_perm[blocks$row_parts[[a]]]
  col_idx = blocks$col_parts[[b]]
  M_target = sim$M[row_idx, col_idx, k, drop = FALSE][, , 1]
  
  set.seed(seed)
  
  x_list = vector("list", n_queries)
  y_list = vector("list", n_queries)
  mu_true = numeric(n_queries)
  
  for (q in seq_len(n_queries)) {
    x = rnorm(blocks$N_part[a])
    y = rnorm(blocks$T_part[b])
    x = x / sqrt(sum(x^2))
    y = y / sqrt(sum(y^2))
    
    x_list[[q]] = x
    y_list[[q]] = y
    mu_true[q] = as.numeric(t(x) %*% M_target %*% y)
  }
  
  list(x = x_list, y = y_list, mu_true = mu_true)
}


run_tau_sensitivity_compare = function(
    sim,
    targets,
    tau_seq,
    n_queries = 100,
    query_seed = 999
) {
  results = vector("list", nrow(targets) * length(tau_seq))
  result_index = 1
  
  for (m in seq_len(nrow(targets))) {
    k = targets$k[m]
    a = targets$a[m]
    b = targets$b[m]
    
    queries = prepare_target_queries(
      sim = sim,
      k = k,
      a = a,
      b = b,
      n_queries = n_queries,
      seed = query_seed + m
    )
    
    for (tau in tau_seq) {
      mu_hat_tensor = numeric(n_queries)
      mu_hat_matrix = numeric(n_queries)
      
      for (q in seq_len(n_queries)) {
        mu_hat_tensor[q] = bilinearTensorStaggered(
          Y = sim$Y_obs,
          k = k,
          a = a,
          b = b,
          r = sim$r,
          x = queries$x[[q]],
          y = queries$y[[q]],
          A = sim$A,
          Omega = sim$Omega,
          tau = tau
        )
        
        mu_hat_matrix[q] = bilinearMatrixStaggered(
          Y_mat = sim$Y_obs[, , k],
          a = a,
          b = b,
          r = sim$r,
          x = queries$x[[q]],
          y = queries$y[[q]],
          A = sim$A[, k],
          Omega = sim$Omega[, , k],
          tau = tau
        )
      }
      
      results[[result_index]] = bind_rows(
        data.frame(
          target_id = m,
          target_k = k,
          target_a = a,
          target_b = b,
          tau = tau,
          query = seq_len(n_queries),
          method = "Tensor-pooled",
          mu_true = queries$mu_true,
          mu_hat = mu_hat_tensor
        ),
        data.frame(
          target_id = m,
          target_k = k,
          target_a = a,
          target_b = b,
          tau = tau,
          query = seq_len(n_queries),
          method = "Matrix-only",
          mu_true = queries$mu_true,
          mu_hat = mu_hat_matrix
        )
      ) %>%
        mutate(
          error = mu_hat - mu_true,
          abs_error = abs(error),
          sq_error = error^2
        )
      
      result_index = result_index + 1
    }
  }
  
  diagnostics = bind_rows(results)
  
  summary = diagnostics %>%
    group_by(target_id, target_k, target_a, target_b, tau, method) %>%
    summarise(
      mean_abs_error = mean(abs_error),
      sd_abs_error = sd(abs_error),
      se_abs_error = sd_abs_error / sqrt(n()),
      n_queries = n(),
      .groups = "drop"
    )
  
  overall_summary = diagnostics %>%
    group_by(method) %>%
    summarise(
      mean_abs_error = mean(abs_error),
      rmse = sqrt(mean(sq_error)),
      bias = mean(error),
      .groups = "drop"
    )
  
  list(
    diagnostics = diagnostics,
    summary = summary,
    overall_summary = overall_summary,
    targets = targets,
    tau_seq = tau_seq
  )
}


experiment_seed = 221198
query_seed = 222198
n_adopt_times = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 5)
tau_grid = 10^seq(-0.5, 0, by = 0.025)
n_queries = 100

sim = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = 10,
  r = 5,
  sigma = 0.03,
  n_adopt_times = n_adopt_times,
  p_never = 0.20,
  p_initial = 0.10,
  seed = experiment_seed
)

targets = make_missing_targets(sim, k = 1, n_targets = 6)

tau_results = run_tau_sensitivity_compare(
  sim = sim,
  targets = targets,
  tau_seq = tau_grid,
  n_queries = n_queries,
  query_seed = query_seed
)

config = data.frame(
  experiment_seed = experiment_seed,
  query_seed = query_seed,
  N = 150,
  Tt = 200,
  K = 10,
  r = 5,
  sigma = 0.02,
  p_never = 0.20,
  p_initial = 0.10,
  n_queries = n_queries,
  tau_min = min(tau_grid),
  tau_max = max(tau_grid),
  tau_count = length(tau_grid)
)

write.csv(
  tau_results$diagnostics,
  file.path(results_dir, "tau_sensitivity_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  tau_results$summary,
  file.path(results_dir, "tau_sensitivity_summary.csv"),
  row.names = FALSE
)
write.csv(
  tau_results$overall_summary,
  file.path(results_dir, "tau_sensitivity_overall_summary.csv"),
  row.names = FALSE
)
write.csv(
  tau_results$targets,
  file.path(results_dir, "tau_sensitivity_targets.csv"),
  row.names = FALSE
)
write.csv(
  config,
  file.path(results_dir, "tau_sensitivity_config.csv"),
  row.names = FALSE
)
saveRDS(
  tau_results,
  file.path(results_dir, "tau_sensitivity_results.rds")
)

message("Results saved to: ", normalizePath(results_dir))
