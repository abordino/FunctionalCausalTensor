#------------------------------------------------------------------------------
# Tucker2 tensor simulation with general staggered-adoption missingness
#------------------------------------------------------------------------------

simulate_general_staggered_tucker2 = function(
    N = 150,
    Tt = 200,
    K = 3,
    r = 5,
    sigma = 1,
    n_adopt_times = c(4, 5, 6),
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
  
  # -------------------------------------------------------------
  # General staggered-adoption missingness
  # -------------------------------------------------------------
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
# Construct target-layer staircase rearrangement
#------------------------------------------------------------------------------

rearrange_by_target_layer = function(sim, k) {
  Y_obs = sim$Y_obs
  Y_full = sim$Y_full
  M = sim$M
  E = sim$E
  Omega = sim$Omega
  A = sim$A
  
  K = dim(Y_obs)[3]
  N = dim(Y_obs)[1]
  Tt = dim(Y_obs)[2]
  
  stopifnot(k >= 1, k <= K)
  
  obs_len_k = ifelse(is.infinite(A[, k]), Tt, pmin(Tt, A[, k] - 1))
  row_perm = order(-obs_len_k, seq_len(N))
  
  Y_obs_perm = Y_obs[row_perm, , , drop = FALSE]
  Y_full_perm = Y_full[row_perm, , , drop = FALSE]
  M_perm = M[row_perm, , , drop = FALSE]
  E_perm = E[row_perm, , , drop = FALSE]
  Omega_perm = Omega[row_perm, , , drop = FALSE]
  A_perm = A[row_perm, , drop = FALSE]
  
  obs_len_perm = obs_len_k[row_perm]
  m_desc = unique(obs_len_perm)
  
  o_k = length(m_desc)
  row_parts_k = lapply(m_desc, function(m) which(obs_len_perm == m))
  
  m_asc = rev(m_desc)
  T_part_k = diff(c(0, m_asc))
  
  ends = cumsum(T_part_k)
  starts = c(1, head(ends, -1) + 1)
  col_parts_k = Map(seq, starts, ends)
  
  out = sim
  
  out$Y_obs = Y_obs_perm
  out$Y_full = Y_full_perm
  out$M = M_perm
  out$E = E_perm
  out$Omega = Omega_perm
  out$A = A_perm
  
  out$target_layer = k
  out$row_perm = row_perm
  out$o_k = o_k
  out$N_part_k = vapply(row_parts_k, length, integer(1))
  out$T_part_k = T_part_k
  out$row_parts_k = row_parts_k
  out$col_parts_k = col_parts_k
  
  class(out) = "target_rearranged_staggered_tucker2_sim"
  
  out
}


#------------------------------------------------------------------------------
# Test whether a layer is staircase in current row ordering
#------------------------------------------------------------------------------

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


# ===============================================================
# Test rearrangement by common target-layer permutation
# ===============================================================

test_target_rearrangement = function(sim, k) {
  sim_perm = rearrange_by_target_layer(sim, k)
  
  K = dim(sim$Omega)[3]
  
  before = logical(K)
  after = logical(K)
  
  for (j in seq_len(K)) {
    before[j] = is_staircase_layer(sim$Omega[, , j])
    after[j] = is_staircase_layer(sim_perm$Omega[, , j])
  }
  
  out = data.frame(
    layer = seq_len(K),
    staircase_before = before,
    staircase_after_target_permutation = after
  )
  
  cat("\n")
  cat("Target layer:", k, "\n")
  cat("Target layer staircase after permutation:",
      after[k], "\n")
  cat("\n")
  
  print(out)
  
  invisible(list(
    summary = out,
    sim_rearranged = sim_perm
  ))
}


# ===============================================================
# Print adoption-time summary
# ===============================================================

print_general_staggered_missingness = function(sim) {
  K = dim(sim$Omega)[3]
  Tt = dim(sim$Omega)[2]
  
  for (k in seq_len(K)) {
    A_k = sim$A[, k]
    
    cat("\n")
    cat("============================================\n")
    cat("Layer k =", k, "\n")
    cat("Distinct adoption times:",
        paste(sort(unique(A_k[is.finite(A_k)])), collapse = " "),
        "\n")
    cat("Never adopters:", sum(is.infinite(A_k)), "\n")
    cat("Observed entries:",
        sum(sim$Omega[, , k]),
        "out of",
        length(sim$Omega[, , k]),
        "\n")
    cat("Staircase in current row ordering:",
        is_staircase_layer(sim$Omega[, , k]),
        "\n")
  }
  
  invisible(NULL)
}


# ===============================================================
# Plot full N x T missingness matrix for one layer
# Blue = observed
# Red  = missing
# ===============================================================

plot_missingness_layer = function(
    sim,
    layer = 1,
    show_adoption_lines = FALSE,
    main = NULL
) {
  Omega_k = sim$Omega[, , layer]
  
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
  
  if (show_adoption_lines && !is.null(sim$row_parts_k) && layer == sim$target_layer) {
    row_cuts = cumsum(sim$N_part_k)
    col_cuts = cumsum(sim$T_part_k)
    
    abline(v = col_cuts + 0.5, col = "black", lwd = 1)
    abline(h = N - row_cuts + 0.5, col = "black", lwd = 1)
  }
  
  invisible(NULL)
}


#------------------------------------------------------------------------------
# Plot all layers
# ===============================================================

plot_all_missingness_layers = function(
    sim,
    show_adoption_lines = FALSE,
    main_prefix = "Layer "
) {
  K = dim(sim$Omega)[3]
  
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  nrow_plot = ceiling(sqrt(K))
  ncol_plot = ceiling(K / nrow_plot)
  
  par(mfrow = c(nrow_plot, ncol_plot), mar = c(4, 4, 3, 1))
  
  for (k in seq_len(K)) {
    plot_missingness_layer(
      sim = sim,
      layer = k,
      show_adoption_lines = show_adoption_lines,
      main = paste0(main_prefix, k)
    )
  }
  
  invisible(NULL)
}


# ===============================================================
# Plot before and after target-layer rearrangement
# ===============================================================

plot_before_after_rearrangement = function(sim, k) {
  sim_perm = rearrange_by_target_layer(sim, k)
  
  K = dim(sim$Omega)[3]
  
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
      sim = sim_perm,
      layer = j,
      show_adoption_lines = (j == k),
      main = paste0("After perm by k=", k, ", layer ", j)
    )
  }
  
  invisible(sim_perm)
}


# ===============================================================
# Plot the observed tensor values for one layer
# ===============================================================

plot_tensor_values_layer = function(
    sim,
    layer = 1,
    main = NULL
) {
  Y_k = sim$Y_obs[, , layer]
  
  N = nrow(Y_k)
  Tt = ncol(Y_k)
  
  if (is.null(main)) {
    main = paste0("Observed tensor values, layer ", layer)
  }
  
  image(
    x = seq_len(Tt),
    y = seq_len(N),
    z = t(Y_k[N:1, ]),
    axes = FALSE,
    xlab = "Time index t",
    ylab = "Unit index i",
    main = main,
    col = heat.colors(100),
    useRaster = TRUE
  )
  
  axis(1)
  axis(2, at = pretty(seq_len(N)), labels = rev(pretty(seq_len(N))))
  box()
  
  invisible(NULL)
}


# ===============================================================
# Example run
# ===============================================================

sim = simulate_general_staggered_tucker2(
  N = 150,
  Tt = 200,
  K = 3,
  r = 5,
  sigma = 1,
  n_adopt_times = c(4, 5, 6),
  p_never = 0.20,
  p_initial = 0.10,
  seed = 123
)

print_general_staggered_missingness(sim)

plot_all_missingness_layers(sim)

test = test_target_rearrangement(sim, k = 1)

sim_perm = test$sim_rearranged

plot_before_after_rearrangement(sim, k = 1)

plot_missingness_layer(
  sim = sim_perm,
  layer = 1,
  show_adoption_lines = TRUE,
  main = "Target layer after common row permutation"
)

plot_tensor_values_layer(sim, layer = 1)