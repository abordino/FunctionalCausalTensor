setwd("~/Documents/phd/projects/causalMatrix/code")

source("bilinearTensorStaggered.R")
source("bilinearMatrixStaggered.R")
source("bilinearTensorStaggeredPsi.R")
source("bilinearMatrixStaggeredPsi.R")
source("pluginPsi_c1.R")


# ===============================================================
# 1. Setup
# ===============================================================

set.seed(221198)

N = 8
Tt = 8
K = 3
target_k = 1

r = 1
tau = 1e-10


# ===============================================================
# 2. Inclusive adoption times
# ===============================================================

A_inclusive = matrix(NA_real_, nrow = N, ncol = K)

A_target_sorted = c(Inf, Inf, 6, 6, 4, 4, 2, 2)

row_shuffle = c(5, 1, 7, 3, 8, 2, 4, 6)

A_inclusive[, 1] = A_target_sorted[row_shuffle]

A_inclusive[, 2] = sample(c(Inf, Inf, 6, 6, 4, 4, 2, 2), size = N, replace = FALSE)
A_inclusive[, 3] = sample(c(Inf, Inf, 6, 6, 4, 4, 2, 2), size = N, replace = FALSE)

cat("\nInclusive adoption-time matrix A_inclusive\n")
cat("Toy convention: observed if t <= A_inclusive[i,k]\n")
print(A_inclusive)


# ===============================================================
# 3. Convert to estimator convention
# ===============================================================

A_exclusive = A_inclusive

finite_idx = is.finite(A_exclusive)
A_exclusive[finite_idx] = A_exclusive[finite_idx] + 1

cat("\nExclusive adoption-time matrix A_exclusive\n")
cat("Estimator convention: observed if t < A_exclusive[i,k]\n")
print(A_exclusive)


# ===============================================================
# 4. Build Omega and noiseless observed Y
# ===============================================================

Omega = array(FALSE, dim = c(N, Tt, K))

for (k in seq_len(K)) {
  for (i in seq_len(N)) {
    if (is.infinite(A_exclusive[i, k])) {
      Omega[i, , k] = TRUE
    } else if (A_exclusive[i, k] > 1) {
      Omega[i, seq_len(A_exclusive[i, k] - 1), k] = TRUE
    }
  }
}


# ===============================================================
# 5. Random rank-one signal, no noise
# ===============================================================

u = rnorm(N)
v = rnorm(Tt)

u = u / sqrt(sum(u^2))
v = v / sqrt(sum(v^2))

theta = c(4, -2, 3)

M = array(0, dim = c(N, Tt, K))

for (k in seq_len(K)) {
  M[, , k] = theta[k] * tcrossprod(u, v)
}

Y_full = M
Y_obs = M
Y_obs[!Omega] = NA_real_


# ===============================================================
# 6. Target-layer row ordering and staircase blocks
# ===============================================================

get_target_plugin_inputs_from_exclusive_A = function(A, Tt, k) {
  A_k = A[, k]
  
  obs_len = ifelse(
    is.infinite(A_k),
    Tt,
    pmin(Tt, A_k - 1)
  )
  
  row_perm = order(-obs_len, seq_along(obs_len))
  obs_len_sorted = obs_len[row_perm]
  
  m_desc = unique(obs_len_sorted)
  
  row_parts = lapply(m_desc, function(m) which(obs_len_sorted == m))
  
  m_asc = rev(m_desc)
  T_part = diff(c(0, m_asc))
  
  make_parts = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  col_parts = make_parts(T_part)
  
  N_part = vapply(row_parts, length, integer(1))
  
  list(
    row_perm = row_perm,
    obs_len = obs_len,
    obs_len_sorted = obs_len_sorted,
    row_parts = row_parts,
    col_parts = col_parts,
    N_part = N_part,
    T_part = T_part,
    N_parts = list(N_part),
    T_parts = list(T_part),
    o_k = length(row_parts)
  )
}

blocks = get_target_plugin_inputs_from_exclusive_A(
  A = A_exclusive,
  Tt = Tt,
  k = target_k
)

row_perm = blocks$row_perm

cat("\nTarget-layer observed prefix lengths before ordering\n")
print(blocks$obs_len)

cat("\nTarget-layer row permutation\n")
print(row_perm)

cat("\nTarget-layer observed prefix lengths after ordering\n")
print(blocks$obs_len_sorted)

cat("\nN_parts\n")
print(blocks$N_part)

cat("\nT_parts\n")
print(blocks$T_part)

cat("\nTarget-layer mask before ordering\n")
print(1 * Omega[, , target_k])

cat("\nTarget-layer mask after ordering\n")
print(1 * Omega[row_perm, , target_k])


# ===============================================================
# 7. Prepare sorted target layer for plug-in
# ===============================================================

M_target_sorted = M[row_perm, , target_k]

N_parts_plugin = blocks$N_parts
T_parts_plugin = blocks$T_parts
k_plugin = 1

eta_original = rep(c(1, -1), length.out = N)

eta_sorted = eta_original[row_perm]

row_index_sorted = blocks$row_parts[[3]][1]
row_index_original = row_perm[row_index_sorted]

cat("\nRowHet eta in original order\n")
print(eta_original)

cat("\nRowHet eta in sorted target-layer order\n")
print(eta_sorted)

cat("\nLocal row\n")
cat("row_index_sorted   =", row_index_sorted, "\n")
cat("row_index_original =", row_index_original, "\n")


# ===============================================================
# 8. Compute plug-in, matrix, and tensor Psi
# ===============================================================

functionals = c("ATE", "RowHet", "Local", "Trend")

toy_compare = data.frame(
  functional = character(0),
  plugin = numeric(0),
  matrix = numeric(0),
  tensor = numeric(0),
  matrix_minus_plugin = numeric(0),
  tensor_minus_plugin = numeric(0),
  stringsAsFactors = FALSE
)

for (h in functionals) {
  
  plugin_val = pluginPsi_c1(
    Y = M_target_sorted,
    k = k_plugin,
    N_parts = N_parts_plugin,
    T_parts = T_parts_plugin,
    functional = h,
    eta = if (h == "RowHet") eta_sorted else NULL,
    row_index = if (h == "Local") row_index_sorted else NULL
  )
  
  matrix_val = switch(
    h,
    ATE = bilinearMatrixStaggeredATE(
      Y_mat = Y_obs[, , target_k],
      r = r,
      tau = tau,
      A = A_exclusive[, target_k],
      Omega = Omega[, , target_k]
    ),
    RowHet = bilinearMatrixStaggeredRowHet(
      Y_mat = Y_obs[, , target_k],
      r = r,
      eta = eta_original,
      tau = tau,
      A = A_exclusive[, target_k],
      Omega = Omega[, , target_k]
    ),
    Local = bilinearMatrixStaggeredLocal(
      Y_mat = Y_obs[, , target_k],
      r = r,
      row_index = row_index_original,
      tau = tau,
      A = A_exclusive[, target_k],
      Omega = Omega[, , target_k]
    ),
    Trend = bilinearMatrixStaggeredTrend(
      Y_mat = Y_obs[, , target_k],
      r = r,
      tau = tau,
      A = A_exclusive[, target_k],
      Omega = Omega[, , target_k]
    )
  )
  
  tensor_val = switch(
    h,
    ATE = bilinearTensorStaggeredATE(
      Y = Y_obs,
      k = target_k,
      r = r,
      tau = tau,
      A = A_exclusive,
      Omega = Omega
    ),
    RowHet = bilinearTensorStaggeredRowHet(
      Y = Y_obs,
      k = target_k,
      r = r,
      eta = eta_original,
      tau = tau,
      A = A_exclusive,
      Omega = Omega
    ),
    Local = bilinearTensorStaggeredLocal(
      Y = Y_obs,
      k = target_k,
      r = r,
      row_index = row_index_original,
      tau = tau,
      A = A_exclusive,
      Omega = Omega
    ),
    Trend = bilinearTensorStaggeredTrend(
      Y = Y_obs,
      k = target_k,
      r = r,
      tau = tau,
      A = A_exclusive,
      Omega = Omega
    )
  )
  
  toy_compare = rbind(
    toy_compare,
    data.frame(
      functional = h,
      plugin = plugin_val,
      matrix = matrix_val,
      tensor = tensor_val,
      matrix_minus_plugin = matrix_val - plugin_val,
      tensor_minus_plugin = tensor_val - plugin_val
    )
  )
}


cat("\nComparison: plug-in vs matrix vs tensor, no noise\n")
print(toy_compare)


# ===============================================================
# 9. Numerical equality checks
# ===============================================================

tol = 1e-6

cat("\nMaximum absolute difference: matrix vs plug-in\n")
print(max(abs(toy_compare$matrix_minus_plugin)))

cat("\nMaximum absolute difference: tensor vs plug-in\n")
print(max(abs(toy_compare$tensor_minus_plugin)))

cat("\nAll matrix values equal plug-in within tolerance?\n")
print(all(abs(toy_compare$matrix_minus_plugin) < tol))

cat("\nAll tensor values equal plug-in within tolerance?\n")
print(all(abs(toy_compare$tensor_minus_plugin) < tol))


# ===============================================================
# 10. Print target matrices to ease numerical checks
# ===============================================================

cat("\n==================================================\n")
cat("LAST PRINT: target matrix before ordering\n")
cat("==================================================\n")
print(round(M[, , target_k], 4))

cat("\nTarget observed matrix before ordering, NA = missing\n")
print(round(Y_obs[, , target_k], 4))

cat("\n==================================================\n")
cat("LAST PRINT: target matrix after ordering\n")
cat("==================================================\n")
print(round(M[row_perm, , target_k], 4))

cat("\nTarget observed matrix after ordering, NA = missing\n")
print(round(Y_obs[row_perm, , target_k], 4))

data.frame(
  functional = c("ATE", "RowHet", "Local", "Trend"),
  value = c(
    -0.3360958,
    -0.2375292,
    -0.3604500,
    -0.4326083
  )
)