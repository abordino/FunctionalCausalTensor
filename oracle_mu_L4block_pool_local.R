# ---------------------------------------------------------------
# oracle_mu_L4block_pool_local : oracle local version for a general query
#                              L = U_L %*% Sigma_L %*% t(V_L)
# ---------------------------------------------------------------
# Args:
#   Y   : either
#         (i)  a list of length K, each element an N x T matrix
#              with NA in the bottom-right block, OR
#         (ii) an array of dim c(N, T, K)
#   k   : target arm in {1, ..., K}
#   U   : N x r oracle left singular vectors
#   V   : T x r oracle right singular vectors
#   R_list : list of length K, each element an r x r matrix
#   U_L : N2^(k) x r_L matrix
#   Sigma_L : r_L x r_L diagonal matrix
#   V_L : T2^(k) x r_L matrix
#   N1, T1 : optional integer vectors of length K giving block sizes;
#            if NULL, detected from the NA pattern
#
# Returns:
#   scalar mu_hat = <L, M^{(d)}_{.,.,k}> plus the local oracle correction term
# ---------------------------------------------------------------

oracle_mu_L4block_pool_local = function(Y, k, U, V, R_list, U_L, Sigma_L, V_L, N1 = NULL, T1 = NULL) {
  if (is.array(Y)) {
    K = dim(Y)[3]
    Y = lapply(seq_len(K), function(j) Y[, , j])
  } else {
    K = length(Y)
  }
  
  N = nrow(Y[[1]])
  Tt = ncol(Y[[1]])
  r = nrow(R_list[[1]])
  
  detect_blocks = function(M) {
    row_has_NA = apply(M, 1, function(z) any(is.na(z)))
    col_has_NA = apply(M, 2, function(z) any(is.na(z)))
    list(
      N1 = if (any(row_has_NA)) which(row_has_NA)[1] - 1 else nrow(M),
      T1 = if (any(col_has_NA)) which(col_has_NA)[1] - 1 else ncol(M)
    )
  }
  
  if (is.null(N1) || is.null(T1)) {
    blk = lapply(Y, detect_blocks)
    N1 = vapply(blk, function(z) z$N1, integer(1))
    T1 = vapply(blk, function(z) z$T1, integer(1))
  }
  
  N1k = N1[k]
  T1k = T1[k]
  
  W_left = do.call(rbind, lapply(seq_len(K), function(j) {
    V1j = V[seq_len(T1[j]), seq_len(r), drop = FALSE]
    t(R_list[[j]] %*% t(V1j))
  }))
  
  W_up = do.call(rbind, lapply(seq_len(K), function(j) {
    U1j = U[seq_len(N1[j]), seq_len(r), drop = FALSE]
    U1j %*% R_list[[j]]
  }))
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  E_up_pool = Y_up_pool - W_up %*% t(V[, seq_len(r), drop = FALSE])
  
  Rk = R_list[[k]]
  U1k = U[seq_len(N1k), seq_len(r), drop = FALSE]
  U2k = U[(N1k + 1):N, seq_len(r), drop = FALSE]
  V2k = V[(T1k + 1):Tt, seq_len(r), drop = FALSE]
  
  alpha = crossprod(U2k, U_L)
  gamma = crossprod(V2k, V_L)
  
  mu_d = sum((alpha %*% Sigma_L) * (Rk %*% gamma))
  
  row_start = if (k == 1) 1 else sum(N1[seq_len(k - 1)]) + 1
  row_end = sum(N1[seq_len(k)])
  
  G_U1_inv = solve(crossprod(U1k))
  E_up_k = E_up_pool[row_start:row_end, , drop = FALSE]
  
  B3 = G_U1_inv %*% crossprod(U1k, E_up_k %*% V[, seq_len(r), drop = FALSE] %*% gamma)
  mu_3 = sum((alpha %*% Sigma_L) * B3)
  
  mu_d + mu_3
}