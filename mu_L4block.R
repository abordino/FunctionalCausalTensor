# ---------------------------------------------------------------
# mu_L4block_pool : implements Algorithm 1 for a general query
#                 L = U_L %*% Sigma_L %*% t(V_L)
# ---------------------------------------------------------------
# Args:
#   Y   : either
#         (i)  a list of length K, each element an N x T matrix
#              with NA in the bottom-right block, OR
#         (ii) an array of dim c(N, T, K)
#   k   : target arm in {1, ..., K}
#   r   : target rank
#   U_L : N2^(k) x r_L matrix
#   Sigma_L : r_L x r_L diagonal matrix
#   V_L : T2^(k) x r_L matrix
#   tau : threshold used in H_k^\dagger
#   N1, T1 : optional integer vectors of length K giving block sizes;
#            if NULL, detected from the NA pattern
#
# Returns:
#   scalar mu_hat = <L, M^{(d)}_{.,.,k}>
# ---------------------------------------------------------------

mu_L4block_pool = function(Y, k, r, U_L, Sigma_L, V_L, tau, N1 = NULL, T1 = NULL) {
  if (is.array(Y)) {
    K = dim(Y)[3]
    N = dim(Y)[1]
    Tt = dim(Y)[2]
    Y_list = lapply(seq_len(K), function(j) Y[, , j])
  } else {
    Y_list = Y
    K = length(Y_list)
    N = nrow(Y_list[[1]])
    Tt = ncol(Y_list[[1]])
  }
  
  detect_blocks = function(M) {
    row_has_NA = apply(M, 1, function(z) any(is.na(z)))
    col_has_NA = apply(M, 2, function(z) any(is.na(z)))
    list(
      N1 = if (any(row_has_NA)) which(row_has_NA)[1] - 1 else nrow(M),
      T1 = if (any(col_has_NA)) which(col_has_NA)[1] - 1 else ncol(M)
    )
  }
  
  if (is.null(N1) || is.null(T1)) {
    blk = lapply(Y_list, detect_blocks)
    N1 = vapply(blk, function(z) z$N1, integer(1))
    T1 = vapply(blk, function(z) z$T1, integer(1))
  }
  
  N1k = N1[k]
  T1k = T1[k]
  r_eff = min(r, N, Tt, sum(T1), sum(N1))
  
  Y_left_pool = do.call(cbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][, seq_len(T1[j]), drop = FALSE]
  }))
  
  s_left = svd(Y_left_pool, nu = r_eff, nv = 0)
  U_left = s_left$u[, seq_len(r_eff), drop = FALSE]
  
  U1k_hat = U_left[seq_len(N1k), , drop = FALSE]
  U2k_hat = U_left[(N1k + 1):N, , drop = FALSE]
  
  H_k = crossprod(U1k_hat)
  eH = eigen(H_k, symmetric = TRUE)
  H_k_dagger = eH$vectors %*%
    diag(1 / pmax(eH$values, tau), r_eff) %*%
    t(eH$vectors)
  
  alpha_hat = crossprod(U2k_hat, U_L)
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  
  s_up = svd(Y_up_pool, nu = r_eff, nv = r_eff)
  U_up = s_up$u[, seq_len(r_eff), drop = FALSE]
  V_up = s_up$v[, seq_len(r_eff), drop = FALSE]
  D_up = s_up$d[seq_len(r_eff)]
  
  s_k = if (k == 1) 0 else sum(N1[seq_len(k - 1)])
  U_up_k = U_up[(s_k + 1):(s_k + N1k), , drop = FALSE]
  V2k_hat = V_up[(T1k + 1):Tt, , drop = FALSE]
  
  T_L_hat = crossprod(V2k_hat, V_L)
  W_L = diag(D_up, r_eff) %*% T_L_hat
  X_L = U_up_k %*% W_L
  
  beta_hat = H_k_dagger %*% crossprod(U1k_hat, X_L)
  
  sum((alpha_hat %*% Sigma_L) * beta_hat)
}