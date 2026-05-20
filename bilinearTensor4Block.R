# ---------------------------------------------------------------
# bilinearTensor4Block
# Implements Algorithm 1 for a bilinear query x' M_k y
# under four-block missingness.
#
# Args:
#   Y      : list of K matrices or N x T x K array with NA in lower-right blocks
#   k      : target slice
#   r      : target rank
#   x      : vector in R^{N2_k}
#   y      : vector in R^{T2_k}
#   tau    : threshold for H_k^\dagger
#   N1,T1  : integer vectors of length K giving four-block observed sizes
#
# Returns:
#   scalar estimate of x' M_k y
# ---------------------------------------------------------------

bilinearTensor4Block = function(Y, k, r, x, y, tau, N1, T1) {
  if (is.array(Y)) {
    K  = dim(Y)[3]
    N  = dim(Y)[1]
    Tt = dim(Y)[2]
    Y_list = lapply(seq_len(K), function(j) Y[, , j])
  } else {
    Y_list = Y
    K  = length(Y_list)
    N  = nrow(Y_list[[1]])
    Tt = ncol(Y_list[[1]])
  }
  
  stopifnot(length(N1) == K, length(T1) == K)
  stopifnot(k >= 1, k <= K)
  
  N2k = N - N1[k]
  T2k = Tt - T1[k]
  
  if (length(x) != N2k) {
    stop("length(x) must equal N2_k = nrow(Y) - N1[k].")
  }
  if (length(y) != T2k) {
    stop("length(y) must equal T2_k = ncol(Y) - T1[k].")
  }
  
  r_eff = min(r, N, Tt, sum(T1), sum(N1))
  
  if (r_eff < 1) {
    stop("Effective rank is zero. Check r, N1, and T1.")
  }
  
  Y_left_pool = do.call(cbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][, seq_len(T1[j]), drop = FALSE]
  }))
  
  s_left = svd(Y_left_pool, nu = r_eff, nv = 0)
  U_left = s_left$u[, seq_len(r_eff), drop = FALSE]
  
  U1k_hat = U_left[seq_len(N1[k]), , drop = FALSE]
  U2k_hat = U_left[(N1[k] + 1):N, , drop = FALSE]
  
  H_k = crossprod(U1k_hat)
  eH = eigen(H_k, symmetric = TRUE)
  
  H_k_dagger = eH$vectors %*%
    diag(1 / pmax(eH$values, tau), r_eff) %*%
    t(eH$vectors)
  
  alpha_hat = crossprod(U2k_hat, x)
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  
  s_up = svd(Y_up_pool, nu = r_eff, nv = r_eff)
  U_up = s_up$u[, seq_len(r_eff), drop = FALSE]
  V_up = s_up$v[, seq_len(r_eff), drop = FALSE]
  D_up = s_up$d[seq_len(r_eff)]
  
  s_k = if (k == 1) 0 else sum(N1[seq_len(k - 1)])
  
  U_up_k = U_up[(s_k + 1):(s_k + N1[k]), , drop = FALSE]
  V2k_hat = V_up[(T1[k] + 1):Tt, , drop = FALSE]
  
  T_y = crossprod(V2k_hat, y)
  W_y = diag(D_up, r_eff) %*% T_y
  X_y = U_up_k %*% W_y
  
  beta_hat = H_k_dagger %*% crossprod(U1k_hat, X_y)
  
  as.numeric(crossprod(alpha_hat, beta_hat))
}