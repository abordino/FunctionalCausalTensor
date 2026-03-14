# ---------------------------------------------------------------
#  oracle_mu4block_pool : oracle pooled version to check leading Gaussian terms
#  (fixed scalar algebra; matches the step-by-step debug)
# ---------------------------------------------------------------
# Args:
#   Y   : list length K of N×T matrices (NA only in bottom-right), OR array dim (N,T,K)
#   k   : target arm in {1,...,K}
#   U   : N×r oracle left subspace (orthonormal cols assumed)
#   V   : T×r oracle right subspace (orthonormal cols assumed)
#   R_list : list length K of r×r matrices R^(k)
#   x   : length N2^(k) vector; if NULL -> ones (normalized if normalize_xy)
#   y   : length T2^(k) vector; if NULL -> ones (normalized if normalize_xy)
#   N1, T1 : optional integer vectors length K; if NULL -> detected from NA pattern
#   normalize_xy : if TRUE normalize x,y to unit norm
#
# Returns:
#   scalar mu_hat
# ---------------------------------------------------------------

oracle_muXY4block = function(Y, k, U, V, R_list,
                                x = NULL, y = NULL,
                                N1 = NULL, T1 = NULL,
                                normalize_xy = TRUE) {
  
  # ---- coerce Y to list (if array) ----
  if (is.array(Y)) {
    K = dim(Y)[3]
    Y = lapply(1:K, function(j) Y[, , j])
  } else {
    K = length(Y)
  }
  
  N  = nrow(Y[[1]])
  Tt = ncol(Y[[1]])
  
  # ---- detect (N1,T1) if not supplied ----
  detect_blocks = function(M) {
    N1j = which(apply(M, 1, function(z) any(is.na(z))))[1] - 1
    T1j = which(apply(M, 2, function(z) any(is.na(z))))[1] - 1
    list(N1 = N1j, T1 = T1j)
  }
  
  if (is.null(N1) || is.null(T1)) {
    bt = lapply(Y, detect_blocks)
    N1 = sapply(bt, `[[`, "N1")
    T1 = sapply(bt, `[[`, "T1")
  }
  
  N1k = N1[k]; T1k = T1[k]
  N2k = N - N1k
  T2k = Tt - T1k
  r   = nrow(R_list[[1]])
  
  # ---- defaults x,y ----
  if (is.null(x)) x = rep(1, N2k)
  if (is.null(y)) y = rep(1, T2k)
  if (normalize_xy) {
    x = x / sqrt(sum(x^2))
    y = y / sqrt(sum(y^2))
  }
  
  # ==============================================================
  # Build pooled W_left, W_up (oracle)
  # ==============================================================
  # W_left^T = [R^(1) (V1^(1))^T | ... | R^(K) (V1^(K))^T]
  W_left = do.call(rbind, lapply(1:K, function(j) {
    V1j = V[1:T1[j], 1:r, drop = FALSE]         # T1j×r
    t(R_list[[j]] %*% t(V1j))                   # T1j×r
  }))
  G_left_inv = solve(crossprod(W_left))         # r×r
  
  # W_up = [U1^(1) R^(1); ...; U1^(K) R^(K)]
  W_up = do.call(rbind, lapply(1:K, function(j) {
    U1j = U[1:N1[j], 1:r, drop = FALSE]         # N1j×r
    U1j %*% R_list[[j]]                         # N1j×r
  }))
  G_up_inv = solve(crossprod(W_up))             # r×r
  
  # ==============================================================
  # Pooled residuals E_left^pool, E_up^pool
  # ==============================================================
  Y_left_pool = do.call(cbind, lapply(1:K, function(j) Y[[j]][, 1:T1[j], drop = FALSE]))
  E_left_pool = Y_left_pool - (U[, 1:r, drop = FALSE] %*% t(W_left))
  
  Y_up_pool = do.call(rbind, lapply(1:K, function(j) Y[[j]][1:N1[j], , drop = FALSE]))
  E_up_pool = Y_up_pool - (W_up %*% t(V[, 1:r, drop = FALSE]))
  
  # ==============================================================
  # Target oracle pieces and mu = x^T \hat M_d^{(k)} y
  # ==============================================================
  Rk  = R_list[[k]]
  U2k = U[(N1k + 1):N, 1:r, drop = FALSE]       # N2k×r
  V2k = V[(T1k + 1):Tt, 1:r, drop = FALSE]      # T2k×r
  by  = crossprod(V2k, y)                        # r×1
  
  # ---- base oracle signal (scalar): (U2^T x)^T [R (V2^T y)]
  ax   = crossprod(U2k, x)                       # r×1
  mu_d = as.numeric(crossprod(ax, Rk %*% by))    # scalar
  
  # ---- left correction (scalar): (E_left_2^T x)^T [W_left G_left^{-1} R (V2^T y)]
  E_left_2 = E_left_pool[(N1k + 1):N, , drop = FALSE]  # N2k×(sum T1)
  s1 = crossprod(E_left_2, x)                            # (sum T1)×1
  u1 = W_left %*% (G_left_inv %*% (Rk %*% by))           # (sum T1)×1
  mu_1 = as.numeric(crossprod(s1, u1))                   # scalar

  # ---- up correction (scalar): x^T U2 R G_up^{-1} W_up^T (E_up_sub y)
  E_up_sub = E_up_pool[, (T1k + 1):Tt, drop = FALSE]     # (sum N1)×T2k
  s2 = E_up_sub %*% y                                     # (sum N1)×1
  t2 = G_up_inv %*% crossprod(W_up, s2)                  # r×1
  mu_2 = as.numeric(crossprod(x, U2k %*% (Rk %*% t2)))   # scalar
  
  mu_d + mu_2 + mu_1
}