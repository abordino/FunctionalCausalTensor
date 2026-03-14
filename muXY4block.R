# ---------------------------------------------------------------
#  mu4block_pool : implements Algorithm 1 (bilinear 4-block panel)
# ---------------------------------------------------------------
# Args:
#   Y   : either
#         (i)  a list of length K, each element an N×T matrix with NA in bottom-right block, OR
#         (ii) an array of dim c(N, T, K)
#   k   : target arm in {1,...,K}
#   r   : target rank
#   x   : length N2^(k) vector (ideally unit); if NULL -> all-ones, normalized
#   y   : length T2^(k) vector (ideally unit); if NULL -> all-ones, normalized
#   N1, T1 : optional integer vectors of length K giving block sizes; if NULL -> detected from NA pattern
#   normalize_xy : if TRUE, normalize x and y to unit norm
#
# Returns:
#   list(mu_hat = scalar, N1 = N1_vec, T1 = T1_vec, r_eff = used_rank)
# ---------------------------------------------------------------

mu4block_pool = function(Y, k, r = 1, x = NULL, y = NULL,
                          N1 = NULL, T1 = NULL,
                          normalize_xy = TRUE) {
  
  ## ---- coerce input to list of matrices ------------------------
  if (is.array(Y)) {
    d = dim(Y)
    if (length(d) != 3) stop("If Y is an array, it must have dim (N, T, K).")
    N = d[1]; Tt = d[2]; K = d[3]
    Y_list = lapply(1:K, function(j) Y[, , j, drop = FALSE])
  } else if (is.list(Y)) {
    Y_list = Y
    K = length(Y_list)
    if (K < 1) stop("Y list is empty.")
    N = nrow(Y_list[[1]])
    Tt = ncol(Y_list[[1]])
  } else {
    stop("Y must be either a list of matrices or an (N,T,K) array.")
  }
  
  if (k < 1 || k > K) stop("k must be between 1 and K.")
  if (r < 1) stop("r must be >= 1.")
  
  ## ---- validate dims -------------------------------------------
  for (j in 1:K) {
    if (!is.matrix(Y_list[[j]])) Y_list[[j]] = as.matrix(Y_list[[j]])
    if (nrow(Y_list[[j]]) != N || ncol(Y_list[[j]]) != Tt) {
      stop("All slices must have the same dimensions N×T.")
    }
  }
  
  ## ---- detect N1^(j), T1^(j) if not provided -------------------
  detect_blocks = function(M) {
    row_has_NA = apply(M, 1, function(z) any(is.na(z)))
    N1j = if (any(row_has_NA)) which(row_has_NA)[1] - 1 else nrow(M)
    if (N1j < 1) stop("Could not find a fully-observed row block (N1).")
    
    col_has_NA = apply(M, 2, function(z) any(is.na(z)))
    T1j = if (any(col_has_NA)) which(col_has_NA)[1] - 1 else ncol(M)
    if (T1j < 1) stop("Could not find a fully-observed column block (T1).")
    
    list(N1 = N1j, T1 = T1j)
  }
  
  if (is.null(N1) || is.null(T1)) {
    N1_vec = integer(K)
    T1_vec = integer(K)
    for (j in 1:K) {
      bt = detect_blocks(Y_list[[j]])
      N1_vec[j] = bt$N1
      T1_vec[j] = bt$T1
    }
  } else {
    if (length(N1) != K || length(T1) != K) stop("N1 and T1 must be length K.")
    N1_vec = as.integer(N1)
    T1_vec = as.integer(T1)
  }
  
  ## ---- target arm block sizes ----------------------------------
  N1k = N1_vec[k]; T1k = T1_vec[k]
  N2k = N - N1k
  T2k = Tt - T1k
  if (N2k < 1 || T2k < 1) stop("Target arm k has no missing bottom-right block (N2 or T2 is 0).")
  
  ## ---- default x,y (all-ones, normalized) ----------------------
  if (is.null(x)) x = rep(1, N2k)
  if (is.null(y)) y = rep(1, T2k)
  
  if (length(x) != N2k) stop(sprintf("x must have length N2^(k)=%d.", N2k))
  if (length(y) != T2k) stop(sprintf("y must have length T2^(k)=%d.", T2k))
  
  if (normalize_xy) {
    nx = sqrt(sum(x^2)); ny = sqrt(sum(y^2))
    if (nx == 0 || ny == 0) stop("x and y must be nonzero.")
    x = x / nx
    y = y / ny
  }
  
  ## ---- choose an effective rank that is feasible ---------------
  sumT1 = sum(T1_vec)
  sumN1 = sum(N1_vec)
  r_eff = min(r, N, Tt, sumT1, sumN1)
  if (r_eff < 1) stop("Effective rank became < 1; check dimensions/block sizes.")
  
  ## ==============================================================
  ##  Step 1–4: pooled left SVD -> alpha_x^(k)
  ## ==============================================================
  Y_left_pool = do.call(cbind, lapply(1:K, function(j) {
    Yj = Y_list[[j]][, 1:T1_vec[j], drop = FALSE]
    if (anyNA(Yj)) stop(sprintf("Found NA in left block for arm %d; expected fully observed.", j))
    Yj
  }))
  
  svd_left = svd(Y_left_pool, nu = r_eff, nv = 0)
  U_left_hat = svd_left$u[, 1:r_eff, drop = FALSE]  # N × r
  
  U1_hat_k = U_left_hat[1:N1k, , drop = FALSE]              # N1k × r
  U2_hat_k = U_left_hat[(N1k + 1):N, , drop = FALSE]        # N2k × r
  
  alpha_hat = crossprod(U2_hat_k, x)  # r×1
  
  ## ---- Step 5: QR of U1_hat_k ----------------------------------
  qrU1 = qr(U1_hat_k)
  rankU1 = qrU1$rank
  
  Qk = qr.Q(qrU1)  # N1k × rankU1 (orthonormal columns)
  Rk = qr.R(qrU1)  # rankU1 × r_eff (upper trapezoidal in general)
  
  ## ==============================================================
  ##  Step 6–14: pooled upper SVD -> beta_y^(k)
  ## ==============================================================
  Y_up_pool = do.call(rbind, lapply(1:K, function(j) {
    Yj = Y_list[[j]][1:N1_vec[j], , drop = FALSE]
    if (anyNA(Yj)) stop(sprintf("Found NA in upper block for arm %d; expected fully observed.", j))
    Yj
  }))
  
  svd_up = svd(Y_up_pool, nu = r_eff, nv = r_eff)
  U_up_hat = svd_up$u[, 1:r_eff, drop = FALSE]  # (sum N1) × r
  V_up_hat = svd_up$v[, 1:r_eff, drop = FALSE]  # T × r
  d_up     = svd_up$d[1:r_eff]                  # length r
  
  s_k = if (k == 1) 0 else sum(N1_vec[1:(k - 1)])
  idx = (s_k + 1):(s_k + N1k)
  U_up_hat_k = U_up_hat[idx, , drop = FALSE]    # N1k × r
  
  V2_hat_k = V_up_hat[(T1k + 1):Tt, , drop = FALSE]  # T2k × r
  
  t_vec = crossprod(V2_hat_k, y)          # r×1
  u_vec = d_up * as.numeric(t_vec)        # r×1 (Sigma_up %*% t)
  v_vec = U_up_hat_k %*% u_vec            # N1k×1
  
  ## ---- Step 15: beta via QR
  if (rankU1 == r_eff && ncol(Rk) == r_eff && nrow(Rk) == r_eff) {
    beta_hat = solve(Rk, crossprod(Qk, v_vec))  # r×1
  } else {
    # fallback: least-squares solution to U1_hat_k * beta ≈ v_vec
    beta_hat = qr.solve(U1_hat_k, v_vec)
  }
  
  ## ---- Step 16: mu_hat = <alpha_hat, beta_hat> -----------------
  mu_hat = as.numeric(crossprod(alpha_hat, beta_hat))
  return(mu_hat)
}
