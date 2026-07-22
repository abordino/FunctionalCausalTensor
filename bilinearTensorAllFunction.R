# ---------------------------------------------------------------
# 1. mu_xy4block_pool : Algorithm 1
# 
# Args:
#   Y  : list of K N x T matrices, or an N x T x K array
#   k  : target arm
#   r  : target rank
#   x  : vector of length N_2k
#   y  : vector of length T_2k
#   tau: threshold used in H_k^\dagger
#   N1, T1: optional block-size vectors; detected from NA patterns if omitted
#
# Returns:
#   scalar mu_hat = x^T \cal{M}^{(d)}_{\bullet, \bullet, k} y
# ---------------------------------------------------------------

mu_xy4block_pool = function(Y, k, r, x, y, tau, N1 = NULL, T1 = NULL) {
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
    row_has_NA = apply(M, 1, anyNA)
    col_has_NA = apply(M, 2, anyNA)
    
    list(
      N1 = if (any(row_has_NA)) which(row_has_NA)[1] - 1L else nrow(M),
      T1 = if (any(col_has_NA)) which(col_has_NA)[1] - 1L else ncol(M)
    )
  }
  
  if (is.null(N1) || is.null(T1)) {
    blocks = lapply(Y_list, detect_blocks)
    N1 = vapply(blocks, `[[`, integer(1), "N1")
    T1 = vapply(blocks, `[[`, integer(1), "T1")
  }
  
  N1k = N1[k]
  T1k = T1[k]
  r_eff = min(r, N, Tt, sum(N1), sum(T1))
  
  Y_left_pool = do.call(cbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][, seq_len(T1[j]), drop = FALSE]
  }))
  
  U_left = svd(Y_left_pool, nu = r_eff, nv = 0)$u[, seq_len(r_eff), drop = FALSE]
  
  U1k_hat = U_left[seq_len(N1k), , drop = FALSE]
  U2k_hat = U_left[seq.int(N1k + 1L, N), , drop = FALSE]
  
  H_k = crossprod(U1k_hat)
  H_eig = eigen(H_k, symmetric = TRUE)
  H_k_dagger = H_eig$vectors %*%
    (t(H_eig$vectors) / pmax(H_eig$values, tau))
  
  alpha_hat = drop(crossprod(U2k_hat, x))
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  
  s_up = svd(Y_up_pool, nu = r_eff, nv = r_eff)
  U_up = s_up$u[, seq_len(r_eff), drop = FALSE]
  V_up = s_up$v[, seq_len(r_eff), drop = FALSE]
  D_up = s_up$d[seq_len(r_eff)]
  
  row_start = sum(N1[seq_len(k - 1L)]) + 1L
  U_up_k = U_up[seq.int(row_start, row_start + N1k - 1L), , drop = FALSE]
  V2k_hat = V_up[seq.int(T1k + 1L, Tt), , drop = FALSE]
  
  gamma_hat = drop(crossprod(V2k_hat, y))
  X_y = drop(U_up_k %*% (D_up * gamma_hat))
  
  beta_hat = drop(
    H_k_dagger %*% crossprod(U1k_hat, X_y)
  )
  
  sum(alpha_hat * beta_hat)
}

# ---------------------------------------------------------------
# 2. oracle_mu_xy4block_pool : oracle pooled version
# 
# Args:
#   Y      : list of K N x T matrices, or an N x T x K array
#   k      : target arm
#   U      : N x r oracle left singular vectors
#   V      : T x r oracle right singular vectors
#   R_list : list of K r x r matrices
#   x      : vector of length N_2k
#   y      : vector of length T_2k
#   N1, T1 : optional block-size vectors; detected from NA patterns if omitted
#
# Returns:
#   scalar oracle estimate for x^T \cal{M}^{(d)}_{\bullet, \bullet, k} y,
#   including the two oracle correction terms
# ---------------------------------------------------------------

oracle_mu_xy4block_pool = function(Y, k, U, V, R_list, x, y,
                                   N1 = NULL, T1 = NULL) {
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
    row_has_NA = apply(M, 1, anyNA)
    col_has_NA = apply(M, 2, anyNA)
    
    list(
      N1 = if (any(row_has_NA)) which(row_has_NA)[1] - 1L else nrow(M),
      T1 = if (any(col_has_NA)) which(col_has_NA)[1] - 1L else ncol(M)
    )
  }
  
  if (is.null(N1) || is.null(T1)) {
    blocks = lapply(Y, detect_blocks)
    N1 = vapply(blocks, `[[`, integer(1), "N1")
    T1 = vapply(blocks, `[[`, integer(1), "T1")
  }
  
  N1k = N1[k]
  T1k = T1[k]
  
  W_left = do.call(rbind, lapply(seq_len(K), function(j) {
    V1j = V[seq_len(T1[j]), seq_len(r), drop = FALSE]
    V1j %*% t(R_list[[j]])
  }))
  
  W_up = do.call(rbind, lapply(seq_len(K), function(j) {
    U1j = U[seq_len(N1[j]), seq_len(r), drop = FALSE]
    U1j %*% R_list[[j]]
  }))
  
  G_left_inv = solve(crossprod(W_left))
  G_up_inv = solve(crossprod(W_up))
  
  Y_left_pool = do.call(cbind, lapply(seq_len(K), function(j) {
    Y[[j]][, seq_len(T1[j]), drop = FALSE]
  }))
  
  E_left_pool =
    Y_left_pool -
    U[, seq_len(r), drop = FALSE] %*% t(W_left)
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  
  E_up_pool =
    Y_up_pool -
    W_up %*% t(V[, seq_len(r), drop = FALSE])
  
  Rk = R_list[[k]]
  U2k = U[seq.int(N1k + 1L, N), seq_len(r), drop = FALSE]
  V2k = V[seq.int(T1k + 1L, Tt), seq_len(r), drop = FALSE]
  
  alpha = drop(crossprod(U2k, x))
  gamma = drop(crossprod(V2k, y))
  
  Rk_gamma = drop(Rk %*% gamma)
  mu_d = sum(alpha * Rk_gamma)
  
  E_left_2 = E_left_pool[seq.int(N1k + 1L, N), , drop = FALSE]
  a1 = drop(crossprod(E_left_2, x))
  b1 = drop(W_left %*% (G_left_inv %*% Rk_gamma))
  mu_1 = sum(a1 * b1)
  
  E_up_2 = E_up_pool[, seq.int(T1k + 1L, Tt), drop = FALSE]
  b2 = drop(
    G_up_inv %*%
      crossprod(W_up, drop(E_up_2 %*% y))
  )
  mu_2 = sum(alpha * drop(Rk %*% b2))
  
  mu_d + mu_1 + mu_2
}

# ---------------------------------------------------------------
# 3. oracle_mu_xy4block_pool_local : oracle local version 
# 
# Args:
#   Y      : list of K N x T matrices, or an N x T x K array
#   k      : target arm
#   U      : N x r oracle left singular vectors
#   V      : T x r oracle right singular vectors
#   R_list : list of K r x r matrices
#   x      : vector of length N_2k
#   y      : vector of length T_2k
#   N1, T1 : optional block-size vectors; detected from NA patterns if omitted
#
# Returns:
#   scalar oracle estimate for x^T \cal{M}^{(d)}_{\bullet, \bullet, k} y,
#   including the local oracle correction term
# ---------------------------------------------------------------

oracle_mu_xy4block_pool_local = function(Y, k, U, V, R_list, x, y,
                                         N1 = NULL, T1 = NULL) {
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
    row_has_NA = apply(M, 1, anyNA)
    col_has_NA = apply(M, 2, anyNA)
    
    list(
      N1 = if (any(row_has_NA)) which(row_has_NA)[1] - 1L else nrow(M),
      T1 = if (any(col_has_NA)) which(col_has_NA)[1] - 1L else ncol(M)
    )
  }
  
  if (is.null(N1) || is.null(T1)) {
    blocks = lapply(Y, detect_blocks)
    N1 = vapply(blocks, `[[`, integer(1), "N1")
    T1 = vapply(blocks, `[[`, integer(1), "T1")
  }
  
  N1k = N1[k]
  T1k = T1[k]
  
  W_up = do.call(rbind, lapply(seq_len(K), function(j) {
    U1j = U[seq_len(N1[j]), seq_len(r), drop = FALSE]
    U1j %*% R_list[[j]]
  }))
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y[[j]][seq_len(N1[j]), , drop = FALSE]
  }))
  
  V_r = V[, seq_len(r), drop = FALSE]
  E_up_pool = Y_up_pool - W_up %*% t(V_r)
  
  Rk = R_list[[k]]
  U1k = U[seq_len(N1k), seq_len(r), drop = FALSE]
  U2k = U[seq.int(N1k + 1L, N), seq_len(r), drop = FALSE]
  V2k = V[seq.int(T1k + 1L, Tt), seq_len(r), drop = FALSE]
  
  alpha = drop(crossprod(U2k, x))
  gamma = drop(crossprod(V2k, y))
  
  mu_d = sum(alpha * drop(Rk %*% gamma))
  
  row_start = sum(N1[seq_len(k - 1L)]) + 1L
  row_end = row_start + N1k - 1L
  
  E_up_k = E_up_pool[row_start:row_end, , drop = FALSE]
  
  beta_local = drop(
    solve(
      crossprod(U1k),
      crossprod(U1k, drop(E_up_k %*% (V_r %*% gamma)))
    )
  )
  
  mu_3 = sum(alpha * beta_local)
  
  mu_d + mu_3
}

# ---------------------------------------------------------------
# 4. bilinearTensorStaggered
#
# Args:
#   Y       : list of K N x T matrices, or N x T x K array
#   k       : target slice
#   a,b     : target staggered-adoption block indices in slice k
#   r       : target rank
#   x       : unit vector in R^{N_{ak}}
#   y       : unit vector in R^{T_{bk}}
#   A       : N x K adoption-time matrix, A[i,j] in {1,...,T+1, Inf}
#             Omega[i,t,j] = 1{t < A[i,j]}; Inf means never adopter
#   Omega   : optional N x T x K logical observation mask
#   tau     : threshold for H_k^\dagger
#
# Returns:
#   scalar estimate of mu_{xy}^{(k,a,b)}
# ---------------------------------------------------------------

bilinearTensorStaggered = function(Y, k, a, b, r, x, y, A = NULL,
                                   Omega = NULL, tau) {
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
  
  stopifnot(k >= 1, k <= K)
  stopifnot(r >= 1)
  stopifnot(tau > 0)
  
  if (is.null(A) && is.null(Omega)) {
    Omega = array(FALSE, dim = c(N, Tt, K))
    for (j in seq_len(K)) {
      Omega[, , j] = !is.na(Y_list[[j]])
    }
    
    A = matrix(Inf, nrow = N, ncol = K)
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        obs = Omega[i, , j]
        
        if (any(diff(as.integer(obs)) > 0)) {
          stop("Rows must have staggered-adoption missingness.")
        }
        
        first_unobs = which(!obs)[1]
        
        if (!is.na(first_unobs)) {
          A[i, j] = first_unobs
        }
      }
    }
  }
  
  if (!is.null(A)) {
    A = as.matrix(A)
    
    if (nrow(A) != N || ncol(A) != K) {
      stop("A must be an N x K adoption-time matrix.")
    }
  }
  
  if (is.null(Omega)) {
    Omega = array(FALSE, dim = c(N, Tt, K))
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        if (is.infinite(A[i, j])) {
          Omega[i, , j] = TRUE
        } else if (A[i, j] > 1) {
          Omega[i, seq_len(A[i, j] - 1), j] = TRUE
        }
      }
    }
  } else {
    if (!is.array(Omega) || !identical(dim(Omega), c(N, Tt, K))) {
      stop("Omega must be an N x T x K logical array.")
    }
    
    Omega = array(as.logical(Omega), dim = c(N, Tt, K))
    
    if (is.null(A)) {
      A = matrix(Inf, nrow = N, ncol = K)
      
      for (j in seq_len(K)) {
        for (i in seq_len(N)) {
          obs = Omega[i, , j]
          
          if (any(diff(as.integer(obs)) > 0)) {
            stop("Rows must have staggered-adoption missingness.")
          }
          
          first_unobs = which(!obs)[1]
          
          if (!is.na(first_unobs)) {
            A[i, j] = first_unobs
          }
        }
      }
    }
  }
  
  # -------------------------------------------------------------
  # Reorder rows by target-layer adoption times.
  # -------------------------------------------------------------
  
  obs_len_k = ifelse(is.infinite(A[, k]), Tt, pmin(Tt, A[, k] - 1))
  row_perm = order(-obs_len_k, seq_len(N))
  
  Y_list = lapply(seq_len(K), function(j) {
    Y_list[[j]][row_perm, , drop = FALSE]
  })
  
  A = A[row_perm, , drop = FALSE]
  Omega = Omega[row_perm, , , drop = FALSE]
  
  obs_len_k = obs_len_k[row_perm]
  
  if (any(diff(obs_len_k) > 0)) {
    stop("Target-layer ordering failed.")
  }
  
  # -------------------------------------------------------------
  # Construct target-layer staircase blocks from adoption times.
  # -------------------------------------------------------------
  
  m_desc = unique(obs_len_k)
  ok = length(m_desc)
  
  if (m_desc[1] != Tt) {
    stop("Target layer must contain at least one fully observed row block.")
  }
  
  if (tail(m_desc, 1) <= 0) {
    stop("Target layer must contain a non-empty initial observed time period.")
  }
  
  row_parts = lapply(m_desc, function(m) which(obs_len_k == m))
  
  m_asc = rev(m_desc)
  T_sizes = diff(c(0, m_asc))
  
  ends = cumsum(T_sizes)
  starts = c(1, head(ends, -1) + 1)
  col_parts = Map(seq, starts, ends)
  
  if (a < 1 || a > ok || b < 1 || b > ok) {
    stop("a and b must be valid block indices for slice k.")
  }
  
  if (a + b <= ok + 1) {
    stop("Target block is observed under staggered missingness. Need a + b > o_k + 1.")
  }
  
  if (length(x) != length(row_parts[[a]])) {
    stop("length(x) must equal N_{ak}.")
  }
  
  if (length(y) != length(col_parts[[b]])) {
    stop("length(y) must equal T_{bk}.")
  }
  
  # -------------------------------------------------------------
  # S^+, S^-, Q^+, Q^- for target block (a,b)
  # -------------------------------------------------------------
  
  safe_seq = function(u, v) {
    if (v < u) integer(0) else seq(u, v)
  }
  
  take_parts = function(parts, ids) {
    ids = ids[ids >= 1 & ids <= length(parts)]
    
    if (length(ids) == 0) {
      integer(0)
    } else {
      unlist(parts[ids], use.names = FALSE)
    }
  }
  
  S_plus_end  = ok + 1 - b
  S_minus_beg = ok + 2 - b
  
  Q_plus_end  = ok + 1 - a
  Q_minus_beg = ok + 2 - a
  
  S_plus = take_parts(row_parts, safe_seq(1, S_plus_end))
  S_minus = take_parts(row_parts, safe_seq(S_minus_beg, a))
  
  Q_plus = take_parts(col_parts, safe_seq(1, Q_plus_end))
  Q_minus = take_parts(col_parts, safe_seq(Q_minus_beg, b))
  
  S_idx = c(S_plus, S_minus)
  Q_idx = c(Q_plus, Q_minus)
  
  n_aux = length(S_idx)
  t_aux = length(Q_idx)
  
  S_plus_local = seq_len(length(S_plus))
  S_minus_local = (length(S_plus) + 1):n_aux
  
  Q_minus_local = (length(Q_plus) + 1):t_aux
  
  # -------------------------------------------------------------
  # Anchor rows and columns, with target layer using S^+ and Q^+.
  # -------------------------------------------------------------
  
  RowAnc = vector("list", K)
  ColAnc = vector("list", K)
  
  for (j in seq_len(K)) {
    if (j == k) {
      RowAnc[[j]] = S_plus
      ColAnc[[j]] = Q_plus
    } else {
      Omega_j = Omega[S_idx, Q_idx, j, drop = FALSE][, , 1]
      
      ColAnc[[j]] = Q_idx[which(colSums(Omega_j) == n_aux)]
      RowAnc[[j]] = S_idx[which(rowSums(Omega_j) == t_aux)]
    }
  }
  
  N1_aux = vapply(RowAnc, length, integer(1))
  T1_aux = vapply(ColAnc, length, integer(1))
  
  r_eff = min(r, n_aux, t_aux, sum(N1_aux), sum(T1_aux))
  
  if (r_eff < 1) {
    stop("Effective rank is zero")
  }
  
  # -------------------------------------------------------------
  # Embed x and y into S^- and Q^-.
  # -------------------------------------------------------------
  
  lx = length(S_minus) - length(row_parts[[a]])
  ly = length(Q_minus) - length(col_parts[[b]])
  
  if (lx < 0 || ly < 0) {
    stop("Invalid embedding dimensions for x or y.")
  }
  
  x_bar = c(rep(0, lx), x)
  y_bar = c(rep(0, ly), y)
  
  # -------------------------------------------------------------
  # Left pooled matrix
  # Y_left_pool = (Y_{S,ColAnc(1),1}, ..., Y_{S,ColAnc(K),K})
  # -------------------------------------------------------------
  
  Y_left_pool = do.call(cbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][S_idx, ColAnc[[j]], drop = FALSE]
  }))
  
  if (anyNA(Y_left_pool)) {
    stop("Y_left_pool contains NA. Anchor columns must be fully observed.")
  }
  
  s_left = svd(Y_left_pool, nu = r_eff, nv = 0)
  U_left = s_left$u[, seq_len(r_eff), drop = FALSE]
  
  U_plus_hat = U_left[S_plus_local, , drop = FALSE]
  U_minus_hat = U_left[S_minus_local, , drop = FALSE]
  
  H_k = crossprod(U_plus_hat)
  eH = eigen(H_k, symmetric = TRUE)
  
  H_k_dagger = eH$vectors %*%
    diag(1 / pmax(eH$values, tau), r_eff) %*%
    t(eH$vectors)
  
  alpha_hat = crossprod(U_minus_hat, x_bar)
  
  # -------------------------------------------------------------
  # Upper pooled matrix
  # Y_up_pool = (Y_{RowAnc(1),Q,1}; ... ; Y_{RowAnc(K),Q,K})
  # -------------------------------------------------------------
  
  Y_up_pool = do.call(rbind, lapply(seq_len(K), function(j) {
    Y_list[[j]][RowAnc[[j]], Q_idx, drop = FALSE]
  }))
  
  if (anyNA(Y_up_pool)) {
    stop("Y_up_pool contains NA. Anchor rows must be fully observed.")
  }
  
  s_up = svd(Y_up_pool, nu = r_eff, nv = r_eff)
  
  U_up = s_up$u[, seq_len(r_eff), drop = FALSE]
  V_up = s_up$v[, seq_len(r_eff), drop = FALSE]
  D_up = s_up$d[seq_len(r_eff)]
  
  s_k = if (k == 1) 0 else sum(N1_aux[seq_len(k - 1)])
  U_up_k = U_up[(s_k + 1):(s_k + length(S_plus)), , drop = FALSE]
  
  V_minus_hat = V_up[Q_minus_local, , drop = FALSE]
  
  T_y = crossprod(V_minus_hat, y_bar)
  W_y = diag(D_up, r_eff) %*% T_y
  X_y = U_up_k %*% W_y
  
  beta_hat = H_k_dagger %*% crossprod(U_plus_hat, X_y)
  
  as.numeric(crossprod(alpha_hat, beta_hat))
}

# ---------------------------------------------------------------
# 5. bilinearMatrixStaggered
#
# Args:
#   Y_mat : N x T matrix with staggered missingness
#   a,b   : target missing block indices after sorting rows by adoption time
#   r     : target rank
#   x     : unit vector in R^{N_a}
#   y     : unit vector in R^{T_b}
#   A     : optional length-N adoption-time vector,
#           Omega[i,t] = 1{t < A[i]}; Inf means never adopter
#   Omega : optional N x T logical observation mask
#   tau   : threshold for H^\dagger
#
# Returns:
#   scalar estimate of mu_{xy}^{(k,a,b)}
# ---------------------------------------------------------------

bilinearMatrixStaggered = function(Y_mat, a, b, r, x, y,
                                   A = NULL, Omega = NULL, tau) {
  
  N = nrow(Y_mat)
  Tt = ncol(Y_mat)
  
  if (is.null(A) && is.null(Omega)) {
    Omega = !is.na(Y_mat)
    
    A = rep(Inf, N)
    
    for (i in seq_len(N)) {
      obs = Omega[i, ]
      
      if (any(diff(as.integer(obs)) > 0)) {
        stop("Rows must have staggered-adoption missingness.")
      }
      
      first_unobs = which(!obs)[1]
      
      if (!is.na(first_unobs)) {
        A[i] = first_unobs
      }
    }
  }
  
  if (!is.null(A)) {
    if (length(A) != N) {
      stop("A must have length nrow(Y_mat).")
    }
  }
  
  if (is.null(Omega)) {
    Omega = matrix(FALSE, nrow = N, ncol = Tt)
    
    for (i in seq_len(N)) {
      if (is.infinite(A[i])) {
        Omega[i, ] = TRUE
      } else if (A[i] > 1) {
        Omega[i, seq_len(A[i] - 1)] = TRUE
      }
    }
  } else {
    if (!is.matrix(Omega) || !identical(dim(Omega), c(N, Tt))) {
      stop("Omega must be an N x T logical matrix.")
    }
    
    Omega = matrix(as.logical(Omega), nrow = N, ncol = Tt)
    
    if (is.null(A)) {
      A = rep(Inf, N)
      
      for (i in seq_len(N)) {
        obs = Omega[i, ]
        
        if (any(diff(as.integer(obs)) > 0)) {
          stop("Rows must have staggered-adoption missingness.")
        }
        
        first_unobs = which(!obs)[1]
        
        if (!is.na(first_unobs)) {
          A[i] = first_unobs
        }
      }
    }
  }
  
  # -------------------------------------------------------------
  # Reorder rows by adoption times in this matrix.
  # -------------------------------------------------------------
  
  obs_len = ifelse(is.infinite(A), Tt, pmin(Tt, A - 1))
  row_perm = order(-obs_len, seq_len(N))
  
  Y_mat = Y_mat[row_perm, , drop = FALSE]
  Omega = Omega[row_perm, , drop = FALSE]
  A = A[row_perm]
  
  obs_len = obs_len[row_perm]
  
  if (any(diff(obs_len) > 0)) {
    stop("Row ordering failed.")
  }
  
  # -------------------------------------------------------------
  # Construct staircase blocks from adoption times.
  # -------------------------------------------------------------
  
  m_desc = unique(obs_len)
  o = length(m_desc)
  
  if (m_desc[1] != Tt) {
    stop("Matrix must contain at least one fully observed row block.")
  }
  
  if (tail(m_desc, 1) <= 0) {
    stop("Matrix must contain a non-empty initial observed time period.")
  }
  
  row_parts = lapply(m_desc, function(m) which(obs_len == m))
  
  m_asc = rev(m_desc)
  T_part = diff(c(0, m_asc))
  
  ends = cumsum(T_part)
  starts = c(1, head(ends, -1) + 1)
  col_parts = Map(seq, starts, ends)
  
  if (a < 1 || a > o || b < 1 || b > o) {
    stop("a and b must be valid staggered block indices.")
  }
  
  if (a + b <= o + 1) {
    stop("Target block is observed under staggered missingness. Need a + b > o + 1.")
  }
  
  if (length(x) != length(row_parts[[a]])) {
    stop("length(x) must equal N_a.")
  }
  
  if (length(y) != length(col_parts[[b]])) {
    stop("length(y) must equal T_b.")
  }
  
  # -------------------------------------------------------------
  # S^+, S^-, Q^+, Q^- for target block (a,b)
  # -------------------------------------------------------------
  
  safe_seq = function(u, v) {
    if (v < u) integer(0) else seq(u, v)
  }
  
  take_parts = function(parts, ids) {
    ids = ids[ids >= 1 & ids <= length(parts)]
    
    if (length(ids) == 0) {
      integer(0)
    } else {
      unlist(parts[ids], use.names = FALSE)
    }
  }
  
  S_plus_end  = o + 1 - b
  S_minus_beg = o + 2 - b
  
  Q_plus_end  = o + 1 - a
  Q_minus_beg = o + 2 - a
  
  S_plus = take_parts(row_parts, safe_seq(1, S_plus_end))
  S_minus = take_parts(row_parts, safe_seq(S_minus_beg, a))
  
  Q_plus = take_parts(col_parts, safe_seq(1, Q_plus_end))
  Q_minus = take_parts(col_parts, safe_seq(Q_minus_beg, b))
  
  S_idx = c(S_plus, S_minus)
  Q_idx = c(Q_plus, Q_minus)
  
  n_aux = length(S_idx)
  t_aux = length(Q_idx)
  
  S_plus_local = seq_len(length(S_plus))
  S_minus_local = (length(S_plus) + 1):n_aux
  
  Q_minus_local = (length(Q_plus) + 1):t_aux
  
  # -------------------------------------------------------------
  # Embed x and y into S^- and Q^-.
  # -------------------------------------------------------------
  
  lx = length(S_minus) - length(row_parts[[a]])
  ly = length(Q_minus) - length(col_parts[[b]])
  
  if (lx < 0 || ly < 0) {
    stop("Invalid embedding dimensions for x or y.")
  }
  
  x_bar = c(rep(0, lx), x)
  y_bar = c(rep(0, ly), y)
  
  r_eff = min(r, n_aux, t_aux, length(S_plus), length(Q_plus))
  
  if (r_eff < 1) {
    stop("Effective rank is zero")
  }
  
  # -------------------------------------------------------------
  # Left block only from target matrix:
  #
  # Y_left = Y_{S,Q^+}
  # -------------------------------------------------------------
  
  Y_left = Y_mat[S_idx, Q_plus, drop = FALSE]
  
  if (anyNA(Y_left)) {
    stop("Y_left contains NA. Q^+ must be observed for all rows in S.")
  }
  
  s_left = svd(Y_left, nu = r_eff, nv = 0)
  U_left = s_left$u[, seq_len(r_eff), drop = FALSE]
  
  U_plus_hat = U_left[S_plus_local, , drop = FALSE]
  U_minus_hat = U_left[S_minus_local, , drop = FALSE]
  
  H = crossprod(U_plus_hat)
  eH = eigen(H, symmetric = TRUE)
  
  H_dagger = eH$vectors %*%
    diag(1 / pmax(eH$values, tau), r_eff) %*%
    t(eH$vectors)
  
  alpha_hat = crossprod(U_minus_hat, x_bar)
  
  # -------------------------------------------------------------
  # Upper block only from target matrix:
  #
  # Y_up = Y_{S^+,Q}
  # -------------------------------------------------------------
  
  Y_up = Y_mat[S_plus, Q_idx, drop = FALSE]
  
  if (anyNA(Y_up)) {
    stop("Y_up contains NA. S^+ must be observed for all columns in Q.")
  }
  
  s_up = svd(Y_up, nu = r_eff, nv = r_eff)
  
  U_up = s_up$u[, seq_len(r_eff), drop = FALSE]
  V_up = s_up$v[, seq_len(r_eff), drop = FALSE]
  D_up = s_up$d[seq_len(r_eff)]
  
  V_minus_hat = V_up[Q_minus_local, , drop = FALSE]
  
  T_y = crossprod(V_minus_hat, y_bar)
  W_y = diag(D_up, r_eff) %*% T_y
  X_y = U_up %*% W_y
  
  beta_hat = H_dagger %*% crossprod(U_plus_hat, X_y)
  
  as.numeric(crossprod(alpha_hat, beta_hat))
}


# ---------------------------------------------------------------
# 6. bilinearTensorStaggeredPsi
#
# Estimates Psi_0^{(h)}(k) over all policy-on / missing target
# blocks in target slice k, for h in {"ATE", "RowHet", "Local", "Trend"}
#
# Args:
#   Y          : list of K N x T matrices, or N x T x K array
#   k          : target slice
#   r          : target rank
#   tau        : threshold for H_k^\dagger
#   functional : one of "ATE", "RowHet", "Local", "Trend"
#   eta        : length-N vector in {+1,-1}; required for RowHet.
#                eta is indexed in the original, unsorted row order.
#   row_index  : original global row index i0; required for Local
#   A          : optional N x K adoption-time matrix,
#                Omega[i,t,j] = 1{t < A[i,j]};
#                Inf means never adopter
#   Omega      : optional N x T x K logical observation mask
#
# Returns:
#   scalar estimate of Psi_0^{(functional)}(k)
# ---------------------------------------------------------------

bilinearTensorStaggeredPsi = function(Y, k, r, tau,
                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                      eta = NULL,
                                      row_index = NULL,
                                      A = NULL,
                                      Omega = NULL) {
  functional = match.arg(functional)
  
  if (is.array(Y)) {
    K  = dim(Y)[3]
    N  = dim(Y)[1]
    Tt = dim(Y)[2]
    Y_list = lapply(seq_len(K), function(j) Y[, , j])
  } else if (is.list(Y)) {
    Y_list = Y
    K  = length(Y_list)
    N  = nrow(Y_list[[1]])
    Tt = ncol(Y_list[[1]])
    
    for (j in seq_len(K)) {
      if (!is.matrix(Y_list[[j]]) ||
          nrow(Y_list[[j]]) != N ||
          ncol(Y_list[[j]]) != Tt) {
        stop("All elements of Y must be N x T matrices.")
      }
    }
  } else {
    stop("Y must be either an N x T x K array or a list of K N x T matrices.")
  }
  
  stopifnot(k >= 1, k <= K)
  stopifnot(r >= 1)
  stopifnot(tau > 0)
  
  # -------------------------------------------------------------
  # Infer A and Omega
  # -------------------------------------------------------------
  
  if (is.null(A) && is.null(Omega)) {
    Omega = array(FALSE, dim = c(N, Tt, K))
    
    for (j in seq_len(K)) {
      Omega[, , j] = !is.na(Y_list[[j]])
    }
    
    A = matrix(Inf, nrow = N, ncol = K)
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        obs = Omega[i, , j]
        
        if (any(diff(as.integer(obs)) > 0)) {
          stop("Rows must have staggered-adoption missingness.")
        }
        
        first_unobs = which(!obs)[1]
        
        if (!is.na(first_unobs)) {
          A[i, j] = first_unobs
        }
      }
    }
  }
  
  if (!is.null(A)) {
    A = as.matrix(A)
    
    if (nrow(A) != N || ncol(A) != K) {
      stop("A must be an N x K adoption-time matrix.")
    }
  }
  
  if (is.null(Omega)) {
    Omega = array(FALSE, dim = c(N, Tt, K))
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        if (is.infinite(A[i, j])) {
          Omega[i, , j] = TRUE
        } else if (A[i, j] > 1) {
          last_obs = min(Tt, A[i, j] - 1)
          Omega[i, seq_len(last_obs), j] = TRUE
        }
      }
    }
  } else {
    if (!is.array(Omega) || !all(dim(Omega) == c(N, Tt, K))) {
      stop("Omega must be an N x T x K logical array.")
    }
    
    Omega = array(as.logical(Omega), dim = c(N, Tt, K))
    
    if (is.null(A)) {
      A = matrix(Inf, nrow = N, ncol = K)
      
      for (j in seq_len(K)) {
        for (i in seq_len(N)) {
          obs = Omega[i, , j]
          
          if (any(diff(as.integer(obs)) > 0)) {
            stop("Rows must have staggered-adoption missingness.")
          }
          
          first_unobs = which(!obs)[1]
          
          if (!is.na(first_unobs)) {
            A[i, j] = first_unobs
          }
        }
      }
    }
  }
  
  # -------------------------------------------------------------
  # Compute the target-layer row permutation.
  # -------------------------------------------------------------
  
  obs_len_k = ifelse(is.infinite(A[, k]), Tt, pmin(Tt, A[, k] - 1))
  
  row_perm = order(-obs_len_k, seq_len(N))
  
  obs_len_k_sorted = obs_len_k[row_perm]
  
  if (any(diff(obs_len_k_sorted) > 0)) {
    stop("Target-layer ordering failed.")
  }
  
  # -------------------------------------------------------------
  # Construct target-layer staircase blocks after sorting.
  # -------------------------------------------------------------
  
  m_desc = unique(obs_len_k_sorted)
  o_k = length(m_desc)
  # print(o_k)
  
  if (m_desc[1] != Tt) {
    stop("Target layer must contain at least one fully observed row block.")
  }
  
  if (tail(m_desc, 1) <= 0) {
    stop("Target layer must contain a non-empty initial observed time period.")
  }
  
  row_parts = lapply(m_desc, function(m) which(obs_len_k_sorted == m))
  
  original_row_parts = lapply(row_parts, function(idx) row_perm[idx])
  
  m_asc = rev(m_desc)
  T_parts = diff(c(0, m_asc))
  
  col_ends = cumsum(T_parts)
  col_starts = c(1, head(col_ends, -1) + 1)
  col_parts = Map(seq, col_starts, col_ends)
  
  if (length(col_parts) != o_k) {
    stop("Row and column partitions have different lengths.")
  }
  
  if (functional == "RowHet") {
    if (is.null(eta)) {
      stop("eta must be supplied for functional = 'RowHet'.")
    }
    
    if (length(eta) != N) {
      stop("eta must have length N.")
    }
    
    if (!all(eta %in% c(-1, 1))) {
      stop("eta should be a vector in {+1, -1}^N for RowHet.")
    }
  }
  
  if (functional == "Local") {
    if (is.null(row_index)) {
      stop("row_index must be supplied for functional = 'Local'.")
    }
    
    if (length(row_index) != 1 ||
        row_index < 1 ||
        row_index > N ||
        row_index != as.integer(row_index)) {
      stop("row_index must be a single valid original row index.")
    }
    
    sorted_position = match(row_index, row_perm)
    
    if (is.na(sorted_position)) {
      stop("Could not locate row_index after target-layer sorting.")
    }
    
    a0 = which(vapply(
      row_parts,
      function(idx) sorted_position %in% idx,
      logical(1)
    ))
    
    if (length(a0) != 1) {
      stop("Could not locate row_index in the target-layer row partition.")
    }
    
    local_pos = match(sorted_position, row_parts[[a0]])
    
    if (is.na(local_pos)) {
      stop("Could not compute the local position of row_index.")
    }
  }
  
  # -------------------------------------------------------------
  # Aggregate over target policy-on / missing blocks
  # -------------------------------------------------------------
  
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      
      if (a + b <= o_k + 1) {
        next
      }
      
      if (functional == "Local" && a != a0) {
        next
      }
      
      Nik = length(row_parts[[a]])
      Tbk = length(col_parts[[b]])
      
      if (functional == "ATE") {
        
        x = rep(1 / sqrt(Nik), Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "RowHet") {
        
        rows_original = original_row_parts[[a]]
        
        x = eta[rows_original] / sqrt(Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "Local") {
        
        x = rep(0, Nik)
        x[local_pos] = 1
        
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Tbk)
        normalizer_increment = Tbk
        
      } else if (functional == "Trend") {
        
        if (Tbk <= 1) {
          next
        }
        
        z = seq_len(Tbk)
        z_centered = z - mean(z)
        
        x = rep(1 / sqrt(Nik), Nik)
        y = z_centered / sqrt(sum(z_centered^2))
        
        weight = 1 / sqrt(Nik * Tbk * (Tbk^2 - 1) / 12)
        
        normalizer_increment = 1
      }
      
      mu_hat = bilinearTensorStaggered(
        Y     = Y,
        k     = k,
        a     = a,
        b     = b,
        r     = r,
        x     = x,
        y     = y,
        A     = A,
        Omega = Omega,
        tau   = tau
      )
      
      weighted_sum = weighted_sum + weight * mu_hat
      normalizer = normalizer + normalizer_increment
    }
  }
  
  if (normalizer == 0) {
    stop("No target blocks were included in the aggregation.")
  }
  
  weighted_sum / normalizer
}


# 
# Wrappers
# 

bilinearTensorStaggeredATE = function(Y, k, r, tau,
                                      A = NULL,
                                      Omega = NULL) {
  bilinearTensorStaggeredPsi(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "ATE",
    A = A,
    Omega = Omega
  )
}


bilinearTensorStaggeredRowHet = function(Y, k, r, eta, tau,
                                         A = NULL,
                                         Omega = NULL) {
  bilinearTensorStaggeredPsi(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "RowHet",
    eta = eta,
    A = A,
    Omega = Omega
  )
}


bilinearTensorStaggeredLocal = function(Y, k, r, row_index, tau,
                                        A = NULL,
                                        Omega = NULL) {
  bilinearTensorStaggeredPsi(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "Local",
    row_index = row_index,
    A = A,
    Omega = Omega
  )
}


bilinearTensorStaggeredTrend = function(Y, k, r, tau,
                                        A = NULL,
                                        Omega = NULL) {
  bilinearTensorStaggeredPsi(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "Trend",
    A = A,
    Omega = Omega
  )
}


# ---------------------------------------------------------------
# 7. bilinearMatrixStaggeredPsi
#
# Estimates Psi_0^{(h)} over all missing staggered blocks using
# one N x T matrix Y_mat.
#
# Args:
#   Y_mat      : N x T matrix with staggered missingness
#   r          : target rank
#   tau        : threshold for H^\dagger
#   functional : one of "ATE", "RowHet", "Local", "Trend"
#   eta        : length-N vector in {+1,-1}; required for RowHet.
#                eta is indexed in the original, unsorted row order.
#   row_index  : original global row index; required for Local
#   A          : optional length-N adoption-time vector,
#                Omega[i,t] = 1{t < A[i]};
#                Inf means never adopter
#   Omega      : optional N x T logical observation mask
#
# Returns:
#   scalar estimate of Psi_0^{(functional)}(k)
# ---------------------------------------------------------------

bilinearMatrixStaggeredPsi = function(Y_mat, r, tau,
                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                      eta = NULL,
                                      row_index = NULL,
                                      A = NULL,
                                      Omega = NULL) {
  functional = match.arg(functional)
  N = nrow(Y_mat)
  Tt = ncol(Y_mat)
  
  # -------------------------------------------------------------
  # Infer A and Omega
  # -------------------------------------------------------------
  
  if (is.null(A) && is.null(Omega)) {
    Omega = !is.na(Y_mat)
    
    A = rep(Inf, N)
    
    for (i in seq_len(N)) {
      obs = Omega[i, ]
      
      if (any(diff(as.integer(obs)) > 0)) {
        stop("Rows must have staggered-adoption missingness.")
      }
      
      first_unobs = which(!obs)[1]
      
      if (!is.na(first_unobs)) {
        A[i] = first_unobs
      }
    }
  }
  
  if (!is.null(A)) {
    if (length(A) != N) {
      stop("A must have length nrow(Y_mat).")
    }
    
    A = as.numeric(A)
  }
  
  if (is.null(Omega)) {
    Omega = matrix(FALSE, nrow = N, ncol = Tt)
    
    for (i in seq_len(N)) {
      if (is.infinite(A[i])) {
        Omega[i, ] = TRUE
      } else if (A[i] > 1) {
        last_obs = min(Tt, A[i] - 1)
        Omega[i, seq_len(last_obs)] = TRUE
      }
    }
  } else {
    if (!is.matrix(Omega) || !identical(dim(Omega), c(N, Tt))) {
      stop("Omega must be an N x T logical matrix.")
    }
    
    Omega = matrix(as.logical(Omega), nrow = N, ncol = Tt)
    
    if (is.null(A)) {
      A = rep(Inf, N)
      
      for (i in seq_len(N)) {
        obs = Omega[i, ]
        
        if (any(diff(as.integer(obs)) > 0)) {
          stop("Rows must have staggered-adoption missingness.")
        }
        
        first_unobs = which(!obs)[1]
        
        if (!is.na(first_unobs)) {
          A[i] = first_unobs
        }
      }
    }
  }
  
  # -------------------------------------------------------------
  # Compute row permutation inducing staircase form.
  # -------------------------------------------------------------
  
  obs_len = ifelse(is.infinite(A), Tt, pmin(Tt, A - 1))
  
  row_perm = order(-obs_len, seq_len(N))
  
  obs_len_sorted = obs_len[row_perm]
  
  if (any(diff(obs_len_sorted) > 0)) {
    stop("Row ordering failed.")
  }
  
  # -------------------------------------------------------------
  # Construct staircase row and column blocks after sorting.
  # -------------------------------------------------------------
  
  m_desc = unique(obs_len_sorted)
  o = length(m_desc)
  # print(o)
  
  if (m_desc[1] != Tt) {
    stop("Matrix must contain at least one fully observed row block.")
  }
  
  if (tail(m_desc, 1) <= 0) {
    stop("Matrix must contain a non-empty initial observed time period.")
  }
  
  row_parts = lapply(m_desc, function(m) which(obs_len_sorted == m))
  
  original_row_parts = lapply(row_parts, function(idx) row_perm[idx])
  
  m_asc = rev(m_desc)
  T_part = diff(c(0, m_asc))
  
  col_ends = cumsum(T_part)
  col_starts = c(1, head(col_ends, -1) + 1)
  col_parts = Map(seq, col_starts, col_ends)
  
  if (length(col_parts) != o) {
    stop("Row and column partitions have different lengths.")
  }
  
  # -------------------------------------------------------------
  # Validate functional-specific inputs
  # -------------------------------------------------------------
  
  if (functional == "RowHet") {
    if (is.null(eta)) {
      stop("eta must be supplied for functional = 'RowHet'.")
    }
    
    if (length(eta) != N) {
      stop("eta must have length N.")
    }
    
    if (!all(eta %in% c(-1, 1))) {
      stop("eta should be a vector in {+1, -1}^N for RowHet.")
    }
  }
  
  if (functional == "Local") {
    if (is.null(row_index)) {
      stop("row_index must be supplied for functional = 'Local'.")
    }
    
    if (length(row_index) != 1 ||
        row_index < 1 ||
        row_index > N ||
        row_index != as.integer(row_index)) {
      stop("row_index must be a single valid original row index.")
    }
    
    sorted_position = match(row_index, row_perm)
    
    if (is.na(sorted_position)) {
      stop("Could not locate row_index after row sorting.")
    }
    
    a0 = which(vapply(
      row_parts,
      function(idx) sorted_position %in% idx,
      logical(1)
    ))
    
    if (length(a0) != 1) {
      stop("Could not locate row_index in the row partition.")
    }
    
    local_pos = match(sorted_position, row_parts[[a0]])
    
    if (is.na(local_pos)) {
      stop("Could not compute the local position of row_index.")
    }
  }
  
  # -------------------------------------------------------------
  # Aggregate over missing staggered blocks
  # -------------------------------------------------------------
  
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o)) {
    for (b in seq_len(o)) {
      
      if (a + b <= o + 1) {
        next
      }
      
      if (functional == "Local" && a != a0) {
        next
      }
      
      Nik = length(row_parts[[a]])
      Tbk = length(col_parts[[b]])
      
      if (functional == "ATE") {
        
        x = rep(1 / sqrt(Nik), Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "RowHet") {
        
        rows_original = original_row_parts[[a]]
        
        x = eta[rows_original] / sqrt(Nik)
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Nik * Tbk)
        normalizer_increment = Nik * Tbk
        
      } else if (functional == "Local") {
        
        x = rep(0, Nik)
        x[local_pos] = 1
        
        y = rep(1 / sqrt(Tbk), Tbk)
        
        weight = sqrt(Tbk)
        normalizer_increment = Tbk
        
      } else if (functional == "Trend") {
        
        if (Tbk <= 1) {
          next
        }
        
        z = seq_len(Tbk)
        z_centered = z - mean(z)
        
        x = rep(1 / sqrt(Nik), Nik)
        y = z_centered / sqrt(sum(z_centered^2))
        
        weight = 1 / sqrt(Nik * Tbk * (Tbk^2 - 1) / 12)
        
        normalizer_increment = 1
      }
      
      mu_hat = bilinearMatrixStaggered(
        Y_mat = Y_mat,
        a     = a,
        b     = b,
        r     = r,
        x     = x,
        y     = y,
        A     = A,
        Omega = Omega,
        tau   = tau
      )
      
      weighted_sum = weighted_sum + weight * mu_hat
      normalizer = normalizer + normalizer_increment
    }
  }
  
  if (normalizer == 0) {
    stop("No target blocks were included in the aggregation.")
  }
  
  weighted_sum / normalizer
}


# 
# Wrappers
# 

bilinearMatrixStaggeredATE = function(Y_mat, r, tau,
                                      A = NULL,
                                      Omega = NULL) {
  bilinearMatrixStaggeredPsi(
    Y_mat = Y_mat,
    r = r,
    tau = tau,
    functional = "ATE",
    A = A,
    Omega = Omega
  )
}


bilinearMatrixStaggeredRowHet = function(Y_mat, r, eta, tau,
                                         A = NULL,
                                         Omega = NULL) {
  bilinearMatrixStaggeredPsi(
    Y_mat = Y_mat,
    r = r,
    tau = tau,
    functional = "RowHet",
    eta = eta,
    A = A,
    Omega = Omega
  )
}


bilinearMatrixStaggeredLocal = function(Y_mat, r, row_index, tau,
                                        A = NULL,
                                        Omega = NULL) {
  bilinearMatrixStaggeredPsi(
    Y_mat = Y_mat,
    r = r,
    tau = tau,
    functional = "Local",
    row_index = row_index,
    A = A,
    Omega = Omega
  )
}


bilinearMatrixStaggeredTrend = function(Y_mat, r, tau,
                                        A = NULL,
                                        Omega = NULL) {
  bilinearMatrixStaggeredPsi(
    Y_mat = Y_mat,
    r = r,
    tau = tau,
    functional = "Trend",
    A = A,
    Omega = Omega
  )
}




# ---------------------------------------------------------------
# 8. pluginPsi_c1
#
# Estimates Psi_1^{(h)} over all c = 1 staggered blocks using
# the observed entries in the kth slice.
#
# Args:
#   Y          : N x T x K tensor
#   k          : slice index
#   N_parts    : list of row-block sizes for each partition
#   T_parts    : list of column-block sizes for each partition
#   functional : one of "ATE", "RowHet", "Local", "Trend"
#   eta        : length-N vector in {+1,-1}; required for RowHet.
#                eta is indexed in the original row order.
#   row_index  : global row index; required for Local
#
# Returns:
#   scalar estimate of Psi_1^{(functional)}
# ---------------------------------------------------------------

pluginPsi_c1 = function(Y, k, N_parts, T_parts,
                        functional = c("ATE", "RowHet", "Local", "Trend"),
                        eta = NULL,
                        row_index = NULL) {
  
  functional = match.arg(functional)
  
  make_parts = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  row_parts = make_parts(N_parts[[k]])
  col_parts = make_parts(T_parts[[k]])
  o_k = length(row_parts)
  
  Yk = if (length(dim(Y)) == 3) Y[, , k] else Y
  
  if (functional == "RowHet") stopifnot(!is.null(eta))
  
  if (functional == "Local") {
    stopifnot(!is.null(row_index))
    a0 = which(vapply(row_parts, function(idx) row_index %in% idx, logical(1)))
    stopifnot(length(a0) == 1)
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      
      ## c = 1 region: policy-on / complementary staggered region
      if (a + b <= o_k + 1) next
      
      if (functional == "Local" && a != a0) next
      
      rows = row_parts[[a]]
      cols = col_parts[[b]]
      
      Nik = length(rows)
      Tbk = length(cols)
      Y_ab = Yk[rows, cols, drop = FALSE]
      
      if (functional == "ATE") {
        
        val = mean(Y_ab, na.rm = TRUE)
        weight = Nik * Tbk
        
      } else if (functional == "RowHet") {
        
        val = sum(eta[rows] * rowSums(Y_ab, na.rm = TRUE)) / (Nik * Tbk)
        weight = Nik * Tbk
        
      } else if (functional == "Local") {
        
        val = mean(Yk[row_index, cols], na.rm = TRUE)
        weight = Tbk
        
      } else if (functional == "Trend") {
        
        if (Tbk <= 1) {
          next
        }
        
        z = seq_len(Tbk)
        zc = z - mean(z)
        
        y = zc / sqrt(sum(zc^2))
        x = rep(1 / sqrt(Nik), Nik)
        
        mu_ab = as.numeric(t(x) %*% Y_ab %*% y)
        val = mu_ab / sqrt(Nik * Tbk * (Tbk^2 - 1) / 12)
        
        weight = 1
      }
      
      weighted_sum = weighted_sum + weight * val
      normalizer = normalizer + weight
    }
  }
  
  weighted_sum / normalizer
}

## 
## Wrappers
## 

pluginPsi1_ATE = function(Y, k, N_parts, T_parts) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "ATE"
  )
}

pluginPsi1_RowHet = function(Y, k, N_parts, T_parts, eta) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "RowHet",
    eta = eta
  )
}

pluginPsi1_Local = function(Y, k, N_parts, T_parts, row_index) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "Local",
    row_index = row_index
  )
}

pluginPsi1_Trend = function(Y, k, N_parts, T_parts) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "Trend"
  )
}

# ---------------------------------------------------------------
# 9. Reduced-Anchor Bilinear Estimation for Staggered-Adoption Tensor Data
#
# Inputs:
#   Y      : N x T x K array, or list of K matrices of size N x T
#   k      : target layer index
#   r      : requested rank
#   tau    : eigenvalue threshold for stable inverse
#   A      : optional N x K adoption-time matrix
#   Omega  : optional N x T x K logical observation mask
#
# Output:
#   scalar estimate of Psi_0^{(functional)}(k)
# ---------------------------------------------------------------

bilinearTensorStaggeredLinearReducedAnchor = function(Y, k, r, tau,
                                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                                      eta = NULL,
                                                      row_index = NULL,
                                                      A = NULL,
                                                      Omega = NULL,
                                                      return_fit = FALSE) {
  functional = match.arg(functional)
  stopifnot(r >= 1, tau > 0)
  
  if (is.array(Y)) {
    N = dim(Y)[1]
    Tt = dim(Y)[2]
    K = dim(Y)[3]
    Y_list = vector("list", K)
    for (j in seq_len(K)) Y_list[[j]] = Y[, , j]
  } else {
    Y_list = Y
    K = length(Y_list)
    N = nrow(Y_list[[1]])
    Tt = ncol(Y_list[[1]])
  }
  
  if (k < 1 || k > K) stop("k must index a layer of Y.")
  
  if (is.null(A) && is.null(Omega)) {
    Omega = array(FALSE, c(N, Tt, K))
    for (j in seq_len(K)) Omega[, , j] = !is.na(Y_list[[j]])
  }
  
  if (!is.null(Omega)) {
    Omega = array(as.logical(Omega), c(N, Tt, K))
  }
  
  if (is.null(A)) {
    A = matrix(Inf, N, K)
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        observed = Omega[i, , j]
        
        if (any(diff(as.integer(observed)) > 0)) {
          stop("Rows must have staggered-adoption missingness.")
        }
        
        first_unobserved = which(!observed)[1]
        if (!is.na(first_unobserved)) A[i, j] = first_unobserved
      }
    }
  } else {
    A = as.matrix(A)
  }
  
  if (is.null(Omega)) {
    Omega = array(FALSE, c(N, Tt, K))
    
    for (j in seq_len(K)) {
      for (i in seq_len(N)) {
        if (is.infinite(A[i, j])) {
          Omega[i, , j] = TRUE
        } else if (A[i, j] > 1) {
          Omega[i, seq_len(min(Tt, A[i, j] - 1)), j] = TRUE
        }
      }
    }
  }
  
  target_observed_length = ifelse(
    is.infinite(A[, k]),
    Tt,
    pmin(Tt, A[, k] - 1)
  )
  
  row_order = order(-target_observed_length, seq_len(N))
  
  Y_sorted = vector("list", K)
  for (j in seq_len(K)) {
    Y_sorted[[j]] = Y_list[[j]][row_order, , drop = FALSE]
  }
  
  A_sorted = A[row_order, , drop = FALSE]
  Omega_sorted = Omega[row_order, , , drop = FALSE]
  observed_length_sorted = target_observed_length[row_order]
  
  row_observed_lengths = unique(observed_length_sorted)
  block_count = length(row_observed_lengths)
  
  if (row_observed_lengths[1] != Tt) {
    stop("Target layer must contain at least one fully observed row block.")
  }
  
  if (tail(row_observed_lengths, 1) <= 0) {
    stop("Target layer must contain a non-empty initial observed time period.")
  }
  
  row_blocks = vector("list", block_count)
  for (a in seq_len(block_count)) {
    row_blocks[[a]] = which(
      observed_length_sorted == row_observed_lengths[a]
    )
  }
  
  time_block_sizes = diff(c(0, rev(row_observed_lengths)))
  time_block_ends = cumsum(time_block_sizes)
  time_block_starts = c(1, head(time_block_ends, -1) + 1)
  
  col_blocks = vector("list", block_count)
  for (b in seq_len(block_count)) {
    col_blocks[[b]] = seq.int(
      time_block_starts[b],
      time_block_ends[b]
    )
  }
  
  blocks = list(
    Y_list = Y_sorted,
    A = A_sorted,
    Omega = Omega_sorted,
    row_perm = row_order,
    obs_len_k = observed_length_sorted,
    row_parts = row_blocks,
    col_parts = col_blocks,
    N_part = lengths(row_blocks),
    T_part = lengths(col_blocks),
    o_k = block_count,
    N = N,
    Tt = Tt,
    K = K
  )
  
  if (functional == "RowHet") {
    if (
      is.null(eta) ||
      length(eta) != N ||
      !all(eta %in% c(-1, 1))
    ) {
      stop("eta must be a vector in {+1, -1}^N.")
    }
  }
  
  local_row_block = NULL
  local_row_position = NULL
  
  if (functional == "Local") {
    if (
      length(row_index) != 1 ||
      row_index < 1 ||
      row_index > N ||
      row_index != as.integer(row_index)
    ) {
      stop("row_index must be a single valid original row index.")
    }
    
    sorted_position = match(row_index, row_order)
    
    for (a in seq_len(block_count)) {
      local_row_position = match(
        sorted_position,
        row_blocks[[a]]
      )
      
      if (!is.na(local_row_position)) {
        local_row_block = a
        break
      }
    }
  }
  
  active = list()
  
  for (a in seq_len(block_count)) {
    for (b in seq_len(block_count)) {
      if (a + b <= block_count + 1) next
      if (functional == "Local" && a != local_row_block) next
      if (functional == "Trend" && length(col_blocks[[b]]) <= 1) next
      
      active[[length(active) + 1]] = c(a = a, b = b)
    }
  }
  
  if (!length(active)) {
    stop("No target blocks were included in the aggregation.")
  }
  
  D = do.call(rbind, active)
  active_row_blocks = sort(unique(D[, "a"]))
  active_col_blocks = sort(unique(D[, "b"]))
  
  first_rows = row_blocks[[1]]
  first_cols = col_blocks[[1]]
  r_eff = r
  
  for (a in active_row_blocks) {
    rows_for_fit = unlist(
      row_blocks[seq_len(a)],
      use.names = FALSE
    )
    
    target_cols_for_fit = unlist(
      col_blocks[seq_len(block_count + 1 - a)],
      use.names = FALSE
    )
    
    available_cols = 0
    
    for (j in seq_len(K)) {
      if (j == k) {
        available_cols =
          available_cols + length(target_cols_for_fit)
      } else {
        observed = matrix(
          as.logical(
            Omega_sorted[
              rows_for_fit,
              first_cols,
              j,
              drop = FALSE
            ]
          ),
          nrow = length(rows_for_fit),
          ncol = length(first_cols)
        )
        
        available_cols =
          available_cols +
          sum(colSums(observed) == length(rows_for_fit))
      }
    }
    
    r_eff = min(
      r_eff,
      length(rows_for_fit),
      available_cols
    )
  }
  
  for (b in active_col_blocks) {
    rows_for_fit = unlist(
      row_blocks[seq_len(block_count + 1 - b)],
      use.names = FALSE
    )
    
    cols_for_fit = unlist(
      col_blocks[seq_len(b)],
      use.names = FALSE
    )
    
    available_rows = 0
    
    for (j in seq_len(K)) {
      if (j == k) {
        available_rows =
          available_rows + length(rows_for_fit)
      } else {
        observed = matrix(
          as.logical(
            Omega_sorted[
              first_rows,
              cols_for_fit,
              j,
              drop = FALSE
            ]
          ),
          nrow = length(first_rows),
          ncol = length(cols_for_fit)
        )
        
        available_rows =
          available_rows +
          sum(rowSums(observed) == length(cols_for_fit))
      }
    }
    
    r_eff = min(
      r_eff,
      available_rows,
      length(cols_for_fit)
    )
  }
  
  if (r_eff < 1) stop("Effective rank is zero.")
  
  left_fits = vector("list", block_count)
  
  for (a in active_row_blocks) {
    rows_for_fit = unlist(
      row_blocks[seq_len(a)],
      use.names = FALSE
    )
    
    target_rows = row_blocks[[a]]
    
    target_layer_cols = unlist(
      col_blocks[seq_len(block_count + 1 - a)],
      use.names = FALSE
    )
    
    column_anchors = vector("list", K)
    
    for (j in seq_len(K)) {
      if (j == k) {
        column_anchors[[j]] = target_layer_cols
      } else {
        observed = matrix(
          as.logical(
            Omega_sorted[
              rows_for_fit,
              first_cols,
              j,
              drop = FALSE
            ]
          ),
          nrow = length(rows_for_fit),
          ncol = length(first_cols)
        )
        
        column_anchors[[j]] = first_cols[
          colSums(observed) == length(rows_for_fit)
        ]
      }
    }
    
    left_mats = list()
    
    for (j in seq_len(K)) {
      cols = column_anchors[[j]]
      
      if (length(cols)) {
        left_mats[[length(left_mats) + 1]] =
          Y_sorted[[j]][
            rows_for_fit,
            cols,
            drop = FALSE
          ]
      }
    }
    
    Y_left = do.call(cbind, left_mats)
    svd_left = svd(Y_left, nu = r_eff, nv = 0)
    
    U_left = svd_left$u[
      ,
      seq_len(r_eff),
      drop = FALSE
    ]
    
    target_rows_local = match(
      target_rows,
      rows_for_fit
    )
    
    n_rows = length(target_rows)
    
    if (functional == "RowHet") {
      x = eta[row_order[target_rows]] / sqrt(n_rows)
    } else if (functional == "Local") {
      x = numeric(n_rows)
      x[local_row_position] = 1
    } else {
      x = rep(1 / sqrt(n_rows), n_rows)
    }
    
    alpha = as.vector(
      crossprod(
        U_left[target_rows_local, , drop = FALSE],
        x
      )
    )
    
    left_fits[[a]] = list(
      U_left = U_left,
      alpha = alpha,
      S_a = rows_for_fit,
      R_a_local = target_rows_local,
      ColAnc = column_anchors,
      Y_left_dim = dim(Y_left)
    )
  }
  
  upper_fits = vector("list", block_count)
  
  for (b in active_col_blocks) {
    target_layer_rows = unlist(
      row_blocks[seq_len(block_count + 1 - b)],
      use.names = FALSE
    )
    
    cols_for_fit = unlist(
      col_blocks[seq_len(b)],
      use.names = FALSE
    )
    
    target_cols = col_blocks[[b]]
    row_anchors = vector("list", K)
    
    for (j in seq_len(K)) {
      if (j == k) {
        row_anchors[[j]] = target_layer_rows
      } else {
        observed = matrix(
          as.logical(
            Omega_sorted[
              first_rows,
              cols_for_fit,
              j,
              drop = FALSE
            ]
          ),
          nrow = length(first_rows),
          ncol = length(cols_for_fit)
        )
        
        row_anchors[[j]] = first_rows[
          rowSums(observed) == length(cols_for_fit)
        ]
      }
    }
    
    upper_mats = list()
    
    for (j in seq_len(K)) {
      rows = row_anchors[[j]]
      
      if (length(rows)) {
        upper_mats[[length(upper_mats) + 1]] =
          Y_sorted[[j]][
            rows,
            cols_for_fit,
            drop = FALSE
          ]
      }
    }
    
    Y_upper = do.call(rbind, upper_mats)
    svd_upper = svd(
      Y_upper,
      nu = r_eff,
      nv = r_eff
    )
    
    U_upper = svd_upper$u[
      ,
      seq_len(r_eff),
      drop = FALSE
    ]
    
    V_upper = svd_upper$v[
      ,
      seq_len(r_eff),
      drop = FALSE
    ]
    
    singular_values = svd_upper$d[seq_len(r_eff)]
    
    target_layer_start = 0
    
    if (k > 1) {
      for (j in seq_len(k - 1)) {
        target_layer_start =
          target_layer_start + length(row_anchors[[j]])
      }
    }
    
    target_layer_indices =
      target_layer_start + seq_along(target_layer_rows)
    
    U_target = U_upper[
      target_layer_indices,
      ,
      drop = FALSE
    ]
    
    target_cols_local = match(
      target_cols,
      cols_for_fit
    )
    
    V_target = V_upper[
      target_cols_local,
      ,
      drop = FALSE
    ]
    
    n_cols = length(target_cols)
    
    if (functional == "Trend") {
      y = seq_len(n_cols)
      y = y - mean(y)
      y = y / sqrt(sum(y^2))
    } else {
      y = rep(1 / sqrt(n_cols), n_cols)
    }
    
    projected_y = crossprod(V_target, y)
    
    X = as.vector(
      U_target %*%
        (singular_values * as.vector(projected_y))
    )
    
    upper_fits[[b]] = list(
      X = X,
      S_plus = target_layer_rows,
      Q_b_all = cols_for_fit,
      C_b_local = target_cols_local,
      RowAnc = row_anchors,
      U_up_k = U_target,
      V_b_hat = V_target,
      D_up = singular_values,
      Y_up_dim = dim(Y_upper)
    )
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (ell in seq_len(nrow(D))) {
    a = D[ell, "a"]
    b = D[ell, "b"]
    
    left_fit = left_fits[[a]]
    upper_fit = upper_fits[[b]]
    
    common_rows = match(
      upper_fit$S_plus,
      left_fit$S_a
    )
    
    U_common = left_fit$U_left[
      common_rows,
      ,
      drop = FALSE
    ]
    
    gram = crossprod(U_common)
    eig = eigen(gram, symmetric = TRUE)
    
    gram_inverse =
      eig$vectors %*%
      diag(
        1 / pmax(eig$values, tau),
        nrow = length(eig$values)
      ) %*%
      t(eig$vectors)
    
    beta = as.vector(
      gram_inverse %*%
        crossprod(U_common, upper_fit$X)
    )
    
    n_rows = length(row_blocks[[a]])
    n_cols = length(col_blocks[[b]])
    
    if (functional %in% c("ATE", "RowHet")) {
      weight = sqrt(n_rows * n_cols)
      normalizer_increment = n_rows * n_cols
    } else if (functional == "Local") {
      weight = sqrt(n_cols)
      normalizer_increment = n_cols
    } else {
      weight =
        1 /
        sqrt(
          n_rows *
            n_cols *
            (n_cols^2 - 1) /
            12
        )
      
      normalizer_increment = 1
    }
    
    weighted_sum =
      weighted_sum +
      weight * sum(left_fit$alpha * beta)
    
    normalizer =
      normalizer +
      normalizer_increment
  }
  
  psi_hat = weighted_sum / normalizer
  
  if (!return_fit) return(psi_hat)
  
  list(
    psi_hat = psi_hat,
    weighted_sum = weighted_sum,
    normalizer = normalizer,
    r_eff = r_eff,
    blocks = blocks,
    active_blocks = D,
    left_fits = left_fits,
    upper_fits = upper_fits,
    functional = functional
  )
}

##
# Wrappers
##

bilinearTensorStaggeredATELinearReducedAnchor = function(Y, k, r, tau,
                                                         A = NULL,
                                                         Omega = NULL,
                                                         return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y, k, r, tau,
    functional = "ATE",
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredRowHetLinearReducedAnchor = function(Y, k, r, eta, tau,
                                                            A = NULL,
                                                            Omega = NULL,
                                                            return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y, k, r, tau,
    functional = "RowHet",
    eta = eta,
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredLocalLinearReducedAnchor = function(Y, k, r, row_index, tau,
                                                           A = NULL,
                                                           Omega = NULL,
                                                           return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y, k, r, tau,
    functional = "Local",
    row_index = row_index,
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredTrendLinearReducedAnchor = function(Y, k, r, tau,
                                                           A = NULL,
                                                           Omega = NULL,
                                                           return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y, k, r, tau,
    functional = "Trend",
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}