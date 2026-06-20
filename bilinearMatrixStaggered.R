# ---------------------------------------------------------------
# bilinearMatrixStaggered
#
# Matrix-only wrapper for staggered adoption missingness.
#
# This estimates
#
#   mu_{xy}^{(a,b)} = x' M^{(a,b)} y
#
# using only one N x T matrix Y_mat.
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
#   scalar estimate of x' M^{(a,b)} y
# ---------------------------------------------------------------

bilinearMatrixStaggered = function(Y_mat, a, b, r, x, y,
                                   A = NULL, Omega = NULL, tau) {
  if (!is.matrix(Y_mat)) {
    stop("Y_mat must be a matrix.")
  }
  
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
    stop("Effective rank is zero. Check r and the target-layer anchor blocks.")
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