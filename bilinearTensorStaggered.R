# ---------------------------------------------------------------
# bilinearTensorStaggered
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