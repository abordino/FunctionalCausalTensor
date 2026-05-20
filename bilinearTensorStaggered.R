# ---------------------------------------------------------------
# bilinearTensorStaggered
#
# Args:
#   Y       : list of K N x T matrices, or N x T x K array
#   k       : target slice
#   i0,t0   : target staggered block indices in slice k
#   r       : target rank
#   x       : unit vector in R^{N_{i0,k}}
#   y       : unit vector in R^{T_{t0,k}}
#   N_parts : list of length K; N_parts[[j]] = c(N_{1j}, ..., N_{o_j j})
#   T_parts : list of length K; T_parts[[j]] = c(T_{1j}, ..., T_{o_j j})
#   tau     : threshold for H_k^\dagger
#
# Returns:
#   scalar estimate of mu_{xy}^{(k,i0,t0)}
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

bilinearTensorStaggered = function(Y, k, i0, t0, r, x, y,
                                      N_parts, T_parts, tau) {
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
  
  stopifnot(length(N_parts) == K, length(T_parts) == K)
  stopifnot(k >= 1, k <= K)
  
  make_partition_indices = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  row_parts = lapply(N_parts, make_partition_indices)
  col_parts = lapply(T_parts, make_partition_indices)
  
  o = vapply(N_parts, length, integer(1))
  
  for (j in seq_len(K)) {
    if (length(N_parts[[j]]) != length(T_parts[[j]])) {
      stop("N_parts[[j]] and T_parts[[j]] must have the same length for every j.")
    }
    if (sum(N_parts[[j]]) != N) {
      stop("Each N_parts[[j]] must sum to N.")
    }
    if (sum(T_parts[[j]]) != Tt) {
      stop("Each T_parts[[j]] must sum to T.")
    }
  }
  
  ok = o[k]
  
  if (i0 < 1 || i0 > ok || t0 < 1 || t0 > ok) {
    stop("i0 and t0 must be valid block indices for slice k.")
  }
  
  if (i0 + t0 <= ok + 1) {
    stop("Target block is observed under staggered missingness. Need i0 + t0 > o_k + 1.")
  }
  
  if (length(x) != N_parts[[k]][i0]) {
    stop("length(x) must equal N_parts[[k]][i0].")
  }
  
  if (length(y) != T_parts[[k]][t0]) {
    stop("length(y) must equal T_parts[[k]][t0].")
  }
  
  a_plus_end  = ok + 1 - t0
  a_minus_beg = ok + 2 - t0
  
  b_plus_end  = ok + 1 - i0
  b_minus_beg = ok + 2 - i0
  
  S_plus = unlist(row_parts[[k]][seq_len(a_plus_end)], use.names = FALSE)
  S_minus = unlist(row_parts[[k]][a_minus_beg:i0], use.names = FALSE)
  
  Q_plus = unlist(col_parts[[k]][seq_len(b_plus_end)], use.names = FALSE)
  Q_minus = unlist(col_parts[[k]][b_minus_beg:t0], use.names = FALSE)
  
  S_idx = c(S_plus, S_minus)
  Q_idx = c(Q_plus, Q_minus)
  
  n_aux = length(S_idx)
  t_aux = length(Q_idx)
  
  S_plus_local = seq_len(length(S_plus))
  S_minus_local = (length(S_plus) + 1):n_aux
  
  Q_plus_local = seq_len(length(Q_plus))
  Q_minus_local = (length(Q_plus) + 1):t_aux
  
  make_omega = function(j) {
    row_block = rep(seq_along(N_parts[[j]]), times = N_parts[[j]])
    col_block = rep(seq_along(T_parts[[j]]), times = T_parts[[j]])
    
    outer(row_block, col_block, function(a, b) {
      a + b <= o[j] + 1
    })
  }
  
  Omega_list = lapply(seq_len(K), make_omega)
  

  Y_aux = array(NA_real_, dim = c(n_aux, t_aux, K))
  
  N1_aux = integer(K)
  T1_aux = integer(K)
  
  M_k = matrix(NA_real_, nrow = n_aux, ncol = t_aux)
  obs_k = matrix(FALSE, nrow = n_aux, ncol = t_aux)
  obs_k[S_plus_local, ] = TRUE
  obs_k[, Q_plus_local] = TRUE
  
  Y_sub_k = Y_list[[k]][S_idx, Q_idx, drop = FALSE]
  M_k[obs_k] = Y_sub_k[obs_k]
  
  Y_aux[, , k] = M_k
  
  N1_aux[k] = length(S_plus)
  T1_aux[k] = length(Q_plus)
  
  # S and Q as original indices
  S_set = S_idx
  Q_set = Q_idx
  
  for (j in setdiff(seq_len(K), k)) {
    Omega_j = Omega_list[[j]]
    
    # ColAnc(j)
    ColAnc_j = which(colSums(Omega_j[S_set, , drop = FALSE]) == length(S_set))
    
    # RowAnc(j)
    RowAnc_j = which(rowSums(Omega_j[, Q_set, drop = FALSE]) == length(Q_set))
    
    row_obs_local = which(S_idx %in% RowAnc_j)
    col_obs_local = which(Q_idx %in% ColAnc_j)
    
    N1_aux[j] = length(row_obs_local)
    T1_aux[j] = length(col_obs_local)
    
    if (N1_aux[j] > 0 &&
        !identical(row_obs_local, seq_len(N1_aux[j]))) {
      stop(
        paste0(
          "For slice ", j,
          ", RowAnc(j) intersect S is not a prefix of S. ",
          "Reordering would be needed before calling the four-block estimator."
        )
      )
    }
    
    if (T1_aux[j] > 0 &&
        !identical(col_obs_local, seq_len(T1_aux[j]))) {
      stop(
        paste0(
          "For slice ", j,
          ", ColAnc(j) intersect Q is not a prefix of Q. ",
          "Reordering would be needed before calling the four-block estimator."
        )
      )
    }
    
    M_j = matrix(NA_real_, nrow = n_aux, ncol = t_aux)
    
    obs_j = matrix(FALSE, nrow = n_aux, ncol = t_aux)
    
    if (N1_aux[j] > 0) {
      obs_j[seq_len(N1_aux[j]), ] = TRUE
    }
    
    if (T1_aux[j] > 0) {
      obs_j[, seq_len(T1_aux[j])] = TRUE
    }
    
    Y_sub_j = Y_list[[j]][S_idx, Q_idx, drop = FALSE]
    M_j[obs_j] = Y_sub_j[obs_j]
    
    Y_aux[, , j] = M_j
  }
  
  
  lx = length(S_minus) - N_parts[[k]][i0]
  ly = length(Q_minus) - T_parts[[k]][t0]
  
  if (lx < 0 || ly < 0) {
    stop("Invalid embedding dimensions for x or y.")
  }
  
  x_bar = c(rep(0, lx), x)
  y_bar = c(rep(0, ly), y)
  
  bilinearTensor4Block(
    Y   = Y_aux,
    k   = k,
    r   = r,
    x   = x_bar,
    y   = y_bar,
    tau = tau,
    N1  = N1_aux,
    T1  = T1_aux
  )
}