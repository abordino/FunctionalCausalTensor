# ---------------------------------------------------------------
# LinearATEStaggered
#
# Fast all-missing ATE estimator for the default/all-ones target:
#
#   average over all missing entries in slice k,
#   i.e. over blocks satisfying i + t > o_k + 1.
#
# It uses one auxiliary four-block call per missing column frontier t,
# rather than one call per missing block (i,t).
#
# Requires:
#   bilinearTensor4Block()
# ---------------------------------------------------------------

LinearATEStaggered = function(Y, k, r, N_parts, T_parts, tau,
                              return_details = FALSE) {
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
  
  o = vapply(N_parts, length, integer(1))
  
  for (j in seq_len(K)) {
    if (length(N_parts[[j]]) != length(T_parts[[j]])) {
      stop("N_parts[[j]] and T_parts[[j]] must have the same length for every j.")
    }
    if (sum(N_parts[[j]]) != N) {
      stop("Each N_parts[[j]] must sum to nrow(Y).")
    }
    if (sum(T_parts[[j]]) != Tt) {
      stop("Each T_parts[[j]] must sum to ncol(Y).")
    }
  }
  
  o_k = o[k]
  
  if (o_k < 2) {
    stop("Need at least two staggered blocks in slice k.")
  }
  
 make_partition_indices = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq.int, starts, ends)
  }
  
  concat_blocks = function(parts, a, b) {
    if (a > b) {
      integer(0)
    } else {
      unlist(parts[a:b], use.names = FALSE)
    }
  }
  
  row_parts = lapply(N_parts, make_partition_indices)
  col_parts = lapply(T_parts, make_partition_indices)
  
  make_omega = function(j) {
    row_block = rep(seq_along(N_parts[[j]]), times = N_parts[[j]])
    col_block = rep(seq_along(T_parts[[j]]), times = T_parts[[j]])
    
    outer(row_block, col_block, function(a, b) {
      a + b <= o[j] + 1
    })
  }
  
  Omega_list = lapply(seq_len(K), make_omega)
  
  numerator = 0
  denominator = 0
  
  column_contributions = numeric(o_k)
  column_denominators  = numeric(o_k)
  
  # Missing columns start at t0 = 2.
  # For fixed t0, the missing row blocks are
  # i0 = o_k + 2 - t0, ..., o_k.
  for (t0 in 2:o_k) {
    i_min = o_k + 2 - t0
    i_max = o_k
    
    S_plus  = concat_blocks(row_parts[[k]], 1, o_k + 1 - t0)
    S_minus = concat_blocks(row_parts[[k]], i_min, i_max)
    
    Q_plus  = concat_blocks(col_parts[[k]], 1, 1)
    Q_minus = concat_blocks(col_parts[[k]], 2, t0)
    
    S_idx = c(S_plus, S_minus)
    Q_idx = c(Q_plus, Q_minus)
    
    n_aux = length(S_idx)
    t_aux = length(Q_idx)
    
    S_plus_local  = seq_len(length(S_plus))
    S_minus_local = (length(S_plus) + 1):n_aux
    
    Q_plus_local  = seq_len(length(Q_plus))
    Q_minus_local = (length(Q_plus) + 1):t_aux
    
    Y_aux = array(NA_real_, dim = c(n_aux, t_aux, K))
    
    N1_aux = integer(K)
    T1_aux = integer(K)
    
    M_k = matrix(NA_real_, nrow = n_aux, ncol = t_aux)
    obs_k = matrix(FALSE, nrow = n_aux, ncol = t_aux)
    
    if (length(S_plus_local) > 0) {
      obs_k[S_plus_local, ] = TRUE
    }
    
    if (length(Q_plus_local) > 0) {
      obs_k[, Q_plus_local] = TRUE
    }
    
    Y_sub_k = Y_list[[k]][S_idx, Q_idx, drop = FALSE]
    M_k[obs_k] = Y_sub_k[obs_k]
    
    Y_aux[, , k] = M_k
    
    N1_aux[k] = length(S_plus)
    T1_aux[k] = length(Q_plus)
    
    S_set = S_idx
    Q_set = Q_idx
    
    for (j in setdiff(seq_len(K), k)) {
      Omega_j = Omega_list[[j]]
      
      ColAnc_j = which(
        colSums(Omega_j[S_set, , drop = FALSE]) == length(S_set)
      )
      
      RowAnc_j = which(
        rowSums(Omega_j[, Q_set, drop = FALSE]) == length(Q_set)
      )
      
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
            "Reordering would be needed before calling bilinearTensor4Block()."
          )
        )
      }
      
      if (T1_aux[j] > 0 &&
          !identical(col_obs_local, seq_len(T1_aux[j]))) {
        stop(
          paste0(
            "For slice ", j,
            ", ColAnc(j) intersect Q is not a prefix of Q. ",
            "Reordering would be needed before calling bilinearTensor4Block()."
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
    
    x_query = rep(1, length(S_minus))
    
    q_before_t0 = concat_blocks(col_parts[[k]], 2, t0 - 1)
    q_t0        = concat_blocks(col_parts[[k]], t0, t0)
    
    y_query = c(
      rep(0, length(q_before_t0)),
      rep(1, length(q_t0))
    )
    
    if (length(x_query) != n_aux - N1_aux[k]) {
      stop("Internal error: x_query has wrong length.")
    }
    
    if (length(y_query) != t_aux - T1_aux[k]) {
      stop("Internal error: y_query has wrong length.")
    }
  
    col_hat = bilinearTensor4Block(
      Y   = Y_aux,
      k   = k,
      r   = r,
      x   = x_query,
      y   = y_query,
      tau = tau,
      N1  = N1_aux,
      T1  = T1_aux
    )
    
    col_denom = sum(N_parts[[k]][i_min:i_max]) * T_parts[[k]][t0]
    
    numerator = numerator + col_hat
    denominator = denominator + col_denom
    
    column_contributions[t0] = col_hat
    column_denominators[t0]  = col_denom
  }
  
  if (denominator == 0) {
    stop("No missing staggered entries found.")
  }
  
  ate_hat = numerator / denominator
  
  if (return_details) {
    return(
      list(
        ate = ate_hat,
        numerator = numerator,
        denominator = denominator,
        column_contributions = column_contributions,
        column_denominators = column_denominators
      )
    )
  }
  
  ate_hat
}