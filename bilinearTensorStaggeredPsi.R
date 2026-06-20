# ---------------------------------------------------------------
# bilinearTensorStaggeredPsi
#
# Estimates Psi_0^{(h)}(k) over all policy-on / missing target
# blocks in target slice k, for h in {"ATE", "RowHet", "Local", "Trend"}.
#
# This wrapper is compatible with
#
#   bilinearTensorStaggered(Y, k, a, b, r, x, y, A = NULL,
#                           Omega = NULL, tau)
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
  # Infer or validate A and Omega
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
    stop("Internal error: row and column partitions have different lengths.")
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


# ---------------------------------------------------------------
# Wrappers
# ---------------------------------------------------------------

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