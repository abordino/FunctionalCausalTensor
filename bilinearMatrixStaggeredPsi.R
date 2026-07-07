# ---------------------------------------------------------------
# bilinearMatrixStaggeredPsi
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
#   scalar estimate of Psi_0^{(functional)}
# ---------------------------------------------------------------

bilinearMatrixStaggeredPsi = function(Y_mat, r, tau,
                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                      eta = NULL,
                                      row_index = NULL,
                                      A = NULL,
                                      Omega = NULL) {
  functional = match.arg(functional)
  
  # -------------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------------
  
  if (!is.matrix(Y_mat)) {
    stop("Y_mat must be a matrix.")
  }
  
  N = nrow(Y_mat)
  Tt = ncol(Y_mat)
  
  stopifnot(r >= 1)
  stopifnot(tau > 0)
  
  # -------------------------------------------------------------
  # Infer or validate A and Omega
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


# ---------------------------------------------------------------
# Wrappers
# ---------------------------------------------------------------

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