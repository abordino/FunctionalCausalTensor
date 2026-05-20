# ---------------------------------------------------------------
# bilinearTensorStaggeredPsi
#
# Estimates overall Psi_0^{(h)}(k) over all missing staggered blocks
# for h in {"ATE", "RowHet", "Local", "Trend"}.
#
# Args:
#   Y          : list of K N x T matrices, or N x T x K array
#   k          : target slice
#   r          : target rank
#   N_parts    : list of length K; N_parts[[j]] = c(N_{1j}, ..., N_{o_j j})
#   T_parts    : list of length K; T_parts[[j]] = c(T_{1j}, ..., T_{o_j j})
#   tau        : threshold
#   functional : one of "ATE", "RowHet", "Local", "Trend"
#   eta        : length-N vector in {+1,-1}; required for RowHet
#   row_index  : global row index i0; required for Local
#
# Returns:
#   scalar estimate of Psi_0^{(functional)}(k)
# ---------------------------------------------------------------

bilinearTensorStaggeredPsi = function(Y, k, r,
                                      N_parts, T_parts, tau,
                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                      eta = NULL,
                                      row_index = NULL) {
  functional = match.arg(functional)
  
  if (is.array(Y)) {
    K  = dim(Y)[3]
    N  = dim(Y)[1]
    Tt = dim(Y)[2]
  } else {
    K  = length(Y)
    N  = nrow(Y[[1]])
    Tt = ncol(Y[[1]])
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
  
  o_k = length(N_parts[[k]])
  
  if (length(T_parts[[k]]) != o_k) {
    stop("N_parts[[k]] and T_parts[[k]] must have the same length.")
  }
  
  if (sum(N_parts[[k]]) != N) {
    stop("N_parts[[k]] must sum to N.")
  }
  
  if (sum(T_parts[[k]]) != Tt) {
    stop("T_parts[[k]] must sum to T.")
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
    
    if (length(row_index) != 1 || row_index < 1 || row_index > N) {
      stop("row_index must be a single valid global row index.")
    }
    
    # Find the unique row block a0 such that row_index in R_{a0,k}
    a0 = which(vapply(row_parts[[k]], function(idx) row_index %in% idx, logical(1)))
    
    if (length(a0) != 1) {
      stop("Could not locate row_index in the row partition for slice k.")
    }
    
    local_pos = match(row_index, row_parts[[k]][[a0]])
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (row_block in seq_len(o_k)) {
    for (col_block in seq_len(o_k)) {
      
      # Only missing staggered blocks
      if (row_block + col_block <= o_k + 1) {
        next
      }
      
      # For Local-i0, only include blocks in the row block containing row_index
      if (functional == "Local" && row_block != a0) {
        next
      }
      
      Nik = N_parts[[k]][row_block]
      Ttk = T_parts[[k]][col_block]
      
      if (functional == "ATE") {
        
        x = rep(1 / sqrt(Nik), Nik)
        y = rep(1 / sqrt(Ttk), Ttk)
        
        weight = sqrt(Nik * Ttk)
        normalizer_increment = Nik * Ttk
        
      } else if (functional == "RowHet") {
        
        rows_this_block = row_parts[[k]][[row_block]]
        
        x = eta[rows_this_block] / sqrt(Nik)
        y = rep(1 / sqrt(Ttk), Ttk)
        
        weight = sqrt(Nik * Ttk)
        normalizer_increment = Nik * Ttk
        
      } else if (functional == "Local") {
        
        x = rep(0, Nik)
        x[local_pos] = 1
        
        y = rep(1 / sqrt(Ttk), Ttk)
        
        weight = sqrt(Ttk)
        normalizer_increment = Ttk
        
      } else if (functional == "Trend") {
        
        ## Skip one-year column blocks, because the centred trend vector is zero.
        if (Ttk <= 1) {
          next
        }
        
        z = seq_len(Ttk)
        z_centered = z - mean(z)
        
        x = rep(1 / sqrt(Nik), Nik)
        y = z_centered / sqrt(sum(z_centered^2))
        
        ## c_ab^{Trend} = {N_ab T_ab (T_ab^2 - 1) / 12}^{-1/2}
        weight = 1 / sqrt(Nik * Ttk * (Ttk^2 - 1) / 12)
        
        ## W_h(k) = number of usable trend blocks
        normalizer_increment = 1
      }
      
      mu_hat = bilinearTensorStaggered(
        Y       = Y,
        k       = k,
        i0      = row_block,
        t0      = col_block,
        r       = r,
        x       = x,
        y       = y,
        N_parts = N_parts,
        T_parts = T_parts,
        tau     = tau
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


bilinearTensorStaggeredATE = function(Y, k, r, N_parts, T_parts, tau) {
  bilinearTensorStaggeredPsi(
    Y = Y, k = k, r = r,
    N_parts = N_parts, T_parts = T_parts, tau = tau,
    functional = "ATE"
  )
}

bilinearTensorStaggeredRowHet = function(Y, k, r, eta, N_parts, T_parts, tau) {
  bilinearTensorStaggeredPsi(
    Y = Y, k = k, r = r,
    N_parts = N_parts, T_parts = T_parts, tau = tau,
    functional = "RowHet",
    eta = eta
  )
}

bilinearTensorStaggeredLocal = function(Y, k, r, row_index, N_parts, T_parts, tau) {
  bilinearTensorStaggeredPsi(
    Y = Y, k = k, r = r,
    N_parts = N_parts, T_parts = T_parts, tau = tau,
    functional = "Local",
    row_index = row_index
  )
}

bilinearTensorStaggeredTrend = function(Y, k, r, N_parts, T_parts, tau) {
  bilinearTensorStaggeredPsi(
    Y = Y, k = k, r = r,
    N_parts = N_parts, T_parts = T_parts, tau = tau,
    functional = "Trend"
  )
}