# ---------------------------------------------------------------
# bilinearMatrixStaggeredPsi
#
# Matrix-only Psi wrapper for staggered adoption missingness.
#
# This estimates Psi_0^{(h)} over all missing staggered blocks
# using only one N x T matrix Y_mat.
#
# This assumes bilinearMatrixStaggered(...) is already defined.
#
# Args:
#   Y_mat      : N x T matrix with staggered missingness
#   r          : target rank
#   N_part     : vector c(N_1, ..., N_o)
#   T_part     : vector c(T_1, ..., T_o)
#   tau        : threshold for H^\dagger
#   functional : one of "ATE", "RowHet", "Local", "Trend"
#   eta        : length-N vector in {+1,-1}; required for RowHet
#   row_index  : global row index; required for Local
#
# Returns:
#   scalar estimate of Psi_0^{(functional)}
# ---------------------------------------------------------------

bilinearMatrixStaggeredPsi = function(Y_mat, r,
                                      N_part, T_part, tau,
                                      functional = c("ATE", "RowHet", "Local", "Trend"),
                                      eta = NULL,
                                      row_index = NULL) {
  functional = match.arg(functional)
  
  if (!is.matrix(Y_mat)) {
    stop("Y_mat must be a matrix.")
  }
  
  N = nrow(Y_mat)
  Tt = ncol(Y_mat)
  
  if (sum(N_part) != N) {
    stop("sum(N_part) must equal nrow(Y_mat).")
  }
  
  if (sum(T_part) != Tt) {
    stop("sum(T_part) must equal ncol(Y_mat).")
  }
  
  if (length(N_part) != length(T_part)) {
    stop("N_part and T_part must have the same length.")
  }
  
  o = length(N_part)
  
  make_partition_indices = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  row_parts = make_partition_indices(N_part)
  col_parts = make_partition_indices(T_part)
  
  
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
    
    # Find the unique row block a0 such that row_index is in R_a0
    a0 = which(vapply(row_parts, function(idx) row_index %in% idx, logical(1)))
    
    if (length(a0) != 1) {
      stop("Could not locate row_index in the row partition.")
    }
    
    local_pos = match(row_index, row_parts[[a0]])
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (row_block in seq_len(o)) {
    for (col_block in seq_len(o)) {
      
      # Only missing staggered blocks
      if (row_block + col_block <= o + 1) {
        next
      }
      
      # For Local, only include blocks in the row block containing row_index
      if (functional == "Local" && row_block != a0) {
        next
      }
      
      Nik = N_part[row_block]
      Ttk = T_part[col_block]
      
    
      if (functional == "ATE") {
        
        x = rep(1 / sqrt(Nik), Nik)
        y = rep(1 / sqrt(Ttk), Ttk)
        
        weight = sqrt(Nik * Ttk)
        normalizer_increment = Nik * Ttk
        
      } else if (functional == "RowHet") {
        
        rows_this_block = row_parts[[row_block]]
        
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
        
        # Skip one-period column blocks, because the centered trend vector is zero.
        if (Ttk <= 1) {
          next
        }
        
        z = seq_len(Ttk)
        z_centered = z - mean(z)
        
        x = rep(1 / sqrt(Nik), Nik)
        y = z_centered / sqrt(sum(z_centered^2))
        
        # c_ab^{Trend} = {N_ab T_ab (T_ab^2 - 1) / 12}^{-1/2}
        weight = 1 / sqrt(Nik * Ttk * (Ttk^2 - 1) / 12)
        
        # W_h = number of usable trend blocks
        normalizer_increment = 1
      }
      
    
      mu_hat = bilinearMatrixStaggered(
        Y_mat  = Y_mat,
        i0     = row_block,
        t0     = col_block,
        r      = r,
        x      = x,
        y      = y,
        N_part = N_part,
        T_part = T_part,
        tau    = tau
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


bilinearMatrixStaggeredATE = function(Y_mat, r, N_part, T_part, tau) {
  bilinearMatrixStaggeredPsi(
    Y_mat      = Y_mat,
    r          = r,
    N_part     = N_part,
    T_part     = T_part,
    tau        = tau,
    functional = "ATE"
  )
}


bilinearMatrixStaggeredRowHet = function(Y_mat, r, eta,
                                         N_part, T_part, tau) {
  bilinearMatrixStaggeredPsi(
    Y_mat      = Y_mat,
    r          = r,
    N_part     = N_part,
    T_part     = T_part,
    tau        = tau,
    functional = "RowHet",
    eta        = eta
  )
}


bilinearMatrixStaggeredLocal = function(Y_mat, r, row_index,
                                        N_part, T_part, tau) {
  bilinearMatrixStaggeredPsi(
    Y_mat      = Y_mat,
    r          = r,
    N_part     = N_part,
    T_part     = T_part,
    tau        = tau,
    functional = "Local",
    row_index  = row_index
  )
}


bilinearMatrixStaggeredTrend = function(Y_mat, r, N_part, T_part, tau) {
  bilinearMatrixStaggeredPsi(
    Y_mat      = Y_mat,
    r          = r,
    N_part     = N_part,
    T_part     = T_part,
    tau        = tau,
    functional = "Trend"
  )
}