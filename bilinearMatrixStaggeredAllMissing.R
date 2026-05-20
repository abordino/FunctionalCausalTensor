# ---------------------------------------------------------------
# bilinearMatrixStaggeredAllMissing
#
# Matrix-only wrapper for staggered adoption missingness.
#
# This estimates the entrywise average over all missing staggered
# blocks:
#
#   [1 / sum_{i+t>o+1} N_i T_t]
#   sum_{i+t>o+1} sqrt(N_i T_t)
#     x_i' M^{(i,t)} y_t
#
# with default choices
#
#   x_i = N_i^(-1/2) 1_{N_i}
#   y_t = T_t^(-1/2) 1_{T_t}
#
# so that the target becomes the average of all missing entries.
#
# Args:
#   Y_mat  : N x T matrix with staggered missingness
#   r      : target rank
#   x      : optional list of row-block vectors
#   y      : optional list of column-block vectors
#   N_part : vector c(N_1, ..., N_o)
#   T_part : vector c(T_1, ..., T_o)
#   tau    : threshold for H^\dagger
#
# Returns:
#   scalar estimate of the entrywise missing-region average
# ---------------------------------------------------------------

bilinearMatrixStaggeredAllMissing = function(Y_mat, r, x = NULL, y = NULL,
                                             N_part, T_part, tau) {
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
  
  if (is.null(x)) {
    x = lapply(N_part, function(Ni) rep(1 / sqrt(Ni), Ni))
  }
  
  if (is.null(y)) {
    y = lapply(T_part, function(Tj) rep(1 / sqrt(Tj), Tj))
  }
  
  if (!is.list(x) || length(x) < o) {
    stop("x must be a list with x[[i0]] for each row block i0.")
  }
  
  if (!is.list(y) || length(y) < o) {
    stop("y must be a list with y[[t0]] for each column block t0.")
  }
  
  Y_list = list(Y_mat)
  N_parts = list(N_part)
  T_parts = list(T_part)
  
  weighted_sum = 0
  normalizer = 0
  
  for (i0 in seq_len(o)) {
    for (t0 in seq_len(o)) {
      
    
      if (i0 + t0 <= o + 1) {
        next
      }
      
      Ni = N_part[i0]
      Tj = T_part[t0]
      
      if (length(x[[i0]]) != Ni) {
        stop(
          paste0(
            "length(x[[", i0, "]]) must equal N_part[", i0, "]."
          )
        )
      }
      
      if (length(y[[t0]]) != Tj) {
        stop(
          paste0(
            "length(y[[", t0, "]]) must equal T_part[", t0, "]."
          )
        )
      }
      
      mu_hat = bilinearTensorStaggered(
        Y       = Y_list,
        k       = 1,
        i0      = i0,
        t0      = t0,
        r       = r,
        x       = x[[i0]],
        y       = y[[t0]],
        N_parts = N_parts,
        T_parts = T_parts,
        tau     = tau
      )
      
      weight = sqrt(Ni * Tj)
      
      weighted_sum = weighted_sum + weight * mu_hat
      normalizer = normalizer + Ni * Tj
    }
  }
  
  if (normalizer == 0) {
    stop("No missing staggered blocks satisfy i0 + t0 > o + 1.")
  }
  
  return(weighted_sum / normalizer)
}