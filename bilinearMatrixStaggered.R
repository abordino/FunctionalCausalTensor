# ---------------------------------------------------------------
# bilinearMatrixStaggered
#
# Matrix-only wrapper for staggered adoption missingness.
#
# This estimates
#
#   mu_{xy}^{(i0,t0)} = x' M^{(i0,t0)} y
#
# using only one N x T matrix Y_k.
#
# Args:
#   Y_mat  : N x T matrix with staggered missingness
#   i0,t0  : target missing block indices
#   r      : target rank
#   x      : unit vector in R^{N_{i0}}
#   y      : unit vector in R^{T_{t0}}
#   N_part : vector c(N_1, ..., N_o)
#   T_part : vector c(T_1, ..., T_o)
#   tau    : threshold for H^\dagger
#
# Returns:
#   scalar estimate of x' M^{(i0,t0)} y
# ---------------------------------------------------------------

bilinearMatrixStaggered = function(Y_mat, i0, t0, r, x, y,
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
  
  if (i0 < 1 || i0 > o || t0 < 1 || t0 > o) {
    stop("i0 and t0 must be valid staggered block indices.")
  }
  
  if (i0 + t0 <= o + 1) {
    stop("Target block is observed under staggered missingness. Need i0 + t0 > o + 1.")
  }
  
  if (length(x) != N_part[i0]) {
    stop("length(x) must equal N_part[i0].")
  }
  
  if (length(y) != T_part[t0]) {
    stop("length(y) must equal T_part[t0].")
  }
  
  # Treat the matrix as a one-layer tensor/list.
  Y_list = list(Y_mat)
  
  N_parts = list(N_part)
  T_parts = list(T_part)
  
  bilinearTensorStaggered(
    Y       = Y_list,
    k       = 1,
    i0      = i0,
    t0      = t0,
    r       = r,
    x       = x,
    y       = y,
    N_parts = N_parts,
    T_parts = T_parts,
    tau     = tau
  )
}