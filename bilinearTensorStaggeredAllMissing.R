bilinearTensorStaggeredAllMissing = function(Y, k, r, x = NULL, y = NULL,
                                             N_parts, T_parts, tau) {
  if (is.array(Y)) {
    K = dim(Y)[3]
  } else {
    K = length(Y)
  }
  
  stopifnot(length(N_parts) == K, length(T_parts) == K)
  stopifnot(k >= 1, k <= K)
  
  o_k = length(N_parts[[k]])
  
  if (length(T_parts[[k]]) != o_k) {
    stop("N_parts[[k]] and T_parts[[k]] must have the same length.")
  }
  
  # Default choices:
  # x_i = N_ik^(-1/2) 1_{N_ik}
  # y_t = T_tk^(-1/2) 1_{T_tk}
  if (is.null(x)) {
    x = lapply(N_parts[[k]], function(Nik) rep(1 / sqrt(Nik), Nik))
  }
  
  if (is.null(y)) {
    y = lapply(T_parts[[k]], function(Ttk) rep(1 / sqrt(Ttk), Ttk))
  }
  
  if (!is.list(x) || length(x) < o_k) {
    stop("x must be a list with x[[i0]] for each row block i0.")
  }
  
  if (!is.list(y) || length(y) < o_k) {
    stop("y must be a list with y[[t0]] for each column block t0.")
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (i0 in seq_len(o_k)) {
    for (t0 in seq_len(o_k)) {
      
      # Only estimate missing staggered blocks: i + t > o_k + 1
      if (i0 + t0 <= o_k + 1) {
        next
      }
      
      Nik = N_parts[[k]][i0]
      Ttk = T_parts[[k]][t0]
      
      if (length(x[[i0]]) != Nik) {
        stop(
          paste0(
            "length(x[[", i0, "]]) must equal N_parts[[k]][", i0, "]."
          )
        )
      }
      
      if (length(y[[t0]]) != Ttk) {
        stop(
          paste0(
            "length(y[[", t0, "]]) must equal T_parts[[k]][", t0, "]."
          )
        )
      }
      
      mu_hat = bilinearTensorStaggered(
        Y       = Y,
        k       = k,
        i0      = i0,
        t0      = t0,
        r       = r,
        x       = x[[i0]],
        y       = y[[t0]],
        N_parts = N_parts,
        T_parts = T_parts,
        tau     = tau
      )
      
      weight = sqrt(Nik * Ttk)
      
      weighted_sum = weighted_sum + weight * mu_hat
      normalizer = normalizer + Nik * Ttk
    }
  }
  
  if (normalizer == 0) {
    stop("No missing staggered blocks satisfy i0 + t0 > o_k + 1.")
  }
  
  return(weighted_sum / normalizer)
}