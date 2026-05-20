pluginPsi_c1 = function(Y, k, N_parts, T_parts,
                        functional = c("ATE", "RowHet", "Local", "Trend"),
                        eta = NULL,
                        row_index = NULL) {
  
  functional = match.arg(functional)
  
  make_parts = function(sizes) {
    ends = cumsum(sizes)
    starts = c(1, head(ends, -1) + 1)
    Map(seq, starts, ends)
  }
  
  row_parts = make_parts(N_parts[[k]])
  col_parts = make_parts(T_parts[[k]])
  o_k = length(row_parts)
  
  Yk = if (length(dim(Y)) == 3) Y[, , k] else Y
  
  if (functional == "RowHet") stopifnot(!is.null(eta))
  
  if (functional == "Local") {
    stopifnot(!is.null(row_index))
    a0 = which(vapply(row_parts, function(idx) row_index %in% idx, logical(1)))
    stopifnot(length(a0) == 1)
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (a in seq_len(o_k)) {
    for (b in seq_len(o_k)) {
      
      ## c = 1 region: policy-on / complementary staggered region
      if (a + b <= o_k + 1) next
      
      if (functional == "Local" && a != a0) next
      
      rows = row_parts[[a]]
      cols = col_parts[[b]]
      
      Nik = length(rows)
      Tbk = length(cols)
      Y_ab = Yk[rows, cols, drop = FALSE]
      
      if (functional == "ATE") {
        
        val = mean(Y_ab, na.rm = TRUE)
        weight = Nik * Tbk
        
      } else if (functional == "RowHet") {
        
        val = sum(eta[rows] * rowSums(Y_ab, na.rm = TRUE)) / (Nik * Tbk)
        weight = Nik * Tbk
        
      } else if (functional == "Local") {
        
        val = mean(Yk[row_index, cols], na.rm = TRUE)
        weight = Tbk
        
      } else if (functional == "Trend") {
        
        ## Skip one-year column blocks, because the centred trend vector is zero.
        if (Tbk <= 1) {
          next
        }
        
        z = seq_len(Tbk)
        zc = z - mean(z)
        
        y = zc / sqrt(sum(zc^2))
        x = rep(1 / sqrt(Nik), Nik)
        
        mu_ab = as.numeric(t(x) %*% Y_ab %*% y)
        val = mu_ab / sqrt(Nik * Tbk * (Tbk^2 - 1) / 12)
        
        weight = 1
      }
      
      weighted_sum = weighted_sum + weight * val
      normalizer = normalizer + weight
    }
  }
  
  weighted_sum / normalizer
}

## ---------------------------------------------------------------
## Wrappers for c = 1 plug-in quantities
## ---------------------------------------------------------------

pluginPsi1_ATE = function(Y, k, N_parts, T_parts) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "ATE"
  )
}

pluginPsi1_RowHet = function(Y, k, N_parts, T_parts, eta) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "RowHet",
    eta = eta
  )
}

pluginPsi1_Local = function(Y, k, N_parts, T_parts, row_index) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "Local",
    row_index = row_index
  )
}

pluginPsi1_Trend = function(Y, k, N_parts, T_parts) {
  pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = "Trend"
  )
}