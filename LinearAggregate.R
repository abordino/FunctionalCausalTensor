# =============================================================================
# Reduced-Anchor Bilinear Estimation for Staggered-Adoption Tensor Data
#
# Inputs:
#   Y      : N x T x K array, or list of K matrices of size N x T
#   k      : target layer index
#   r      : requested rank
#   tau    : eigenvalue threshold for stable inverse
#   A      : optional N x K adoption-time matrix
#   Omega  : optional N x T x K logical observation mask
#
# Output:
#   psi_hat
# =============================================================================


sym_threshold_inverse = function(A, tau) {
  eig = eigen(A, symmetric = TRUE)
  eig$vectors %*%
    diag(1 / pmax(eig$values, tau), nrow = length(eig$values)) %*%
    t(eig$vectors)
}

safe_seq = function(from, to) {
  if (to < from) integer(0) else seq(from, to)
}

take_blocks = function(blocks, ids) {
  ids = ids[ids >= 1 & ids <= length(blocks)]
  
  if (length(ids) == 0) {
    integer(0)
  } else {
    unlist(blocks[ids], use.names = FALSE)
  }
}

omega_slice = function(Omega, rows, cols, layer) {
  if (length(rows) == 0 || length(cols) == 0) {
    return(matrix(FALSE, nrow = length(rows), ncol = length(cols)))
  }
  
  out = Omega[rows, cols, layer, drop = FALSE]
  matrix(as.logical(out), nrow = length(rows), ncol = length(cols))
}

standardize_Y = function(Y) {
  if (is.array(Y)) {
    K = dim(Y)[3]
    return(list(
      Y_list = lapply(seq_len(K), function(j) Y[, , j]),
      N = dim(Y)[1],
      Tt = dim(Y)[2],
      K = K
    ))
  }
  
  if (!is.list(Y)) {
    stop("Y must be either an N x T x K array or a list of K N x T matrices.")
  }
  
  K = length(Y)
  N = nrow(Y[[1]])
  Tt = ncol(Y[[1]])
  
  for (j in seq_len(K)) {
    if (!is.matrix(Y[[j]]) || nrow(Y[[j]]) != N || ncol(Y[[j]]) != Tt) {
      stop("All elements of Y must be N x T matrices.")
    }
  }
  
  list(Y_list = Y, N = N, Tt = Tt, K = K)
}

infer_adoption_times = function(Omega) {
  N = dim(Omega)[1]
  Tt = dim(Omega)[2]
  K = dim(Omega)[3]
  A = matrix(Inf, nrow = N, ncol = K)
  
  for (j in seq_len(K)) {
    for (i in seq_len(N)) {
      observed = Omega[i, , j]
      
      if (any(diff(as.integer(observed)) > 0)) {
        stop("Rows must have staggered-adoption missingness.")
      }
      
      first_unobserved = which(!observed)[1]
      if (!is.na(first_unobserved)) A[i, j] = first_unobserved
    }
  }
  
  A
}

build_omega_from_A = function(A, N, Tt, K) {
  Omega = array(FALSE, dim = c(N, Tt, K))
  
  for (j in seq_len(K)) {
    for (i in seq_len(N)) {
      if (is.infinite(A[i, j])) {
        Omega[i, , j] = TRUE
      } else if (A[i, j] > 1) {
        last_observed = min(Tt, A[i, j] - 1)
        Omega[i, seq_len(last_observed), j] = TRUE
      }
    }
  }
  
  Omega
}

makeTargetStaircaseBlocks = function(Y, k, A = NULL, Omega = NULL) {
  data = standardize_Y(Y)
  
  Y_list = data$Y_list
  N = data$N
  Tt = data$Tt
  K = data$K
  
  stopifnot(k >= 1, k <= K)
  
  if (is.null(A) && is.null(Omega)) {
    Omega = array(FALSE, dim = c(N, Tt, K))
    
    for (j in seq_len(K)) {
      Omega[, , j] = !is.na(Y_list[[j]])
    }
    
    A = infer_adoption_times(Omega)
  }
  
  if (!is.null(A)) {
    A = as.matrix(A)
    
    if (nrow(A) != N || ncol(A) != K) {
      stop("A must be an N x K adoption-time matrix.")
    }
  }
  
  if (is.null(Omega)) {
    Omega = build_omega_from_A(A, N, Tt, K)
  } else {
    if (!is.array(Omega) || !all(dim(Omega) == c(N, Tt, K))) {
      stop("Omega must be an N x T x K logical array.")
    }
    
    Omega = array(as.logical(Omega), dim = c(N, Tt, K))
    
    if (is.null(A)) {
      A = infer_adoption_times(Omega)
    }
  }
  
  target_observed_length = ifelse(
    is.infinite(A[, k]),
    Tt,
    pmin(Tt, A[, k] - 1)
  )
  
  row_order = order(-target_observed_length, seq_len(N))
  
  Y_sorted = lapply(seq_len(K), function(j) {
    Y_list[[j]][row_order, , drop = FALSE]
  })
  
  A_sorted = A[row_order, , drop = FALSE]
  Omega_sorted = Omega[row_order, , , drop = FALSE]
  observed_length_sorted = target_observed_length[row_order]
  
  if (any(diff(observed_length_sorted) > 0)) {
    stop("Target-layer ordering failed.")
  }
  
  row_observed_lengths = unique(observed_length_sorted)
  block_count = length(row_observed_lengths)
  
  if (row_observed_lengths[1] != Tt) {
    stop("Target layer must contain at least one fully observed row block.")
  }
  
  if (tail(row_observed_lengths, 1) <= 0) {
    stop("Target layer must contain a non-empty initial observed time period.")
  }
  
  row_blocks = lapply(row_observed_lengths, function(m) {
    which(observed_length_sorted == m)
  })
  
  ascending_lengths = rev(row_observed_lengths)
  time_block_sizes = diff(c(0, ascending_lengths))
  
  time_block_ends = cumsum(time_block_sizes)
  time_block_starts = c(1, head(time_block_ends, -1) + 1)
  
  col_blocks = Map(seq, time_block_starts, time_block_ends)
  
  if (length(row_blocks) != length(col_blocks)) {
    stop("Row and column partitions have different lengths.")
  }
  
  list(
    Y_list = Y_sorted,
    A = A_sorted,
    Omega = Omega_sorted,
    row_perm = row_order,
    obs_len_k = observed_length_sorted,
    row_parts = row_blocks,
    col_parts = col_blocks,
    N_part = vapply(row_blocks, length, integer(1)),
    T_part = vapply(col_blocks, length, integer(1)),
    o_k = block_count,
    N = N,
    Tt = Tt,
    K = K
  )
}

linear_active_blocks = function(blocks, functional, eta = NULL, row_index = NULL) {
  functional = match.arg(functional, c("ATE", "RowHet", "Local", "Trend"))
  
  row_blocks = blocks$row_parts
  col_blocks = blocks$col_parts
  row_order = blocks$row_perm
  block_count = blocks$o_k
  N = blocks$N
  
  if (functional == "RowHet") {
    if (is.null(eta)) stop("eta must be supplied for functional = 'RowHet'.")
    if (length(eta) != N) stop("eta must have length N.")
    if (!all(eta %in% c(-1, 1))) stop("eta should be a vector in {+1, -1}^N.")
  }
  
  local_row_block = NULL
  local_row_position = NULL
  
  if (functional == "Local") {
    if (is.null(row_index)) {
      stop("row_index must be supplied for functional = 'Local'.")
    }
    
    if (
      length(row_index) != 1 ||
      row_index < 1 ||
      row_index > N ||
      row_index != as.integer(row_index)
    ) {
      stop("row_index must be a single valid original row index.")
    }
    
    sorted_position = match(row_index, row_order)
    
    if (is.na(sorted_position)) {
      stop("Could not locate row_index after target-layer sorting.")
    }
    
    local_row_block = which(vapply(
      row_blocks,
      function(idx) sorted_position %in% idx,
      logical(1)
    ))
    
    if (length(local_row_block) != 1) {
      stop("Could not locate row_index in the row partition.")
    }
    
    local_row_position = match(sorted_position, row_blocks[[local_row_block]])
    
    if (is.na(local_row_position)) {
      stop("Could not compute the local position of row_index.")
    }
  }
  
  active = list()
  
  for (a in seq_len(block_count)) {
    for (b in seq_len(block_count)) {
      if (a + b <= block_count + 1) next
      if (functional == "Local" && a != local_row_block) next
      if (functional == "Trend" && length(col_blocks[[b]]) <= 1) next
      
      active[[length(active) + 1]] = c(a = a, b = b)
    }
  }
  
  if (length(active) == 0) {
    stop("No target blocks were included in the aggregation.")
  }
  
  active_matrix = do.call(rbind, active)
  
  list(
    D = active_matrix,
    A_h = sort(unique(active_matrix[, "a"])),
    B_h = sort(unique(active_matrix[, "b"])),
    a0 = local_row_block,
    local_pos = local_row_position
  )
}

make_x = function(blocks, a, functional, eta = NULL, local_pos = NULL) {
  rows = blocks$row_parts[[a]]
  n_rows = length(rows)
  
  if (functional == "ATE") {
    return(rep(1 / sqrt(n_rows), n_rows))
  }
  
  if (functional == "RowHet") {
    original_rows = blocks$row_perm[rows]
    return(eta[original_rows] / sqrt(n_rows))
  }
  
  if (functional == "Local") {
    x = rep(0, n_rows)
    x[local_pos] = 1
    return(x)
  }
  
  rep(1 / sqrt(n_rows), n_rows)
}

make_y = function(blocks, b, functional) {
  n_cols = length(blocks$col_parts[[b]])
  
  if (functional %in% c("ATE", "RowHet", "Local")) {
    return(rep(1 / sqrt(n_cols), n_cols))
  }
  
  z = seq_len(n_cols)
  z_centered = z - mean(z)
  
  z_centered / sqrt(sum(z_centered^2))
}

block_weight = function(blocks, a, b, functional) {
  n_rows = length(blocks$row_parts[[a]])
  n_cols = length(blocks$col_parts[[b]])
  
  if (functional %in% c("ATE", "RowHet")) {
    return(list(
      weight = sqrt(n_rows * n_cols),
      normalizer_increment = n_rows * n_cols
    ))
  }
  
  if (functional == "Local") {
    return(list(
      weight = sqrt(n_cols),
      normalizer_increment = n_cols
    ))
  }
  
  list(
    weight = 1 / sqrt(n_rows * n_cols * (n_cols^2 - 1) / 12),
    normalizer_increment = 1
  )
}

determine_linear_rank = function(blocks, k, r, active_row_blocks, active_col_blocks) {
  Omega = blocks$Omega
  K = blocks$K
  
  row_blocks = blocks$row_parts
  col_blocks = blocks$col_parts
  block_count = blocks$o_k
  
  first_cols = col_blocks[[1]]
  first_rows = row_blocks[[1]]
  
  r_eff = r
  
  for (a in active_row_blocks) {
    rows_for_left_fit = take_blocks(row_blocks, seq_len(a))
    target_cols_for_left_fit = take_blocks(
      col_blocks,
      safe_seq(1, block_count + 1 - a)
    )
    
    available_cols = 0
    
    for (j in seq_len(K)) {
      if (j == k) {
        available_cols = available_cols + length(target_cols_for_left_fit)
      } else {
        observed = omega_slice(Omega, rows_for_left_fit, first_cols, j)
        available_cols = available_cols + sum(colSums(observed) == length(rows_for_left_fit))
      }
    }
    
    r_eff = min(r_eff, length(rows_for_left_fit), available_cols)
  }
  
  for (b in active_col_blocks) {
    rows_for_upper_fit = take_blocks(
      row_blocks,
      safe_seq(1, block_count + 1 - b)
    )
    
    cols_for_upper_fit = take_blocks(col_blocks, seq_len(b))
    
    available_rows = 0
    
    for (j in seq_len(K)) {
      if (j == k) {
        available_rows = available_rows + length(rows_for_upper_fit)
      } else {
        observed = omega_slice(Omega, first_rows, cols_for_upper_fit, j)
        available_rows = available_rows + sum(rowSums(observed) == length(cols_for_upper_fit))
      }
    }
    
    r_eff = min(r_eff, available_rows, length(cols_for_upper_fit))
  }
  
  if (r_eff < 1) {
    stop("Effective rank is zero.")
  }
  
  r_eff
}

build_left_fit = function(blocks, k, a, r_eff, functional, eta = NULL, local_pos = NULL) {
  Y_list = blocks$Y_list
  Omega = blocks$Omega
  K = blocks$K
  
  block_count = blocks$o_k
  row_blocks = blocks$row_parts
  col_blocks = blocks$col_parts
  
  rows_for_fit = take_blocks(row_blocks, seq_len(a))
  target_rows = row_blocks[[a]]
  
  first_cols = col_blocks[[1]]
  target_layer_cols = take_blocks(
    col_blocks,
    safe_seq(1, block_count + 1 - a)
  )
  
  column_anchors = vector("list", K)
  
  for (j in seq_len(K)) {
    if (j == k) {
      column_anchors[[j]] = target_layer_cols
    } else {
      observed = omega_slice(Omega, rows_for_fit, first_cols, j)
      good_cols = which(colSums(observed) == length(rows_for_fit))
      column_anchors[[j]] = first_cols[good_cols]
    }
  }
  
  if (sum(vapply(column_anchors, length, integer(1))) < r_eff) {
    stop(sprintf("Too few reduced column anchors for left fit at a = %s.", a))
  }
  
  left_mats = lapply(seq_len(K), function(j) {
    cols = column_anchors[[j]]
    if (length(cols) == 0) return(NULL)
    
    Y_list[[j]][rows_for_fit, cols, drop = FALSE]
  })
  
  Y_left = do.call(cbind, Filter(Negate(is.null), left_mats))
  
  if (anyNA(Y_left)) {
    stop(sprintf("Y_left contains NA for a = %s.", a))
  }
  
  svd_left = svd(Y_left, nu = r_eff, nv = 0)
  U_left = svd_left$u[, seq_len(r_eff), drop = FALSE]
  
  target_rows_local = match(target_rows, rows_for_fit)
  
  if (anyNA(target_rows_local)) {
    stop(sprintf("Target rows are not contained in rows_for_fit for a = %s.", a))
  }
  
  x = make_x(blocks, a, functional, eta = eta, local_pos = local_pos)
  alpha = as.vector(crossprod(U_left[target_rows_local, , drop = FALSE], x))
  
  list(
    U_left = U_left,
    alpha = alpha,
    S_a = rows_for_fit,
    R_a_local = target_rows_local,
    ColAnc = column_anchors,
    Y_left_dim = dim(Y_left)
  )
}

build_upper_fit = function(blocks, k, b, r_eff, functional) {
  Y_list = blocks$Y_list
  Omega = blocks$Omega
  K = blocks$K
  
  block_count = blocks$o_k
  row_blocks = blocks$row_parts
  col_blocks = blocks$col_parts
  
  first_rows = row_blocks[[1]]
  
  target_layer_rows = take_blocks(
    row_blocks,
    safe_seq(1, block_count + 1 - b)
  )
  
  cols_for_fit = take_blocks(col_blocks, seq_len(b))
  target_cols = col_blocks[[b]]
  
  row_anchors = vector("list", K)
  
  for (j in seq_len(K)) {
    if (j == k) {
      row_anchors[[j]] = target_layer_rows
    } else {
      observed = omega_slice(Omega, first_rows, cols_for_fit, j)
      good_rows = which(rowSums(observed) == length(cols_for_fit))
      row_anchors[[j]] = first_rows[good_rows]
    }
  }
  
  if (sum(vapply(row_anchors, length, integer(1))) < r_eff) {
    stop(sprintf("Too few reduced row anchors for upper fit at b = %s.", b))
  }
  
  upper_mats = lapply(seq_len(K), function(j) {
    rows = row_anchors[[j]]
    if (length(rows) == 0) return(NULL)
    
    Y_list[[j]][rows, cols_for_fit, drop = FALSE]
  })
  
  Y_upper = do.call(rbind, Filter(Negate(is.null), upper_mats))
  
  if (anyNA(Y_upper)) {
    stop(sprintf("Y_upper contains NA for b = %s.", b))
  }
  
  svd_upper = svd(Y_upper, nu = r_eff, nv = r_eff)
  
  U_upper = svd_upper$u[, seq_len(r_eff), drop = FALSE]
  V_upper = svd_upper$v[, seq_len(r_eff), drop = FALSE]
  singular_values = svd_upper$d[seq_len(r_eff)]
  
  target_layer_start = if (k == 1) {
    0
  } else {
    sum(vapply(row_anchors[seq_len(k - 1)], length, integer(1)))
  }
  
  target_layer_indices = (target_layer_start + 1):(target_layer_start + length(target_layer_rows))
  U_target = U_upper[target_layer_indices, , drop = FALSE]
  
  target_cols_local = match(target_cols, cols_for_fit)
  
  if (anyNA(target_cols_local)) {
    stop(sprintf("Target columns are not contained in cols_for_fit for b = %s.", b))
  }
  
  V_target = V_upper[target_cols_local, , drop = FALSE]
  y = make_y(blocks, b, functional)
  
  projected_y = crossprod(V_target, y)
  weighted_projected_y = as.vector(singular_values) * as.vector(projected_y)
  
  X = as.vector(U_target %*% weighted_projected_y)
  
  list(
    X = X,
    S_plus = target_layer_rows,
    Q_b_all = cols_for_fit,
    C_b_local = target_cols_local,
    RowAnc = row_anchors,
    U_up_k = U_target,
    V_b_hat = V_target,
    D_up = singular_values,
    Y_up_dim = dim(Y_upper)
  )
}

bilinearTensorStaggeredLinearReducedAnchor = function(Y, k, r, tau,
                                                       functional = c("ATE", "RowHet", "Local", "Trend"),
                                                       eta = NULL,
                                                       row_index = NULL,
                                                       A = NULL,
                                                       Omega = NULL,
                                                       return_fit = FALSE) {
  functional = match.arg(functional)
  
  stopifnot(r >= 1)
  stopifnot(tau > 0)
  
  blocks = makeTargetStaircaseBlocks(
    Y = Y,
    k = k,
    A = A,
    Omega = Omega
  )
  
  active = linear_active_blocks(
    blocks = blocks,
    functional = functional,
    eta = eta,
    row_index = row_index
  )
  
  r_eff = determine_linear_rank(
    blocks = blocks,
    k = k,
    r = r,
    active_row_blocks = active$A_h,
    active_col_blocks = active$B_h
  )
  
  left_fits = vector("list", blocks$o_k)
  
  for (a in active$A_h) {
    left_fits[[a]] = build_left_fit(
      blocks = blocks,
      k = k,
      a = a,
      r_eff = r_eff,
      functional = functional,
      eta = eta,
      local_pos = active$local_pos
    )
  }
  
  upper_fits = vector("list", blocks$o_k)
  
  for (b in active$B_h) {
    upper_fits[[b]] = build_upper_fit(
      blocks = blocks,
      k = k,
      b = b,
      r_eff = r_eff,
      functional = functional
    )
  }
  
  weighted_sum = 0
  normalizer = 0
  
  for (ell in seq_len(nrow(active$D))) {
    a = active$D[ell, "a"]
    b = active$D[ell, "b"]
    
    left_fit = left_fits[[a]]
    upper_fit = upper_fits[[b]]
    
    upper_rows_in_left_fit = match(upper_fit$S_plus, left_fit$S_a)
    
    if (anyNA(upper_rows_in_left_fit)) {
      stop(sprintf(
        "Upper rows for b = %s are not contained in left rows for a = %s.",
        b, a
      ))
    }
    
    U_common = left_fit$U_left[upper_rows_in_left_fit, , drop = FALSE]
    
    gram = crossprod(U_common)
    gram_inverse = sym_threshold_inverse(gram, tau)
    
    beta = as.vector(gram_inverse %*% crossprod(U_common, upper_fit$X))
    
    weight = block_weight(blocks, a, b, functional)
    
    weighted_sum = weighted_sum + weight$weight * sum(left_fit$alpha * beta)
    normalizer = normalizer + weight$normalizer_increment
  }
  
  if (normalizer == 0) {
    stop("No target blocks were included in the aggregation.")
  }
  
  psi_hat = weighted_sum / normalizer
  
  if (!return_fit) {
    return(psi_hat)
  }
  
  list(
    psi_hat = psi_hat,
    weighted_sum = weighted_sum,
    normalizer = normalizer,
    r_eff = r_eff,
    blocks = blocks,
    active_blocks = active$D,
    left_fits = left_fits,
    upper_fits = upper_fits,
    functional = functional
  )
}

bilinearTensorStaggeredATELinearReducedAnchor = function(Y, k, r, tau,
                                                          A = NULL,
                                                          Omega = NULL,
                                                          return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "ATE",
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredRowHetLinearReducedAnchor = function(Y, k, r, eta, tau,
                                                             A = NULL,
                                                             Omega = NULL,
                                                             return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "RowHet",
    eta = eta,
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredLocalLinearReducedAnchor = function(Y, k, r, row_index, tau,
                                                            A = NULL,
                                                            Omega = NULL,
                                                            return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "Local",
    row_index = row_index,
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}

bilinearTensorStaggeredTrendLinearReducedAnchor = function(Y, k, r, tau,
                                                            A = NULL,
                                                            Omega = NULL,
                                                            return_fit = FALSE) {
  bilinearTensorStaggeredLinearReducedAnchor(
    Y = Y,
    k = k,
    r = r,
    tau = tau,
    functional = "Trend",
    A = A,
    Omega = Omega,
    return_fit = return_fit
  )
}
