# Choi-Yuan matrix completion for a staggered-adoption panel.
#
# The numerical steps match the MATLAB code in 
# https://www.tandfonline.com/doi/suppl/10.1080/01621459.2024.2380105?scroll=top

cy_as_numeric_matrix = function(x, name) {
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  storage.mode(x) = "double"
  
  if (length(dim(x)) != 2L) {
    stop(name, " must be a two-dimensional numeric matrix.")
  }
  
  x
}

cy_spectral_norm = function(x) {
  x = cy_as_numeric_matrix(x, "x")
  
  if (any(!is.finite(x))) {
    stop("x must contain only finite values.")
  }
  
  if (length(x) == 0L) {
    return(0)
  }
  
  as.numeric(svd(x, nu = 0L, nv = 0L)$d[1L])
}

cy_reconstruct_svd = function(svd_fit, singular_values) {
  singular_values = as.numeric(singular_values)
  keep = which(singular_values != 0)
  
  if (length(keep) == 0L) {
    return(matrix(0, nrow = nrow(svd_fit$u), ncol = nrow(svd_fit$v)))
  }
  
  u = svd_fit$u[, keep, drop = FALSE]
  v = svd_fit$v[, keep, drop = FALSE]
  d = singular_values[keep]
  
  sweep(u, 2L, d, `*`) %*% t(v)
}

cy_truncated_svd = function(x, rank = 2L) {
  x = cy_as_numeric_matrix(x, "x")
  rank = as.integer(rank)
  
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be a positive integer.")
  }
  
  if (rank > min(dim(x))) {
    stop("rank cannot exceed min(nrow(x), ncol(x)).")
  }
  
  fit = svd(x, nu = rank, nv = rank)
  sweep(fit$u[, seq_len(rank), drop = FALSE],
        2L,
        fit$d[seq_len(rank)],
        `*`) %*%
    t(fit$v[, seq_len(rank), drop = FALSE])
}

cy_split_groups = function(index, group_size = 1L) {
  index = as.integer(index)
  group_size = as.integer(group_size)
  
  if (length(group_size) != 1L || is.na(group_size) || group_size < 1L) {
    stop("group_size must be a positive integer.")
  }
  
  if (length(index) == 0L) {
    return(list())
  }
  
  split(index, ceiling(seq_along(index) / group_size))
}

cy_validate_observed_submatrix = function(y_observed, omega) {
  y_observed = cy_as_numeric_matrix(y_observed, "y_observed")
  omega = cy_as_numeric_matrix(omega, "omega")
  
  if (!identical(dim(y_observed), dim(omega))) {
    stop("y_observed and omega must have identical dimensions.")
  }
  
  if (any(!(omega %in% c(0, 1)))) {
    stop("omega must contain only 0 and 1.")
  }
  
  if (any(!is.finite(y_observed))) {
    stop("y_observed must be finite; use zero where omega equals zero.")
  }
  
  if (any(y_observed[omega == 0] != 0)) {
    stop("Entries of y_observed with omega == 0 must equal zero.")
  }
  
  if (sum(omega) == 0) {
    stop("At least one entry must be observed.")
  }
  
  list(y_observed = y_observed, omega = omega)
}

# Lambda selection

cy_matlab_quantile = function(x, probability = 0.95) {
  x = as.numeric(x)
  probability = as.numeric(probability)
  
  if (length(x) == 0L || any(!is.finite(x))) {
    stop("x must be a non-empty finite numeric vector.")
  }
  
  if (length(probability) != 1L || !is.finite(probability) ||
      probability < 0 || probability > 1) {
    stop("probability must lie in [0, 1].")
  }
  
  as.numeric(stats::quantile(
    x,
    probs = probability,
    type = 5,
    names = FALSE
  ))
}

cy_lambda_matlab = function(omega,
                            sigma_squared,
                            monte_carlo_draws = 100L,
                            probability = 0.95,
                            lambda_multiplier = 1.14) {
  omega = cy_as_numeric_matrix(omega, "omega")
  sigma_squared = as.numeric(sigma_squared)
  monte_carlo_draws = as.integer(monte_carlo_draws)
  
  if (any(!(omega %in% c(0, 1)))) {
    stop("omega must contain only 0 and 1.")
  }
  
  if (length(sigma_squared) != 1L || !is.finite(sigma_squared) ||
      sigma_squared <= 0) {
    stop("sigma_squared must be a positive finite scalar.")
  }
  
  if (length(monte_carlo_draws) != 1L || is.na(monte_carlo_draws) ||
      monte_carlo_draws < 1L) {
    stop("monte_carlo_draws must be a positive integer.")
  }
  
  simulated_norms = numeric(monte_carlo_draws)
  standard_deviation = sqrt(sigma_squared)
  
  for (draw in seq_len(monte_carlo_draws)) {
    generated_error = matrix(
      stats::rnorm(length(omega), mean = 0, sd = standard_deviation),
      nrow = nrow(omega),
      ncol = ncol(omega)
    )
    
    simulated_norms[draw] = cy_spectral_norm(generated_error * omega)
  }
  
  noise_quantile = cy_matlab_quantile(simulated_norms, probability)
  
  lambda = (1 / 2) * noise_quantile * lambda_multiplier * 2
  
  list(
    lambda = as.numeric(lambda),
    simulated_norms = simulated_norms,
    noise_quantile = noise_quantile,
    sigma_squared = sigma_squared
  )
}

cy_lambda = function(omega,
                     sigma,
                     monte_carlo_draws = 100L,
                     probability = 0.95,
                     lambda_multiplier = 1.14) {
  result = cy_lambda_matlab(
    omega = omega,
    sigma_squared = as.numeric(sigma)^2,
    monte_carlo_draws = monte_carlo_draws,
    probability = probability,
    lambda_multiplier = lambda_multiplier
  )
  
  result$lambda
}

# Nuclear-norm estimator

cy_nuclear_proximal_matlab = function(y_observed,
                                      omega,
                                      lambda,
                                      tau = 0.5,
                                      tolerance = 1e-7,
                                      max_iterations = Inf) {
  checked = cy_validate_observed_submatrix(y_observed, omega)
  y_observed = checked$y_observed
  omega = checked$omega
  
  lambda = as.numeric(lambda)
  tau = as.numeric(tau)
  tolerance = as.numeric(tolerance)
  max_iterations = as.numeric(max_iterations)
  
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a non-negative finite scalar.")
  }
  
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be a positive finite scalar.")
  }
  
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a non-negative finite scalar.")
  }
  
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      max_iterations <= 0) {
    stop("max_iterations must be positive or Inf.")
  }
  
  estimate = y_observed
  difference = 1
  iteration = 0L
  
  while (difference > tolerance) {
    iteration = iteration + 1L
    
    if (is.finite(max_iterations) && iteration > max_iterations) {
      stop("Nuclear-norm iteration did not converge within max_iterations.")
    }
    
    old = estimate
    gradient_step = old - tau * (omega * (old - y_observed))
    
    fit = svd(
      gradient_step,
      nu = min(dim(gradient_step)),
      nv = min(dim(gradient_step))
    )
    
    threshold = lambda * tau
    thresholded_values = (fit$d - threshold) * (fit$d > threshold)
    estimate = cy_reconstruct_svd(fit, thresholded_values)
    
    difference = sum((estimate - old)^2) /
      (nrow(y_observed) * ncol(y_observed))
  }
  
  list(
    estimate = estimate,
    iterations = iteration,
    convergence_measure = difference
  )
}

# Error variance

cy_estimate_sigma_squared_matlab = function(
    y_observed,
    omega,
    tau = 0.5,
    nuclear_tolerance = 1e-7,
    sigma_tolerance = 1e-2,
    monte_carlo_draws = 100L,
    lambda_probability = 0.95,
    lambda_multiplier = 1.14,
    max_nuclear_iterations = Inf,
    max_sigma_iterations = Inf) {
  checked = cy_validate_observed_submatrix(y_observed, omega)
  y_observed = checked$y_observed
  omega = checked$omega
  
  sigma_tolerance = as.numeric(sigma_tolerance)
  max_sigma_iterations = as.numeric(max_sigma_iterations)
  
  if (length(sigma_tolerance) != 1L || !is.finite(sigma_tolerance) ||
      sigma_tolerance < 0) {
    stop("sigma_tolerance must be a non-negative finite scalar.")
  }
  
  if (length(max_sigma_iterations) != 1L || is.na(max_sigma_iterations) ||
      max_sigma_iterations <= 0) {
    stop("max_sigma_iterations must be positive or Inf.")
  }
  
  n_rows = nrow(y_observed)
  n_columns = ncol(y_observed)
  
  initial_mean_matrix = matrix(0, nrow = n_rows, ncol = n_columns)
  
  for (column in seq_len(n_columns)) {
    observed_count = sum(omega[, column])
    
    if (observed_count <= 0) {
      stop("Every submatrix column must contain at least one observed entry.")
    }
    
    initial_mean_matrix[, column] =
      sum(y_observed[, column]) / observed_count
  }
  
  initial_residual = y_observed - omega * initial_mean_matrix
  initial_sigma_squared = sum(initial_residual^2) / sum(omega)
  
  if (!is.finite(initial_sigma_squared) || initial_sigma_squared <= 0) {
    stop("The initial sigma^2 estimate is not positive.")
  }
  
  estimated_sigma_squared = 1 / initial_sigma_squared
  sigma_difference = 1
  sigma_iteration = 0L
  last_weighted_fit = NULL
  last_lambda_initial = NA_real_
  last_nuclear_iterations = NA_integer_
  
  while (sigma_difference > sigma_tolerance) {
    sigma_iteration = sigma_iteration + 1L
    
    if (is.finite(max_sigma_iterations) &&
        sigma_iteration > max_sigma_iterations) {
      stop("The iteration did not converge.")
    }
    
    old_sigma_squared = estimated_sigma_squared
    
    lambda_draw = cy_lambda_matlab(
      omega = omega,
      sigma_squared = old_sigma_squared,
      monte_carlo_draws = monte_carlo_draws,
      probability = lambda_probability,
      lambda_multiplier = lambda_multiplier
    )
    last_lambda_initial = lambda_draw$lambda
    
    weighted_fit = cy_nuclear_proximal_matlab(
      y_observed = y_observed,
      omega = omega,
      lambda = last_lambda_initial,
      tau = tau,
      tolerance = nuclear_tolerance,
      max_iterations = max_nuclear_iterations
    )
    
    last_weighted_fit = weighted_fit$estimate
    last_nuclear_iterations = weighted_fit$iterations
    
    estimated_residual = y_observed - omega * last_weighted_fit
    estimated_sigma_squared = sum(estimated_residual^2) / sum(omega)
    
    if (!is.finite(estimated_sigma_squared) || estimated_sigma_squared <= 0) {
      stop("The iteration produced a non-positive value.")
    }
    
    sigma_difference = abs(old_sigma_squared - estimated_sigma_squared)
  }
  
  list(
    sigma_squared = estimated_sigma_squared,
    sigma = sqrt(estimated_sigma_squared),
    sigma_iterations = sigma_iteration,
    sigma_convergence_measure = sigma_difference,
    initial_sigma_squared = initial_sigma_squared,
    last_initial_lambda = last_lambda_initial,
    last_weighted_fit = last_weighted_fit,
    last_weighted_fit_iterations = last_nuclear_iterations
  )
}

# Complete one small submatrix

cy_complete_submatrix = function(y_observed,
                                 omega,
                                 rank = 2L,
                                 sigma = NULL,
                                 tau = 0.5,
                                 tolerance = 1e-7,
                                 sigma_tolerance = 1e-2,
                                 max_iterations = Inf,
                                 max_sigma_iterations = Inf,
                                 monte_carlo_draws = 100L,
                                 lambda_probability = 0.95,
                                 lambda_multiplier = 1.14) {
  checked = cy_validate_observed_submatrix(y_observed, omega)
  y_observed = checked$y_observed
  omega = checked$omega
  
  rank = as.integer(rank)
  
  if (length(rank) != 1L || is.na(rank) || rank < 1L ||
      rank > min(dim(y_observed))) {
    stop("rank must be between 1 and min(nrow(y_observed), ncol(y_observed)).")
  }
  
  sigma_estimation = NULL
  
  if (is.null(sigma)) {
    sigma_estimation = cy_estimate_sigma_squared_matlab(
      y_observed = y_observed,
      omega = omega,
      tau = tau,
      nuclear_tolerance = tolerance,
      sigma_tolerance = sigma_tolerance,
      monte_carlo_draws = monte_carlo_draws,
      lambda_probability = lambda_probability,
      lambda_multiplier = lambda_multiplier,
      max_nuclear_iterations = max_iterations,
      max_sigma_iterations = max_sigma_iterations
    )
    sigma_squared = sigma_estimation$sigma_squared
  } else {
    sigma = as.numeric(sigma)
    
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
      stop("sigma must be NULL or a positive finite scalar.")
    }
    
    sigma_squared = sigma^2
  }
  
  final_lambda_draw = cy_lambda_matlab(
    omega = omega,
    sigma_squared = sigma_squared,
    monte_carlo_draws = monte_carlo_draws,
    probability = lambda_probability,
    lambda_multiplier = lambda_multiplier
  )
  
  nuclear_fit = cy_nuclear_proximal_matlab(
    y_observed = y_observed,
    omega = omega,
    lambda = final_lambda_draw$lambda,
    tau = tau,
    tolerance = tolerance,
    max_iterations = max_iterations
  )
  
  m_nuclear = nuclear_fit$estimate
  
  filled_for_debiasing = (1 - omega) * m_nuclear + y_observed
  m_debiased = cy_truncated_svd(filled_for_debiasing, rank = rank)
  
  list(
    nuclear = m_nuclear,
    filled_for_debiasing = filled_for_debiasing,
    debiased = m_debiased,
    lambda = final_lambda_draw$lambda,
    lambda_simulated_norms = final_lambda_draw$simulated_norms,
    sigma = sqrt(sigma_squared),
    sigma_squared = sigma_squared,
    sigma_was_estimated = is.null(sigma),
    sigma_estimation = sigma_estimation,
    iterations = nuclear_fit$iterations,
    convergence_measure = nuclear_fit$convergence_measure,
    settings = list(
      rank = rank,
      tau = tau,
      tolerance = tolerance,
      sigma_tolerance = sigma_tolerance,
      monte_carlo_draws = monte_carlo_draws,
      lambda_probability = lambda_probability,
      lambda_multiplier = lambda_multiplier
    )
  )
}

# Whole staggered panel completion

cy_build_matlab_tasks = function(omega,
                                 adoption_time,
                                 time,
                                 subgroup_size = 1L,
                                 processing_order = c("matlab", "cohort")) {
  processing_order = match.arg(processing_order)
  subgroup_size = as.integer(subgroup_size)
  
  finite_adoption_times = unique(adoption_time[is.finite(adoption_time)])
  
  if (length(finite_adoption_times) == 0L) {
    return(list())
  }
  
  tasks = list()
  
  if (processing_order == "matlab") {
    block_starts = sort(finite_adoption_times, decreasing = TRUE)
    
    for (block_index in seq_along(block_starts)) {
      block_start = block_starts[block_index]
      upper_bound = if (block_index == 1L) {
        Inf
      } else {
        block_starts[block_index - 1L]
      }
      
      target_columns = which(time >= block_start & time < upper_bound)
      
      if (length(target_columns) == 0L) {
        next
      }
      
      eligible_rows = which(
        is.finite(adoption_time) & adoption_time <= block_start
      )
      
      cohort_order = unique(adoption_time[eligible_rows])
      row_groups = list()
      row_group_cohort = numeric(0)
      
      for (cohort_time in cohort_order) {
        cohort_rows = eligible_rows[adoption_time[eligible_rows] == cohort_time]
        cohort_groups = cy_split_groups(cohort_rows, subgroup_size)
        
        for (group in cohort_groups) {
          row_groups[[length(row_groups) + 1L]] = group
          row_group_cohort = c(row_group_cohort, cohort_time)
        }
      }
      
      if (length(row_groups) > 1L) {
        group_order = order(vapply(row_groups, min, integer(1L)))
        row_groups = row_groups[group_order]
        row_group_cohort = row_group_cohort[group_order]
      }
      
      for (group_index in seq_along(row_groups)) {
        target_rows = row_groups[[group_index]]
        cohort_time = row_group_cohort[group_index]
        
        for (target_column in target_columns) {
          if (all(omega[target_rows, target_column] == 0)) {
            tasks[[length(tasks) + 1L]] = list(
              cohort_time = cohort_time,
              target_rows = as.integer(target_rows),
              target_column = as.integer(target_column)
            )
          }
        }
      }
    }
  } else {
    cohort_order = sort(finite_adoption_times)
    
    for (cohort_time in cohort_order) {
      cohort_rows = which(adoption_time == cohort_time)
      target_columns = which(time >= cohort_time)
      row_groups = cy_split_groups(cohort_rows, subgroup_size)
      
      for (target_column in target_columns) {
        for (target_rows in row_groups) {
          missing_rows = target_rows[omega[target_rows, target_column] == 0]
          
          if (length(missing_rows) > 0L) {
            tasks[[length(tasks) + 1L]] = list(
              cohort_time = cohort_time,
              target_rows = as.integer(missing_rows),
              target_column = as.integer(target_column)
            )
          }
        }
      }
    }
  }
  
  tasks
}

cy_complete_staggered_layer = function(y,
                                       omega,
                                       adoption_time,
                                       time,
                                       rank = 2L,
                                       sigma = NULL,
                                       subgroup_size = 1L,
                                       tau = 0.5,
                                       tolerance = 1e-7,
                                       sigma_tolerance = 1e-2,
                                       max_iterations = Inf,
                                       max_sigma_iterations = Inf,
                                       monte_carlo_draws = 100L,
                                       lambda_probability = 0.95,
                                       lambda_multiplier = 1.14,
                                       processing_order = c("matlab", "cohort"),
                                       seed = NULL,
                                       verbose = TRUE) {
  y = cy_as_numeric_matrix(y, "y")
  omega = cy_as_numeric_matrix(omega, "omega")
  processing_order = match.arg(processing_order)
  
  if (!identical(dim(y), dim(omega))) {
    stop("y and omega must have identical dimensions.")
  }
  
  if (any(!(omega %in% c(0, 1)))) {
    stop("omega must contain only 0 and 1.")
  }
  
  adoption_time = as.numeric(adoption_time)
  time = as.numeric(time)
  
  if (length(adoption_time) != nrow(y)) {
    stop("adoption_time must have one value per row of y.")
  }
  
  if (length(time) != ncol(y)) {
    stop("time must have one value per column of y.")
  }
  
  if (any(!is.finite(y[omega == 1]))) {
    stop("Observed entries of y must be finite.")
  }
  
  rank = as.integer(rank)
  subgroup_size = as.integer(subgroup_size)
  
  if (length(rank) != 1L || is.na(rank) || rank < 1L) {
    stop("rank must be a positive integer.")
  }
  
  if (length(subgroup_size) != 1L || is.na(subgroup_size) ||
      subgroup_size < 1L) {
    stop("subgroup_size must be a positive integer.")
  }
  
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  
  completed_debiased = y
  completed_nuclear = y
  completed_debiased[omega == 0] = NA_real_
  completed_nuclear[omega == 0] = NA_real_
  
  tasks = cy_build_matlab_tasks(
    omega = omega,
    adoption_time = adoption_time,
    time = time,
    subgroup_size = subgroup_size,
    processing_order = processing_order
  )
  
  diagnostics = vector("list", length(tasks))
  
  for (task_index in seq_along(tasks)) {
    task = tasks[[task_index]]
    cohort_time = task$cohort_time
    target_rows = task$target_rows
    target_column = task$target_column
    
    pre_treatment_columns = which(time < cohort_time)
    
    if (length(pre_treatment_columns) < rank) {
      stop(
        "Cohort ", cohort_time,
        " has fewer than rank pre-treatment periods."
      )
    }
    
    anchor_rows = which(omega[, target_column] == 1)
    
    if (length(anchor_rows) < rank) {
      stop(
        "Target period ", time[target_column],
        " has fewer than rank untreated anchor rows."
      )
    }
    
    submatrix_rows = c(anchor_rows, target_rows)
    submatrix_columns = c(pre_treatment_columns, target_column)
    
    omega_sub = omega[
      submatrix_rows,
      submatrix_columns,
      drop = FALSE
    ]
    
    y_sub = y[
      submatrix_rows,
      submatrix_columns,
      drop = FALSE
    ]
    y_sub[omega_sub == 0] = 0
    
    n_anchor = length(anchor_rows)
    target_positions = n_anchor + seq_along(target_rows)
    target_column_position = ncol(y_sub)
    
    if (any(omega_sub[target_positions, target_column_position] != 0)) {
      stop("The requested target block is not missing.")
    }
    
    if (any(omega_sub[seq_len(n_anchor), target_column_position] != 1)) {
      stop("An anchor-row target-period entry is not observed.")
    }
    
    expected_missing = matrix(FALSE, nrow = nrow(omega_sub), ncol = ncol(omega_sub))
    expected_missing[target_positions, target_column_position] = TRUE
    
    if (any((omega_sub == 0) != expected_missing)) {
      stop(
        "The selected panel is not compatible with the staggered-adoption ",
        "submatrix construction."
      )
    }
    
    fit = cy_complete_submatrix(
      y_observed = y_sub,
      omega = omega_sub,
      rank = rank,
      sigma = sigma,
      tau = tau,
      tolerance = tolerance,
      sigma_tolerance = sigma_tolerance,
      max_iterations = max_iterations,
      max_sigma_iterations = max_sigma_iterations,
      monte_carlo_draws = monte_carlo_draws,
      lambda_probability = lambda_probability,
      lambda_multiplier = lambda_multiplier
    )
    
    completed_debiased[target_rows, target_column] = fit$debiased[
      target_positions,
      target_column_position
    ]
    
    completed_nuclear[target_rows, target_column] = fit$nuclear[
      target_positions,
      target_column_position
    ]
    
    diagnostics[[task_index]] = data.frame(
      task = task_index,
      cohort_time = cohort_time,
      target_time = time[target_column],
      target_rows = paste(target_rows, collapse = ","),
      subgroup_size = length(target_rows),
      anchor_rows = length(anchor_rows),
      pre_treatment_periods = length(pre_treatment_columns),
      lambda = fit$lambda,
      sigma = fit$sigma,
      sigma_squared = fit$sigma_squared,
      sigma_was_estimated = fit$sigma_was_estimated,
      sigma_iterations = if (is.null(fit$sigma_estimation)) {
        NA_integer_
      } else {
        fit$sigma_estimation$sigma_iterations
      },
      iterations = fit$iterations,
      convergence_measure = fit$convergence_measure,
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(verbose)) {
      message(
        "task=", task_index, "/", length(tasks),
        ", cohort=", cohort_time,
        ", period=", time[target_column],
        ", rows=", paste(target_rows, collapse = ","),
        ", anchors=", length(anchor_rows),
        ", sigma2=", signif(fit$sigma_squared, 6),
        ", lambda=", signif(fit$lambda, 6)
      )
    }
  }
  
  if (anyNA(completed_debiased[omega == 0])) {
    stop("At least one target entry was not reconstructed by the debiased estimator.")
  }
  
  if (anyNA(completed_nuclear[omega == 0])) {
    stop("At least one target entry was not reconstructed by the nuclear estimator.")
  }
  
  completed_debiased[omega == 1] = y[omega == 1]
  completed_nuclear[omega == 1] = y[omega == 1]
  
  counterfactual_debiased = matrix(
    NA_real_,
    nrow = nrow(y),
    ncol = ncol(y),
    dimnames = dimnames(y)
  )
  counterfactual_nuclear = counterfactual_debiased
  
  counterfactual_debiased[omega == 0] = completed_debiased[omega == 0]
  counterfactual_nuclear[omega == 0] = completed_nuclear[omega == 0]
  
  diagnostic_table = if (length(diagnostics) == 0L) {
    data.frame()
  } else {
    do.call(rbind, diagnostics)
  }
  
  list(
    completed_debiased = completed_debiased,
    completed_nuclear = completed_nuclear,
    counterfactual_debiased = counterfactual_debiased,
    counterfactual_nuclear = counterfactual_nuclear,
    
    completed = completed_debiased,
    counterfactual = counterfactual_debiased,
    
    observed = y,
    omega = omega,
    adoption_time = adoption_time,
    time = time,
    diagnostics = diagnostic_table,
    settings = list(
      rank = rank,
      sigma = sigma,
      subgroup_size = subgroup_size,
      tau = tau,
      tolerance = tolerance,
      sigma_tolerance = sigma_tolerance,
      monte_carlo_draws = monte_carlo_draws,
      lambda_probability = lambda_probability,
      lambda_multiplier = lambda_multiplier,
      processing_order = processing_order,
      seed = seed
    )
  )
}

# Data helpers

cy_first_existing_name = function(data, candidates, label) {
  found = candidates[candidates %in% names(data)]
  
  if (length(found) == 0L) {
    stop(
      "Could not identify ", label, ". Tried: ",
      paste(candidates, collapse = ", "), "."
    )
  }
  
  found[1L]
}

cy_panel_matrix = function(data,
                           value_column,
                           unit_column,
                           time_column,
                           unit_order,
                           time_order) {
  result = matrix(
    NA_real_,
    nrow = length(unit_order),
    ncol = length(time_order),
    dimnames = list(unit_order, as.character(time_order))
  )
  
  row_index = match(as.character(data[[unit_column]]), unit_order)
  column_index = match(data[[time_column]], time_order)
  valid = !is.na(row_index) & !is.na(column_index)
  
  cell_id = paste(row_index[valid], column_index[valid], sep = ":")
  if (anyDuplicated(cell_id)) {
    stop("The data contain duplicate unit-time observations.")
  }
  
  result[cbind(row_index[valid], column_index[valid])] =
    as.numeric(data[[value_column]][valid])
  
  result
}
