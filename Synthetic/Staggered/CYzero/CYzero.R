# Synthetic noiseless staggered-panel diagnostic for CYNuclear.R

setwd("~/Desktop/code")

source("CYNuclear.R")

make_noiseless_staggered_panel = function() {
  number_units = 30L
  number_periods = 20L
  true_rank = 2L
  
  unit = seq_len(number_units)
  time = seq_len(number_periods)
  
  unit_factor_1 = 10 * (1 + 0.20 * sin(2 * pi * unit / number_units))
  unit_factor_2 = 0.50 * cos(2 * pi * unit / number_units)
  
  time_factor_1 = 1 + 0.10 * time / number_periods
  time_factor_2 = sin(2 * pi * time / number_periods)
  
  truth = tcrossprod(unit_factor_1, time_factor_1) +
    tcrossprod(unit_factor_2, time_factor_2)
  
  rownames(truth) = paste0("unit_", unit)
  colnames(truth) = paste0("time_", time)
  
  adoption_time = c(
    rep(Inf, 15L),
    rep(11, 5L),
    rep(15, 5L),
    rep(18, 5L)
  )
  
  omega = outer(
    adoption_time,
    time,
    FUN = function(first_treated_period, current_period) {
      as.numeric(current_period < first_treated_period)
    }
  )
  
  storage.mode(omega) = "double"
  dimnames(omega) = dimnames(truth)
  
  observed = truth
  observed[omega == 0] = 0
  
  singular_values = svd(truth, nu = 0L, nv = 0L)$d
  rank_tolerance = max(dim(truth)) * max(singular_values) * .Machine$double.eps
  numerical_rank = sum(singular_values > rank_tolerance)
  
  if (numerical_rank != true_rank) {
    stop(
      "Internal test error: generated matrix has numerical rank ",
      numerical_rank,
      ", not ",
      true_rank,
      "."
    )
  }
  
  list(
    truth = truth,
    observed = observed,
    omega = omega,
    adoption_time = adoption_time,
    time = time,
    true_rank = true_rank,
    singular_values = singular_values
  )
}

run_cy_noiseless_staggered_test = function(
    seed = 20260803L,
    verbose = FALSE) {
  source(method_path, local = .GlobalEnv)
  
  if (!exists("cy_complete_staggered_layer", mode = "function")) {
    stop("CYNuclear.R did not define cy_complete_staggered_layer().")
  }
  
  panel = make_noiseless_staggered_panel()
  missing = panel$omega == 0
  observed = panel$omega == 1
  
  cat("\n============================================================\n")
  cat("CY noiseless staggered-panel diagnostic\n")
  cat("============================================================\n")
  cat("Method file:       ", method_path, "\n", sep = "")
  cat("Seed:              ", seed, "\n", sep = "")
  cat("Matrix dimensions: ", nrow(panel$truth), " x ", ncol(panel$truth), "\n", sep = "")
  cat("Numerical rank:    ", panel$true_rank, "\n", sep = "")
  cat("Observed cells:    ", sum(observed), "\n", sep = "")
  cat("Missing cells:     ", sum(missing), "\n", sep = "")
  
  cohort_label = ifelse(
    is.finite(panel$adoption_time),
    as.character(panel$adoption_time),
    "never"
  )
  
  cohort_table = as.data.frame(table(cohort_label), stringsAsFactors = FALSE)
  names(cohort_table) = c("first_treated_period", "number_units")
  
  cat("\nAdoption cohorts:\n")
  print(cohort_table, row.names = FALSE)
  
  cat("\nRunning cy_complete_staggered_layer()...\n")
  
  fit = cy_complete_staggered_layer(
    y = panel$observed,
    omega = panel$omega,
    adoption_time = panel$adoption_time,
    time = panel$time,
    rank = panel$true_rank,
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
    processing_order = "matlab",
    seed = seed,
    verbose = verbose
  )
  
  truth_missing = panel$truth[missing]
  debiased_missing = fit$completed_debiased[missing]
  nuclear_missing = fit$completed_nuclear[missing]
  
  zero_fill_rmse = sqrt(mean(truth_missing^2))
  debiased_rmse = sqrt(mean((debiased_missing - truth_missing)^2))
  nuclear_rmse = sqrt(mean((nuclear_missing - truth_missing)^2))
  
  debiased_max_error = max(abs(debiased_missing - truth_missing))
  nuclear_max_error = max(abs(nuclear_missing - truth_missing))
  
  observed_debiased_error = max(abs(
    fit$completed_debiased[observed] - panel$truth[observed]
  ))
  observed_nuclear_error = max(abs(
    fit$completed_nuclear[observed] - panel$truth[observed]
  ))
  
  metrics = data.frame(
    estimator = c("Zero-fill baseline", "CY nuclear", "CY debiased"),
    RMSE = c(zero_fill_rmse, nuclear_rmse, debiased_rmse),
    relative_RMSE = c(
      1,
      nuclear_rmse / zero_fill_rmse,
      debiased_rmse / zero_fill_rmse
    ),
    maximum_absolute_error = c(
      max(abs(truth_missing)),
      nuclear_max_error,
      debiased_max_error
    ),
    stringsAsFactors = FALSE
  )
  
  cat("\nCompletion metrics on genuinely missing cells:\n")
  print(metrics, row.names = FALSE, digits = 7)
  
  cat("\nAssembly checks:\n")
  cat("Number of local CY tasks:     ", nrow(fit$diagnostics), "\n", sep = "")
  cat("Number of missing cells:      ", sum(missing), "\n", sep = "")
  cat("One task per missing cell:    ",
      identical(as.integer(nrow(fit$diagnostics)), as.integer(sum(missing))),
      "\n", sep = "")
  cat("Max observed-cell error, CY:  ",
      format(observed_debiased_error, scientific = TRUE), "\n", sep = "")
  cat("Max observed-cell error, NNM: ",
      format(observed_nuclear_error, scientific = TRUE), "\n", sep = "")
  
  missing_index = which(missing, arr.ind = TRUE)
  number_to_print = min(12L, nrow(missing_index))
  selected = seq_len(number_to_print)
  
  comparison = data.frame(
    unit = rownames(panel$truth)[missing_index[selected, "row"]],
    period = panel$time[missing_index[selected, "col"]],
    truth = truth_missing[selected],
    CY_debiased = debiased_missing[selected],
    error = debiased_missing[selected] - truth_missing[selected],
    stringsAsFactors = FALSE
  )
  
  cat("\nFirst reconstructed missing cells:\n")
  print(comparison, row.names = FALSE, digits = 7)
  
  cat("\nTuning and convergence diagnostics:\n")
  print(summary(
    fit$diagnostics[
      c("lambda", "sigma_squared", "sigma_iterations", "iterations")
    ]
  ))
  
  relative_debiased_rmse = debiased_rmse / zero_fill_rmse
  task_check = nrow(fit$diagnostics) == sum(missing)
  observed_check = observed_debiased_error <= 100 * .Machine$double.eps
  
  status = if (
    is.finite(relative_debiased_rmse) &&
    relative_debiased_rmse < 0.10 &&
    task_check &&
    observed_check
  ) {
    "PASS"
  } else {
    "CHECK"
  }
  
  cat("\nDiagnostic status: ", status, "\n", sep = "")
  cat(
    "Criterion: debiased relative RMSE < 0.10, one local task per ",
    "missing cell, and unchanged observed cells.\n",
    sep = ""
  )
  cat("============================================================\n\n")
  
  invisible(NULL)
}

run_cy_noiseless_staggered_test()
