setwd("~/Desktop/code")

source("bilinearTensorAllFunction.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
})

# Configuration -----------------------------------------------------------

tensor_file = "Real/CovidOx/data/Omega_Y_until_2020-04-05_delay_28.rds"
results_dir = "Results/oxford_deaths_target_layer_bootstrap"

n_outcomes_to_keep = 2
target_outcome_pattern = "death"
r_grid = 3
tau = 1e-4
functionals = c("ATE", "Trend")

B = 500
bootstrap_seed = 123
stratify_target_bootstrap = TRUE

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Data preparation --------------------------------------------------------

`%||%` = function(x, y) if (is.null(x)) y else x

layer_names = function(x, n, prefix) {
  if (is.null(x)) paste0(prefix, "_", seq_len(n)) else as.character(x[seq_len(n)])
}

first_true = function(x) {
  index = which(x)
  if (length(index)) index[1] else Inf
}

fill_time_neighbours = function(Y, max_passes = 5) {
  for (pass in seq_len(max_passes)) {
    missing = which(is.na(Y), arr.ind = TRUE)
    if (!nrow(missing)) break
    
    missing_before = nrow(missing)
    
    for (row in seq_len(nrow(missing))) {
      i = missing[row, 1]
      t = missing[row, 2]
      k = missing[row, 3]
      
      neighbours = c(
        if (t > 1) Y[i, t - 1, k] else NA_real_,
        if (t < dim(Y)[2]) Y[i, t + 1, k] else NA_real_
      )
      neighbours = neighbours[!is.na(neighbours)]
      
      if (length(neighbours)) Y[i, t, k] = mean(neighbours)
    }
    
    if (sum(is.na(Y)) == missing_before) break
  }
  
  Y
}

align_layers = function(data, target_pattern, n_keep) {
  Y = data$Y
  policy = data$Omega
  
  storage.mode(Y) = "numeric"
  storage.mode(policy) = "logical"
  
  countries = data$country_levels %||% paste0("country_", seq_len(dim(Y)[1]))
  times = as.character(data$date_names %||% seq_len(dim(Y)[2]))
  outcomes = layer_names(data$outcomes, dim(Y)[3], "outcome")
  policies = layer_names(data$policies, dim(policy)[3], "policy")
  
  dimnames(Y) = list(country = countries, time = times, layer = outcomes)
  dimnames(policy) = list(country = countries, time = times, layer = policies)
  
  map = data$policy_outcome_map
  if (is.null(map) || !all(c("policy", "outcome") %in% names(map))) {
    map = tibble(outcome = outcomes, policy = policies, panel_label = outcomes)
  } else {
    map = as_tibble(map) %>%
      transmute(
        outcome = as.character(outcome),
        policy = as.character(policy),
        panel_label = if ("panel_label" %in% names(map)) {
          as.character(panel_label)
        } else {
          as.character(outcome)
        }
      ) %>%
      filter(outcome %in% outcomes, policy %in% policies)
  }
  
  map = map %>%
    mutate(
      outcome_index = match(outcome, outcomes),
      policy_index = match(policy, policies)
    ) %>%
    arrange(outcome_index)
  
  target = map %>%
    filter(
      str_detect(str_to_lower(outcome), target_pattern) |
        str_detect(str_to_lower(panel_label), target_pattern)
    ) %>%
    mutate(exact = !str_to_lower(outcome) %in% c("death", "deaths")) %>%
    arrange(exact, outcome_index) %>%
    slice(1)
  
  if (!nrow(target)) stop("No target outcome matched: ", target_pattern)
  
  selected = if (nrow(map) == n_keep) {
    map
  } else {
    bind_rows(
      target,
      map %>% filter(outcome != target$outcome[[1]]) %>% slice_head(n = n_keep - 1)
    ) %>%
      distinct(outcome, .keep_all = TRUE) %>%
      arrange(outcome_index)
  }
  
  if (nrow(selected) != n_keep) stop("Could not select the requested outcome layers.")
  
  Y = Y[, , selected$outcome_index, drop = FALSE]
  policy = policy[, , selected$policy_index, drop = FALSE]
  dimnames(Y)[[3]] = selected$outcome
  dimnames(policy)[[3]] = selected$outcome
  
  layer_map = selected %>%
    transmute(
      layer_index = row_number(),
      outcome,
      policy,
      panel_label,
      original_outcome_index = outcome_index,
      original_policy_index = policy_index
    )
  
  target_k = match(target$outcome[[1]], selected$outcome)
  
  list(
    Y = fill_time_neighbours(Y),
    policy = policy,
    layer_map = layer_map,
    target_k = target_k,
    target_outcome = selected$outcome[[target_k]],
    target_policy = selected$policy[[target_k]]
  )
}

# Estimation --------------------------------------------------------------

make_staircase = function(policy, target_k) {
  target_policy = policy[, , target_k]
  treatment_time = apply(target_policy, 1, first_true)
  observed_length = ifelse(is.infinite(treatment_time), ncol(target_policy), treatment_time - 1)
  row_names = dimnames(policy)[[1]] %||% as.character(seq_len(nrow(target_policy)))
  
  permutation = order(-observed_length, row_names, seq_along(row_names))
  observed_sorted = observed_length[permutation]
  observed_blocks = unique(observed_sorted)
  row_blocks = lapply(observed_blocks, function(x) which(observed_sorted == x))
  
  list(
    permutation = permutation,
    N_parts = as.integer(lengths(row_blocks)),
    T_parts = as.integer(diff(c(0, rev(observed_blocks))))
  )
}

make_A = function(policy) {
  A = apply(policy, c(1, 3), first_true)
  dimnames(A) = list(row = dimnames(policy)[[1]], layer = dimnames(policy)[[3]])
  A
}

estimate_one = function(Y, policy, target_k, r, functional, setting, layer_map) {
  staircase = make_staircase(policy, target_k)
  permutation = staircase$permutation
  
  Y = Y[permutation, , , drop = FALSE]
  policy = policy[permutation, , , drop = FALSE]
  observed = !policy
  
  Y_observed = Y
  Y_observed[!observed] = NA_real_
  
  K = dim(Y)[3]
  N_parts = replicate(K, staircase$N_parts, simplify = FALSE)
  T_parts = replicate(K, staircase$T_parts, simplify = FALSE)
  A = make_A(policy)
  
  psi0_tensor = bilinearTensorStaggeredPsi(
    Y = Y_observed,
    k = target_k,
    r = r,
    tau = tau,
    functional = functional,
    eta = NULL,
    row_index = NULL,
    A = A,
    Omega = observed
  )
  
  psi0_matrix = bilinearMatrixStaggeredPsi(
    Y_mat = Y_observed[, , target_k],
    r = r,
    tau = tau,
    functional = functional,
    eta = NULL,
    row_index = NULL,
    A = A[, target_k],
    Omega = observed[, , target_k]
  )
  
  psi1 = pluginPsi_c1(
    Y = Y,
    k = target_k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = functional,
    eta = NULL,
    row_index = NULL
  )
  
  target_outcome = dimnames(Y)[[3]][target_k]
  target_policy = layer_map$policy[match(target_outcome, layer_map$outcome)]
  
  tibble(
    setting = setting,
    target_outcome = target_outcome,
    target_policy = target_policy,
    target_k = target_k,
    r = r,
    functional = functional,
    n_staircase_blocks = length(staircase$N_parts),
    N_parts = paste(staircase$N_parts, collapse = ","),
    T_parts = paste(staircase$T_parts, collapse = ","),
    Psi0 = as.numeric(psi0_tensor),
    Psi1 = as.numeric(psi1),
    Delta_h = as.numeric(psi1 - psi0_tensor),
    Psi0_matrix = as.numeric(psi0_matrix),
    Delta_h_matrix = as.numeric(psi1 - psi0_matrix),
    tensor_minus_matrix = as.numeric(psi0_tensor - psi0_matrix)
  )
}

run_estimator = function(Y, policy, target_k, layer_map) {
  expand_grid(r = r_grid, functional = functionals) %>%
    pmap_dfr(
      ~ estimate_one(
        Y = Y,
        policy = policy,
        target_k = target_k,
        r = ..1,
        functional = ..2,
        setting = "plotY_slice",
        layer_map = layer_map
      )
    )
}

# Bootstrap ---------------------------------------------------------------

result_columns = c(
  "Psi0",
  "Psi1",
  "Delta_h",
  "Psi0_matrix",
  "Delta_h_matrix",
  "tensor_minus_matrix"
)

row_strata = function(Y, policy, target_k) {
  target_Y = Y[, , target_k]
  target_policy = policy[, , target_k]
  
  fully_observed = which(
    rowSums(target_policy == TRUE, na.rm = TRUE) == 0 &
      rowSums(is.na(target_Y)) == 0
  )
  
  list(
    target_layer_index = target_k,
    fully_observed_idx = fully_observed,
    non_fully_observed_idx = setdiff(seq_len(dim(Y)[1]), fully_observed)
  )
}

resample_rows = function(Y, policy, target_k) {
  strata = row_strata(Y, policy, target_k)
  full = strata$fully_observed_idx
  nonfull = strata$non_fully_observed_idx
  
  if (stratify_target_bootstrap && length(full) && length(nonfull)) {
    index = seq_len(dim(Y)[1])
    index[full] = sample(full, length(full), replace = TRUE)
    index[nonfull] = sample(nonfull, length(nonfull), replace = TRUE)
    index
  } else {
    sample(seq_len(dim(Y)[1]), dim(Y)[1], replace = TRUE)
  }
}

bootstrap_once = function(b, Y, policy, target_k, layer_map, point_results) {
  set.seed(bootstrap_seed + b)
  index = resample_rows(Y, policy, target_k)
  
  Y_boot = Y
  Y_boot[, , target_k] = Y[index, , target_k]
  
  tryCatch(
    run_estimator(Y_boot, policy, target_k, layer_map) %>%
      mutate(
        bootstrap_id = b,
        bootstrap_type = "target_layer_only_fixed_anchors",
        bootstrap_error = NA_character_
      ),
    error = function(error) {
      point_results %>%
        mutate(
          across(all_of(result_columns), ~ NA_real_),
          bootstrap_id = b,
          bootstrap_type = "target_layer_only_fixed_anchors",
          bootstrap_error = conditionMessage(error)
        )
    }
  )
}

make_ci = function(point_results, bootstrap_results) {
  bootstrap_summary = bootstrap_results %>%
    select(
      bootstrap_id,
      setting,
      target_outcome,
      target_policy,
      target_k,
      r,
      functional,
      all_of(result_columns)
    ) %>%
    pivot_longer(all_of(result_columns), names_to = "quantity", values_to = "boot_value") %>%
    group_by(setting, target_outcome, target_policy, target_k, r, functional, quantity) %>%
    summarise(
      boot_mean = mean(boot_value, na.rm = TRUE),
      boot_se = sd(boot_value, na.rm = TRUE),
      n_boot_nonmissing = sum(!is.na(boot_value)),
      .groups = "drop"
    )
  
  point_results %>%
    select(
      setting,
      target_outcome,
      target_policy,
      target_k,
      r,
      functional,
      all_of(result_columns)
    ) %>%
    pivot_longer(all_of(result_columns), names_to = "quantity", values_to = "point_estimate") %>%
    left_join(
      bootstrap_summary,
      by = c(
        "setting",
        "target_outcome",
        "target_policy",
        "target_k",
        "r",
        "functional",
        "quantity"
      )
    ) %>%
    mutate(
      ci_low = if_else(
        n_boot_nonmissing >= 2,
        point_estimate - 1.96 * boot_se,
        NA_real_
      ),
      ci_high = if_else(
        n_boot_nonmissing >= 2,
        point_estimate + 1.96 * boot_se,
        NA_real_
      ),
      ci_type = "bootstrap_se_centered",
      across(c(point_estimate, boot_mean, boot_se, ci_low, ci_high), ~ round(.x, 6)),
      estimate_ci = if_else(
        is.na(ci_low) | is.na(ci_high),
        sprintf("%.4f [NA, NA]", point_estimate),
        sprintf("%.4f [%.4f, %.4f]", point_estimate, ci_low, ci_high)
      )
    ) %>%
    arrange(setting, target_outcome, r, functional, quantity)
}

# Run ---------------------------------------------------------------------

tensor_data = readRDS(tensor_file)
prepared = align_layers(tensor_data, target_outcome_pattern, n_outcomes_to_keep)

Y = prepared$Y
policy = prepared$policy
layer_map = prepared$layer_map
target_k = prepared$target_k

if (anyNA(Y)) stop("Y still contains missing values after neighbour filling.")

set.seed(bootstrap_seed)
point_results = run_estimator(Y, policy, target_k, layer_map) %>%
  mutate(across(all_of(result_columns), ~ round(.x, 6)))

bootstrap_results = map_dfr(
  seq_len(B),
  bootstrap_once,
  Y = Y,
  policy = policy,
  target_k = target_k,
  layer_map = layer_map,
  point_results = point_results
)

results_with_ci = make_ci(point_results, bootstrap_results)

main_summary = point_results %>%
  filter(r == r_grid[1]) %>%
  select(
    setting,
    target_outcome,
    target_policy,
    functional,
    Psi0,
    Psi0_matrix,
    Psi1,
    Delta_h,
    Delta_h_matrix,
    tensor_minus_matrix,
    N_parts,
    T_parts
  ) %>%
  arrange(setting, target_outcome, functional)

slice_label = paste0(
  "until_",
  tensor_data$time_horizon %||% "unknown",
  "_delay_",
  tensor_data$delay_days %||% "unknown"
)

analysis_results = list(
  metadata = list(
    tensor_file = tensor_file,
    slice_label = slice_label,
    selected_layer_map = layer_map,
    target_outcome = prepared$target_outcome,
    target_policy = prepared$target_policy,
    target_k = target_k,
    r_grid = r_grid,
    tau = tau,
    functionals = functionals,
    B = B,
    bootstrap_seed = bootstrap_seed,
    rng_kind = RNGkind(),
    stratify_target_bootstrap = stratify_target_bootstrap,
    ci_vars = result_columns,
    project_dir = getwd()
  ),
  layer_map = layer_map,
  point_results = point_results,
  bootstrap_results = bootstrap_results,
  results_with_ci = results_with_ci,
  main_summary = main_summary,
  target_row_strata = row_strata(Y, policy, target_k),
  Y = Y,
  Omega_policy = policy
)

saveRDS(
  analysis_results,
  file.path(results_dir, paste0("oxford_deaths_target_layer_bootstrap_", slice_label, ".rds"))
)

outputs = list(
  layer_map = layer_map,
  point_results = point_results,
  bootstrap_results = bootstrap_results,
  results_with_ci = results_with_ci,
  main_summary = main_summary
)

filenames = c(
  layer_map = paste0("oxford_deaths_layer_map_", slice_label, ".csv"),
  point_results = paste0("oxford_deaths_point_results_", slice_label, ".csv"),
  bootstrap_results = paste0("oxford_deaths_bootstrap_results_", slice_label, ".csv"),
  results_with_ci = paste0("oxford_deaths_results_with_ci_", slice_label, ".csv"),
  main_summary = paste0("oxford_deaths_main_summary_r", r_grid[1], "_", slice_label, ".csv")
)

walk2(outputs, filenames, ~ write_csv(.x, file.path(results_dir, .y)))
