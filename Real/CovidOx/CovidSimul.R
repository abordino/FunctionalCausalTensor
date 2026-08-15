setwd("~/Desktop/code")

source("bilinearTensorAllFunction.R")
source("CYNuclear.R")

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

n_outcomes_to_keep = 2
target_outcome_pattern = "death"
r_grid = 1 # tensor method uses rank = 3 later on
tau = 1e-4
cy_seed = 123
functionals = c("ATE", "Trend")

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
      map %>%
        filter(outcome != target$outcome[[1]]) %>%
        slice_head(n = n_keep - 1)
    ) %>%
      distinct(outcome, .keep_all = TRUE) %>%
      arrange(outcome_index)
  }
  
  if (nrow(selected) != n_keep) {
    stop("Could not select the requested outcome layers.")
  }
  
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
  observed_length = ifelse(
    is.infinite(treatment_time),
    ncol(target_policy),
    treatment_time - 1
  )
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
  dimnames(A) = list(
    row = dimnames(policy)[[1]],
    layer = dimnames(policy)[[3]]
  )
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
  
  cy_rank = r
  cy_fit = NULL
  
  while (is.null(cy_fit)) {
    cy_error = NULL
    
    cy_fit = tryCatch(
      cy_complete_staggered_layer(
        y = Y_observed[, , target_k],
        omega = observed[, , target_k],
        adoption_time = A[, target_k],
        time = seq_len(dim(Y)[2]),
        rank = cy_rank,
        seed = cy_seed,
        verbose = FALSE
      ),
      error = function(e) {
        cy_error <<- e
        NULL
      }
    )
    
    if (is.null(cy_fit)) {
      if (!grepl(
        "fewer than rank",
        conditionMessage(cy_error),
        fixed = TRUE
      )) {
        stop(cy_error)
      }
      
      cy_rank = cy_rank - 1L
      
      if (cy_rank < 1L) {
        stop(cy_error)
      }
    }
  }
  
  psi0_tensor = bilinearTensorStaggeredPsi(
    Y = Y_observed,
    k = target_k,
    r = 3,
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
  
  psi0_CY_nuclear = pluginPsi_c1(
    Y = cy_fit$completed_nuclear,
    k = target_k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = functional,
    eta = NULL,
    row_index = NULL
  )
  
  psi0_CY_debiased = pluginPsi_c1(
    Y = cy_fit$completed_debiased,
    k = target_k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = functional,
    eta = NULL,
    row_index = NULL
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
  target_policy = layer_map$policy[
    match(target_outcome, layer_map$outcome)
  ]
  
  tibble(
    setting = setting,
    target_outcome = target_outcome,
    target_policy = target_policy,
    target_k = target_k,
    r = r,
    r_CY = cy_rank,
    functional = functional,
    n_staircase_blocks = length(staircase$N_parts),
    N_parts = paste(staircase$N_parts, collapse = ","),
    T_parts = paste(staircase$T_parts, collapse = ","),
    Psi0 = as.numeric(psi0_tensor),
    Psi1 = as.numeric(psi1),
    Delta_h = as.numeric(psi1 - psi0_tensor),
    Psi0_matrix = as.numeric(psi0_matrix),
    Delta_h_matrix = as.numeric(psi1 - psi0_matrix),
    tensor_minus_matrix = as.numeric(psi0_tensor - psi0_matrix),
    Psi0_CY_nuclear = as.numeric(psi0_CY_nuclear),
    Delta_h_CY_nuclear = as.numeric(psi1 - psi0_CY_nuclear),
    Psi0_CY_debiased = as.numeric(psi0_CY_debiased),
    Delta_h_CY_debiased = as.numeric(psi1 - psi0_CY_debiased)
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

# Run ---------------------------------------------------------------------

tensor_data = readRDS(tensor_file)

prepared = align_layers(
  tensor_data,
  target_outcome_pattern,
  n_outcomes_to_keep
)

Y = prepared$Y
policy = prepared$policy
layer_map = prepared$layer_map
target_k = prepared$target_k

if (anyNA(Y)) stop("Y still contains missing values after neighbour filling.")

result_columns = c(
  "Psi0",
  "Psi1",
  "Delta_h",
  "Psi0_matrix",
  "Delta_h_matrix",
  "tensor_minus_matrix",
  "Psi0_CY_nuclear",
  "Delta_h_CY_nuclear",
  "Psi0_CY_debiased",
  "Delta_h_CY_debiased"
)

point_results = run_estimator(Y, policy, target_k, layer_map) %>%
  mutate(across(all_of(result_columns), ~ round(.x, 6)))

main_summary = point_results %>%
  filter(r == r_grid[1]) %>%
  select(
    setting,
    target_outcome,
    target_policy,
    functional,
    r,
    r_CY,
    Psi0,
    Psi0_matrix,
    Psi0_CY_nuclear,
    Psi0_CY_debiased,
    Psi1,
    Delta_h,
    Delta_h_matrix,
    Delta_h_CY_nuclear,
    Delta_h_CY_debiased,
    tensor_minus_matrix,
    N_parts,
    T_parts
  ) %>%
  arrange(setting, target_outcome, functional)

print(main_summary, n = Inf, width = Inf)