setwd("~/Desktop/code")

suppressPackageStartupMessages(library(tidyverse))

source("bilinearTensorAllFunction.R")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

results_dir = "results/masked_nonmasked_rank3_bootstrap"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

target_layer = "l_robbery"
target_name = "robbery"

n_observed_years_for_other_rows = 3
n_fully_observed_grid = c(5, 10, 15)
rank_value = n_observed_years_for_other_rows

B = 500
bootstrap_seed = 231198
tau = 0.01

RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

keep_local_target_rows = FALSE
stratify_target_bootstrap = TRUE
local_states = c("Montana", "Texas", "Florida")

functional_order = c(
  "ATE",
  "Local-Florida",
  "Local-Montana",
  "Local-Texas",
  "RowHet",
  "Trend"
)

ci_vars = c(
  "Psi0",
  "Psi1",
  "Delta_h",
  "Psi0_matrix",
  "Delta_h_matrix",
  "tensor_minus_matrix"
)

crime_vars = c("l_motor", "l_robbery", "l_assault", "l_homicide")
crime_labels = c(
  l_motor = "Motor theft rate, log",
  l_robbery = "Robbery rate, log",
  l_assault = "Aggravated assault rate, log",
  l_homicide = "Murder rate, log"
)

# -----------------------------------------------------------------------------
# Load and shape data
# -----------------------------------------------------------------------------

castle = readRDS(
  url("https://github.com/guerramarcelino/PolicyEval/raw/main/Datasets/castle.RDS")
) %>%
  as_tibble()

state_col = if ("state" %in% names(castle)) "state" else "sid"
treat_col = if ("cdl" %in% names(castle)) "cdl" else "post"

castle_long = castle %>%
  transmute(
    state_id = as.character(.data[[state_col]]),
    year = as.integer(year),
    treatment_on = as.numeric(.data[[treat_col]]) > 0,
    across(all_of(crime_vars), as.numeric)
  ) %>%
  pivot_longer(
    cols = all_of(crime_vars),
    names_to = "crime",
    values_to = "outcome_value"
  ) %>%
  mutate(crime = factor(crime, levels = crime_vars)) %>%
  group_by(state_id, year, crime) %>%
  summarize(
    treatment_value = as.numeric(any(treatment_on, na.rm = TRUE)),
    outcome_value = mean(outcome_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(outcome_value = if_else(is.nan(outcome_value), NA_real_, outcome_value))

state_order = castle_long %>%
  distinct(state_id, year, treatment_value) %>%
  group_by(state_id) %>%
  summarize(
    ever_treated = any(treatment_value == 1),
    first_treat_year = if (ever_treated) min(year[treatment_value == 1]) else NA_integer_,
    .groups = "drop"
  ) %>%
  arrange(ever_treated, desc(first_treat_year), state_id) %>%
  mutate(row_index = row_number())

state_order_vec = state_order$state_id
year_order = sort(unique(castle_long$year))
crime_order = crime_vars

make_layer_matrix = function(data, crime_name, value_col) {
  data %>%
    filter(crime == crime_name) %>%
    mutate(
      state_id = factor(state_id, levels = state_order_vec),
      year = factor(year, levels = year_order)
    ) %>%
    arrange(state_id, year) %>%
    transmute(state_id, year, value = .data[[value_col]]) %>%
    pivot_wider(names_from = year, values_from = value) %>%
    arrange(state_id) %>%
    select(-state_id) %>%
    as.matrix()
}

array_dims = c(length(state_order_vec), length(year_order), length(crime_order))
array_names = list(
  state = state_order_vec,
  year = as.character(year_order),
  crime = crime_order
)

Y = array(NA_real_, dim = array_dims, dimnames = array_names)
D_original = array(NA_real_, dim = array_dims, dimnames = array_names)

for (k in seq_along(crime_order)) {
  Y[, , k] = make_layer_matrix(castle_long, crime_order[k], "outcome_value")
  D_original[, , k] = make_layer_matrix(castle_long, crime_order[k], "treatment_value")
}

storage.mode(Y) = "numeric"
storage.mode(D_original) = "numeric"

rn = dimnames(Y)$state
target_label = unname(crime_labels[target_layer])

# -----------------------------------------------------------------------------
# Design and estimator helpers
# -----------------------------------------------------------------------------

make_A_Omega_from_D = function(D) {
  N = dim(D)[1]
  Tt = dim(D)[2]
  K = dim(D)[3]
  
  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(FALSE, dim = c(N, Tt, K))
  
  for (k in seq_len(K)) {
    for (i in seq_len(N)) {
      first_on = which(D[i, , k] == 1)[1]
      
      if (is.na(first_on)) {
        Omega[i, , k] = TRUE
      } else {
        A[i, k] = first_on
        if (first_on > 1) Omega[i, seq_len(first_on - 1), k] = TRUE
      }
    }
  }
  
  list(A = A, Omega = Omega)
}

make_staggered_parts_from_D = function(D, layer) {
  Dmat = D[, , layer]
  years = as.integer(colnames(Dmat))
  
  first_on = apply(Dmat, 1, function(row) {
    treated_years = years[row == 1]
    if (length(treated_years) == 0) Inf else min(treated_years)
  })
  
  adopt_years = sort(unique(first_on[is.finite(first_on)]))
  m = length(adopt_years)
  
  if (m == 0) {
    return(list(
      N_sizes = nrow(Dmat),
      T_sizes = ncol(Dmat)
    ))
  }
  
  row_groups = vector("list", m + 1)
  row_groups[[1]] = which(!is.finite(first_on))
  
  for (group in seq_len(m)) {
    row_groups[[group + 1]] = which(first_on == rev(adopt_years)[group])
  }
  
  col_groups = vector("list", m + 1)
  col_groups[[1]] = which(years < adopt_years[1])
  
  if (m >= 2) {
    for (group in 2:m) {
      col_groups[[group]] = which(
        years >= adopt_years[group - 1] & years < adopt_years[group]
      )
    }
  }
  
  col_groups[[m + 1]] = which(years >= adopt_years[m])
  
  list(
    N_sizes = lengths(row_groups),
    T_sizes = lengths(col_groups)
  )
}

make_masked_target_D = function(
    D_base,
    target_layer,
    n_fully_observed_rows,
    n_observed_years_for_other_rows
) {
  D = D_base
  target_idx = match(target_layer, dimnames(D)$crime)
  
  full_rows = seq_len(min(n_fully_observed_rows, dim(D)[1]))
  short_rows = setdiff(seq_len(dim(D)[1]), full_rows)
  early_years = seq_len(min(n_observed_years_for_other_rows, dim(D)[2]))
  late_years = setdiff(seq_len(dim(D)[2]), early_years)
  
  D[, , target_idx] = 0
  D[short_rows, late_years, target_idx] = 1
  D
}

make_design_D = function(design_row) {
  if (design_row$design_type[[1]] == "nonmasked") return(D_original)
  
  make_masked_target_D(
    D_base = D_original,
    target_layer = target_layer,
    n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
    n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]]
  )
}

bush_2000_states = c(
  "Alabama", "Alaska", "Arizona", "Arkansas", "Colorado",
  "Florida", "Georgia", "Idaho", "Indiana", "Kansas",
  "Kentucky", "Louisiana", "Mississippi", "Missouri", "Montana",
  "Nebraska", "Nevada", "New Hampshire", "North Carolina",
  "North Dakota", "Ohio", "Oklahoma", "South Carolina",
  "South Dakota", "Tennessee", "Texas", "Utah", "Virginia",
  "West Virginia", "Wyoming"
)

eta = ifelse(rn %in% bush_2000_states, 1, -1)
names(eta) = rn

run_estimator = function(Y_in, D_in) {
  Y0 = Y_in
  Y0[D_in == 1] = NA_real_
  
  AO = make_A_Omega_from_D(D_in)
  A = AO$A
  Omega = AO$Omega & !is.na(Y0)
  
  target_idx = match(target_layer, dimnames(Y_in)$crime)
  parts = make_staggered_parts_from_D(D_in, target_idx)
  N_parts = replicate(dim(Y_in)[3], parts$N_sizes, simplify = FALSE)
  T_parts = replicate(dim(Y_in)[3], parts$T_sizes, simplify = FALSE)
  
  estimate_one = function(functional, local_state = NA_character_) {
    row_index = if (functional == "Local") match(local_state, rn) else NULL
    functional_label = if (functional == "Local") {
      paste0("Local-", local_state)
    } else {
      functional
    }
    
    if (functional == "Local" && is.na(row_index)) {
      return(tibble(
        r = rank_value,
        crime = target_name,
        layer = target_layer,
        functional = functional_label,
        local_state = local_state,
        Psi0 = NA_real_,
        Psi1 = NA_real_,
        Delta_h = NA_real_,
        Psi0_matrix = NA_real_,
        Delta_h_matrix = NA_real_,
        tensor_minus_matrix = NA_real_
      ))
    }
    
    eta_arg = if (functional == "RowHet") eta else NULL
    
    psi0_tensor = bilinearTensorStaggeredPsi(
      Y = Y0,
      k = target_idx,
      r = rank_value,
      tau = tau,
      functional = functional,
      eta = eta_arg,
      row_index = row_index,
      A = A,
      Omega = Omega
    )
    
    psi0_matrix = bilinearMatrixStaggeredPsi(
      Y_mat = Y0[, , target_idx],
      r = rank_value,
      tau = tau,
      functional = functional,
      eta = eta_arg,
      row_index = row_index,
      A = A[, target_idx],
      Omega = Omega[, , target_idx]
    )
    
    psi1 = pluginPsi_c1(
      Y = Y_in,
      k = target_idx,
      N_parts = N_parts,
      T_parts = T_parts,
      functional = functional,
      eta = eta_arg,
      row_index = row_index
    )
    
    tibble(
      r = rank_value,
      crime = target_name,
      layer = target_layer,
      functional = functional_label,
      local_state = if (functional == "Local") local_state else NA_character_,
      Psi0 = as.numeric(psi0_tensor),
      Psi1 = as.numeric(psi1),
      Delta_h = as.numeric(psi1 - psi0_tensor),
      Psi0_matrix = as.numeric(psi0_matrix),
      Delta_h_matrix = as.numeric(psi1 - psi0_matrix),
      tensor_minus_matrix = as.numeric(psi0_tensor - psi0_matrix)
    )
  }
  
  specs = bind_rows(
    tibble(
      functional = c("ATE", "RowHet", "Trend"),
      local_state = NA_character_
    ),
    tibble(
      functional = "Local",
      local_state = local_states
    )
  )
  
  pmap_dfr(specs, estimate_one) %>%
    mutate(functional = factor(functional, levels = functional_order)) %>%
    arrange(functional) %>%
    mutate(functional = as.character(functional))
}

# -----------------------------------------------------------------------------
# Bootstrap helpers
# -----------------------------------------------------------------------------

get_target_row_strata = function(Y, D) {
  target_idx = match(target_layer, dimnames(D)$crime)
  Y_target = Y[, , target_idx]
  D_target = D[, , target_idx]
  
  fully_observed = which(
    rowSums(D_target == 1, na.rm = TRUE) == 0 &
      rowSums(is.na(Y_target)) == 0
  )
  
  list(
    fully_observed = fully_observed,
    other = setdiff(seq_len(dim(D)[1]), fully_observed)
  )
}

make_target_resample_index = function(Y, D) {
  strata = get_target_row_strata(Y, D)
  N = dim(D)[1]
  
  if (
    stratify_target_bootstrap &&
    length(strata$fully_observed) > 0 &&
    length(strata$other) > 0
  ) {
    index = seq_len(N)
    index[strata$fully_observed] = sample(
      strata$fully_observed,
      length(strata$fully_observed),
      replace = TRUE
    )
    index[strata$other] = sample(
      strata$other,
      length(strata$other),
      replace = TRUE
    )
  } else {
    index = sample(seq_len(N), N, replace = TRUE)
  }
  
  if (keep_local_target_rows) {
    local_index = match(local_states, rn)
    local_index = local_index[!is.na(local_index)]
    index[local_index] = local_index
  }
  
  index
}

run_bootstrap_replication = function(b, D_design) {
  set.seed(bootstrap_seed + b)
  
  target_idx = match(target_layer, dimnames(Y)$crime)
  index = make_target_resample_index(Y, D_design)
  
  Y_boot = Y
  Y_boot[, , target_idx] = Y[index, , target_idx]
  
  run_estimator(Y_boot, D_design) %>%
    mutate(
      bootstrap_id = b,
      bootstrap_type = "target_layer_only_fixed_anchors",
      bootstrap_error = NA_character_
    )
}

failed_bootstrap_rows = function(point_results, b, message) {
  point_results %>%
    mutate(
      across(all_of(ci_vars), ~ NA_real_),
      bootstrap_id = b,
      bootstrap_type = "target_layer_only_fixed_anchors",
      bootstrap_error = message
    )
}

make_results_with_ci = function(point_results, bootstrap_results) {
  bootstrap_summary = bootstrap_results %>%
    select(
      bootstrap_id,
      r,
      crime,
      layer,
      functional,
      local_state,
      all_of(ci_vars)
    ) %>%
    pivot_longer(
      cols = all_of(ci_vars),
      names_to = "quantity",
      values_to = "boot_value"
    ) %>%
    group_by(r, crime, layer, functional, local_state, quantity) %>%
    summarize(
      boot_mean = mean(boot_value, na.rm = TRUE),
      boot_se = sd(boot_value, na.rm = TRUE),
      n_boot_nonmissing = sum(!is.na(boot_value)),
      .groups = "drop"
    )
  
  point_results %>%
    select(
      r,
      crime,
      layer,
      functional,
      local_state,
      all_of(ci_vars)
    ) %>%
    pivot_longer(
      cols = all_of(ci_vars),
      names_to = "quantity",
      values_to = "point_estimate"
    ) %>%
    left_join(
      bootstrap_summary,
      by = c("r", "crime", "layer", "functional", "local_state", "quantity")
    ) %>%
    mutate(
      ci_low = if_else(
        !is.na(point_estimate) & !is.na(boot_se) & n_boot_nonmissing >= 2,
        point_estimate - 1.96 * boot_se,
        NA_real_
      ),
      ci_high = if_else(
        !is.na(point_estimate) & !is.na(boot_se) & n_boot_nonmissing >= 2,
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
    )
}

add_design_columns = function(results, design_row) {
  results %>%
    mutate(
      design_id = design_row$design_id[[1]],
      design_label = design_row$design_label[[1]],
      design_type = design_row$design_type[[1]],
      n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
      n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]],
      .before = 1
    )
}

# -----------------------------------------------------------------------------
# Designs and simulation
# -----------------------------------------------------------------------------

design_settings = bind_rows(
  tibble(
    design_id = paste0(
      "masked_full",
      n_fully_observed_grid,
      "_obs",
      n_observed_years_for_other_rows
    ),
    design_label = paste0(
      "Masked: ",
      n_fully_observed_grid,
      " fully observed rows, ",
      n_observed_years_for_other_rows,
      " observed years elsewhere"
    ),
    design_type = "masked",
    n_fully_observed_rows = n_fully_observed_grid,
    n_observed_years_for_other_rows = n_observed_years_for_other_rows,
    rank_value = rank_value,
    B = B
  ),
  tibble(
    design_id = "nonmasked_original",
    design_label = "Non-masked: original Castle staggered design",
    design_type = "nonmasked",
    n_fully_observed_rows = NA_integer_,
    n_observed_years_for_other_rows = NA_integer_,
    rank_value = rank_value,
    B = B
  )
)

run_one_design = function(design_row) {
  design_id = design_row$design_id[[1]]
  D_design = make_design_D(design_row)
  
  point_results = run_estimator(Y, D_design) %>%
    add_design_columns(design_row) %>%
    mutate(across(all_of(ci_vars), ~ round(.x, 6)))
  
  bootstrap_results = map_dfr(seq_len(design_row$B[[1]]), function(b) {
    tryCatch(
      run_bootstrap_replication(b, D_design),
      error = function(e) failed_bootstrap_rows(point_results, b, conditionMessage(e))
    )
  }) %>%
    add_design_columns(design_row)
  
  results_with_ci = make_results_with_ci(point_results, bootstrap_results) %>%
    add_design_columns(design_row)
  
  result = list(
    design_settings = design_row,
    point_results = point_results,
    bootstrap_results = bootstrap_results,
    results_with_ci = results_with_ci,
    D_design = D_design
  )
  
  saveRDS(
    result,
    file.path(results_dir, paste0("design_", target_name, "_", design_id, ".rds"))
  )
  write_csv(
    point_results,
    file.path(results_dir, paste0("point_results_", target_name, "_", design_id, ".csv"))
  )
  write_csv(
    bootstrap_results,
    file.path(results_dir, paste0("bootstrap_results_", target_name, "_", design_id, ".csv"))
  )
  write_csv(
    results_with_ci,
    file.path(results_dir, paste0("results_with_ci_", target_name, "_", design_id, ".csv"))
  )
  
  result
}

design_objects = map(seq_len(nrow(design_settings)), ~ run_one_design(design_settings[.x, ]))
names(design_objects) = design_settings$design_id

point_results_all = map_dfr(design_objects, "point_results")
bootstrap_results_all = map_dfr(design_objects, "bootstrap_results")
results_with_ci_all = map_dfr(design_objects, "results_with_ci")

analysis_results = list(
  metadata = list(
    target_layer = target_layer,
    target_name = target_name,
    target_label = target_label,
    n_observed_years_for_other_rows = n_observed_years_for_other_rows,
    n_fully_observed_grid = n_fully_observed_grid,
    rank_value = rank_value,
    B = B,
    bootstrap_seed = bootstrap_seed,
    tau = tau,
    keep_local_target_rows = keep_local_target_rows,
    stratify_target_bootstrap = stratify_target_bootstrap,
    local_states = local_states,
    functional_order = functional_order,
    ci_vars = ci_vars,
    crime_vars = crime_vars,
    crime_labels = crime_labels
  ),
  design_settings = design_settings,
  point_results = point_results_all,
  bootstrap_results = bootstrap_results_all,
  results_with_ci = results_with_ci_all,
  design_objects = design_objects,
  Y = Y,
  D_original = D_original,
  state_order = state_order,
  year_order = year_order,
  crime_order = crime_order
)

saveRDS(
  analysis_results,
  file.path(
    results_dir,
    paste0("castle_masked_nonmasked_bootstrap_results_", target_name, ".rds")
  )
)

write_csv(
  design_settings,
  file.path(results_dir, paste0("design_settings_", target_name, ".csv"))
)
write_csv(
  point_results_all,
  file.path(results_dir, paste0("point_results_all_designs_", target_name, ".csv"))
)
write_csv(
  bootstrap_results_all,
  file.path(results_dir, paste0("bootstrap_results_all_designs_", target_name, ".csv"))
)
write_csv(
  results_with_ci_all,
  file.path(results_dir, paste0("results_with_ci_all_designs_", target_name, ".csv"))
)
