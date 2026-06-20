PROJECT_DIR = "~/Documents/phd/projects/causalMatrix/code/real-world/CastleDoctrine"
if (dir.exists(path.expand(PROJECT_DIR))) {
  setwd(path.expand(PROJECT_DIR))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ===============================================================
# 0. Set par
# ===============================================================

results_dir = "results/masked_nonmasked_rank3_bootstrap"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# "l_motor", "l_robbery", "l_assault", "l_homicide"
target_layer = "l_robbery"
target_name = "robbery"
target_file_stub = target_name

n_observed_years_for_other_rows = 3
n_fully_observed_grid = c(5, 10, 15)

rank_value = n_observed_years_for_other_rows
B = 500
bootstrap_seed = 231198
tau = 0.01

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

# ===============================================================
# 1. Load Castle Doctrine data
# ===============================================================

castle_url = "https://github.com/guerramarcelino/PolicyEval/raw/main/Datasets/castle.RDS"

castle = readRDS(url(castle_url)) %>%
  as_tibble()

crime_vars = c("l_motor", "l_robbery", "l_assault", "l_homicide")

crime_labels = c(
  l_motor    = "Motor theft rate, log",
  l_robbery  = "Robbery rate, log",
  l_assault  = "Aggravated assault rate, log",
  l_homicide = "Murder rate, log"
)

state_col = if ("state" %in% names(castle)) "state" else "sid"
treat_col = if ("cdl" %in% names(castle)) "cdl" else "post"

# ===============================================================
# 2. Long data
# ===============================================================

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
  mutate(
    crime = factor(crime, levels = crime_vars),
    crime_label = factor(
      unname(crime_labels[as.character(crime)]),
      levels = unname(crime_labels)
    )
  ) %>%
  group_by(state_id, year, crime, crime_label) %>%
  summarize(
    treatment_on = any(treatment_on, na.rm = TRUE),
    treatment_value = as.numeric(treatment_on),
    treatment_status = if_else(treatment_on, "On", "Off"),
    outcome_value = mean(outcome_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    outcome_value = if_else(is.nan(outcome_value), NA_real_, outcome_value)
  )

# ===============================================================
# 3. Common staircase row ordering
# ===============================================================

state_order = castle_long %>%
  distinct(state_id, year, treatment_on) %>%
  group_by(state_id) %>%
  summarize(
    ever_treated = any(treatment_on),
    first_treat_year = if_else(
      ever_treated,
      min(year[treatment_on]),
      NA_integer_
    ),
    .groups = "drop"
  ) %>%
  arrange(
    ever_treated,
    desc(first_treat_year),
    state_id
  ) %>%
  mutate(row_index = row_number())

state_order_vec = state_order$state_id
year_order = sort(unique(castle_long$year))
crime_order = crime_vars

n_states = length(state_order_vec)
n_years = length(year_order)
n_crimes = length(crime_order)

state_order = state_order %>%
  mutate(y_plot = n_states - row_index + 1)

castle_tensor_long = castle_long %>%
  left_join(
    state_order %>%
      select(state_id, ever_treated, first_treat_year, row_index, y_plot),
    by = "state_id"
  ) %>%
  arrange(row_index, year, crime)

# ===============================================================
# 4. Build Y and original D tensors
# ===============================================================

Y = array(
  NA_real_,
  dim = c(n_states, n_years, n_crimes),
  dimnames = list(
    state = state_order_vec,
    year = as.character(year_order),
    crime = crime_order
  )
)

D_original = Y

for (cc in seq_along(crime_order)) {

  this_crime = crime_order[cc]

  tmp = castle_tensor_long %>%
    filter(as.character(crime) == this_crime) %>%
    mutate(
      state_id = factor(state_id, levels = state_order_vec),
      year = factor(year, levels = year_order)
    ) %>%
    arrange(state_id, year)

  Y[, , cc] = tmp %>%
    select(state_id, year, outcome_value) %>%
    pivot_wider(names_from = year, values_from = outcome_value) %>%
    arrange(state_id) %>%
    select(-state_id) %>%
    as.matrix()

  D_original[, , cc] = tmp %>%
    select(state_id, year, treatment_value) %>%
    pivot_wider(names_from = year, values_from = treatment_value) %>%
    arrange(state_id) %>%
    select(-state_id) %>%
    as.matrix()
}

storage.mode(Y) = "numeric"
storage.mode(D_original) = "numeric"

rn = dimnames(Y)$state
cn = dimnames(Y)$year
ln = dimnames(Y)$crime

target_layer_index = match(target_layer, ln)
if (is.na(target_layer_index)) {
  stop("target_layer is not present in dimnames(Y)$crime.")
}

target_crimes = setNames(target_layer, target_name)
target_label = unname(crime_labels[target_layer])

cat("\nCreated tensors:\n")
cat("Y:", paste(dim(Y), collapse = " x "), "\n")
cat("D_original:", paste(dim(D_original), collapse = " x "), "\n")
cat("Dimension order: state x year x crime\n")
cat("Target layer:", target_layer, "\n")
cat("Rank:", rank_value, "\n")
cat("Bootstrap B:", B, "\n")

# ===============================================================
# 5. Source functions
# ===============================================================

source("../../bilinearMatrixStaggered.R")
source("../../bilinearTensorStaggered.R")
source("../../bilinearMatrixStaggeredPsi.R")
source("../../bilinearTensorStaggeredPsi.R")
source("../../pluginPsi_c1.R")

# ===============================================================
# 6. Helper functions for masks and staircase parts
# ===============================================================

make_A_Omega_from_D = function(D) {

  N = dim(D)[1]
  Tt = dim(D)[2]
  K = dim(D)[3]

  A = matrix(Inf, nrow = N, ncol = K)
  Omega = array(FALSE, dim = c(N, Tt, K))

  for (k in seq_len(K)) {
    for (i in seq_len(N)) {

      on_idx = which(D[i, , k] == 1)

      if (length(on_idx) == 0) {
        A[i, k] = Inf
        Omega[i, , k] = TRUE
      } else {
        A[i, k] = min(on_idx)

        if (A[i, k] > 1) {
          Omega[i, seq_len(A[i, k] - 1), k] = TRUE
        }
      }
    }
  }

  list(A = A, Omega = Omega)
}

make_staggered_parts_from_D = function(D, layer = 1) {

  Dmat = D[, , layer]
  years = as.integer(colnames(Dmat))

  first_on = apply(Dmat, 1, function(drow) {
    on_years = years[drow == 1]
    if (length(on_years) == 0) Inf else min(on_years)
  })

  adopt_years = sort(unique(first_on[is.finite(first_on)]))
  m = length(adopt_years)

  if (m == 0) {
    return(
      list(
        N_sizes = nrow(Dmat),
        T_sizes = ncol(Dmat),
        row_groups = list(seq_len(nrow(Dmat))),
        col_groups = list(seq_len(ncol(Dmat))),
        adopt_years = adopt_years,
        first_on = first_on
      )
    )
  }

  o = m + 1

  row_groups = vector("list", o)
  row_groups[[1]] = which(!is.finite(first_on))

  for (ell in seq_len(m)) {
    row_groups[[ell + 1]] = which(first_on == rev(adopt_years)[ell])
  }

  col_groups = vector("list", o)
  col_groups[[1]] = which(years < adopt_years[1])

  if (m >= 2) {
    for (ell in 2:m) {
      col_groups[[ell]] = which(
        years >= adopt_years[ell - 1] &
          years < adopt_years[ell]
      )
    }
  }

  col_groups[[o]] = which(years >= adopt_years[m])

  list(
    N_sizes = lengths(row_groups),
    T_sizes = lengths(col_groups),
    row_groups = row_groups,
    col_groups = col_groups,
    adopt_years = adopt_years,
    first_on = first_on
  )
}

make_masked_target_D = function(
    D_base,
    target_layer,
    n_fully_observed_rows,
    n_observed_years_for_other_rows
) {

  D_out = D_base

  target_idx = match(target_layer, dimnames(D_out)$crime)
  if (is.na(target_idx)) {
    stop("target_layer is not present in dimnames(D_base)$crime.")
  }

  N = dim(D_out)[1]
  Tt = dim(D_out)[2]

  full_rows = seq_len(min(n_fully_observed_rows, N))
  short_rows = setdiff(seq_len(N), full_rows)

  early_years = seq_len(min(n_observed_years_for_other_rows, Tt))
  late_years = setdiff(seq_len(Tt), early_years)

  # Overwrite only the target layer's analysis mask
  D_out[, , target_idx] = 0

  if (length(short_rows) > 0 && length(late_years) > 0) {
    D_out[short_rows, late_years, target_idx] = 1
  }

  attr(D_out, "mask_info") = list(
    design_type = "masked",
    target_layer = target_layer,
    target_layer_index = target_idx,
    n_fully_observed_rows = n_fully_observed_rows,
    n_observed_years_for_other_rows = n_observed_years_for_other_rows,
    full_rows = full_rows,
    short_rows = short_rows,
    early_years = early_years,
    late_years = late_years
  )

  D_out
}

make_design_D = function(design_row, D_original, target_layer) {

  design_type = design_row$design_type[[1]]

  if (identical(design_type, "masked")) {
    return(
      make_masked_target_D(
        D_base = D_original,
        target_layer = target_layer,
        n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
        n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]]
      )
    )
  }

  if (identical(design_type, "nonmasked")) {
    D_out = D_original
    attr(D_out, "mask_info") = list(
      design_type = "nonmasked",
      target_layer = target_layer,
      target_layer_index = match(target_layer, dimnames(D_original)$crime),
      n_fully_observed_rows = NA_integer_,
      n_observed_years_for_other_rows = NA_integer_,
      full_rows = integer(),
      short_rows = integer(),
      early_years = integer(),
      late_years = integer()
    )
    return(D_out)
  }

  stop("Unknown design_type: ", design_type)
}

is_staircase_layer = function(Omega_k) {

  N = nrow(Omega_k)
  Tt = ncol(Omega_k)

  obs_len = integer(N)

  for (i in seq_len(N)) {
    obs = Omega_k[i, ]

    if (any(diff(as.integer(obs)) > 0)) {
      return(FALSE)
    }

    first_unobs = which(!obs)[1]

    if (is.na(first_unobs)) {
      obs_len[i] = Tt
    } else {
      obs_len[i] = first_unobs - 1
    }
  }

  all(diff(obs_len) <= 0)
}

# ===============================================================
# 7. Estimation function, reusable for bootstrap samples
# ===============================================================

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

run_estimator_given_arrays = function(
    Y_in,
    D_in,
    rn_in,
    eta_in,
    rank_value,
    tau,
    target_crimes,
    local_states
) {

  Y0_in = Y_in
  Y0_in[D_in == 1] = NA_real_

  AO_in = make_A_Omega_from_D(D_in)

  A_in = AO_in$A
  Omega_in = AO_in$Omega & !is.na(Y0_in)

  staggered_layer_in = match(target_crimes[[1]], dimnames(Y_in)$crime)

  if (is.na(staggered_layer_in)) {
    stop("The requested target layer is not present in dimnames(Y_in)$crime.")
  }

  staggered_in = make_staggered_parts_from_D(
    D = D_in,
    layer = staggered_layer_in
  )

  K_in = dim(Y_in)[3]

  N_parts_in = replicate(K_in, staggered_in$N_sizes, simplify = FALSE)
  T_parts_in = replicate(K_in, staggered_in$T_sizes, simplify = FALSE)

  estimate_one_inner = function(
      crime_name,
      functional,
      local_state = NA_character_
  ) {

    k = match(target_crimes[[crime_name]], dimnames(Y_in)$crime)

    row_index = if (functional == "Local") {
      match(local_state, rn_in)
    } else {
      NULL
    }

    functional_label = if (functional == "Local") {
      paste0("Local-", local_state)
    } else {
      functional
    }

    if (functional == "Local" && is.na(row_index)) {
      return(
        tibble(
          r = rank_value,
          crime = crime_name,
          layer = target_crimes[[crime_name]],
          functional = functional_label,
          local_state = local_state,
          Psi0 = NA_real_,
          Psi1 = NA_real_,
          Delta_h = NA_real_,
          Psi0_matrix = NA_real_,
          Delta_h_matrix = NA_real_,
          tensor_minus_matrix = NA_real_
        )
      )
    }

    psi0_tensor = bilinearTensorStaggeredPsi(
      Y = Y0_in,
      k = k,
      r = rank_value,
      tau = tau,
      functional = functional,
      eta = if (functional == "RowHet") eta_in else NULL,
      row_index = row_index,
      A = A_in,
      Omega = Omega_in
    )

    psi0_matrix = bilinearMatrixStaggeredPsi(
      Y_mat = Y0_in[, , k],
      r = rank_value,
      tau = tau,
      functional = functional,
      eta = if (functional == "RowHet") eta_in else NULL,
      row_index = row_index,
      A = A_in[, k],
      Omega = Omega_in[, , k]
    )

    psi1 = pluginPsi_c1(
      Y = Y_in,
      k = k,
      N_parts = N_parts_in,
      T_parts = T_parts_in,
      functional = functional,
      eta = if (functional == "RowHet") eta_in else NULL,
      row_index = row_index
    )

    tibble(
      r = rank_value,
      crime = crime_name,
      layer = target_crimes[[crime_name]],
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

  nonlocal_results = map_dfr(
    names(target_crimes),
    function(crime_name) {
      map_dfr(
        c("ATE", "RowHet", "Trend"),
        function(f) {
          estimate_one_inner(
            crime_name = crime_name,
            functional = f
          )
        }
      )
    }
  )

  local_results = map_dfr(
    names(target_crimes),
    function(crime_name) {
      map_dfr(
        local_states,
        function(st) {
          estimate_one_inner(
            crime_name = crime_name,
            functional = "Local",
            local_state = st
          )
        }
      )
    }
  )

  bind_rows(nonlocal_results, local_results) %>%
    mutate(
      functional = factor(functional, levels = functional_order)
    ) %>%
    arrange(functional, crime) %>%
    mutate(functional = as.character(functional))
}

# ===============================================================
# 8. Target-layer-only bootstrap helpers
# ===============================================================

get_target_row_strata = function(Y, D, target_layer) {

  target_idx = match(target_layer, dimnames(D)$crime)

  if (is.na(target_idx)) {
    stop("target_layer is not present in dimnames(D)$crime.")
  }

  Y_target = Y[, , target_idx]
  D_target = D[, , target_idx]

  fully_observed_idx = which(
    rowSums(D_target == 1, na.rm = TRUE) == 0 &
      rowSums(is.na(Y_target)) == 0
  )

  non_fully_observed_idx = setdiff(
    seq_len(dim(D)[1]),
    fully_observed_idx
  )

  list(
    target_layer_index = target_idx,
    fully_observed_idx = fully_observed_idx,
    non_fully_observed_idx = non_fully_observed_idx
  )
}

make_target_layer_row_resample_index = function(
    rn,
    Y,
    D,
    target_layer,
    local_states = c("Montana", "Texas", "Florida"),
    keep_local_target_rows = FALSE,
    stratify_target_bootstrap = TRUE
) {

  strata = get_target_row_strata(
    Y = Y,
    D = D,
    target_layer = target_layer
  )
  
  # we sample fully observed rows from fully observed rows
  # and non-fully observed rows from non-fully observed rows
  full_idx = strata$fully_observed_idx
  nonfull_idx = strata$non_fully_observed_idx

  N = length(rn)
  idx = seq_len(N)

  if (
    isTRUE(stratify_target_bootstrap) &&
      length(full_idx) > 0 &&
      length(nonfull_idx) > 0
  ) {

    idx[full_idx] = sample(
      full_idx,
      size = length(full_idx),
      replace = TRUE
    )

    idx[nonfull_idx] = sample(
      nonfull_idx,
      size = length(nonfull_idx),
      replace = TRUE
    )

  } else {

    idx = sample(
      seq_len(N),
      size = N,
      replace = TRUE
    )
  }

  if (keep_local_target_rows) {
    local_idx = match(local_states, rn)
    local_idx = local_idx[!is.na(local_idx)]
    idx[local_idx] = local_idx
  }

  idx
}

make_failed_bootstrap_rows = function(point_results, b, bootstrap_type, error_message) {

  point_results %>%
    mutate(
      across(all_of(ci_vars), ~ NA_real_),
      bootstrap_id = b,
      bootstrap_type = bootstrap_type,
      bootstrap_error = error_message
    )
}

run_target_layer_only_bootstrap_one = function(
    b,
    Y_base,
    D_design,
    seed_base,
    point_results
) {

  set.seed(seed_base + b)

  target_idx = match(target_layer, dimnames(Y_base)$crime)

  if (is.na(target_idx)) {
    stop("target_layer is not present in dimnames(Y_base)$crime.")
  }

  other_layers = setdiff(seq_len(dim(Y_base)[3]), target_idx)

  idx = make_target_layer_row_resample_index(
    rn = rn,
    Y = Y_base,
    D = D_design,
    target_layer = target_layer,
    local_states = local_states,
    keep_local_target_rows = keep_local_target_rows,
    stratify_target_bootstrap = stratify_target_bootstrap
  )

  Y_b = Y_base
  D_b = D_design

  # Resample only the target layer.
  Y_b[, , target_idx] = Y_base[idx, , target_idx]

  # Keep anchor layers fixed.
  if (length(other_layers) > 0) {
    Y_b[, , other_layers] = Y_base[, , other_layers]
  }

  dimnames(Y_b) = dimnames(Y_base)
  dimnames(D_b) = dimnames(D_design)

  if (!identical(dimnames(Y_b), dimnames(Y_base))) {
    stop("Y_b dimnames changed.")
  }

  if (!identical(D_b, D_design)) {
    stop("D_b changed, but this bootstrap should keep D fixed.")
  }

  if (length(other_layers) > 0) {
    anchors_equal = isTRUE(
      all.equal(
        Y_b[, , other_layers, drop = FALSE],
        Y_base[, , other_layers, drop = FALSE],
        check.attributes = TRUE
      )
    )

    if (!anchors_equal) {
      stop("Anchor layers changed, but they should be fixed.")
    }
  }

  run_estimator_given_arrays(
    Y_in = Y_b,
    D_in = D_b,
    rn_in = rn,
    eta_in = eta,
    rank_value = rank_value,
    tau = tau,
    target_crimes = target_crimes,
    local_states = local_states
  ) %>%
    mutate(
      bootstrap_id = b,
      bootstrap_type = "target_layer_only_fixed_anchors",
      bootstrap_error = NA_character_
    )
}

# ===============================================================
# 9. Confidence interval helper
#   CI = point_estimate +/- 1.96 * bootstrap_se
# ===============================================================

make_results_with_ci = function(point_results, bootstrap_results) {
  
  bootstrap_se_long = bootstrap_results %>%
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
  
  point_long = point_results %>%
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
    )
  
  point_long %>%
    left_join(
      bootstrap_se_long,
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
      point_estimate = round(point_estimate, 6),
      boot_mean = round(boot_mean, 6),
      boot_se = round(boot_se, 6),
      ci_low = round(ci_low, 6),
      ci_high = round(ci_high, 6),
      estimate_ci = if_else(
        is.na(ci_low) | is.na(ci_high),
        sprintf("%.4f [NA, NA]", point_estimate),
        sprintf("%.4f [%.4f, %.4f]", point_estimate, ci_low, ci_high)
      )
    )
}

# ===============================================================
# 10. Design grid
# ===============================================================

design_settings = bind_rows(
  tibble(
    design_id = paste0("masked_full", n_fully_observed_grid, "_obs", n_observed_years_for_other_rows),
    design_label = paste0(
      "Masked: ", n_fully_observed_grid,
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

print(design_settings)

# ===============================================================
# 11. Run one design
# ===============================================================

run_one_design = function(design_row) {

  design_id = design_row$design_id[[1]]
  design_label = design_row$design_label[[1]]

  cat("\n============================================================\n")
  cat("Running design:", design_id, "\n")
  cat(design_label, "\n")
  cat("============================================================\n")

  D_design = make_design_D(
    design_row = design_row,
    D_original = D_original,
    target_layer = target_layer
  )

  Y0_design = Y
  Y0_design[D_design == 1] = NA_real_

  AO_design = make_A_Omega_from_D(D_design)
  Omega_design = AO_design$Omega & !is.na(Y0_design)

  staircase_tbl = tibble(
    layer = dimnames(Y)$crime,
    staircase = map_lgl(seq_len(dim(Y)[3]), ~ is_staircase_layer(Omega_design[, , .x]))
  )

  cat("\nStaircase checks for Omega by layer:\n")
  print(staircase_tbl)

  strata = get_target_row_strata(
    Y = Y,
    D = D_design,
    target_layer = target_layer
  )

  cat("\nTarget-layer bootstrap strata:\n")
  cat("Fully observed target-layer rows:", length(strata$fully_observed_idx), "\n")
  cat("Non-fully-observed target-layer rows:", length(strata$non_fully_observed_idx), "\n")

  point_results = run_estimator_given_arrays(
    Y_in = Y,
    D_in = D_design,
    rn_in = rn,
    eta_in = eta,
    rank_value = rank_value,
    tau = tau,
    target_crimes = target_crimes,
    local_states = local_states
  ) %>%
    mutate(
      design_id = design_id,
      design_label = design_label,
      design_type = design_row$design_type[[1]],
      n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
      n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]],
      .before = 1
    ) %>%
    mutate(across(all_of(ci_vars), ~ round(.x, 6)))

  cat("\nPoint estimates:\n")
  print(point_results, n = Inf)

  cat("\nRunning target-layer-only bootstrap for design ", design_id, "\n", sep = "")
  cat("B =", B, "\n")

  bootstrap_results = map_dfr(
    seq_len(B),
    function(b) {

      if (b %% 25 == 0) {
        cat("Design", design_id, "| bootstrap replication", b, "of", B, "\n")
      }

      tryCatch(
        {
          run_target_layer_only_bootstrap_one(
            b = b,
            Y_base = Y,
            D_design = D_design,
            seed_base = bootstrap_seed,
            point_results = point_results
          )
        },
        error = function(e) {
          warning(
            paste0(
              "Bootstrap replication ", b,
              " failed for design ", design_id,
              ": ", conditionMessage(e)
            )
          )

          make_failed_bootstrap_rows(
            point_results = point_results,
            b = b,
            bootstrap_type = "target_layer_only_fixed_anchors",
            error_message = conditionMessage(e)
          )
        }
      )
    }
  ) %>%
    mutate(
      design_id = design_id,
      design_label = design_label,
      design_type = design_row$design_type[[1]],
      n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
      n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]],
      .before = 1
    )

  results_with_ci = make_results_with_ci(
    point_results = point_results,
    bootstrap_results = bootstrap_results
  ) %>%
    mutate(
      design_id = design_id,
      design_label = design_label,
      design_type = design_row$design_type[[1]],
      n_fully_observed_rows = design_row$n_fully_observed_rows[[1]],
      n_observed_years_for_other_rows = design_row$n_observed_years_for_other_rows[[1]],
      .before = 1
    )

  design_object = list(
    design_settings = design_row,
    point_results = point_results,
    bootstrap_results = bootstrap_results,
    results_with_ci = results_with_ci,
    staircase_checks = staircase_tbl,
    target_row_strata = strata,
    D_design = D_design
  )

  saveRDS(
    design_object,
    file = file.path(results_dir, paste0("design_", target_file_stub, "_", design_id, ".rds"))
  )
  
  write_csv(
    point_results,
    file = file.path(results_dir, paste0("point_results_", target_file_stub, "_", design_id, ".csv"))
  )
  
  write_csv(
    bootstrap_results,
    file = file.path(results_dir, paste0("bootstrap_results_", target_file_stub, "_", design_id, ".csv"))
  )
  
  write_csv(
    results_with_ci,
    file = file.path(results_dir, paste0("results_with_ci_", target_file_stub, "_", design_id, ".csv"))
  )

  cat("\nFinished design:", design_id, "\n")

  design_object
}

# ===============================================================
# 12. Run all four designs and save combined results
# ===============================================================

design_objects = map(
  seq_len(nrow(design_settings)),
  function(ii) run_one_design(design_settings[ii, ])
)

names(design_objects) = design_settings$design_id

point_results_all = map_dfr(design_objects, "point_results")
bootstrap_results_all = map_dfr(design_objects, "bootstrap_results")
results_with_ci_all = map_dfr(design_objects, "results_with_ci")
staircase_checks_all = map2_dfr(
  design_objects,
  names(design_objects),
  function(obj, did) {
    obj$staircase_checks %>%
      mutate(design_id = did, .before = 1)
  }
)

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
    crime_labels = crime_labels,
    project_dir = getwd()
  ),
  design_settings = design_settings,
  point_results = point_results_all,
  bootstrap_results = bootstrap_results_all,
  results_with_ci = results_with_ci_all,
  staircase_checks = staircase_checks_all,
  design_objects = design_objects,
  Y = Y,
  D_original = D_original,
  state_order = state_order,
  year_order = year_order,
  crime_order = crime_order
)

combined_rds = file.path(
  results_dir,
  paste0("castle_masked_nonmasked_bootstrap_results_", target_file_stub, ".rds")
)

saveRDS(analysis_results, file = combined_rds)

write_csv(
  design_settings,
  file.path(results_dir, paste0("design_settings_", target_file_stub, ".csv"))
)

write_csv(
  point_results_all,
  file.path(results_dir, paste0("point_results_all_designs_", target_file_stub, ".csv"))
)

write_csv(
  bootstrap_results_all,
  file.path(results_dir, paste0("bootstrap_results_all_designs_", target_file_stub, ".csv"))
)

write_csv(
  results_with_ci_all,
  file.path(results_dir, paste0("results_with_ci_all_designs_", target_file_stub, ".csv"))
)

write_csv(
  staircase_checks_all,
  file.path(results_dir, paste0("staircase_checks_all_designs_", target_file_stub, ".csv"))
)
cat("\n============================================================\n")
cat("All designs finished.\n")
cat("Saved combined RDS to:\n")
cat(combined_rds, "\n")
cat("Saved CSV files to:\n")
cat(results_dir, "\n")
cat("============================================================\n")
