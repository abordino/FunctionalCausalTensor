setwd("~/Desktop/code")

suppressPackageStartupMessages(library(tidyverse))

source("bilinearTensorAllFunction.R")
source("CYNuclear.R")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

target_layer = "l_motor"
target_name = "motor"

rank_value = 3
tau = 0.01

cy_seed = 123

local_states = c("Montana", "Texas", "Florida")

functional_order = c(
  "ATE",
  "Local-Florida",
  "Local-Montana",
  "Local-Texas",
  "RowHet",
  "Trend"
)

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
  D_original[, , k] = make_layer_matrix(
    castle_long,
    crime_order[k],
    "treatment_value"
  )
}

storage.mode(Y) = "numeric"
storage.mode(D_original) = "numeric"

rn = dimnames(Y)$state
target_label = unname(crime_labels[target_layer])

# -----------------------------------------------------------------------------
# Estimator helpers
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
        
        if (first_on > 1) {
          Omega[i, seq_len(first_on - 1), k] = TRUE
        }
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
    
    if (length(treated_years) == 0) {
      Inf
    } else {
      min(treated_years)
    }
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
    row_groups[[group + 1]] = which(
      first_on == rev(adopt_years)[group]
    )
  }
  
  col_groups = vector("list", m + 1)
  col_groups[[1]] = which(years < adopt_years[1])
  
  if (m >= 2) {
    for (group in 2:m) {
      col_groups[[group]] = which(
        years >= adopt_years[group - 1] &
          years < adopt_years[group]
      )
    }
  }
  
  col_groups[[m + 1]] = which(years >= adopt_years[m])
  
  list(
    N_sizes = lengths(row_groups),
    T_sizes = lengths(col_groups)
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
  
  parts = make_staggered_parts_from_D(
    D = D_in,
    layer = target_idx
  )
  
  N_parts = replicate(
    dim(Y_in)[3],
    parts$N_sizes,
    simplify = FALSE
  )
  
  T_parts = replicate(
    dim(Y_in)[3],
    parts$T_sizes,
    simplify = FALSE
  )
  
  cy_fit = cy_complete_staggered_layer(
    y = Y0[, , target_idx],
    omega = Omega[, , target_idx],
    adoption_time = A[, target_idx],
    time = seq_len(dim(Y_in)[2]),
    rank = 1, # rank_value,
    seed = cy_seed,
    verbose = FALSE
  )
  
  estimate_one = function(functional, local_state = NA_character_) {
    row_index = if (functional == "Local") {
      match(local_state, rn)
    } else {
      NULL
    }
    
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
        tensor_minus_matrix = NA_real_,
        Psi0_CY_nuclear = NA_real_,
        Delta_h_CY_nuclear = NA_real_,
        Psi0_CY_debiased = NA_real_,
        Delta_h_CY_debiased = NA_real_
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
    
    psi0_CY_nuclear = pluginPsi_c1(
      Y = cy_fit$completed_nuclear,
      k = target_idx,
      N_parts = N_parts,
      T_parts = T_parts,
      functional = functional,
      eta = eta_arg,
      row_index = row_index
    )
    
    psi0_CY_debiased = pluginPsi_c1(
      Y = cy_fit$completed_debiased,
      k = target_idx,
      N_parts = N_parts,
      T_parts = T_parts,
      functional = functional,
      eta = eta_arg,
      row_index = row_index
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
      local_state = if (functional == "Local") {
        local_state
      } else {
        NA_character_
      },
      Psi0 = as.numeric(psi0_tensor),
      Psi1 = as.numeric(psi1),
      Delta_h = as.numeric(psi1 - psi0_tensor),
      Psi0_matrix = as.numeric(psi0_matrix),
      Delta_h_matrix = as.numeric(psi1 - psi0_matrix),
      tensor_minus_matrix = as.numeric(
        psi0_tensor - psi0_matrix
      ),
      Psi0_CY_nuclear = as.numeric(psi0_CY_nuclear),
      Delta_h_CY_nuclear = as.numeric(
        psi1 - psi0_CY_nuclear
      ),
      Psi0_CY_debiased = as.numeric(psi0_CY_debiased),
      Delta_h_CY_debiased = as.numeric(
        psi1 - psi0_CY_debiased
      )
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
    mutate(
      functional = factor(
        functional,
        levels = functional_order
      )
    ) %>%
    arrange(functional) %>%
    mutate(functional = as.character(functional))
}

# -----------------------------------------------------------------------------
# Point estimation
# -----------------------------------------------------------------------------

point_results = run_estimator(
  Y_in = Y,
  D_in = D_original
) %>%
  mutate(
    across(
      all_of(result_columns),
      ~ round(.x, 6)
    )
  )

# -----------------------------------------------------------------------------
# Final table
# -----------------------------------------------------------------------------

main_summary = point_results %>%
  select(
    r,
    crime,
    layer,
    functional,
    local_state,
    Psi0,
    Psi0_matrix,
    Psi0_CY_nuclear,
    Psi0_CY_debiased,
    Psi1,
    Delta_h,
    Delta_h_matrix,
    Delta_h_CY_nuclear,
    Delta_h_CY_debiased,
    tensor_minus_matrix
  ) %>%
  arrange(functional)

print(main_summary, n = Inf, width = Inf)