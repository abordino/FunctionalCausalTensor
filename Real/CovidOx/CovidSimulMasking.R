# Setup -----------------------------------------------------------------------

setwd("~/Desktop/code")

suppressPackageStartupMessages(library(tidyverse))

# Settings --------------------------------------------------------------------

RUN_CY = TRUE

source("bilinearTensorAllFunction.R")
if (RUN_CY) {
  source("CYNuclear.R")
}

JUMP = 30L
STEP_SIZE = 1L
REP = 10L
REP_SEED = 20260805L

ESTIMATION_RANKS = c(1L, 2L, 3L, 4L)

tensor_file =
  "Real/CovidOx/data/Omega_Y_until_2020-04-05_delay_28.rds"

n_outcomes_to_keep = 2L
target_outcome_pattern = "death"
target_name = "deaths"

tau = 1e-4
cy_seed = 123L

N_FULL_ROWS = 2L 
N_FULL_COLUMNS = 3L # 3L

SAVE_CSV_OUTPUTS = FALSE

DATASET_TAG = "oxford_deaths"
RESULTS_STUDY_TAG = "oxford_deaths_target_layer_bootstrap"

cy_file_tag = if (RUN_CY) "" else "_noCY"
rep_file_tag = if (REP == 1L) "" else paste0("_rep", REP)

# COVID-19 Oxford data --------------------------------------------------------

`%||%` = function(x, y) {
  if (is.null(x)) y else x
}

layer_names = function(x, n, prefix) {
  if (is.null(x)) {
    paste0(prefix, "_", seq_len(n))
  } else {
    as.character(x[seq_len(n)])
  }
}

first_true = function(x) {
  index = which(x)
  if (length(index)) index[1] else Inf
}

fill_time_neighbours = function(Y, max_passes = 5L) {
  for (pass in seq_len(max_passes)) {
    missing = which(is.na(Y), arr.ind = TRUE)
    
    if (!nrow(missing)) {
      break
    }
    
    missing_before = nrow(missing)
    
    for (row in seq_len(nrow(missing))) {
      i = missing[row, 1]
      t = missing[row, 2]
      k = missing[row, 3]
      
      neighbours = c(
        if (t > 1L) Y[i, t - 1L, k] else NA_real_,
        if (t < dim(Y)[2]) Y[i, t + 1L, k] else NA_real_
      )
      
      neighbours = neighbours[!is.na(neighbours)]
      
      if (length(neighbours)) {
        Y[i, t, k] = mean(neighbours)
      }
    }
    
    if (sum(is.na(Y)) == missing_before) {
      break
    }
  }
  
  Y
}

align_layers = function(data, target_pattern, n_keep) {
  Y = data$Y
  policy = data$Omega
  
  storage.mode(Y) = "numeric"
  storage.mode(policy) = "logical"
  
  countries =
    data$country_levels %||%
    paste0("country_", seq_len(dim(Y)[1]))
  
  times =
    as.character(
      data$date_names %||%
        seq_len(dim(Y)[2])
    )
  
  outcomes = layer_names(
    data$outcomes,
    dim(Y)[3],
    "outcome"
  )
  
  policies = layer_names(
    data$policies,
    dim(policy)[3],
    "policy"
  )
  
  dimnames(Y) = list(
    country = countries,
    time = times,
    layer = outcomes
  )
  
  dimnames(policy) = list(
    country = countries,
    time = times,
    layer = policies
  )
  
  map = data$policy_outcome_map
  
  if (
    is.null(map) ||
    !all(c("policy", "outcome") %in% names(map))
  ) {
    map = tibble(
      outcome = outcomes,
      policy = policies,
      panel_label = outcomes
    )
  } else {
    map = as_tibble(map) %>%
      transmute(
        outcome = as.character(outcome),
        policy = as.character(policy),
        panel_label = if (
          "panel_label" %in% names(map)
        ) {
          as.character(panel_label)
        } else {
          as.character(outcome)
        }
      ) %>%
      filter(
        outcome %in% outcomes,
        policy %in% policies
      )
  }
  
  map = map %>%
    mutate(
      outcome_index = match(outcome, outcomes),
      policy_index = match(policy, policies)
    ) %>%
    arrange(outcome_index)
  
  target = map %>%
    filter(
      str_detect(
        str_to_lower(outcome),
        target_pattern
      ) |
        str_detect(
          str_to_lower(panel_label),
          target_pattern
        )
    ) %>%
    mutate(
      exact =
        !str_to_lower(outcome) %in%
        c("death", "deaths")
    ) %>%
    arrange(exact, outcome_index) %>%
    slice(1)
  
  if (!nrow(target)) {
    stop(
      "No target outcome matched: ",
      target_pattern
    )
  }
  
  selected = if (nrow(map) == n_keep) {
    map
  } else {
    bind_rows(
      target,
      map %>%
        filter(
          outcome != target$outcome[[1]]
        ) %>%
        slice_head(n = n_keep - 1L)
    ) %>%
      distinct(outcome, .keep_all = TRUE) %>%
      arrange(outcome_index)
  }
  
  if (nrow(selected) != n_keep) {
    stop(
      "Could not select the requested outcome layers."
    )
  }
  
  Y = Y[
    ,
    ,
    selected$outcome_index,
    drop = FALSE
  ]
  
  policy = policy[
    ,
    ,
    selected$policy_index,
    drop = FALSE
  ]
  
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
  
  target_idx = match(
    target$outcome[[1]],
    selected$outcome
  )
  
  if (is.na(target_idx)) {
    stop(
      "The target outcome was not retained."
    )
  }
  
  Y = fill_time_neighbours(Y)
  
  if (anyNA(Y)) {
    stop(
      "Y still contains missing values after neighbour filling."
    )
  }
  
  list(
    Y = Y,
    policy = policy,
    layer_map = layer_map,
    target_idx = target_idx,
    target_outcome =
      selected$outcome[[target_idx]],
    target_policy =
      selected$policy[[target_idx]]
  )
}

prepare_covid_data = function(tensor_file) {
  tensor_data = readRDS(tensor_file)
  
  prepared = align_layers(
    data = tensor_data,
    target_pattern = target_outcome_pattern,
    n_keep = n_outcomes_to_keep
  )
  
  Y = prepared$Y
  policy = prepared$policy
  layer_map = prepared$layer_map
  target_idx = prepared$target_idx
  
  D_original = array(
    as.numeric(policy),
    dim = dim(policy),
    dimnames = dimnames(policy)
  )
  
  storage.mode(D_original) = "numeric"
  
  country_order = dimnames(Y)[[1]]
  time_order = dimnames(Y)[[2]]
  layer_order = dimnames(Y)[[3]]
  
  target_policy_matrix =
    policy[, , target_idx, drop = TRUE]
  
  first_treat_index = apply(
    target_policy_matrix,
    1,
    first_true
  )
  
  first_treat_date = rep(
    NA_character_,
    length(first_treat_index)
  )
  
  treated_rows = which(
    is.finite(first_treat_index)
  )
  
  if (length(treated_rows)) {
    first_treat_date[treated_rows] =
      time_order[
        as.integer(
          first_treat_index[treated_rows]
        )
      ]
  }
  
  state_order = tibble(
    country_id = country_order,
    ever_treated = is.finite(first_treat_index),
    first_treat_index = if_else(
      is.finite(first_treat_index),
      as.integer(first_treat_index),
      NA_integer_
    ),
    first_treat_date = first_treat_date,
    row_index = seq_along(country_order)
  )
  
  layer_labels = setNames(
    layer_map$panel_label,
    layer_map$outcome
  )
  
  target_label = unname(
    layer_labels[prepared$target_outcome]
  )
  
  if (
    length(target_label) == 0L ||
    is.na(target_label)
  ) {
    target_label = prepared$target_outcome
  }
  
  slice_label = paste0(
    "until_",
    tensor_data$time_horizon %||% "unknown",
    "_delay_",
    tensor_data$delay_days %||% "unknown"
  )
  
  slice_label = str_replace_all(
    slice_label,
    "[^A-Za-z0-9._-]+",
    "_"
  )
  
  list(
    tensor_data = tensor_data,
    Y = Y,
    policy = policy,
    D_original = D_original,
    layer_map = layer_map,
    state_order = state_order,
    state_order_vec = country_order,
    year_order = time_order,
    crime_order = layer_order,
    crime_vars = layer_order,
    crime_labels = layer_labels,
    target_idx = target_idx,
    target_layer = prepared$target_outcome,
    target_policy = prepared$target_policy,
    target_label = target_label,
    slice_label = slice_label
  )
}

covid_data = prepare_covid_data(
  tensor_file
)

target_layer = covid_data$target_layer
slice_label = covid_data$slice_label
slice_file_tag = paste0("_", slice_label)

N_UNITS = dim(covid_data$Y)[1]
N_PERIODS = dim(covid_data$Y)[2]
N_LAYERS = dim(covid_data$Y)[3]

SIMULATION_GRID = crossing(
  N_LAYERS = as.integer(N_LAYERS),
  rank_value = ESTIMATION_RANKS
)

# Simulation ------------------------------------------------------------------

run_simulation = function(N_LAYERS, rank_value) {
  N_LAYERS = as.integer(N_LAYERS)
  rank_value = as.integer(rank_value)
  
  results_dir = file.path(
    "Results",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
  )
  
  dir.create(
    results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  Y = covid_data$Y
  policy = covid_data$policy
  D_original = covid_data$D_original
  layer_map = covid_data$layer_map
  state_order = covid_data$state_order
  year_order = covid_data$year_order
  crime_order = covid_data$crime_order
  crime_vars = covid_data$crime_vars
  crime_labels = covid_data$crime_labels
  target_idx = covid_data$target_idx
  target_label = covid_data$target_label
  target_policy = covid_data$target_policy
  
  if (N_LAYERS != dim(Y)[3]) {
    stop(
      "N_LAYERS does not match the selected COVID tensor."
    )
  }
  
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
          
          if (first_on > 1L) {
            Omega[i, seq_len(first_on - 1L), k] = TRUE
          }
        }
      }
    }
    
    list(A = A, Omega = Omega)
  }
  
  observed_lengths_from_mask = function(M) {
    apply(M, 1, function(row) {
      first_missing = which(row)[1]
      if (is.na(first_missing)) ncol(M) else first_missing - 1L
    })
  }
  
  mask_from_observed_lengths = function(observed_lengths, Tt, dimnames_value) {
    M = outer(
      observed_lengths,
      seq_len(Tt),
      function(observed_length, time_index) time_index > observed_length
    )
    dimnames(M) = dimnames_value
    M
  }
  
  make_mask_sequence = function(M0, step_size, jump, rep_id, rep_seed) {
    N = nrow(M0)
    Tt = ncol(M0)
    
    original_observed_lengths =
      unname(observed_lengths_from_mask(M0))
    
    if (
      N_FULL_COLUMNS < 1L ||
      N_FULL_COLUMNS > Tt
    ) {
      stop(
        "N_FULL_COLUMNS must be between 1 and the number of periods."
      )
    }
    
    if (
      any(
        original_observed_lengths <
        N_FULL_COLUMNS
      )
    ) {
      stop(
        "N_FULL_COLUMNS exceeds the original observed length for at least one COVID country."
      )
    }
    
    full_rows = head(
      which(original_observed_lengths == Tt),
      N_FULL_ROWS
    )
    
    if (length(full_rows) < N_FULL_ROWS) {
      stop("Not enough fully observed rows.")
    }
    
    terminal_observed_lengths = rep.int(N_FULL_COLUMNS, N)
    terminal_observed_lengths[full_rows] = Tt
    
    n_target_operations = sum(
      original_observed_lengths - terminal_observed_lengths
    )
    
    if (jump > n_target_operations) {
      stop("JUMP exceeds the number of maskable cells.")
    }
    
    terminal_mask = mask_from_observed_lengths(
      terminal_observed_lengths,
      Tt,
      dimnames(M0)
    )
    
    set.seed(rep_seed)
    
    current_observed_lengths = original_observed_lengths
    operation_rows = integer(n_target_operations)
    operation_columns = integer(n_target_operations)
    
    for (operation_index in seq_len(n_target_operations)) {
      eligible_rows = which(
        current_observed_lengths > terminal_observed_lengths
      )
      
      selected_row = eligible_rows[
        sample.int(length(eligible_rows), size = 1L)
      ]
      
      operation_rows[operation_index] = selected_row
      operation_columns[operation_index] =
        current_observed_lengths[selected_row]
      
      current_observed_lengths[selected_row] =
        current_observed_lengths[selected_row] - 1L
    }
    
    operations = tibble(
      operation_index = seq_len(n_target_operations),
      row_index = operation_rows,
      column_index = operation_columns
    )
    
    operation_counts = seq.int(
      from = jump,
      to = n_target_operations,
      by = step_size
    )
    
    if (tail(operation_counts, 1L) != n_target_operations) {
      operation_counts = c(operation_counts, n_target_operations)
    }
    
    mask_sequence = map2(
      operation_counts,
      seq_along(operation_counts),
      function(n_operations, sequence_index) {
        row_reductions = tabulate(
          operation_rows[seq_len(n_operations)],
          nbins = N
        )
        
        observed_lengths =
          original_observed_lengths - row_reductions
        
        M = mask_from_observed_lengths(
          observed_lengths,
          Tt,
          dimnames(M0)
        )
        
        previous_operation_count = if (sequence_index == 1L) {
          0L
        } else {
          operation_counts[sequence_index - 1L]
        }
        
        list(
          rep_id = as.integer(rep_id),
          rep_seed = as.integer(rep_seed),
          mask_step = sequence_index - 1L,
          M = M,
          observed_lengths = observed_lengths,
          full_rows = full_rows,
          full_row_names = rownames(M0)[full_rows],
          full_column_names = colnames(M0)[seq_len(N_FULL_COLUMNS)],
          n_fully_observed_rows = sum(rowSums(M) == 0L),
          n_fully_observed_columns = sum(colSums(M) == 0L),
          n_artificial_masked = n_operations,
          n_added_from_previous =
            n_operations - previous_operation_count,
          n_remaining_to_terminal =
            n_target_operations - n_operations,
          is_terminal = n_operations == n_target_operations,
          n_target_operations = n_target_operations
        )
      }
    )
    
    names(mask_sequence) = sprintf(
      "M_%04d",
      seq.int(0L, length(mask_sequence) - 1L)
    )
    
    initial_signature = paste(
      which(mask_sequence[[1]]$M & !M0),
      collapse = ","
    )
    
    terminal_increment = n_target_operations - (
      if (length(operation_counts) == 1L) {
        0L
      } else {
        operation_counts[length(operation_counts) - 1L]
      }
    )
    
    structure(
      mask_sequence,
      rep_id = as.integer(rep_id),
      rep_seed = as.integer(rep_seed),
      full_rows = full_rows,
      full_row_names = rownames(M0)[full_rows],
      full_column_names = colnames(M0)[seq_len(N_FULL_COLUMNS)],
      original_observed_lengths = original_observed_lengths,
      terminal_observed_lengths = terminal_observed_lengths,
      terminal_mask = terminal_mask,
      operations = operations,
      initial_signature = initial_signature,
      n_target_operations = n_target_operations,
      jump = jump,
      terminal_increment = terminal_increment
    )
  }
  
  make_replicate_mask_sequences = function(
    M0,
    step_size,
    jump,
    n_rep
  ) {
    replicate_sequences = vector("list", n_rep)
    used_initial_signatures = character()
    
    for (rep_id in seq_len(n_rep)) {
      attempt = 0L
      
      repeat {
        rep_seed = as.integer(
          (REP_SEED + 100000L * rep_id + attempt) %%
            .Machine$integer.max
        )
        
        candidate = make_mask_sequence(
          M0 = M0,
          step_size = step_size,
          jump = jump,
          rep_id = rep_id,
          rep_seed = rep_seed
        )
        
        signature = attr(candidate, "initial_signature")
        
        if (!(signature %in% used_initial_signatures)) {
          break
        }
        
        attempt = attempt + 1L
        
        if (attempt > 10000L) {
          stop("Could not create distinct initial masks.")
        }
      }
      
      replicate_sequences[[rep_id]] = candidate
      used_initial_signatures = c(
        used_initial_signatures,
        signature
      )
    }
    
    names(replicate_sequences) = sprintf(
      "rep_%03d",
      seq_len(n_rep)
    )
    
    replicate_sequences
  }
  
  make_design_D = function(M, mask_all_layers) {
    D = D_original
    artificial_region = M & !M0
    layers = if (mask_all_layers) seq_len(dim(D)[3]) else target_idx
    
    for (k in layers) {
      D_layer = D[, , k]
      D_layer[artificial_region] = 1
      D[, , k] = D_layer
    }
    
    D
  }
  
  make_staggered_geometry = function(A_target, Tt) {
    observed_lengths = ifelse(
      is.infinite(A_target),
      Tt,
      pmin(Tt, A_target - 1L)
    )
    
    row_perm = order(-observed_lengths, seq_along(observed_lengths))
    observed_lengths_sorted = observed_lengths[row_perm]
    m_desc = unique(observed_lengths_sorted)
    
    row_parts = lapply(
      m_desc,
      function(value) which(observed_lengths_sorted == value)
    )
    
    T_sizes = diff(c(0L, rev(m_desc)))
    col_ends = cumsum(T_sizes)
    col_starts = c(1L, head(col_ends, -1L) + 1L)
    col_parts = Map(seq, col_starts, col_ends)
    
    list(
      row_perm = row_perm,
      row_parts = row_parts,
      col_parts = col_parts,
      n_parts = length(row_parts)
    )
  }
  
  estimate_masked_region_ate = function(
    estimator,
    Y0,
    region,
    A,
    Omega
  ) {
    estimator = match.arg(estimator, c("tensor", "matrix"))
    
    geometry = make_staggered_geometry(
      A[, target_idx],
      dim(Y0)[2]
    )
    
    weighted_sum = 0
    normalizer = 0
    
    for (a in seq_len(geometry$n_parts)) {
      for (b in seq_len(geometry$n_parts)) {
        if (a + b <= geometry$n_parts + 1L) {
          next
        }
        
        row_part = geometry$row_parts[[a]]
        column_part = geometry$col_parts[[b]]
        rows_original = geometry$row_perm[row_part]
        region_block =
          region[rows_original, column_part, drop = FALSE]
        
        row_counts = rowSums(region_block)
        
        if (!any(row_counts > 0L)) {
          next
        }
        
        patterns = apply(
          region_block,
          1,
          function(row) paste(as.integer(row), collapse = "")
        )
        
        valid_patterns = unique(patterns[row_counts > 0L])
        
        for (pattern in valid_patterns) {
          rows_local = which(patterns == pattern)
          columns_local =
            which(region_block[rows_local[1], ])
          
          x = rep(0, length(row_part))
          y = rep(0, length(column_part))
          
          x[rows_local] = 1 / sqrt(length(rows_local))
          y[columns_local] = 1 / sqrt(length(columns_local))
          
          mu_hat = if (estimator == "tensor") {
            bilinearTensorStaggered(
              Y = Y0,
              k = target_idx,
              a = a,
              b = b,
              r = rank_value,
              x = x,
              y = y,
              A = A,
              Omega = Omega,
              tau = tau
            )
          } else {
            bilinearMatrixStaggered(
              Y_mat = Y0[, , target_idx],
              a = a,
              b = b,
              r = rank_value,
              x = x,
              y = y,
              A = A[, target_idx],
              Omega = Omega[, , target_idx],
              tau = tau
            )
          }
          
          W = sqrt(length(rows_local) * length(columns_local))
          
          weighted_sum = weighted_sum + W * mu_hat
          normalizer =
            normalizer +
            length(rows_local) * length(columns_local)
        }
      }
    }
    
    if (normalizer == 0) {
      stop("No masked-region blocks were available.")
    }
    
    weighted_sum / normalizer
  }
  
  try_error_message = function(result) {
    condition = attr(result, "condition")
    if (is.null(condition)) {
      as.character(result)
    } else {
      conditionMessage(condition)
    }
  }
  
  run_cy_with_rank_fallback = function(
    Y0,
    A,
    Omega,
    mask_step,
    rep_id
  ) {
    rank_try = rank_value
    minimum_rank = max(1L, rank_value - 1L)
    fit = NULL
    last_error = NA_character_
    
    while (is.null(fit) && rank_try >= minimum_rank) {
      fit_attempt = try(
        cy_complete_staggered_layer(
          y = Y0[, , target_idx],
          omega = Omega[, , target_idx],
          adoption_time = A[, target_idx],
          time = seq_len(dim(Y0)[2]),
          rank = rank_try,
          seed = cy_seed +
            mask_step +
            100000L * (rep_id - 1L),
          verbose = FALSE
        ),
        silent = TRUE
      )
      
      if (inherits(fit_attempt, "try-error")) {
        last_error = try_error_message(fit_attempt)
        rank_try = rank_try - 1L
      } else {
        fit = fit_attempt
      }
    }
    
    list(
      fit = fit,
      rank = if (is.null(fit)) NA_integer_ else rank_try,
      error = if (is.null(fit)) last_error else NA_character_
    )
  }
  
  M0 = D_original[, , target_idx] == 1
  
  replicate_mask_sequences = make_replicate_mask_sequences(
    M0 = M0,
    step_size = STEP_SIZE,
    jump = JUMP,
    n_rep = REP
  )
  
  n_target_operations = attr(
    replicate_mask_sequences[[1]],
    "n_target_operations"
  )
  
  mask_settings = map_dfr(
    replicate_mask_sequences,
    function(mask_sequence) {
      map_dfr(mask_sequence, function(mask_object) {
        tibble(
          rep_id = mask_object$rep_id,
          rep_seed = mask_object$rep_seed,
          mask_step = mask_object$mask_step,
          jump = JUMP,
          step_size = STEP_SIZE,
          n_artificial_masked =
            mask_object$n_artificial_masked,
          n_added_from_previous =
            mask_object$n_added_from_previous,
          n_fully_observed_rows =
            mask_object$n_fully_observed_rows,
          n_fully_observed_columns =
            mask_object$n_fully_observed_columns,
          terminal_full_rows = paste(
            mask_object$full_row_names,
            collapse = ","
          ),
          terminal_full_columns = paste(
            mask_object$full_column_names,
            collapse = ","
          ),
          n_remaining_to_terminal =
            mask_object$n_remaining_to_terminal,
          is_terminal = mask_object$is_terminal,
          masking_fraction =
            mask_object$n_artificial_masked /
            mask_object$n_target_operations
        )
      })
    }
  )
  
  variant_settings = tribble(
    ~variant_id, ~variant_label, ~mask_all_layers,
    "mask_all_layers",
    "Artificial mask on target and auxiliary layers",
    TRUE,
    "mask_target_only",
    "Artificial mask on target layer only",
    FALSE
  )
  
  run_one_mask = function(mask_object, variant_row) {
    run_shared_methods =
      variant_row$mask_all_layers[[1]] &&
      rank_value %in% c(1L, 2L)
    
    M = mask_object$M
    
    D_design = make_design_D(
      M,
      variant_row$mask_all_layers[[1]]
    )
    
    Y0 = Y
    Y0[D_design == 1] = NA_real_
    
    AO = make_A_Omega_from_D(D_design)
    A = AO$A
    Omega = AO$Omega & !is.na(Y0)
    
    region = M & !M0 & !is.na(Y[, , target_idx])
    
    if (!any(region)) {
      stop("The evaluation region is empty.")
    }
    
    ground_truth = mean(Y[, , target_idx][region])
    
    tensor_attempt = try(
      estimate_masked_region_ate(
        estimator = "tensor",
        Y0 = Y0,
        region = region,
        A = A,
        Omega = Omega
      ),
      silent = TRUE
    )
    
    matrix_attempt = if (run_shared_methods) {
      try(
        estimate_masked_region_ate(
          estimator = "matrix",
          Y0 = Y0,
          region = region,
          A = A,
          Omega = Omega
        ),
        silent = TRUE
      )
    } else {
      NULL
    }
    
    estimate_tensor = if (
      inherits(tensor_attempt, "try-error")
    ) {
      NA_real_
    } else {
      as.numeric(tensor_attempt)
    }
    
    estimate_matrix = if (
      is.null(matrix_attempt) ||
      inherits(matrix_attempt, "try-error")
    ) {
      NA_real_
    } else {
      as.numeric(matrix_attempt)
    }
    
    CY = if (run_shared_methods && RUN_CY) {
      run_cy_with_rank_fallback(
        Y0,
        A,
        Omega,
        mask_object$mask_step,
        mask_object$rep_id
      )
    } else {
      list(
        fit = NULL,
        rank = NA_integer_,
        error = NA_character_
      )
    }
    
    estimate_CY_biased = if (is.null(CY$fit)) {
      NA_real_
    } else {
      mean(CY$fit$completed_nuclear[region])
    }
    
    estimate_CY_debiased = if (is.null(CY$fit)) {
      NA_real_
    } else {
      mean(CY$fit$completed_debiased[region])
    }
    
    tibble(
      N_LAYERS = N_LAYERS,
      rank_value = rank_value,
      rep_id = mask_object$rep_id,
      rep_seed = mask_object$rep_seed,
      variant_id = variant_row$variant_id[[1]],
      variant_label = variant_row$variant_label[[1]],
      mask_all_layers =
        variant_row$mask_all_layers[[1]],
      mask_step = mask_object$mask_step,
      jump = JUMP,
      step_size = STEP_SIZE,
      n_artificial_masked =
        mask_object$n_artificial_masked,
      n_added_from_previous =
        mask_object$n_added_from_previous,
      n_fully_observed_rows =
        mask_object$n_fully_observed_rows,
      n_fully_observed_columns =
        mask_object$n_fully_observed_columns,
      terminal_full_rows = paste(
        mask_object$full_row_names,
        collapse = ","
      ),
      terminal_full_columns = paste(
        mask_object$full_column_names,
        collapse = ","
      ),
      n_evaluation_entries = sum(region),
      n_remaining_to_terminal =
        mask_object$n_remaining_to_terminal,
      is_terminal = mask_object$is_terminal,
      masking_fraction =
        mask_object$n_artificial_masked /
        mask_object$n_target_operations,
      rank_tensor = rank_value,
      rank_matrix = if (run_shared_methods) {
        rank_value
      } else {
        NA_integer_
      },
      rank_CY_biased = CY$rank,
      rank_CY_debiased = CY$rank,
      ground_truth = ground_truth,
      estimate_tensor = estimate_tensor,
      estimate_matrix = estimate_matrix,
      estimate_CY_biased = estimate_CY_biased,
      estimate_CY_debiased = estimate_CY_debiased,
      squared_error_tensor =
        (estimate_tensor - ground_truth)^2,
      squared_error_matrix =
        (estimate_matrix - ground_truth)^2,
      squared_error_CY_biased =
        (estimate_CY_biased - ground_truth)^2,
      squared_error_CY_debiased =
        (estimate_CY_debiased - ground_truth)^2,
      tensor_error = if (
        inherits(tensor_attempt, "try-error")
      ) {
        try_error_message(tensor_attempt)
      } else {
        NA_character_
      },
      matrix_error = if (is.null(matrix_attempt)) {
        NA_character_
      } else if (
        inherits(matrix_attempt, "try-error")
      ) {
        try_error_message(matrix_attempt)
      } else {
        NA_character_
      },
      CY_error = CY$error
    )
  }
  
  simulation_results = map_dfr(
    replicate_mask_sequences,
    function(mask_sequence) {
      map_dfr(
        seq_len(nrow(variant_settings)),
        function(variant_index) {
          variant_row =
            variant_settings[variant_index, ]
          
          map_dfr(
            mask_sequence,
            ~ run_one_mask(.x, variant_row)
          )
        }
      )
    }
  )
  
  if (!RUN_CY) {
    simulation_results = simulation_results %>%
      select(
        -rank_CY_biased,
        -rank_CY_debiased,
        -estimate_CY_biased,
        -estimate_CY_debiased,
        -squared_error_CY_biased,
        -squared_error_CY_debiased,
        -CY_error
      )
  }
  
  tensor_squared_errors = simulation_results %>%
    filter(
      variant_id %in%
        c("mask_all_layers", "mask_target_only")
    ) %>%
    transmute(
      N_LAYERS,
      rank_value,
      rep_id,
      rep_seed,
      mask_step,
      jump,
      step_size,
      n_artificial_masked,
      masking_fraction,
      is_terminal,
      method = if_else(
        variant_id == "mask_all_layers",
        "Tensor: all layers masked",
        "Tensor: target layer only"
      ),
      squared_error = squared_error_tensor
    )
  
  matrix_squared_errors = if (
    rank_value %in% c(1L, 2L)
  ) {
    simulation_results %>%
      filter(variant_id == "mask_all_layers") %>%
      transmute(
        N_LAYERS,
        rank_value,
        rep_id,
        rep_seed,
        mask_step,
        jump,
        step_size,
        n_artificial_masked,
        masking_fraction,
        is_terminal,
        method = "Matrix",
        squared_error = squared_error_matrix
      )
  } else {
    tibble()
  }
  
  squared_error_by_rep_step = bind_rows(
    tensor_squared_errors,
    matrix_squared_errors
  )
  
  if (RUN_CY && rank_value %in% c(1L, 2L)) {
    CY_squared_errors = simulation_results %>%
      filter(variant_id == "mask_all_layers") %>%
      select(
        N_LAYERS,
        rank_value,
        rep_id,
        rep_seed,
        mask_step,
        jump,
        step_size,
        n_artificial_masked,
        masking_fraction,
        is_terminal,
        squared_error_CY_biased,
        squared_error_CY_debiased
      ) %>%
      pivot_longer(
        cols = starts_with("squared_error_CY_"),
        names_to = "method",
        values_to = "squared_error"
      ) %>%
      mutate(
        method = recode(
          method,
          squared_error_CY_biased = "CY biased",
          squared_error_CY_debiased = "CY debiased"
        )
      )
    
    squared_error_by_rep_step = bind_rows(
      squared_error_by_rep_step,
      CY_squared_errors
    )
  }
  
  squared_error_by_rep_step =
    squared_error_by_rep_step %>%
    arrange(method, rep_id, mask_step)
  
  average_squared_error_by_step =
    squared_error_by_rep_step %>%
    group_by(
      N_LAYERS,
      rank_value,
      mask_step,
      jump,
      step_size,
      n_artificial_masked,
      masking_fraction,
      is_terminal,
      method
    ) %>%
    summarise(
      n_rep_available =
        sum(!is.na(squared_error)),
      average_squared_error = if (
        n_rep_available > 0L
      ) {
        mean(squared_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      sd_squared_error = if (
        n_rep_available > 1L
      ) {
        sd(squared_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      se_squared_error = if (
        n_rep_available > 1L
      ) {
        sd_squared_error / sqrt(n_rep_available)
      } else {
        NA_real_
      },
      min_squared_error = if (
        n_rep_available > 0L
      ) {
        min(squared_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      max_squared_error = if (
        n_rep_available > 0L
      ) {
        max(squared_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    arrange(method, mask_step)
  
  average_squared_error_by_rep =
    squared_error_by_rep_step %>%
    group_by(
      N_LAYERS,
      rank_value,
      rep_id,
      rep_seed,
      jump,
      step_size,
      method
    ) %>%
    summarise(
      n_steps_available =
        sum(!is.na(squared_error)),
      average_squared_error = if (
        n_steps_available > 0L
      ) {
        mean(squared_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    arrange(method, rep_id)
  
  analysis_results = list(
    metadata = list(
      target_layer = target_layer,
      target_name = target_name,
      target_label = target_label,
      JUMP = JUMP,
      STEP_SIZE = STEP_SIZE,
      REP = REP,
      REP_SEED = REP_SEED,
      MASK_REPLICATION_SCHEME = paste(
        REP,
        "distinct initial masks start with",
        JUMP,
        "additional missing cells.",
        "Each later mask adds",
        STEP_SIZE,
        "cells and all replicates end at the same terminal mask."
      ),
      RUN_CY = RUN_CY,
      N_UNITS = N_UNITS,
      N_PERIODS = N_PERIODS,
      N_LAYERS = N_LAYERS,
      rank_value = rank_value,
      tau = tau,
      cy_seed = cy_seed,
      DATASET_TAG = DATASET_TAG,
      RESULTS_STUDY_TAG = RESULTS_STUDY_TAG,
      tensor_dimensions = dim(Y),
      tensor_file = tensor_file,
      slice_label = slice_label,
      target_outcome = target_layer,
      target_policy = target_policy,
      selected_layer_map = layer_map,
      start_date =
        covid_data$tensor_data$start.date %||%
        NA_character_,
      time_horizon =
        covid_data$tensor_data$time_horizon %||%
        NA_character_,
      delay_days =
        covid_data$tensor_data$delay_days %||%
        NA_integer_,
      N_NEVER_TREATED = sum(
        !state_order$ever_treated
      ),
      N_FULL_ROWS = N_FULL_ROWS,
      N_FULL_COLUMNS = N_FULL_COLUMNS,
      FIRST_BLOCK_SIZE = N_FULL_COLUMNS,
      n_target_operations = n_target_operations,
      terminal_increment = attr(
        replicate_mask_sequences[[1]],
        "terminal_increment"
      ),
      terminal_full_rows_by_rep = map(
        replicate_mask_sequences,
        ~ attr(.x, "full_row_names")
      ),
      terminal_full_columns = attr(
        replicate_mask_sequences[[1]],
        "full_column_names"
      ),
      crime_vars = crime_vars,
      crime_labels = crime_labels
    ),
    variant_settings = variant_settings,
    mask_settings = mask_settings,
    simulation_results = simulation_results,
    squared_error_by_rep_step = squared_error_by_rep_step,
    average_squared_error_by_step =
      average_squared_error_by_step,
    average_squared_error_by_rep =
      average_squared_error_by_rep,
    mask_sequence = replicate_mask_sequences[[1]],
    mask_sequences = replicate_mask_sequences,
    Y = Y,
    Omega_policy = policy,
    D_original = D_original,
    layer_map = layer_map,
    state_order = state_order,
    year_order = year_order,
    crime_order = crime_order
  )
  
  result_stem = paste0(
    DATASET_TAG,
    "_staggered_masking_results_",
    target_name,
    "_rank",
    rank_value,
    "_jump",
    JUMP,
    "_step",
    STEP_SIZE,
    rep_file_tag,
    cy_file_tag,
    slice_file_tag
  )
  
  result_file = file.path(
    results_dir,
    paste0(result_stem, ".rds")
  )
  
  saveRDS(analysis_results, result_file)
  
  mask_settings_file = file.path(
    results_dir,
    paste0(
      DATASET_TAG,
      "_mask_settings_",
      target_name,
      "_rank",
      rank_value,
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      slice_file_tag,
      ".csv"
    )
  )
  
  simulation_results_file = file.path(
    results_dir,
    paste0(
      DATASET_TAG,
      "_simulation_results_",
      target_name,
      "_rank",
      rank_value,
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      slice_file_tag,
      ".csv"
    )
  )
  
  squared_error_by_rep_step_file = file.path(
    results_dir,
    paste0(
      result_stem,
      "_squared_error_by_rep_step.csv"
    )
  )
  
  average_squared_error_by_step_file = file.path(
    results_dir,
    paste0(
      result_stem,
      "_average_squared_error_by_step.csv"
    )
  )
  
  average_squared_error_by_rep_file = file.path(
    results_dir,
    paste0(
      result_stem,
      "_average_squared_error_by_rep.csv"
    )
  )
  
  if (SAVE_CSV_OUTPUTS) {
    write_csv(mask_settings, mask_settings_file)
    write_csv(
      simulation_results,
      simulation_results_file
    )
    write_csv(
      squared_error_by_rep_step,
      squared_error_by_rep_step_file
    )
    write_csv(
      average_squared_error_by_step,
      average_squared_error_by_step_file
    )
    write_csv(
      average_squared_error_by_rep,
      average_squared_error_by_rep_file
    )
    
    walk(variant_settings$variant_id, function(variant_id) {
      write_csv(
        filter(
          simulation_results,
          .data$variant_id == variant_id
        ),
        file.path(
          results_dir,
          paste0(
            DATASET_TAG,
            "_simulation_results_",
            target_name,
            "_rank",
            rank_value,
            "_",
            variant_id,
            "_jump",
            JUMP,
            "_step",
            STEP_SIZE,
            rep_file_tag,
            cy_file_tag,
            slice_file_tag,
            ".csv"
          )
        )
      )
    })
  }
  
  tibble(
    N_LAYERS = N_LAYERS,
    rank_value = rank_value,
    REP = REP,
    result_file = result_file,
    simulation_results_file = if (
      SAVE_CSV_OUTPUTS
    ) {
      simulation_results_file
    } else {
      NA_character_
    },
    mask_settings_file = if (
      SAVE_CSV_OUTPUTS
    ) {
      mask_settings_file
    } else {
      NA_character_
    },
    squared_error_by_rep_step_file = if (
      SAVE_CSV_OUTPUTS
    ) {
      squared_error_by_rep_step_file
    } else {
      NA_character_
    },
    average_squared_error_by_step_file = if (
      SAVE_CSV_OUTPUTS
    ) {
      average_squared_error_by_step_file
    } else {
      NA_character_
    },
    average_squared_error_by_rep_file = if (
      SAVE_CSV_OUTPUTS
    ) {
      average_squared_error_by_rep_file
    } else {
      NA_character_
    },
    mask_plots_dir = NA_character_
  )
}

# Run all ranks ---------------------------------------------------------------

run_summary = pmap_dfr(
  SIMULATION_GRID,
  ~ run_simulation(
    N_LAYERS = ..1,
    rank_value = ..2
  )
)

if (SAVE_CSV_OUTPUTS) {
  combined_results = map(
    run_summary$result_file,
    readRDS
  )
  
  combined_squared_error_by_rep_step = map_dfr(
    combined_results,
    "squared_error_by_rep_step"
  )
  
  combined_average_squared_error_by_step = map_dfr(
    combined_results,
    "average_squared_error_by_step"
  )
  
  combined_average_squared_error_by_rep = map_dfr(
    combined_results,
    "average_squared_error_by_rep"
  )
  
  combined_summary_stem = paste0(
    DATASET_TAG,
    "_all_configurations_staggered_masking_",
    target_name,
    "_jump",
    JUMP,
    "_step",
    STEP_SIZE,
    rep_file_tag,
    cy_file_tag,
    slice_file_tag
  )
  
  combined_results_dir = file.path(
    "Results",
    RESULTS_STUDY_TAG,
    paste0("step", STEP_SIZE)
  )
  
  dir.create(
    combined_results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write_csv(
    combined_squared_error_by_rep_step,
    file.path(
      combined_results_dir,
      paste0(
        combined_summary_stem,
        "_squared_error_by_rep_step.csv"
      )
    )
  )
  
  write_csv(
    combined_average_squared_error_by_step,
    file.path(
      combined_results_dir,
      paste0(
        combined_summary_stem,
        "_average_squared_error_by_step.csv"
      )
    )
  )
  
  write_csv(
    combined_average_squared_error_by_rep,
    file.path(
      combined_results_dir,
      paste0(
        combined_summary_stem,
        "_average_squared_error_by_rep.csv"
      )
    )
  )
  
  write_csv(
    run_summary,
    file.path(
      combined_results_dir,
      paste0(
        combined_summary_stem,
        "_run_summary.csv"
      )
    )
  )
}

print(
  run_summary,
  n = Inf,
  width = Inf
)
