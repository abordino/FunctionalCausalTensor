setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CovidOx")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# ===============================================================
# 1. Set par
# ===============================================================

tensor_file = "data/Omega_Y_until_2020-04-05_delay_28.rds"

n_outcomes_to_keep = 2
target_outcome_pattern = "death"

r_grid = c(3)
tau = 0.0001
functionals = c("ATE", "Trend")

B = 500
bootstrap_seed = 123
stratify_target_bootstrap = TRUE

save_diagnostic_plots = TRUE
results_dir = "results/oxford_deaths_target_layer_bootstrap"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# 2. Source estimator files
# ===============================================================

source_files = c(
  "../../bilinearMatrixStaggered.R",
  "../../bilinearTensorStaggered.R",
  "../../bilinearMatrixStaggeredPsi.R",
  "../../bilinearTensorStaggeredPsi.R",
  "../../pluginPsi_c1.R"
)

walk(source_files, function(f) {
  if (!file.exists(f)) {
    stop("Required source file not found: ", f)
  }
  source(f)
})

# ===============================================================
# 3. Some helpers
# ===============================================================

`%||%` = function(x, y) {
  if (is.null(x)) y else x
}

as_layer_names = function(x, K, prefix) {
  if (is.null(x)) {
    paste0(prefix, "_", seq_len(K))
  } else {
    as.character(x[seq_len(K)])
  }
}

first_true_index = function(x) {
  idx = which(x)
  if (length(idx) == 0) Inf else idx[1]
}

check_policy_is_staggered = function(Omega_policy, label = "Omega_policy") {
  if (!is.array(Omega_policy) || length(dim(Omega_policy)) != 3) {
    stop(label, " must be an N x T x K logical array.")
  }
  
  Omega_policy = array(
    as.logical(Omega_policy),
    dim = dim(Omega_policy),
    dimnames = dimnames(Omega_policy)
  )
  
  for (k in seq_len(dim(Omega_policy)[3])) {
    for (i in seq_len(dim(Omega_policy)[1])) {
      row_i = Omega_policy[i, , k]
      
      if (any(diff(as.integer(row_i)) < 0)) {
        stop(
          label,
          " is not staggered at row ", i,
          ", layer ", k,
          ". Policy-on must not switch back to off."
        )
      }
    }
  }
  
  invisible(TRUE)
}

make_observation_mask = function(Omega_policy) {
  Omega_obs = !Omega_policy
  dimnames(Omega_obs) = dimnames(Omega_policy)
  storage.mode(Omega_obs) = "logical"
  Omega_obs
}

make_A_from_policy = function(Omega_policy) {
  A = apply(Omega_policy, c(1, 3), first_true_index)
  
  dimnames(A) = list(
    row = dimnames(Omega_policy)[[1]],
    layer = dimnames(Omega_policy)[[3]]
  )
  
  A
}

# ===============================================================
# 4. Align plotY outcomes with their policy slices
# ===============================================================

align_plotY_layers = function(
    tensor_data,
    Y,
    Omega_policy,
    target_outcome_pattern = "death",
    n_outcomes_to_keep = 2
) {
  
  K_y = dim(Y)[3]
  K_o = dim(Omega_policy)[3]
  
  country_names = tensor_data$country_levels %||% paste0("country_", seq_len(dim(Y)[1]))
  time_names = as.character(tensor_data$date_names %||% seq_len(dim(Y)[2]))
  outcome_names = as_layer_names(tensor_data$outcomes, K_y, "outcome")
  policy_names = as_layer_names(tensor_data$policies, K_o, "policy")
  
  if (length(country_names) != dim(Y)[1]) {
    stop("country_levels length does not match dim(Y)[1].")
  }
  if (length(time_names) != dim(Y)[2]) {
    stop("date_names length does not match dim(Y)[2].")
  }
  
  dimnames(Y) = list(
    country = country_names,
    time = time_names,
    layer = outcome_names
  )
  
  dimnames(Omega_policy) = list(
    country = country_names,
    time = time_names,
    layer = policy_names
  )
  
  policy_outcome_map = tensor_data$policy_outcome_map
  
  if (
    is.null(policy_outcome_map) ||
    !all(c("policy", "outcome") %in% names(policy_outcome_map))
  ) {
    
    layer_map = tibble(
      outcome = outcome_names,
      policy = policy_names,
      panel_label = outcome_names
    )
  } else {
    layer_map = as_tibble(policy_outcome_map) %>%
      mutate(
        outcome = as.character(outcome),
        policy = as.character(policy),
        panel_label = if ("panel_label" %in% names(.)) {
          as.character(panel_label)
        } else {
          as.character(outcome)
        }
      ) %>%
      filter(
        outcome %in% outcome_names,
        policy %in% policy_names
      )
  }
  
  layer_map = layer_map %>%
    mutate(
      outcome_index = match(outcome, outcome_names),
      policy_index = match(policy, policy_names)
    ) %>%
    arrange(outcome_index)
  
  target_candidates = layer_map %>%
    filter(
      str_detect(str_to_lower(outcome), target_outcome_pattern) |
        str_detect(str_to_lower(panel_label), target_outcome_pattern)
    )
  
  target_row = target_candidates %>%
    mutate(
      exact_priority = if_else(str_to_lower(outcome) %in% c("death", "deaths"), 0L, 1L)
    ) %>%
    arrange(exact_priority, outcome_index) %>%
    slice(1)
  
  if (nrow(layer_map) == n_outcomes_to_keep) {
    selected_map = layer_map
  } else {
    selected_map = bind_rows(
      target_row,
      layer_map %>%
        filter(outcome != target_row$outcome[[1]]) %>%
        slice_head(n = n_outcomes_to_keep - 1)
    ) %>%
      distinct(outcome, .keep_all = TRUE) %>%
      arrange(outcome_index)
  }
  
  if (nrow(selected_map) != n_outcomes_to_keep) {
    stop(
      "Expected ", n_outcomes_to_keep,
      " selected outcomes but found ", nrow(selected_map), "."
    )
  }
  
  Y_selected = Y[, , selected_map$outcome_index, drop = FALSE]
  
  Omega_selected = array(
    FALSE,
    dim = c(dim(Omega_policy)[1], dim(Omega_policy)[2], nrow(selected_map)),
    dimnames = list(
      country = country_names,
      time = time_names,
      layer = selected_map$outcome
    )
  )
  
  for (jj in seq_len(nrow(selected_map))) {
    Omega_selected[, , jj] = Omega_policy[, , selected_map$policy_index[[jj]]]
  }
  
  dimnames(Y_selected) = list(
    country = country_names,
    time = time_names,
    layer = selected_map$outcome
  )
  
  storage.mode(Y_selected) = "numeric"
  storage.mode(Omega_selected) = "logical"
  
  target_k = match(target_row$outcome[[1]], selected_map$outcome)
  
  list(
    Y = Y_selected,
    Omega_policy = Omega_selected,
    layer_map = selected_map %>%
      transmute(
        layer_index = row_number(),
        outcome,
        policy,
        panel_label,
        original_outcome_index = outcome_index,
        original_policy_index = policy_index
      ),
    target_k = target_k,
    target_outcome = selected_map$outcome[[target_k]],
    target_policy = selected_map$policy[[target_k]]
  )
}

# ===============================================================
# 5. Target-layer staircase helpers
# ===============================================================

make_target_staircase = function(Omega_policy, k, tie_names = NULL) {
  if (k < 1 || k > dim(Omega_policy)[3]) {
    stop("k is outside the layer range.")
  }
  
  Omega_k = Omega_policy[, , k, drop = FALSE][, , 1]
  N = nrow(Omega_k)
  Tt = ncol(Omega_k)
  
  A_k = apply(Omega_k, 1, first_true_index)
  obs_len = ifelse(is.infinite(A_k), Tt, A_k - 1)
  
  if (is.null(tie_names)) {
    tie_names = dimnames(Omega_policy)[[1]] %||% as.character(seq_len(N))
  }
  
  perm = tibble(
    row_original = seq_len(N),
    tie_name = tie_names,
    observed_prefix_length = obs_len,
    ever_treated = is.finite(A_k)
  ) %>%
    arrange(
      desc(observed_prefix_length),
      tie_name,
      row_original
    ) %>%
    pull(row_original)
  
  obs_len_sorted = obs_len[perm]
  A_sorted = A_k[perm]
  
  m_desc = unique(obs_len_sorted)
  row_parts = lapply(m_desc, function(m) which(obs_len_sorted == m))
  
  m_asc = rev(m_desc)
  T_sizes = as.integer(diff(c(0, m_asc)))
  N_sizes = as.integer(lengths(row_parts))
  
  if (length(N_sizes) != length(T_sizes)) {
    stop("N_parts and T_parts have different lengths.")
  }
  
  if (any(N_sizes <= 0)) {
    stop("Target-layer row partition contains an empty block.")
  }
  
  if (any(T_sizes <= 0)) {
    stop("Target-layer column partition contains an empty block.")
  }
  
  list(
    permutation = perm,
    A_original = A_k,
    A_sorted = A_sorted,
    observed_prefix_original = obs_len,
    observed_prefix_sorted = obs_len_sorted,
    N_sizes = N_sizes,
    T_sizes = T_sizes,
    n_blocks = length(N_sizes)
  )
}

reorder_to_target_staircase = function(Y, Omega_policy, k) {
  row_names = dimnames(Y)[[1]] %||% dimnames(Omega_policy)[[1]]
  
  stair = make_target_staircase(
    Omega_policy = Omega_policy,
    k = k,
    tie_names = row_names
  )
  
  perm = stair$permutation
  
  Y_re = Y[perm, , , drop = FALSE]
  Omega_policy_re = Omega_policy[perm, , , drop = FALSE]
  
  dimnames(Y_re)[[1]] = row_names[perm]
  dimnames(Omega_policy_re)[[1]] = row_names[perm]
  
  stair_re = make_target_staircase(
    Omega_policy = Omega_policy_re,
    k = k,
    tie_names = dimnames(Y_re)[[1]]
  )
  
  if (!identical(stair_re$permutation, seq_along(stair_re$permutation))) {
    stop("Reordered target layer is not already in staircase order.")
  }
  
  list(
    Y = Y_re,
    Omega_policy = Omega_policy_re,
    staircase = stair_re,
    original_permutation = perm
  )
}

make_plugin_parts_list = function(K, N_sizes, T_sizes) {
  list(
    N_parts = replicate(K, N_sizes, simplify = FALSE),
    T_parts = replicate(K, T_sizes, simplify = FALSE)
  )
}

is_staircase_layer = function(Omega_obs_k) {
  N = nrow(Omega_obs_k)
  Tt = ncol(Omega_obs_k)
  
  obs_len = integer(N)
  
  for (i in seq_len(N)) {
    obs = Omega_obs_k[i, ]
    
    if (any(diff(as.integer(obs)) > 0)) {
      return(FALSE)
    }
    
    first_unobs = which(!obs)[1]
    obs_len[i] = if (is.na(first_unobs)) Tt else first_unobs - 1
  }
  
  all(diff(obs_len) <= 0)
}

# ===============================================================
# 6. Y NA handling and diagnostics
# ===============================================================

fill_Y_na_by_neighbourhood = function(Y, max_passes = 5) {
  Y_filled = Y
  
  for (pass in seq_len(max_passes)) {
    na_positions = which(is.na(Y_filled), arr.ind = TRUE)
    
    if (nrow(na_positions) == 0) {
      return(Y_filled)
    }
    
    n_before = nrow(na_positions)
    
    for (idx in seq_len(nrow(na_positions))) {
      i = na_positions[idx, 1]
      t = na_positions[idx, 2]
      k = na_positions[idx, 3]
      
      neighbours = numeric(0)
      
      if (t > 1) {
        neighbours = c(neighbours, Y_filled[i, t - 1, k])
      }
      
      if (t < dim(Y_filled)[2]) {
        neighbours = c(neighbours, Y_filled[i, t + 1, k])
      }
      
      neighbours = neighbours[!is.na(neighbours)]
      
      if (length(neighbours) > 0) {
        Y_filled[i, t, k] = mean(neighbours)
      }
    }
    
    n_after = sum(is.na(Y_filled))
    
    if (n_after == n_before) {
      break
    }
  }
  
  Y_filled
}

plot_Y_na_check = function(Y, title = "Y tensor NA check") {
  Y_na_df = as.data.frame.table(
    is.na(Y),
    responseName = "is_na",
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  names(Y_na_df)[1:3] = c("country", "time", "layer")
  
  Y_na_df = Y_na_df %>%
    mutate(
      country = factor(country, levels = rev(dimnames(Y)[[1]])),
      time_index = as.integer(factor(time, levels = dimnames(Y)[[2]])),
      missing_status = if_else(is_na, "NA", "not NA")
    )
  
  ggplot(
    Y_na_df,
    aes(x = time_index, y = country, fill = missing_status)
  ) +
    geom_tile(width = 0.95, height = 0.95) +
    facet_wrap(~ layer, ncol = 2) +
    scale_fill_manual(
      values = c(
        "not NA" = "white",
        "NA" = "black"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(dimnames(Y)[[2]]),
      labels = dimnames(Y)[[2]],
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      subtitle = "White = not NA; black = NA",
      x = "Date",
      y = "Country",
      fill = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

# ===============================================================
# 7. Load plotY tensor data and keep deaths + one fixed anchor outcome
# ===============================================================

tensor_data = readRDS(tensor_file)

Y_raw = tensor_data$Y
Omega_raw = tensor_data$Omega

if (!is.array(Y_raw) || length(dim(Y_raw)) != 3) {
  stop("tensor_data$Y must be an N x T x K array.")
}

if (!is.array(Omega_raw) || length(dim(Omega_raw)) != 3) {
  stop("tensor_data$Omega must be an N x T x K array.")
}

storage.mode(Y_raw) = "numeric"
storage.mode(Omega_raw) = "logical"

aligned = align_plotY_layers(
  tensor_data = tensor_data,
  Y = Y_raw,
  Omega_policy = Omega_raw,
  target_outcome_pattern = target_outcome_pattern,
  n_outcomes_to_keep = n_outcomes_to_keep
)

Y = aligned$Y
Omega_policy = aligned$Omega_policy
layer_map = aligned$layer_map
target_k = aligned$target_k
target_outcome = aligned$target_outcome
target_policy = aligned$target_policy

check_policy_is_staggered(Omega_policy)

message(
  "Loaded plotY tensors: Y = ",
  paste(dim(Y), collapse = " x "),
  "; Omega_policy = ",
  paste(dim(Omega_policy), collapse = " x ")
)

message(
  "Selected outcomes: ",
  paste(layer_map$outcome, collapse = ", ")
)

message(
  "Aligned policy slices: ",
  paste(layer_map$policy, collapse = ", ")
)

message(
  "Target outcome layer: ", target_outcome,
  " | target policy slice: ", target_policy,
  " | target_k = ", target_k
)

message("Raw selected Y NA count before filling: ", sum(is.na(Y)))
Y = fill_Y_na_by_neighbourhood(Y, max_passes = 5)
message("Raw selected Y NA count after filling: ", sum(is.na(Y)))

if (save_diagnostic_plots) {
  g_original = plot_Y_na_check(
    Y,
    title = "Selected plotY outcomes NA check after filling"
  )
  
  ggsave(
    filename = file.path(results_dir, "Y_selected_NA_check.png"),
    plot = g_original,
    width = 11,
    height = 8,
    dpi = 200
  )
}

tensor_settings = list(
  plotY_slice = list(
    Y = Y,
    Omega_policy = Omega_policy
  )
)

# ===============================================================
# 8. Estimation
# ===============================================================

estimate_one_from_arrays = function(
    Y_in,
    Omega_policy_in,
    setting,
    k,
    r,
    functional,
    tau,
    layer_map
) {
  functional = match.arg(functional, choices = c("ATE", "Trend"))
  
  if (anyNA(Y_in)) {
    stop("Y_in contains NA before estimation.")
  }
  
  check_policy_is_staggered(Omega_policy_in)
  
  reordered = reorder_to_target_staircase(
    Y = Y_in,
    Omega_policy = Omega_policy_in,
    k = k
  )
  
  Y_re = reordered$Y
  Omega_policy_re = reordered$Omega_policy
  stair = reordered$staircase
  
  Omega_obs_re = make_observation_mask(Omega_policy_re)
  
  if (!is_staircase_layer(Omega_obs_re[, , k])) {
    stop("Target observation mask is not in staircase form after reordering.")
  }
  
  Y0_re = Y_re
  Y0_re[!Omega_obs_re] = NA_real_
  
  A_re = make_A_from_policy(Omega_policy_re)
  
  parts = make_plugin_parts_list(
    K = dim(Y_re)[3],
    N_sizes = stair$N_sizes,
    T_sizes = stair$T_sizes
  )
  
  psi0_tensor = bilinearTensorStaggeredPsi(
    Y = Y0_re,
    k = k,
    r = r,
    tau = tau,
    functional = functional,
    eta = NULL,
    row_index = NULL,
    A = A_re,
    Omega = Omega_obs_re
  )
  
  psi0_matrix = bilinearMatrixStaggeredPsi(
    Y_mat = Y0_re[, , k],
    r = r,
    tau = tau,
    functional = functional,
    eta = NULL,
    row_index = NULL,
    A = A_re[, k],
    Omega = Omega_obs_re[, , k]
  )
  
  psi1 = pluginPsi_c1(
    Y = Y_re,
    k = k,
    N_parts = parts$N_parts,
    T_parts = parts$T_parts,
    functional = functional,
    eta = NULL,
    row_index = NULL
  )
  
  target_outcome_here = dimnames(Y_re)[[3]][[k]]
  target_policy_here = layer_map$policy[match(target_outcome_here, layer_map$outcome)]
  
  tibble(
    setting = setting,
    target_outcome = target_outcome_here,
    target_policy = target_policy_here,
    target_k = k,
    r = r,
    functional = functional,
    n_staircase_blocks = stair$n_blocks,
    N_parts = paste(stair$N_sizes, collapse = ","),
    T_parts = paste(stair$T_sizes, collapse = ","),
    Psi0 = as.numeric(psi0_tensor),
    Psi1 = as.numeric(psi1),
    Delta_h = as.numeric(psi1 - psi0_tensor),
    Psi0_matrix = as.numeric(psi0_matrix),
    Delta_h_matrix = as.numeric(psi1 - psi0_matrix),
    tensor_minus_matrix = as.numeric(psi0_tensor - psi0_matrix)
  )
}

run_estimator_given_arrays = function(
    Y_in,
    Omega_policy_in,
    setting,
    target_k,
    r_grid,
    functionals,
    tau,
    layer_map
) {
  expand_grid(
    r = r_grid,
    functional = functionals
  ) %>%
    pmap_dfr(function(r, functional) {
      estimate_one_from_arrays(
        Y_in = Y_in,
        Omega_policy_in = Omega_policy_in,
        setting = setting,
        k = target_k,
        r = r,
        functional = functional,
        tau = tau,
        layer_map = layer_map
      )
    })
}

# ===============================================================
# 9. Point estimates
# ===============================================================

point_results = run_estimator_given_arrays(
  Y_in = tensor_settings$plotY_slice$Y,
  Omega_policy_in = tensor_settings$plotY_slice$Omega_policy,
  setting = "plotY_slice",
  target_k = target_k,
  r_grid = r_grid,
  functionals = functionals,
  tau = tau,
  layer_map = layer_map
) %>%
  mutate(across(
    c(Psi0, Psi1, Delta_h, Psi0_matrix, Delta_h_matrix, tensor_minus_matrix),
    ~ round(.x, 6)
  ))

cat("\n============================================================\n")
cat("Point estimates: deaths target layer only\n")
cat("============================================================\n")
print(point_results, n = Inf)

# ===============================================================
# 10. Target-layer-only bootstrap helpers
# ===============================================================

ci_vars = c(
  "Psi0",
  "Psi1",
  "Delta_h",
  "Psi0_matrix",
  "Delta_h_matrix",
  "tensor_minus_matrix"
)

get_target_row_strata = function(Y, Omega_policy, target_k) {
  
  if (target_k < 1 || target_k > dim(Y)[3]) {
    stop("target_k is outside the Y layer range.")
  }
  
  Y_target = Y[, , target_k]
  policy_target = Omega_policy[, , target_k]
  
  fully_observed_idx = which(
    rowSums(policy_target == TRUE, na.rm = TRUE) == 0 &
      rowSums(is.na(Y_target)) == 0
  )
  
  non_fully_observed_idx = setdiff(
    seq_len(dim(Y)[1]),
    fully_observed_idx
  )
  
  list(
    target_layer_index = target_k,
    fully_observed_idx = fully_observed_idx,
    non_fully_observed_idx = non_fully_observed_idx
  )
}

make_target_layer_row_resample_index = function(
    Y,
    Omega_policy,
    target_k,
    stratify_target_bootstrap = TRUE
) {
  
  strata = get_target_row_strata(
    Y = Y,
    Omega_policy = Omega_policy,
    target_k = target_k
  )
  
  full_idx = strata$fully_observed_idx
  nonfull_idx = strata$non_fully_observed_idx
  
  N = dim(Y)[1]
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
    Omega_policy_base,
    seed_base,
    point_results
) {
  
  set.seed(seed_base + b)
  
  idx = make_target_layer_row_resample_index(
    Y = Y_base,
    Omega_policy = Omega_policy_base,
    target_k = target_k,
    stratify_target_bootstrap = stratify_target_bootstrap
  )
  
  Y_b = Y_base
  Omega_policy_b = Omega_policy_base
  
  other_layers = setdiff(seq_len(dim(Y_base)[3]), target_k)
  
  Y_b[, , target_k] = Y_base[idx, , target_k]
  
  if (length(other_layers) > 0) {
    Y_b[, , other_layers] = Y_base[, , other_layers]
  }
  
  dimnames(Y_b) = dimnames(Y_base)
  dimnames(Omega_policy_b) = dimnames(Omega_policy_base)
  
  if (!identical(dimnames(Y_b), dimnames(Y_base))) {
    stop("Y_b dimnames changed.")
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
      stop("Anchor outcome layer changed, but it should be fixed.")
    }
  }
  
  run_estimator_given_arrays(
    Y_in = Y_b,
    Omega_policy_in = Omega_policy_b,
    setting = "plotY_slice",
    target_k = target_k,
    r_grid = r_grid,
    functionals = functionals,
    tau = tau,
    layer_map = layer_map
  ) %>%
    mutate(
      bootstrap_id = b,
      bootstrap_type = "target_layer_only_fixed_anchors",
      bootstrap_error = NA_character_
    )
}

# ===============================================================
# 11. Bootstrap
# ===============================================================

strata = get_target_row_strata(
  Y = Y,
  Omega_policy = Omega_policy,
  target_k = target_k
)

cat("\nTarget-layer bootstrap strata:\n")
cat("Fully observed target-layer rows:", length(strata$fully_observed_idx), "\n")
cat("Non-fully-observed target-layer rows:", length(strata$non_fully_observed_idx), "\n")

staircase_tbl = tibble(
  layer = dimnames(Y)[[3]],
  policy_slice = layer_map$policy,
  staircase = map_lgl(seq_len(dim(Y)[3]), function(kk) {
    is_staircase_layer(make_observation_mask(Omega_policy)[, , kk])
  })
)

cat("\nStaircase checks for Omega observation mask by aligned outcome layer:\n")
print(staircase_tbl)

cat("\nRunning target-layer-only bootstrap for deaths target layer\n")
cat("B =", B, "\n")

bootstrap_results = map_dfr(
  seq_len(B),
  function(b) {
    
    if (b %% 25 == 0) {
      cat("Bootstrap replication", b, "of", B, "\n")
    }
    
    tryCatch(
      {
        run_target_layer_only_bootstrap_one(
          b = b,
          Y_base = Y,
          Omega_policy_base = Omega_policy,
          seed_base = bootstrap_seed,
          point_results = point_results
        )
      },
      error = function(e) {
        warning(
          paste0(
            "Bootstrap replication ", b,
            " failed: ", conditionMessage(e)
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
)

# ===============================================================
# 12. Bootstrap confidence intervals
#       point estimate +/- 1.96 * bootstrap SE
# ===============================================================

make_results_with_ci = function(point_results, bootstrap_results) {
  
  bootstrap_se_long = bootstrap_results %>%
    select(
      bootstrap_id,
      setting,
      target_outcome,
      target_policy,
      target_k,
      r,
      functional,
      all_of(ci_vars)
    ) %>%
    pivot_longer(
      cols = all_of(ci_vars),
      names_to = "quantity",
      values_to = "boot_value"
    ) %>%
    group_by(setting, target_outcome, target_policy, target_k, r, functional, quantity) %>%
    summarize(
      boot_mean = mean(boot_value, na.rm = TRUE),
      boot_se = sd(boot_value, na.rm = TRUE),
      n_boot_nonmissing = sum(!is.na(boot_value)),
      .groups = "drop"
    )
  
  point_long = point_results %>%
    select(
      setting,
      target_outcome,
      target_policy,
      target_k,
      r,
      functional,
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
      by = c("setting", "target_outcome", "target_policy", "target_k", "r", "functional", "quantity")
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
    ) %>%
    arrange(setting, target_outcome, r, functional, quantity)
}

results_with_ci = make_results_with_ci(
  point_results = point_results,
  bootstrap_results = bootstrap_results
)

cat("\n============================================================\n")
cat("Point estimates with bootstrap 95% confidence intervals\n")
cat("============================================================\n")
print(results_with_ci, n = Inf)

# ===============================================================
# 13. Summary table
# ===============================================================

make_summary_table = function(point_results, results_with_ci, main_r = r_grid[1]) {
  point_results %>%
    filter(r == main_r) %>%
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
    mutate(
      across(
        c(Psi0, Psi0_matrix, Psi1, Delta_h, Delta_h_matrix, tensor_minus_matrix),
        ~ round(.x, 6)
      )
    ) %>%
    arrange(setting, target_outcome, functional)
}

main_summary = make_summary_table(
  point_results = point_results,
  results_with_ci = results_with_ci,
  main_r = r_grid[1]
)

cat("\n============================================================\n")
cat("Main summary for r = ", r_grid[1], "\n", sep = "")
cat("============================================================\n")
print(main_summary, n = Inf)

# ===============================================================
# 14. Save results
# ===============================================================

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
    target_outcome = target_outcome,
    target_policy = target_policy,
    target_k = target_k,
    r_grid = r_grid,
    tau = tau,
    functionals = functionals,
    B = B,
    bootstrap_seed = bootstrap_seed,
    stratify_target_bootstrap = stratify_target_bootstrap,
    ci_vars = ci_vars,
    project_dir = getwd()
  ),
  layer_map = layer_map,
  point_results = point_results,
  bootstrap_results = bootstrap_results,
  results_with_ci = results_with_ci,
  main_summary = main_summary,
  staircase_checks = staircase_tbl,
  target_row_strata = strata,
  Y = Y,
  Omega_policy = Omega_policy
)

saveRDS(
  analysis_results,
  file.path(results_dir, paste0("oxford_deaths_target_layer_bootstrap_", slice_label, ".rds"))
)

write_csv(
  layer_map,
  file.path(results_dir, paste0("oxford_deaths_layer_map_", slice_label, ".csv"))
)

write_csv(
  point_results,
  file.path(results_dir, paste0("oxford_deaths_point_results_", slice_label, ".csv"))
)

write_csv(
  bootstrap_results,
  file.path(results_dir, paste0("oxford_deaths_bootstrap_results_", slice_label, ".csv"))
)

write_csv(
  results_with_ci,
  file.path(results_dir, paste0("oxford_deaths_results_with_ci_", slice_label, ".csv"))
)

write_csv(
  main_summary,
  file.path(results_dir, paste0("oxford_deaths_main_summary_r", r_grid[1], "_", slice_label, ".csv"))
)

write_csv(
  staircase_tbl,
  file.path(results_dir, paste0("oxford_deaths_staircase_checks_", slice_label, ".csv"))
)

cat("\n============================================================\n")
cat("Done. Results saved in:\n")
cat(results_dir, "\n")
cat("============================================================\n")


print(main_summary, n = Inf, width = Inf)
print(results_with_ci, n = Inf, width = Inf)


# ===============================================================
# Plot
# ===============================================================

method_colors = c(
  "Tensor" = "#1f77b4", 
  "Matrix" = "#ff7f0e"  
)

get_symmetric_ylim = function(dat, pad_mult = 1.08) {
  y_vals = c(dat$ci_low, dat$ci_high, dat$point_estimate)
  y_vals = y_vals[is.finite(y_vals)]
  
  if (length(y_vals) == 0) {
    return(c(-1, 1))
  }
  
  y_abs_max = max(abs(y_vals), na.rm = TRUE)
  
  if (!is.finite(y_abs_max) || y_abs_max == 0) {
    return(c(-1, 1))
  }
  
  y_abs_max = y_abs_max * pad_mult
  c(-y_abs_max, y_abs_max)
}

main_r = r_grid[1]

psi_delta_plot_data = results_with_ci %>%
  filter(
    r == main_r,
    functional %in% c("ATE", "Trend"),
    quantity %in% c("Psi0", "Psi0_matrix", "Delta_h", "Delta_h_matrix")
  ) %>%
  mutate(
    estimand = case_when(
      quantity %in% c("Psi0", "Psi0_matrix") ~ "Psi0",
      quantity %in% c("Delta_h", "Delta_h_matrix") ~ "Delta"
    ),
    method = case_when(
      quantity %in% c("Psi0", "Delta_h") ~ "Tensor",
      quantity %in% c("Psi0_matrix", "Delta_h_matrix") ~ "Matrix"
    ),
    method = factor(method, levels = c("Tensor", "Matrix")),
    x_group = paste(functional, estimand, sep = "\n"),
    x_group = factor(
      x_group,
      levels = c(
        "ATE\nPsi0",
        "ATE\nDelta",
        "Trend\nPsi0",
        "Trend\nDelta"
      )
    ),
    x_num = as.numeric(x_group),
    method_offset = if_else(method == "Tensor", -0.13, 0.13),
    x_pos = x_num + method_offset
  )

y_limits = get_symmetric_ylim(psi_delta_plot_data)

psi_delta_plot = ggplot(psi_delta_plot_data) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = x_pos,
      xend = x_pos,
      y = ci_low,
      yend = ci_high,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = x_pos - 0.055,
      xend = x_pos + 0.055,
      y = ci_low,
      yend = ci_low,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_segment(
    aes(
      x = x_pos - 0.055,
      xend = x_pos + 0.055,
      y = ci_high,
      yend = ci_high,
      color = method
    ),
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      x = x_pos,
      y = point_estimate,
      color = method,
      shape = method
    ),
    size = 2.6,
    na.rm = TRUE
  ) +
  coord_cartesian(ylim = y_limits) +
  scale_x_continuous(
    breaks = seq_along(levels(psi_delta_plot_data$x_group)),
    labels = levels(psi_delta_plot_data$x_group),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 5),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c("Tensor" = 16, "Matrix" = 17)) +
  labs(
    title = paste0("CovidOx | Deaths target layer: Psi0 and Delta"),
    subtitle = paste0(
      "Rank r = ", main_r,
      "; B = ", B,
      "; intervals are point estimate +/- 1.96 x bootstrap SE"
    ),
    x = "",
    y = "Estimated quantity with bootstrap-SE 95% CI",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9)
  )

print(psi_delta_plot)

ggsave(
  filename = file.path(
    results_dir,
    paste0("oxford_deaths_psi0_delta_ATE_Trend_r", main_r, ".png")
  ),
  plot = psi_delta_plot,
  width = 11,
  height = 7.5,
  dpi = 320
)