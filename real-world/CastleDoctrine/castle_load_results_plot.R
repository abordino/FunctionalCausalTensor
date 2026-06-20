PROJECT_DIR = "~/Documents/phd/projects/causalMatrix/code/real-world/CastleDoctrine"
if (dir.exists(path.expand(PROJECT_DIR))) {
  setwd(path.expand(PROJECT_DIR))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ===============================================================
# 0. Sttings par
# ===============================================================

results_dir = "results/masked_nonmasked_rank3_bootstrap"

target_specs = tibble::tribble(
  ~target_name, ~target_label,         ~result_file,
  "robbery",    "Robbery rate, log",   file.path(results_dir, "castle_masked_nonmasked_bootstrap_results_robbery.rds"),
  "assault",      "Aggravated assault rate, log", file.path(results_dir, "castle_masked_nonmasked_bootstrap_results_assault.rds")
)

if (nrow(target_specs) != 2) {
  stop("target_specs must contain exactly two target rows.")
}

contrast_name = paste0(target_specs$target_name[1], "_minus_", target_specs$target_name[2])
contrast_label = paste0(target_specs$target_label[1], " - ", target_specs$target_label[2])

method_colors = c(
  "Tensor" = "#1f77b4",  
  "Matrix" = "#ff7f0e"  
)

output_dir = file.path(results_dir, "tables_and_figures")
fig_dir = file.path(output_dir, "figures")
table_dir = file.path(output_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

missing_files = target_specs$result_file[!file.exists(target_specs$result_file)]
if (length(missing_files) > 0) {
  stop(
    "Could not find saved result file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\nRun the bootstrap script once for each target and save with target-specific names."
  )
}

# ===============================================================
# 1. Load results
# ===============================================================

analysis_objects = target_specs %>%
  mutate(obj = map(result_file, readRDS))

metadata_1 = analysis_objects$obj[[1]]$metadata
functional_order = metadata_1$functional_order
rank_value = metadata_1$rank_value
B = metadata_1$B

cat("\nLoaded target result files:\n")
print(target_specs)

cat("\nTarget metadata in saved files:\n")
analysis_objects %>%
  transmute(
    target_name_user = target_name,
    target_label_user = target_label,
    target_name_saved = map_chr(obj, ~ .x$metadata$target_name),
    target_layer_saved = map_chr(obj, ~ .x$metadata$target_layer),
    target_label_saved = map_chr(obj, ~ .x$metadata$target_label),
    rank_value = map_dbl(obj, ~ .x$metadata$rank_value),
    B = map_dbl(obj, ~ .x$metadata$B)
  ) %>%
  print()

if (length(unique(map_dbl(analysis_objects$obj, ~ .x$metadata$rank_value))) != 1) {
  warning("The two result files have different rank values.")
}

if (length(unique(map_dbl(analysis_objects$obj, ~ .x$metadata$B))) != 1) {
  warning("The two result files have different bootstrap sizes B.")
}

# ===============================================================
# 2. Helpers
# ===============================================================

method_from_quantity = function(quantity) {
  dplyr::recode(
    quantity,
    Delta_h = "Tensor",
    Delta_h_matrix = "Matrix",
    C_Delta_h = "Tensor",
    C_Delta_h_matrix = "Matrix",
    Psi0 = "Tensor",
    Psi0_matrix = "Matrix",
    C_Psi0 = "Tensor",
    C_Psi0_matrix = "Matrix",
    Psi1 = "Tensor",
    Psi1_matrix = "Matrix",
    C_Psi1 = "Tensor",
    C_Psi1_matrix = "Matrix",
    .default = quantity
  )
}

quantity_from_method = function(method, quantity_prefix) {
  method_chr = as.character(method)
  
  if (quantity_prefix == "Delta_h") {
    return(
      dplyr::recode(
        method_chr,
        Tensor = "C_Delta_h",
        Matrix = "C_Delta_h_matrix",
        .default = method_chr
      )
    )
  }
  
  if (quantity_prefix == "Psi0") {
    return(
      dplyr::recode(
        method_chr,
        Tensor = "C_Psi0",
        Matrix = "C_Psi0_matrix",
        .default = method_chr
      )
    )
  }
  
  if (quantity_prefix == "Psi1") {
    return(
      dplyr::recode(
        method_chr,
        Tensor = "C_Psi1",
        Matrix = "C_Psi1_matrix",
        .default = method_chr
      )
    )
  }
  
  stop("Unknown quantity_prefix: ", quantity_prefix)
}

make_design_short = function(design_id, design_label) {
  case_when(
    str_detect(design_id, "masked_full5_") ~ "Masked: 5 full rows",
    str_detect(design_id, "masked_full10_") ~ "Masked: 10 full rows",
    str_detect(design_id, "masked_full15_") ~ "Masked: 15 full rows",
    design_id == "nonmasked_original" ~ "Non-masked",
    TRUE ~ design_label
  )
}

make_slug = function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("(^_+|_+$)", "")
}

make_local_join = function(local_state) {
  local_state_chr = as.character(local_state)
  if_else(is.na(local_state_chr), "__NA__", local_state_chr)
}

join_keys = c(
  "design_id",
  "r",
  "functional",
  "local_state_join",
  "method"
)

quantity_specs = tibble::tribble(
  ~estimand_type, ~quantity_prefix, ~tensor_quantity, ~matrix_quantity, ~effect_prefix, ~contrast_prefix,
  "Delta",       "Delta_h",       "Delta_h",       "Delta_h_matrix", "Delta",     "Contrast",
  "Psi0",        "Psi0",          "Psi0",          "Psi0_matrix",    "Psi0",      "Contrast Psi0",
  "Psi1",        "Psi1",          "Psi1",          "Psi1",           "Psi1",      "Contrast Psi1"
)

extract_point_quantity = function(obj, target_name, tensor_quantity, matrix_quantity) {
  
  base_dat = obj$point_results %>%
    filter(crime == target_name) %>%
    select(
      design_id,
      design_label,
      r,
      crime,
      layer,
      functional,
      local_state,
      all_of(unique(c(tensor_quantity, matrix_quantity)))
    )
  
  if (!identical(tensor_quantity, matrix_quantity)) {
    return(
      base_dat %>%
        pivot_longer(
          cols = all_of(c(tensor_quantity, matrix_quantity)),
          names_to = "quantity",
          values_to = "point_estimate"
        ) %>%
        mutate(
          method = method_from_quantity(quantity),
          method = factor(method, levels = c("Tensor", "Matrix")),
          local_state_join = make_local_join(local_state)
        )
    )
  }
  
  bind_rows(
    base_dat %>%
      transmute(
        design_id,
        design_label,
        r,
        crime,
        layer,
        functional,
        local_state,
        quantity = tensor_quantity,
        point_estimate = .data[[tensor_quantity]],
        method = "Tensor"
      ),
    base_dat %>%
      transmute(
        design_id,
        design_label,
        r,
        crime,
        layer,
        functional,
        local_state,
        quantity = paste0(tensor_quantity, "_matrix"),
        point_estimate = .data[[tensor_quantity]],
        method = "Matrix"
      )
  ) %>%
    mutate(
      method = factor(method, levels = c("Tensor", "Matrix")),
      local_state_join = make_local_join(local_state)
    )
}

extract_boot_quantity = function(obj, target_name, tensor_quantity, matrix_quantity) {
  
  base_dat = obj$bootstrap_results %>%
    filter(crime == target_name) %>%
    select(
      design_id,
      design_label,
      bootstrap_id,
      r,
      crime,
      layer,
      functional,
      local_state,
      all_of(unique(c(tensor_quantity, matrix_quantity)))
    )
  

  if (!identical(tensor_quantity, matrix_quantity)) {
    return(
      base_dat %>%
        pivot_longer(
          cols = all_of(c(tensor_quantity, matrix_quantity)),
          names_to = "quantity",
          values_to = "boot_value"
        ) %>%
        mutate(
          method = method_from_quantity(quantity),
          method = factor(method, levels = c("Tensor", "Matrix")),
          local_state_join = make_local_join(local_state)
        )
    )
  }
  
  bind_rows(
    base_dat %>%
      transmute(
        design_id,
        design_label,
        bootstrap_id,
        r,
        crime,
        layer,
        functional,
        local_state,
        quantity = tensor_quantity,
        boot_value = .data[[tensor_quantity]],
        method = "Tensor"
      ),
    base_dat %>%
      transmute(
        design_id,
        design_label,
        bootstrap_id,
        r,
        crime,
        layer,
        functional,
        local_state,
        quantity = paste0(tensor_quantity, "_matrix"),
        boot_value = .data[[tensor_quantity]],
        method = "Matrix"
      )
  ) %>%
    mutate(
      method = factor(method, levels = c("Tensor", "Matrix")),
      local_state_join = make_local_join(local_state)
    )
}

extract_individual_plot_data = function(
    obj,
    target_name,
    target_label,
    estimand_type,
    tensor_quantity,
    matrix_quantity,
    effect_prefix
) {
  dat = obj$results_with_ci %>%
    filter(
      crime == target_name,
      quantity %in% unique(c(tensor_quantity, matrix_quantity))
    )
  
  if (!identical(tensor_quantity, matrix_quantity)) {
    dat = dat %>%
      mutate(
        estimand_type = estimand_type,
        estimand_id = target_name,
        estimand_label = target_label,
        method = method_from_quantity(quantity),
        method = factor(method, levels = c("Tensor", "Matrix")),
        design_short = make_design_short(design_id, design_label),
        effect_label = paste0(effect_prefix, ": ", target_label)
      )
  } else {
    dat = bind_rows(
      dat %>%
        mutate(
          quantity = tensor_quantity,
          method = "Tensor"
        ),
      dat %>%
        mutate(
          quantity = paste0(tensor_quantity, "_matrix"),
          method = "Matrix"
        )
    ) %>%
      mutate(
        estimand_type = estimand_type,
        estimand_id = target_name,
        estimand_label = target_label,
        method = factor(method, levels = c("Tensor", "Matrix")),
        design_short = make_design_short(design_id, design_label),
        effect_label = paste0(effect_prefix, ": ", target_label)
      )
  }
  
  if ("ci_type" %in% names(dat)) {
    dat$ci_type = as.character(dat$ci_type)
  } else {
    dat$ci_type = "from_saved_results"
  }
  
  dat
}

make_contrast_results_with_ci = function(
    obj_1,
    obj_2,
    name_1,
    name_2,
    label_1,
    label_2,
    estimand_type,
    quantity_prefix,
    tensor_quantity,
    matrix_quantity,
    contrast_prefix
) {
  
  point_1 = extract_point_quantity(
    obj = obj_1,
    target_name = name_1,
    tensor_quantity = tensor_quantity,
    matrix_quantity = matrix_quantity
  ) %>%
    rename(
      point_1 = point_estimate,
      crime_1 = crime,
      layer_1 = layer,
      design_label_1 = design_label
    )
  
  point_2 = extract_point_quantity(
    obj = obj_2,
    target_name = name_2,
    tensor_quantity = tensor_quantity,
    matrix_quantity = matrix_quantity
  ) %>%
    rename(
      point_2 = point_estimate,
      crime_2 = crime,
      layer_2 = layer,
      design_label_2 = design_label
    )
  
  contrast_point = point_1 %>%
    inner_join(
      point_2,
      by = join_keys,
      suffix = c("_1", "_2")
    ) %>%
    mutate(
      point_estimate = point_1 - point_2,
      quantity = quantity_from_method(method, quantity_prefix = quantity_prefix),
      design_label = coalesce(design_label_1, design_label_2),
      local_state = na_if(local_state_join, "__NA__"),
      crime = contrast_name,
      layer = paste0(layer_1, " - ", layer_2)
    ) %>%
    select(
      design_id,
      design_label,
      r,
      crime,
      layer,
      functional,
      local_state,
      quantity,
      method,
      point_estimate
    )
  
  boot_1 = extract_boot_quantity(
    obj = obj_1,
    target_name = name_1,
    tensor_quantity = tensor_quantity,
    matrix_quantity = matrix_quantity
  ) %>%
    rename(
      boot_1 = boot_value,
      crime_1 = crime,
      layer_1 = layer,
      design_label_1 = design_label
    )
  
  boot_2 = extract_boot_quantity(
    obj = obj_2,
    target_name = name_2,
    tensor_quantity = tensor_quantity,
    matrix_quantity = matrix_quantity
  ) %>%
    rename(
      boot_2 = boot_value,
      crime_2 = crime,
      layer_2 = layer,
      design_label_2 = design_label
    )
  
  contrast_boot = boot_1 %>%
    inner_join(
      boot_2,
      by = c("bootstrap_id", join_keys),
      suffix = c("_1", "_2")
    ) %>%
    mutate(
      boot_value = boot_1 - boot_2,
      quantity = quantity_from_method(method, quantity_prefix = quantity_prefix),
      design_label = coalesce(design_label_1, design_label_2),
      local_state = na_if(local_state_join, "__NA__"),
      crime = contrast_name,
      layer = paste0(layer_1, " - ", layer_2)
    ) %>%
    select(
      design_id,
      design_label,
      bootstrap_id,
      r,
      crime,
      layer,
      functional,
      local_state,
      quantity,
      method,
      boot_value
    )
  
  contrast_se = contrast_boot %>%
    group_by(
      design_id,
      design_label,
      r,
      crime,
      layer,
      functional,
      local_state,
      quantity,
      method
    ) %>%
    summarize(
      boot_mean = mean(boot_value, na.rm = TRUE),
      boot_se = sd(boot_value, na.rm = TRUE),
      n_boot_nonmissing = sum(!is.na(boot_value)),
      .groups = "drop"
    )
  
  contrast_point %>%
    left_join(
      contrast_se,
      by = c(
        "design_id",
        "design_label",
        "r",
        "crime",
        "layer",
        "functional",
        "local_state",
        "quantity",
        "method"
      )
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
      estimate_ci = if_else(
        is.na(ci_low) | is.na(ci_high),
        sprintf("%.4f [NA, NA]", point_estimate),
        sprintf("%.4f [%.4f, %.4f]", point_estimate, ci_low, ci_high)
      ),
      estimand_type = estimand_type,
      estimand_id = contrast_name,
      estimand_label = contrast_label,
      effect_label = paste0(contrast_prefix, ": ", contrast_label),
      design_short = make_design_short(design_id, design_label),
      method = factor(method, levels = c("Tensor", "Matrix"))
    )
}

get_symmetric_xlim = function(dat, pad_mult = 1.08) {
  x_vals = c(dat$ci_low, dat$ci_high, dat$point_estimate)
  x_vals = x_vals[is.finite(x_vals)]
  
  if (length(x_vals) == 0) {
    return(c(-1, 1))
  }
  
  x_abs_max = max(abs(x_vals), na.rm = TRUE)
  
  if (!is.finite(x_abs_max) || x_abs_max == 0) {
    return(c(-1, 1))
  }
  
  x_abs_max = x_abs_max * pad_mult
  c(-x_abs_max, x_abs_max)
}

print_tsv_to_console = function(dat) {
  dat_print = dat %>%
    mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ replace_na(.x, "")))
  
  cat(paste(names(dat_print), collapse = "\t"), "\n", sep = "")
  purrr::pwalk(dat_print, function(...) {
    cat(paste(c(...), collapse = "\t"), "\n", sep = "")
  })
}

# ===============================================================
# 3. Individual target-layer effects for Delta, Psi0, and Psi1
# ===============================================================

individual_plot_data = pmap_dfr(
  list(
    obj = rep(analysis_objects$obj, each = nrow(quantity_specs)),
    target_name = rep(analysis_objects$target_name, each = nrow(quantity_specs)),
    target_label = rep(analysis_objects$target_label, each = nrow(quantity_specs)),
    estimand_type = rep(quantity_specs$estimand_type, times = nrow(analysis_objects)),
    tensor_quantity = rep(quantity_specs$tensor_quantity, times = nrow(analysis_objects)),
    matrix_quantity = rep(quantity_specs$matrix_quantity, times = nrow(analysis_objects)),
    effect_prefix = rep(quantity_specs$effect_prefix, times = nrow(analysis_objects))
  ),
  extract_individual_plot_data
)

# ===============================================================
# 4. Contrast effects for Delta, Psi0, and Psi1
# ===============================================================

obj_1 = analysis_objects$obj[[1]]
obj_2 = analysis_objects$obj[[2]]

name_1 = analysis_objects$target_name[[1]]
name_2 = analysis_objects$target_name[[2]]
label_1 = analysis_objects$target_label[[1]]
label_2 = analysis_objects$target_label[[2]]

contrast_results_with_ci = pmap_dfr(
  list(
    estimand_type = quantity_specs$estimand_type,
    quantity_prefix = quantity_specs$quantity_prefix,
    tensor_quantity = quantity_specs$tensor_quantity,
    matrix_quantity = quantity_specs$matrix_quantity,
    contrast_prefix = quantity_specs$contrast_prefix
  ),
  function(estimand_type, quantity_prefix, tensor_quantity, matrix_quantity, contrast_prefix) {
    make_contrast_results_with_ci(
      obj_1 = obj_1,
      obj_2 = obj_2,
      name_1 = name_1,
      name_2 = name_2,
      label_1 = label_1,
      label_2 = label_2,
      estimand_type = estimand_type,
      quantity_prefix = quantity_prefix,
      tensor_quantity = tensor_quantity,
      matrix_quantity = matrix_quantity,
      contrast_prefix = contrast_prefix
    )
  }
)

# ===============================================================
# 5. Combine individual effects and contrasts
# ===============================================================

plot_results_data = bind_rows(
  individual_plot_data %>%
    select(
      design_id,
      design_label,
      design_short,
      r,
      crime,
      layer,
      functional,
      local_state,
      quantity,
      method,
      point_estimate,
      boot_mean,
      boot_se,
      ci_low,
      ci_high,
      n_boot_nonmissing,
      ci_type,
      estimand_type,
      estimand_id,
      estimand_label,
      effect_label,
      estimate_ci
    ),
  contrast_results_with_ci %>%
    select(
      design_id,
      design_label,
      design_short,
      r,
      crime,
      layer,
      functional,
      local_state,
      quantity,
      method,
      point_estimate,
      boot_mean,
      boot_se,
      ci_low,
      ci_high,
      n_boot_nonmissing,
      ci_type,
      estimand_type,
      estimand_id,
      estimand_label,
      effect_label,
      estimate_ci
    )
) %>%
  mutate(
    design_short = factor(
      design_short,
      levels = c(
        "Masked: 5 full rows",
        "Masked: 10 full rows",
        "Masked: 15 full rows",
        "Non-masked"
      )
    ),
    effect_label = factor(
      effect_label,
      levels = c(
        paste0("Delta: ", label_1),
        paste0("Delta: ", label_2),
        paste0("Contrast: ", contrast_label),
        paste0("Psi0: ", label_1),
        paste0("Psi0: ", label_2),
        paste0("Contrast Psi0: ", contrast_label),
        paste0("Psi1: ", label_1),
        paste0("Psi1: ", label_2),
        paste0("Contrast Psi1: ", contrast_label)
      )
    ),
    functional = factor(functional, levels = functional_order),
    functional_num = as.numeric(factor(functional, levels = rev(functional_order))),
    method_offset = if_else(method == "Tensor", -0.13, 0.13),
    y_pos = functional_num + method_offset
  )

plot_delta_data = plot_results_data %>%
  filter(estimand_type == "Delta")

plot_psi0_data = plot_results_data %>%
  filter(estimand_type == "Psi0")

plot_psi1_data = plot_results_data %>%
  filter(estimand_type == "Psi1")

# ===============================================================
# 6. Plots for Delta, Psi0, and Psi1
# ===============================================================

functional_breaks = seq_along(rev(functional_order))
functional_labels = rev(functional_order)

make_one_effect_plot = function(effect_name, plot_data) {
  dat_one = plot_data %>%
    filter(as.character(effect_label) == as.character(effect_name))
  
  if (nrow(dat_one) == 0) {
    warning("No plot data found for effect: ", effect_name)
    return(tibble(effect_label = as.character(effect_name), plot_png = NA_character_))
  }
  
  effect_stub = make_slug(as.character(effect_name))
  
  x_limits = get_symmetric_xlim(dat_one)
  
  fig = ggplot(dat_one) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_segment(
      aes(x = ci_low, xend = ci_high, y = y_pos, yend = y_pos, color = method),
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    geom_segment(
      aes(x = ci_low, xend = ci_low, y = y_pos - 0.055, yend = y_pos + 0.055, color = method),
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    geom_segment(
      aes(x = ci_high, xend = ci_high, y = y_pos - 0.055, yend = y_pos + 0.055, color = method),
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    geom_point(
      aes(x = point_estimate, y = y_pos, color = method, shape = method),
      size = 2.6,
      na.rm = TRUE
    ) +
    facet_wrap(~ design_short, ncol = 2, scales = "fixed") +
    coord_cartesian(xlim = x_limits) +
    scale_x_continuous(
      breaks = scales::breaks_extended(n = 5),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_y_continuous(
      breaks = functional_breaks,
      labels = functional_labels,
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = c("Tensor" = 16, "Matrix" = 17)) +
    labs(
      title = paste0("Castle Doctrine | ", as.character(effect_name)),
      subtitle = paste0(
        "Rank r = ", rank_value,
        "; B = ", B,
        "; intervals are point estimate +/- 1.96 x bootstrap SE"
      ),
      x = expression("Estimated quantity with bootstrap-SE 95% CI"),
      y = "Functional type",
      color = "Method",
      shape = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  plot_png = file.path(
    fig_dir,
    paste0("castle_", effect_stub, "_", name_1, "_", name_2, ".png")
  )
  
  ggsave(plot_png, fig, width = 11, height = 7.5, dpi = 320)
  
  tibble(
    effect_label = as.character(effect_name),
    plot_png = plot_png
  )
}

plot_file_index = map_dfr(
  levels(plot_results_data$effect_label),
  ~ make_one_effect_plot(.x, plot_results_data)
)

print(plot_file_index, n = Inf)

# ===============================================================
# 7. Result tables
# ===============================================================

results_table_long = plot_results_data %>%
  transmute(
    estimand_type = estimand_type,
    estimand = as.character(effect_label),
    design_id = design_id,
    design = as.character(design_short),
    functional = as.character(functional),
    local_state = local_state,
    method = as.character(method),
    quantity = quantity,
    point_estimate = round(point_estimate, 6),
    ci_low = round(ci_low, 6),
    ci_high = round(ci_high, 6),
    estimate_ci = estimate_ci,
    boot_se = round(boot_se, 6),
    n_boot_nonmissing = n_boot_nonmissing,
    ci_type = ci_type
  ) %>%
  arrange(
    factor(estimand_type, levels = c("Delta", "Psi0", "Psi1")),
    factor(estimand, levels = levels(plot_results_data$effect_label)),
    factor(design, levels = levels(plot_results_data$design_short)),
    factor(functional, levels = functional_order),
    factor(method, levels = c("Tensor", "Matrix"))
  )

results_table_wide = results_table_long %>%
  select(
    estimand_type,
    estimand,
    design_id,
    design,
    functional,
    local_state,
    method,
    point_estimate,
    ci_low,
    ci_high,
    estimate_ci,
    boot_se,
    n_boot_nonmissing
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = c(
      point_estimate,
      ci_low,
      ci_high,
      estimate_ci,
      boot_se,
      n_boot_nonmissing
    ),
    names_glue = "{method}_{.value}"
  ) %>%
  arrange(
    factor(estimand_type, levels = c("Delta", "Psi0", "Psi1")),
    factor(estimand, levels = levels(plot_results_data$effect_label)),
    factor(design, levels = levels(plot_results_data$design_short)),
    factor(functional, levels = functional_order)
  )

write_csv(
  results_table_long,
  file.path(table_dir, paste0("castle_pointwise_ci_long_", name_1, "_", name_2, ".csv"))
)

write_tsv(
  results_table_long,
  file.path(table_dir, paste0("castle_pointwise_ci_long_", name_1, "_", name_2, ".tsv"))
)

write_csv(
  results_table_wide,
  file.path(table_dir, paste0("castle_pointwise_ci_wide_", name_1, "_", name_2, ".csv"))
)

write_tsv(
  results_table_wide,
  file.path(table_dir, paste0("castle_pointwise_ci_wide_", name_1, "_", name_2, ".tsv"))
)

cat("\n============================================================\n")
cat("Copy/paste table: wide format, one row per estimand x design x functional\n")
cat("============================================================\n")
print_tsv_to_console(results_table_wide)

cat("\n============================================================\n")
cat("Also saved tables to:\n")
cat(table_dir, "\n")
cat("Files:\n")
cat("- ", paste0("castle_pointwise_ci_long_", name_1, "_", name_2, ".csv"), "\n", sep = "")
cat("- ", paste0("castle_pointwise_ci_long_", name_1, "_", name_2, ".tsv"), "\n", sep = "")
cat("- ", paste0("castle_pointwise_ci_wide_", name_1, "_", name_2, ".csv"), "\n", sep = "")
cat("- ", paste0("castle_pointwise_ci_wide_", name_1, "_", name_2, ".tsv"), "\n", sep = "")
cat("============================================================\n")
