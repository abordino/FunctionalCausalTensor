setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CastleDoctrine")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

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
  mutate(
    row_index = row_number()
  )

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

Y = array(
  NA_real_,
  dim = c(n_states, n_years, n_crimes),
  dimnames = list(
    state = state_order_vec,
    year = as.character(year_order),
    crime = crime_order
  )
)

D = Y

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
  
  D[, , cc] = tmp %>%
    select(state_id, year, treatment_value) %>%
    pivot_wider(names_from = year, values_from = treatment_value) %>%
    arrange(state_id) %>%
    select(-state_id) %>%
    as.matrix()
}

storage.mode(Y) = "numeric"
storage.mode(D) = "numeric"

policy_on = D == 1

rn = dimnames(Y)$state
cn = dimnames(Y)$year
ln = dimnames(Y)$crime

cat("\nCreated tensors:\n")
cat("Y:", paste(dim(Y), collapse = " x "), "\n")
cat("D:", paste(dim(D), collapse = " x "), "\n")
cat("Dimension order: state x year x crime\n")

make_status_integer_outcome_fill = function(d) {
  d %>%
    mutate(
      outcome_integer = if_else(
        is.na(outcome_value),
        NA_integer_,
        as.integer(round(outcome_value))
      ),
      fill_key = case_when(
        is.na(outcome_integer) ~ "Missing",
        treatment_status == "Off" ~ paste0("Off_", outcome_integer),
        treatment_status == "On"  ~ paste0("On_", outcome_integer)
      )
    )
}

make_status_integer_fill_values = function(d) {
  
  vals = d %>%
    filter(!is.na(outcome_value)) %>%
    mutate(outcome_integer = as.integer(round(outcome_value))) %>%
    distinct(outcome_integer) %>%
    arrange(outcome_integer) %>%
    pull(outcome_integer)
  
  off_cols = grDevices::colorRampPalette(c("#deebf7", "#08519c"))(length(vals))
  on_cols  = grDevices::colorRampPalette(c("#fee0d2", "#a50f15"))(length(vals))
  
  c(
    setNames(off_cols, paste0("Off_", vals)),
    setNames(on_cols,  paste0("On_", vals)),
    "Missing" = "grey85"
  )
}

pick_breaks = function(y_vals, n = 25) {
  ys = sort(unique(y_vals))
  if (length(ys) <= n) ys else ys[round(seq(1, length(ys), length.out = n))]
}

d_plot = castle_tensor_long %>%
  make_status_integer_outcome_fill()

brks = pick_breaks(d_plot$y_plot, 25)

labs_tbl = d_plot %>%
  distinct(y_plot, state_id) %>%
  filter(y_plot %in% brks) %>%
  arrange(y_plot)

fig_castle_tensor = ggplot(
  d_plot,
  aes(x = year, y = y_plot, fill = fill_key)
) +
  geom_tile(width = 0.95, height = 0.95) +
  facet_wrap(~ crime_label, ncol = 2) +
  scale_fill_manual(
    values = make_status_integer_fill_values(d_plot),
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = labs_tbl$y_plot,
    labels = labs_tbl$state_id,
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = sort(unique(d_plot$year)),
    labels = sort(unique(d_plot$year)),
    expand = c(0, 0)
  ) +
  labs(
    title = "Castle Doctrine | state-year-crime tensor and policy status",
    subtitle = paste0(
      "Each panel is one layer of Y[state, year, crime]. ",
      "Never adopters are shown on top. ",
      "Blue = policy Off, red = policy On. ",
      "Darker shade = higher rounded log crime rate."
    ),
    x = "Year",
    y = "State, never adopters on top"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(10, 20, 10, 10)
  )

print(fig_castle_tensor)

source("../../bilinearMatrixStaggered.R")
source("../../bilinearMatrixStaggeredPsi.R")
source("../../pluginPsi_c1.R")

r_grid = c(1, 2, 3, 4, 5, 7, 10)
tau = 0.01

target_crimes = c(
  robbery = "l_robbery",
  assault = "l_assault"
)

local_states = c("Montana", "Texas", "Florida")

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

make_staggered_parts_from_D = function(D, layer = 1) {
  
  Dmat = D[, , layer]
  years = as.integer(colnames(Dmat))
  
  first_on = apply(Dmat, 1, function(drow) {
    on_years = years[drow == 1]
    if (length(on_years) == 0) Inf else min(on_years)
  })
  
  adopt_years = sort(unique(first_on[is.finite(first_on)]))
  m = length(adopt_years)
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

staggered = make_staggered_parts_from_D(D, layer = 1)

N_part = staggered$N_sizes
T_part = staggered$T_sizes

N_parts_matrix = list(N_part)
T_parts_matrix = list(T_part)

Y0 = Y
Y0[D == 1] = NA_real_

cat("\nAdoption years:\n")
print(staggered$adopt_years)

cat("\nMatrix staggered parts:\n")
cat("N_part:", paste(N_part, collapse = ", "), "\n")
cat("T_part:", paste(T_part, collapse = ", "), "\n")

estimate_one = function(r_value, crime_name, functional, local_state = NA_character_) {
  
  k = match(target_crimes[[crime_name]], dimnames(Y)$crime)
  
  Y_mat0 = Y0[, , k]
  Y_mat1 = Y[, , k]
  
  rownames(Y_mat0) = rn
  colnames(Y_mat0) = cn
  
  rownames(Y_mat1) = rn
  colnames(Y_mat1) = cn
  
  row_index = if (functional == "Local") {
    match(local_state, rn)
  } else {
    NULL
  }
  
  psi0 = bilinearMatrixStaggeredPsi(
    Y_mat = Y_mat0,
    r = r_value,
    N_part = N_part,
    T_part = T_part,
    tau = tau,
    functional = functional,
    eta = if (functional == "RowHet") eta else NULL,
    row_index = row_index
  )
  
  psi1 = pluginPsi_c1(
    Y = Y_mat1,
    k = 1,
    N_parts = N_parts_matrix,
    T_parts = T_parts_matrix,
    functional = functional,
    eta = if (functional == "RowHet") eta else NULL,
    row_index = row_index
  )
  
  tibble(
    r = r_value,
    crime = crime_name,
    layer = target_crimes[[crime_name]],
    functional = if_else(
      functional == "Local",
      paste0("Local-", local_state),
      functional
    ),
    local_state = if_else(functional == "Local", local_state, NA_character_),
    Psi0 = as.numeric(psi0),
    Psi1 = as.numeric(psi1),
    Delta_h = as.numeric(psi1 - psi0)
  )
}

run_one_r = function(r_value) {
  
  nonlocal_results = map_dfr(
    names(target_crimes),
    function(crime_name) {
      map_dfr(
        c("ATE", "RowHet", "Trend"),
        function(f) estimate_one(r_value, crime_name, f)
      )
    }
  )
  
  local_results = map_dfr(
    names(target_crimes),
    function(crime_name) {
      map_dfr(
        local_states,
        function(st) estimate_one(r_value, crime_name, "Local", local_state = st)
      )
    }
  )
  
  results = bind_rows(nonlocal_results, local_results) %>%
    group_by(functional, local_state) %>%
    mutate(
      C_0_h = Psi0[crime == "robbery"] - Psi0[crime == "assault"],
      C_1_h = Psi1[crime == "robbery"] - Psi1[crime == "assault"],
      C_Delta_h = Delta_h[crime == "robbery"] - Delta_h[crime == "assault"]
    ) %>%
    ungroup() %>%
    mutate(
      across(
        c(Psi0, Psi1, Delta_h, C_0_h, C_1_h, C_Delta_h),
        ~ round(.x, 6)
      )
    ) %>%
    arrange(functional, crime)
  
  cat("\n============================================================\n")
  cat("Results for r =", r_value, "\n")
  cat("============================================================\n")
  print(results, n = Inf)
  
  results
}

all_results_by_r = map(r_grid, run_one_r)
names(all_results_by_r) = paste0("r_", r_grid)

all_results = bind_rows(all_results_by_r)

cat("\n============================================================\n")
cat("Combined results across all r values\n")
cat("============================================================\n")
print(all_results, n = Inf)

View(D[, , "l_homicide"])
View(Y[, , "l_homicide"])