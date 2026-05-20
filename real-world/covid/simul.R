suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
})

project_dir = "~/Documents/phd/projects/causalMatrix/code/real-world/covid"
setwd(project_dir)

source_bilinear_tensor = "../../bilinearTensorStaggered.R"
source_bilinear_psi    = "../../bilinearTensorStaggeredPsi.R"
source_plugin_psi      = "../../pluginPsi_c1.R"

policy_col = "stringency_index"
policy_name_base = "Stringency index >= 50"
policy_code = "stringency_high"
policy_threshold = 50

lag_days = 14
use_log1p_deaths = FALSE

main_outcome_raw = paste0(
  "new_deaths_smoothed_per_million_lag",
  lag_days,
  "d"
)

main_outcome = if (isTRUE(use_log1p_deaths)) {
  paste0(
    "log1p_new_deaths_smoothed_per_million_lag",
    lag_days,
    "d"
  )
} else {
  main_outcome_raw
}

policy_name = paste0(
  policy_name_base,
  "; smoothed deaths per million lagged ",
  lag_days,
  " days"
)

make_plots = TRUE

r = 2
tau = 0.01
functionals_to_run = c("ATE", "Trend")

period_dates = tibble(
  period_index = 1:3,
  period_name = paste0("period_", period_index),
  window_start = as.Date(c(
    "2020-07-14",
    "2020-11-17",
    "2021-10-23"
  )),
  window_end = as.Date(c(
    "2020-10-11",
    "2021-02-14",
    "2022-01-20"
  ))
) %>%
  mutate(
    tensor_k = period_index,
    period_length_days = as.integer(window_end - window_start) + 1L,
    outcome_window_start = window_start + days(lag_days),
    outcome_window_end = window_end + days(lag_days),
    period_label = paste0(
      "Period ",
      period_index,
      ": ",
      format(window_start, "%Y-%m-%d"),
      " to ",
      format(window_end, "%Y-%m-%d"),
      " | Y = new_deaths_smoothed_per_million",
      " | lag = ",
      lag_days,
      " days"
    )
  ) %>%
  select(
    tensor_k,
    period_index,
    period_name,
    period_label,
    window_start,
    window_end,
    outcome_window_start,
    outcome_window_end,
    period_length_days
  )

if (length(unique(period_dates$period_length_days)) != 1L) {
  stop("Selected periods do not all have the same number of days.")
}

period_length_days = unique(period_dates$period_length_days)
n_periods = nrow(period_dates)
K_countries_per_period = 15

period_country_wide = tibble(
  period_index = 1:3,
  
  row_01 = c(
    "Japan",
    "Sudan",
    "Uzbekistan"
  ),
  
  row_02 = c(
    "Nicaragua",
    "Nicaragua",
    "Sudan"
  ),
  
  row_03 = c(
    "Senegal",
    "Afghanistan",
    "Burkina Faso"
  ),
  
  row_04 = c(
    "Austria",
    "Burkina Faso",
    "Mali"
  ),
  
  row_05 = c(
    "Switzerland",
    "Ghana",
    "Ethiopia"
  ),
  
  row_06 = c(
    "Yemen",
    "Japan",
    "Cameroon"
  ),
  
  row_07 = c(
    "Finland",
    "Haiti",
    "Saudi Arabia"
  ),
  
  row_08 = c(
    "Poland",
    "Niger",
    "Pakistan"
  ),
  
  row_09 = c(
    "Belarus",
    "Madagascar",
    "Tunisia"
  ),
  
  row_10 = c(
    "France",
    "Belarus",
    "Bangladesh"
  ),
  
  row_11 = c(
    "Croatia",
    "Democratic Republic of Congo",
    "Nepal"
  ),
  
  row_12 = c(
    "Bulgaria",
    "Senegal",
    "Ecuador"
  ),
  
  row_13 = c(
    "Tunisia",
    "United Arab Emirates",
    "Kenya"
  ),
  
  row_14 = c(
    "Jordan",
    "Sri Lanka",
    "Ghana"
  ),
  
  row_15 = c(
    "Ukraine",
    "Finland",
    "Netherlands"
  )
) %>%
  left_join(
    period_dates %>%
      select(
        period_index,
        tensor_k,
        period_name,
        period_label,
        window_start,
        window_end,
        outcome_window_start,
        outcome_window_end,
        period_length_days
      ),
    by = "period_index"
  ) %>%
  arrange(tensor_k)

selected_period_country_map = period_country_wide %>%
  pivot_longer(
    cols = starts_with("row_"),
    names_to = "row",
    values_to = "country"
  ) %>%
  mutate(
    country_slot = as.integer(str_remove(row, "row_")),
    tensor_row = country_slot
  ) %>%
  arrange(
    tensor_k,
    country_slot
  ) %>%
  mutate(
    row_state_name = paste0(
      country,
      " | row_",
      str_pad(country_slot, 2, pad = "0")
    )
  )

candidate_countries = sort(unique(selected_period_country_map$country))

cat("\n============================================================\n")
cat("Standalone fixed-period setup\n")
cat("============================================================\n")
cat("lag_days =", lag_days, "\n")
cat("Policy/stringency date = t\n")
cat("Death outcome date     = t + lag_days\n")
cat("Outcome used by estimator:", main_outcome, "\n")
cat("use_log1p_deaths =", use_log1p_deaths, "\n")

cat("\nPeriod map:\n")
print(
  period_dates %>%
    select(
      tensor_k,
      period_name,
      period_label,
      window_start,
      window_end,
      outcome_window_start,
      outcome_window_end
    )
)

cat("\nCountries used in selected fixed subset:\n")
print(candidate_countries)

cat("\nRows per selected period:\n")
print(
  selected_period_country_map %>%
    count(
      tensor_k,
      period_name,
      name = "n_countries"
    )
)

if (any(selected_period_country_map %>% count(tensor_k) %>% pull(n) != K_countries_per_period)) {
  stop("At least one selected period does not have exactly K_countries_per_period rows.")
}

owid_url =
  "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"

owid_raw = readr::read_csv(
  owid_url,
  show_col_types = FALSE,
  guess_max = 100000
)

missing_countries = setdiff(
  candidate_countries,
  unique(owid_raw$location)
)

cat("\nCountries requested but not found in OWID:\n")
print(missing_countries)

if (length(missing_countries) > 0L) {
  stop("Some requested countries are not present in OWID. Check country names.")
}

owid_base = owid_raw %>%
  mutate(
    date = as.Date(date),
    location = as.character(location)
  ) %>%
  filter(
    location %in% candidate_countries,
    date >= min(period_dates$window_start),
    date <= max(period_dates$window_end) + days(lag_days),
    !str_starts(iso_code, "OWID_")
  ) %>%
  transmute(
    country = location,
    iso_code,
    date,
    population = as.numeric(population),
    stringency_index = as.numeric(stringency_index),
    new_deaths_smoothed_per_million =
      as.numeric(new_deaths_smoothed_per_million)
  )

cat("\nLoaded fixed-subset OWID data.\n")
cat("Countries:", n_distinct(owid_base$country), "\n")
cat(
  "OWID loaded date range:",
  format(min(owid_base$date, na.rm = TRUE)),
  "to",
  format(max(owid_base$date, na.rm = TRUE)),
  "\n"
)

policy_data = owid_base %>%
  transmute(
    country,
    iso_code,
    policy_date = date,
    population,
    stringency_index,
    policy_value = stringency_index,
    policy_on = as.integer(
      !is.na(stringency_index) &
        stringency_index >= policy_threshold
    ),
    policy_status = case_when(
      is.na(stringency_index) ~ "Missing",
      policy_on == 1L ~ "On",
      TRUE ~ "Off"
    ),
    outcome_date = policy_date + days(lag_days)
  )

outcome_data = owid_base %>%
  transmute(
    country,
    outcome_date = date,
    deaths_raw_lagged = new_deaths_smoothed_per_million
  )

selected_matrix_data = selected_period_country_map %>%
  select(
    tensor_k,
    period_index,
    period_name,
    period_label,
    window_start,
    window_end,
    outcome_window_start,
    outcome_window_end,
    period_length_days,
    country,
    country_slot,
    tensor_row,
    row,
    row_state_name
  ) %>%
  inner_join(
    policy_data,
    by = "country"
  ) %>%
  filter(
    policy_date >= window_start,
    policy_date <= window_end
  ) %>%
  left_join(
    outcome_data,
    by = c("country", "outcome_date")
  ) %>%
  mutate(
    date = policy_date,
    relative_day = as.integer(policy_date - window_start) + 1L,
    outcome_relative_day = relative_day + lag_days,
    deaths_raw = deaths_raw_lagged,
    deaths_nonnegative = pmax(deaths_raw, 0),
    deaths_log1p = log1p(deaths_nonnegative),
    deaths_value = if (isTRUE(use_log1p_deaths)) {
      deaths_log1p
    } else {
      deaths_raw
    }
  ) %>%
  arrange(
    tensor_k,
    country_slot,
    relative_day
  )

bad_lag_alignment = selected_matrix_data %>%
  filter(outcome_date != policy_date + days(lag_days))

if (nrow(bad_lag_alignment) > 0L) {
  cat("\nBad lag-alignment rows:\n")
  print(bad_lag_alignment, n = Inf)
  stop("Some rows do not satisfy outcome_date = policy_date + lag_days.")
}

coverage_check = selected_matrix_data %>%
  group_by(
    tensor_k,
    period_name,
    period_index,
    country_slot,
    country
  ) %>%
  summarize(
    n_policy_days = n_distinct(policy_date),
    n_outcome_days = n_distinct(outcome_date),
    policy_date_min = min(policy_date),
    policy_date_max = max(policy_date),
    outcome_date_min = min(outcome_date),
    outcome_date_max = max(outcome_date),
    n_missing_stringency = sum(is.na(stringency_index)),
    n_missing_deaths = sum(is.na(deaths_raw)),
    .groups = "drop"
  )

bad_coverage = coverage_check %>%
  filter(
    n_policy_days != period_length_days |
      n_outcome_days != period_length_days
  )

if (nrow(bad_coverage) > 0L) {
  cat("\nBad coverage rows:\n")
  print(bad_coverage, n = Inf)
  stop("At least one selected period-country does not have the required number of policy/outcome days.")
}

bad_lagged_deaths = coverage_check %>%
  filter(n_missing_deaths > 0)

if (nrow(bad_lagged_deaths) > 0L) {
  cat("\nRows with missing lagged smoothed deaths per million:\n")
  print(bad_lagged_deaths, n = Inf)
  stop("At least one selected period-country has missing lagged smoothed death per million values.")
}

cat("\n============================================================\n")
cat("Lag alignment and coverage checks passed\n")
cat("============================================================\n")
cat("Every outcome_date equals policy_date +", lag_days, "days.\n")
print(coverage_check, n = Inf)

cat("\nCoverage row count by selected tensor slice:\n")
print(
  coverage_check %>%
    count(
      tensor_k,
      period_name,
      period_index,
      name = "n_rows"
    )
)

make_fixed_owid_tensor = function(
    selected_data,
    period_dates,
    use_log1p_deaths = FALSE
) {
  
  n_country = K_countries_per_period
  n_day = period_length_days
  n_period = nrow(period_dates)
  
  country_dim_names = paste0("country_", seq_len(n_country))
  day_dim_names = paste0("day_", seq_len(n_day))
  period_dim_names = period_dates$period_name
  
  deaths_raw_array = array(
    NA_real_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  deaths_log1p_array = array(
    NA_real_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  deaths_array = array(
    NA_real_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  policy_array = array(
    NA_integer_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  stringency_array = array(
    NA_real_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  policy_date_array = array(
    NA_character_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  outcome_date_array = array(
    NA_character_,
    dim = c(n_country, n_day, n_period),
    dimnames = list(
      country_top_to_bottom = country_dim_names,
      relative_day = day_dim_names,
      period = period_dim_names
    )
  )
  
  country_order = vector("list", n_period)
  date_order = vector("list", n_period)
  
  names(country_order) = period_dim_names
  names(date_order) = period_dim_names
  
  for (k in seq_len(n_period)) {
    
    tensor_k_i = period_dates$tensor_k[k]
    period_name_i = period_dates$period_name[k]
    
    for (rr in seq_len(n_country)) {
      
      d_pr = selected_data %>%
        filter(
          tensor_k == tensor_k_i,
          country_slot == rr
        ) %>%
        arrange(relative_day)
      
      if (nrow(d_pr) != n_day) {
        stop(
          "Tensor construction failed for tensor k ",
          tensor_k_i,
          " (", period_name_i, "), row ",
          rr,
          ": expected ",
          n_day,
          " days but found ",
          nrow(d_pr),
          "."
        )
      }
      
      if (!all(d_pr$outcome_date == d_pr$policy_date + days(lag_days))) {
        stop(
          "Lag alignment failed inside tensor construction for tensor k ",
          tensor_k_i,
          ", row ",
          rr,
          "."
        )
      }
      
      deaths_raw_array[rr, , k] = d_pr$deaths_raw
      deaths_log1p_array[rr, , k] = d_pr$deaths_log1p
      
      deaths_array[rr, , k] = if (isTRUE(use_log1p_deaths)) {
        d_pr$deaths_log1p
      } else {
        d_pr$deaths_raw
      }
      
      policy_array[rr, , k] = d_pr$policy_on
      stringency_array[rr, , k] = d_pr$stringency_index
      
      policy_date_array[rr, , k] = as.character(d_pr$policy_date)
      outcome_date_array[rr, , k] = as.character(d_pr$outcome_date)
    }
    
    country_order[[k]] = selected_data %>%
      filter(tensor_k == tensor_k_i) %>%
      distinct(
        tensor_k,
        period_name,
        period_label,
        period_index,
        country_slot,
        tensor_row,
        row,
        country,
        row_state_name,
        window_start,
        window_end,
        outcome_window_start,
        outcome_window_end
      ) %>%
      arrange(country_slot)
    
    date_order[[k]] = selected_data %>%
      filter(tensor_k == tensor_k_i) %>%
      distinct(
        tensor_k,
        period_name,
        period_index,
        relative_day,
        policy_date,
        outcome_date
      ) %>%
      arrange(relative_day)
  }
  
  if (isTRUE(use_log1p_deaths)) {
    outcome_name = paste0(
      "log1p_new_deaths_smoothed_per_million_lag",
      lag_days,
      "d"
    )
    outcome_transform = paste0(
      "log1p(pmax(new_deaths_smoothed_per_million at policy date + ",
      lag_days,
      " days, 0))"
    )
  } else {
    outcome_name = paste0(
      "new_deaths_smoothed_per_million_lag",
      lag_days,
      "d"
    )
    outcome_transform = "none"
  }
  
  list(
    policy_group = policy_code,
    policy_name = policy_name,
    policy_col = policy_col,
    policy_threshold = policy_threshold,
    lag_days = lag_days,
    timing = paste0(
      "policy_on uses stringency_index at policy date t; deaths use new_deaths_smoothed_per_million at t + ",
      lag_days,
      " days."
    ),
    outcome = outcome_name,
    outcome_raw = paste0(
      "new_deaths_smoothed_per_million_lag",
      lag_days,
      "d"
    ),
    outcome_transform = outcome_transform,
    use_log1p_deaths = use_log1p_deaths,
    structure = "country_top_to_bottom x relative_day x selected_period",
    row_meaning = "Row 1 is row_01 and is the top country in the plot; row 15 is row_15 and is the bottom country.",
    
    period_map = period_dates %>%
      select(
        tensor_k,
        period_name,
        period_label,
        window_start,
        window_end,
        outcome_window_start,
        outcome_window_end
      ),
    
    deaths_raw = deaths_raw_array,
    deaths_log1p = deaths_log1p_array,
    deaths = deaths_array,
    policy_on = policy_array,
    stringency_index = stringency_array,
    policy_date = policy_date_array,
    outcome_date = outcome_date_array,
    
    country_order = country_order,
    date_order = date_order,
    period_dates = period_dates,
    selected_matrix_data = selected_data,
    
    global_death_min_raw = floor(min(deaths_raw_array, na.rm = TRUE)),
    global_death_max_raw = ceiling(max(deaths_raw_array, na.rm = TRUE)),
    global_death_min = floor(min(deaths_array, na.rm = TRUE)),
    global_death_max = ceiling(max(deaths_array, na.rm = TRUE)),
    
    fill_values = c(
      "Off" = "blue",
      "On" = "red",
      "Missing" = "green"
    )
  )
}

A = make_fixed_owid_tensor(
  selected_data = selected_matrix_data,
  period_dates = period_dates,
  use_log1p_deaths = use_log1p_deaths
)

cat("\n============================================================\n")
cat("Tensor dimensions\n")
cat("============================================================\n")
cat("A$deaths dimensions: country_top_to_bottom x relative_day x selected_period\n")
print(dim(A$deaths))

cat("\nA$policy_on dimensions: country_top_to_bottom x relative_day x selected_period\n")
print(dim(A$policy_on))

cat("\nTiming convention:\n")
cat(A$timing, "\n")

cat("\nOutcome used by estimator:\n")
cat(A$outcome, "\n")
cat("Transform:", A$outcome_transform, "\n")

cat("\nPeriod map:\n")
print(A$period_map)

cat("\n============================================================\n")
cat("Country order by selected tensor slice\n")
cat("============================================================\n")

for (pname in names(A$country_order)) {
  cat("\n", pname, "\n", sep = "")
  print(
    A$country_order[[pname]] %>%
      select(
        tensor_k,
        period_name,
        period_index,
        country_slot,
        row,
        country,
        period_label,
        window_start,
        window_end,
        outcome_window_start,
        outcome_window_end
      ),
    n = Inf
  )
}

cat("\n============================================================\n")
cat("Lagged smoothed deaths per million summary by selected tensor slice and country\n")
cat("============================================================\n")

selected_matrix_data %>%
  group_by(
    tensor_k,
    period_name,
    period_index,
    country_slot,
    country
  ) %>%
  summarize(
    n_days = n(),
    policy_date_min = min(policy_date),
    policy_date_max = max(policy_date),
    outcome_date_min = min(outcome_date),
    outcome_date_max = max(outcome_date),
    n_missing_deaths = sum(is.na(deaths_raw)),
    n_positive_deaths = sum(!is.na(deaths_raw) & deaths_raw > 0),
    share_positive_deaths = round(n_positive_deaths / n_days, 2),
    mean_deaths_per_million = round(mean(deaths_raw, na.rm = TRUE), 4),
    median_positive_deaths_per_million = {
      x = deaths_raw[!is.na(deaths_raw) & deaths_raw > 0]
      if (length(x) == 0L) NA_real_ else round(median(x), 4)
    },
    p95_deaths_per_million = round(
      as.numeric(quantile(deaths_raw, 0.95, na.rm = TRUE, names = FALSE)),
      4
    ),
    max_deaths_per_million = round(max(deaths_raw, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(tensor_k, country_slot) %>%
  print(n = Inf)

plot_policy_matrix_from_A = function(A, tensor_k_i) {
  
  period_name_i = dimnames(A$policy_on)$period[tensor_k_i]
  
  country_labs = A$country_order[[period_name_i]] %>%
    arrange(country_slot)
  
  date_labs = A$date_order[[period_name_i]] %>%
    arrange(relative_day)
  
  d_plot = expand_grid(
    country_slot = seq_len(dim(A$policy_on)[1]),
    relative_day = seq_len(dim(A$policy_on)[2])
  ) %>%
    mutate(
      country = country_labs$country[country_slot],
      policy_date = date_labs$policy_date[relative_day],
      outcome_date = date_labs$outcome_date[relative_day],
      policy_on = as.integer(
        A$policy_on[cbind(country_slot, relative_day, tensor_k_i)]
      ),
      stringency_index = as.numeric(
        A$stringency_index[cbind(country_slot, relative_day, tensor_k_i)]
      ),
      policy_status = case_when(
        is.na(stringency_index) ~ "Missing",
        policy_on == 1L ~ "On",
        TRUE ~ "Off"
      )
    )
  
  ggplot(
    d_plot,
    aes(
      x = relative_day,
      y = country_slot,
      fill = policy_status
    )
  ) +
    geom_tile(width = 1, height = 0.95) +
    scale_fill_manual(
      values = c(
        "Off" = "blue",
        "On" = "red",
        "Missing" = "green"
      ),
      drop = FALSE
    ) +
    scale_y_reverse(
      breaks = country_labs$country_slot,
      labels = country_labs$country
    ) +
    scale_x_continuous(
      breaks = pretty(d_plot$relative_day),
      limits = c(1, dim(A$policy_on)[2]),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Policy status: ", period_name_i),
      subtitle = paste0(
        "Stringency >= ",
        A$policy_threshold,
        "; policy dates ",
        format(min(d_plot$policy_date)),
        " to ",
        format(max(d_plot$policy_date))
      ),
      x = "Relative day",
      y = NULL,
      fill = "Policy status"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

safe_quantile = function(x, prob) {
  x = x[!is.na(x)]
  if (length(x) == 0L) {
    NA_real_
  } else {
    as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
  }
}

plot_deaths_matrix_from_A = function(A, tensor_k_i) {
  
  period_name_i = dimnames(A$deaths)$period[tensor_k_i]
  
  country_labs = A$country_order[[period_name_i]] %>%
    arrange(country_slot)
  
  date_labs = A$date_order[[period_name_i]] %>%
    arrange(relative_day)
  
  d_plot = expand_grid(
    country_slot = seq_len(dim(A$deaths)[1]),
    relative_day = seq_len(dim(A$deaths)[2])
  ) %>%
    mutate(
      country = country_labs$country[country_slot],
      policy_date = date_labs$policy_date[relative_day],
      outcome_date = date_labs$outcome_date[relative_day],
      deaths_value = as.numeric(
        A$deaths[cbind(country_slot, relative_day, tensor_k_i)]
      ),
      policy_on = as.integer(
        A$policy_on[cbind(country_slot, relative_day, tensor_k_i)]
      ),
      policy_status = if_else(policy_on == 1L, "On", "Off")
    )
  
  scale_max = safe_quantile(d_plot$deaths_value, 0.95)
  
  if (is.na(scale_max) || is.infinite(scale_max) || scale_max <= 0) {
    scale_max = max(d_plot$deaths_value, na.rm = TRUE)
  }
  
  if (is.na(scale_max) || is.infinite(scale_max) || scale_max <= 0) {
    scale_max = 1
  }
  
  d_plot = d_plot %>%
    mutate(
      Y_scaled = case_when(
        is.na(deaths_value) ~ NA_real_,
        TRUE ~ pmin(pmax(deaths_value / scale_max, 0), 1)
      ),
      fill_color = "green"
    )
  
  valid_off = !is.na(d_plot$Y_scaled) &
    d_plot$policy_status == "Off"
  
  d_plot$fill_color[valid_off] = rgb(
    red = 1 - d_plot$Y_scaled[valid_off],
    green = 1 - d_plot$Y_scaled[valid_off],
    blue = 1,
    maxColorValue = 1
  )
  
  valid_on = !is.na(d_plot$Y_scaled) &
    d_plot$policy_status == "On"
  
  d_plot$fill_color[valid_on] = rgb(
    red = 1,
    green = 1 - d_plot$Y_scaled[valid_on],
    blue = 1 - d_plot$Y_scaled[valid_on],
    maxColorValue = 1
  )
  
  ggplot(
    d_plot,
    aes(
      x = relative_day,
      y = country_slot,
      fill = fill_color
    )
  ) +
    geom_tile(width = 1, height = 0.95) +
    scale_fill_identity() +
    scale_y_reverse(
      breaks = country_labs$country_slot,
      labels = country_labs$country
    ) +
    scale_x_continuous(
      breaks = pretty(d_plot$relative_day),
      limits = c(1, dim(A$deaths)[2]),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0(
        "Smoothed deaths per million: ",
        format(min(d_plot$outcome_date)),
        " to ",
        format(max(d_plot$outcome_date))
      ),
      subtitle = paste0(
        "Value observed ",
        A$lag_days,
        " days after policy date."
      ),
      x = "Relative day",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

if (isTRUE(make_plots)) {
  for (k in seq_len(dim(A$policy_on)[3])) {
    print(plot_policy_matrix_from_A(A, k))
  }
  
  for (k in seq_len(dim(A$deaths)[3])) {
    print(plot_deaths_matrix_from_A(A, k))
  }
}

cat("\nDone constructing A.\n")
cat("Objects available in memory:\n")
cat("  A\n")
cat("  selected_matrix_data\n")
cat("  selected_period_country_map\n")
cat("  period_dates\n")

dir.create("data_clean", showWarnings = FALSE, recursive = TRUE)

period_tag = paste0("periods_", paste(period_dates$period_index, collapse = "_"))
lag_tag = paste0("lag", lag_days, "d")
transform_tag = if_else(use_log1p_deaths, "log1p", "raw")

saveRDS(
  A,
  paste0(
    "data_clean/OWID_stringency_tensor_15rows_",
    period_tag,
    "_",
    lag_tag,
    "_",
    transform_tag,
    ".rds"
  )
)

cat("\n============================================================\n")
cat("Single tensor object checks\n")
cat("============================================================\n")

cat("\nA$deaths dimensions:\n")
print(dim(A$deaths))

cat("\nA$policy_on dimensions:\n")
print(dim(A$policy_on))

cat("\nA$period_map:\n")
print(A$period_map)

cat("\nTiming convention:\n")
cat(A$timing, "\n")

cat("\nOutcome used by estimator:\n")
cat(A$outcome, "\n")

cat("\nCountry order, first tensor slice:\n")
print(A$country_order[[1]])

cat("\nDate order, first tensor slice:\n")
print(A$date_order[[1]])

cat("\nRange of A$deaths:\n")
print(range(A$deaths, na.rm = TRUE))

cat("\nRange of A$deaths_raw:\n")
print(range(A$deaths_raw, na.rm = TRUE))

if ("deaths_log1p" %in% names(A)) {
  cat("\nRange of A$deaths_log1p:\n")
  print(range(A$deaths_log1p, na.rm = TRUE))
}

if (!file.exists(source_bilinear_tensor)) {
  stop("Cannot find source file: ", source_bilinear_tensor)
}

if (!file.exists(source_bilinear_psi)) {
  stop("Cannot find source file: ", source_bilinear_psi)
}

if (!file.exists(source_plugin_psi)) {
  stop("Cannot find source file: ", source_plugin_psi)
}

source(source_bilinear_tensor)
source(source_bilinear_psi)
source(source_plugin_psi)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

k_vec = seq_len(dim(A$deaths)[3])

policy_label = if ("policy_group" %in% names(A)) {
  A$policy_group
} else {
  "A"
}

outcome_label = if ("outcome" %in% names(A)) {
  A$outcome
} else {
  "deaths"
}

make_staggered_parts_from_D = function(D, layer = NULL) {
  
  Dmat = if (length(dim(D)) == 3) {
    if (is.null(layer)) stop("If D is 3D, provide layer.")
    D[, , layer]
  } else {
    D
  }
  
  Dmat = as.matrix(Dmat)
  
  days = seq_len(ncol(Dmat))
  
  first_on = apply(Dmat, 1, function(drow) {
    on_days = days[which(drow == 1)]
    if (length(on_days) == 0) Inf else min(on_days)
  })
  
  adopt_days = sort(unique(first_on[is.finite(first_on)]))
  m = length(adopt_days)
  
  if (m == 0) {
    stop("No treated countries found in this policy slice.")
  }
  
  row_groups = vector("list", m + 1)
  row_groups[[1]] = which(!is.finite(first_on))
  
  for (ell in seq_len(m)) {
    day = rev(adopt_days)[ell]
    row_groups[[ell + 1]] = which(first_on == day)
  }
  
  col_groups = vector("list", m + 1)
  col_groups[[1]] = which(days < adopt_days[1])
  
  if (m >= 2) {
    for (ell in 2:m) {
      col_groups[[ell]] = which(
        days >= adopt_days[ell - 1] &
          days < adopt_days[ell]
      )
    }
  }
  
  col_groups[[m + 1]] = which(days >= adopt_days[m])
  
  N_sizes = lengths(row_groups)
  T_sizes = lengths(col_groups)
  
  if (any(N_sizes == 0)) {
    stop("At least one row block is empty: ", paste(N_sizes, collapse = ", "))
  }
  
  if (any(T_sizes == 0)) {
    stop("At least one column block is empty: ", paste(T_sizes, collapse = ", "))
  }
  
  o = m + 1
  
  row_block_id = integer(nrow(Dmat))
  col_block_id = integer(ncol(Dmat))
  
  for (a in seq_len(o)) {
    row_block_id[row_groups[[a]]] = a
  }
  
  for (b in seq_len(o)) {
    col_block_id[col_groups[[b]]] = b
  }
  
  D_pred = outer(
    row_block_id,
    col_block_id,
    FUN = function(a, b) as.numeric(a + b > o + 1)
  )
  
  if (!all(Dmat == D_pred, na.rm = TRUE)) {
    bad = which(Dmat != D_pred, arr.ind = TRUE)
    
    stop(
      "D does not exactly match inferred staggered block structure. ",
      "First mismatch: row = ", bad[1, 1],
      ", day = ", bad[1, 2],
      ", D = ", Dmat[bad[1, 1], bad[1, 2]],
      ", predicted = ", D_pred[bad[1, 1], bad[1, 2]]
    )
  }
  
  list(
    N_sizes = N_sizes,
    T_sizes = T_sizes,
    row_groups = row_groups,
    col_groups = col_groups,
    adopt_days = adopt_days,
    first_on = first_on
  )
}

make_parts_for_policy_tensor = function(A) {
  
  K = dim(A$deaths)[3]
  
  parts = map(seq_len(K), function(k) {
    make_staggered_parts_from_D(A$policy_on[, , k])
  })
  
  list(
    N_parts = map(parts, "N_sizes"),
    T_parts = map(parts, "T_sizes"),
    parts = parts
  )
}

estimate_one_policy = function(
    A,
    policy_label,
    k,
    functional,
    r,
    tau,
    N_parts,
    T_parts
) {
  
  Y = A$deaths
  D = A$policy_on
  
  Y0 = Y
  Y0[D == 1] = NA_real_
  
  period_name = dimnames(Y)$period[k]
  
  period_row = A$period_map %>%
    filter(tensor_k == k)
  
  period_label_i = period_row$period_label[1]
  
  psi0 = bilinearTensorStaggeredPsi(
    Y = Y0,
    k = k,
    r = r,
    N_parts = N_parts,
    T_parts = T_parts,
    tau = tau,
    functional = functional,
    eta = NULL,
    row_index = NULL
  )
  
  psi1 = pluginPsi_c1(
    Y = Y,
    k = k,
    N_parts = N_parts,
    T_parts = T_parts,
    functional = functional,
    eta = NULL,
    row_index = NULL
  )
  
  tibble(
    policy = policy_label,
    outcome = outcome_label,
    lag_days = A$lag_days,
    tensor_k = k,
    period = period_name,
    period_label = period_label_i,
    functional = functional,
    Psi0 = as.numeric(psi0),
    Psi1 = as.numeric(psi1),
    Delta_h = as.numeric(psi1 - psi0)
  )
}

make_multislice_policy_table = function(
    A,
    policy_label,
    k_vec,
    functionals = c("ATE", "Trend"),
    r = 1,
    tau = 0.01
) {
  
  parts_obj = make_parts_for_policy_tensor(A)
  N_parts = parts_obj$N_parts
  T_parts = parts_obj$T_parts
  
  out = map_dfr(k_vec, function(k) {
    map_dfr(functionals, function(functional) {
      estimate_one_policy(
        A = A,
        policy_label = policy_label,
        k = k,
        functional = functional,
        r = r,
        tau = tau,
        N_parts = N_parts,
        T_parts = T_parts
      )
    })
  }) %>%
    mutate(
      across(
        c(Psi0, Psi1, Delta_h),
        ~ round(.x, 6)
      )
    ) %>%
    arrange(tensor_k, functional) %>%
    select(
      policy,
      outcome,
      lag_days,
      tensor_k,
      period,
      period_label,
      functional,
      Psi0,
      Psi1,
      Delta_h
    )
  
  cat("\n============================================================\n")
  cat(policy_label, ": tensor k =", paste(k_vec, collapse = ", "), "\n")
  cat("Outcome =", outcome_label, "\n")
  cat("Lag days =", A$lag_days, "\n")
  cat("Timing =", A$timing, "\n")
  cat("No additional outcome transformation is applied here.\n")
  cat("Rank r =", r, "\n")
  cat("Tau =", tau, "\n")
  cat("Functionals:", paste(functionals, collapse = ", "), "\n")
  cat("============================================================\n")
  
  print(out, n = Inf)
  
  return(out)
}

A_table = make_multislice_policy_table(
  A = A,
  policy_label = policy_label,
  k_vec = k_vec,
  functionals = functionals_to_run,
  r = r,
  tau = tau
)

all_policy_tables = list(
  A = A_table
)

saveRDS(
  A_table,
  paste0(
    "data_clean/A_policy_table_15rows_",
    period_tag,
    "_",
    lag_tag,
    "_",
    transform_tag,
    "_ATE_Trend.rds"
  )
)

saveRDS(
  all_policy_tables,
  paste0(
    "data_clean/all_policy_tables_A_only_15rows_",
    period_tag,
    "_",
    lag_tag,
    "_",
    transform_tag,
    "_ATE_Trend.rds"
  )
)

cat("\nDone.\n")
cat("Objects available in memory:\n")
cat("  A\n")
cat("  A_table\n")
cat("  all_policy_tables\n")
cat("  selected_matrix_data\n")
cat("  selected_period_country_map\n")
cat("  period_dates\n")

if (interactive()) {
  View(A_table)
}