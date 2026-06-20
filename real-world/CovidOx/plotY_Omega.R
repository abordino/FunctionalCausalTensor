setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CovidOx")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)

# -----------------------------
# Set par
# -----------------------------

data_file = "data/oxcgrt_owid_merged_daily_2020_2022.csv"

time_horizon = as.Date("2020-04-05")
start.date = "2020-03-05"
start_date = as.Date(start.date)

delay_days = 28

policy_string = "
C6_stay_at_home_requirements,
H3_contact_tracing
"

policies = str_split(policy_string, "[,;\\s]+")[[1]]
policies = policies[policies != ""]

print(policies)

outcomes = c(
  "new_deaths_smoothed_per_million",
  "new_cases_smoothed_per_million"
)

outcome_caps = c(
  new_deaths_smoothed_per_million = 1,
  new_cases_smoothed_per_million = 20
)

policy_outcome_map = tibble(
  policy = policies,
  outcome = outcomes,
  policy_code = c("C6", "H3")
) %>%
  mutate(
    outcome_cap = unname(outcome_caps[outcome]),
    panel_label = paste0(
      outcome,
      "\nPolicy: ", policy_code,
      " | cap: ", outcome_cap
    )
  )

# -----------------------------
# Country restriction
# -----------------------------

reliable_europe_iso_codes = c(
  "AUT", "BEL", "CHE", "DEU", "DNK", "FIN", "FRA", "GBR",
  "IRL", "ISL", "LIE", "LUX", "NLD", "NOR", "SWE",
  
  "AND", "CYP", "ESP", "GRC", "ITA", "MLT", "PRT", "SMR",
  
  "CZE", "EST", "HRV", "HUN", "LTU", "LVA", "POL", "SVK", "SVN",
  
  "ALB", "BIH", "BGR", "MDA", "MNE", "MKD", "SRB",
  
  "TUR", "UKR"
)

reliable_global_iso_codes = c(
  "AUS",
  "CAN",
  "CHL",
  "CRI",
  "ISR",
  "JPN",
  "KOR",
  "NZL",
  "SGP",
  "TWN",
  "USA"
)

reliable_iso_codes = c(
  reliable_europe_iso_codes,
  reliable_global_iso_codes
)

# -----------------------------
# Load data
# -----------------------------

df_raw = read_csv(data_file, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

missing_policies = setdiff(policies, names(df_raw))
missing_outcomes = setdiff(outcomes, names(df_raw))

if (length(missing_policies) > 0) {
  stop(
    paste0(
      "These policies are not found in the dataset: ",
      paste(missing_policies, collapse = ", ")
    )
  )
}

if (length(missing_outcomes) > 0) {
  stop(
    paste0(
      "These outcomes are not found in the dataset: ",
      paste(missing_outcomes, collapse = ", ")
    )
  )
}

policy_df = df_raw %>%
  filter(date >= start_date, date <= time_horizon) %>%
  select(
    iso_code,
    country,
    date,
    all_of(policies)
  )

outcome_df = df_raw %>%
  filter(
    date >= start_date + days(delay_days),
    date <= time_horizon + days(delay_days)
  ) %>%
  mutate(date = date - days(delay_days)) %>%
  select(
    iso_code,
    country,
    date,
    all_of(outcomes)
  )

df = policy_df %>%
  left_join(
    outcome_df,
    by = c("iso_code", "country", "date")
  )


all_dates = seq.Date(start_date, time_horizon, by = "day")
date_names = as.character(all_dates)

country_grid = df %>%
  distinct(iso_code, country) %>%
  crossing(date = all_dates)

df_complete = country_grid %>%
  left_join(df, by = c("iso_code", "country", "date"))

# -----------------------------
# Convert policies to binary
# OFF = 0 if original value is 0
# ON  = 1 if original value > 0
# NA countries are excluded later
# -----------------------------

binary_df = df_complete %>%
  mutate(
    across(
      all_of(policies),
      ~ case_when(
        is.na(.x) ~ NA_real_,
        .x == 0 ~ 0,
        .x > 0 ~ 1
      )
    )
  )

# -----------------------------
# Country-selection rule
# -----------------------------

valid_policy_path = function(x) {
  if (any(is.na(x))) {
    return(FALSE)
  }
  
  if (length(x) == 0) {
    return(FALSE)
  }
  
  if (x[1] != 0) {
    return(FALSE)
  }
  
  changes = diff(x)
  
  no_reversion = all(changes >= 0)
  
  return(no_reversion)
}

valid_countries = binary_df %>%
  arrange(iso_code, date) %>%
  group_by(iso_code, country) %>%
  summarise(
    across(
      all_of(policies),
      valid_policy_path,
      .names = "valid_{.col}"
    ),
    .groups = "drop"
  ) %>%
  filter(
    if_all(starts_with("valid_"), ~ .x)
  )

# -----------------------------
# Final selected data
# -----------------------------

final_df = binary_df %>%
  semi_join(valid_countries, by = c("iso_code", "country")) %>%
  filter(iso_code %in% reliable_iso_codes) %>%
  arrange(country, date)

selected_countries = final_df %>%
  distinct(iso_code, country) %>%
  arrange(country)

print(selected_countries)

cat("\nNumber of selected countries:", nrow(selected_countries), "\n")

missing_reliable_countries = tibble(iso_code = reliable_iso_codes) %>%
  anti_join(selected_countries, by = "iso_code")

cat("\nReliable whitelist countries not selected, either because they are absent from the dataset or failed the valid-path rule:\n")
print(missing_reliable_countries)

if (nrow(final_df) == 0) {
  stop("No countries remain after the policy-path and reliable-whitelist selection rules.")
}

# -----------------------------
# Tensor dimensions
# -----------------------------

country_levels = final_df %>%
  distinct(country) %>%
  arrange(country) %>%
  pull(country)

# -----------------------------
# Create Omega tensor
# -----------------------------

Omega = array(
  NA_real_,
  dim = c(
    length(country_levels),
    length(date_names),
    length(policies)
  ),
  dimnames = list(
    country = country_levels,
    date = date_names,
    policy = policies
  )
)

for (k in seq_along(policies)) {
  
  policy_k = policies[k]
  
  mat_k = final_df %>%
    select(
      country,
      date,
      value = all_of(policy_k)
    ) %>%
    mutate(date_chr = as.character(date)) %>%
    select(country, date_chr, value) %>%
    pivot_wider(
      names_from = date_chr,
      values_from = value
    ) %>%
    right_join(
      tibble(country = country_levels),
      by = "country"
    ) %>%
    arrange(match(country, country_levels))
  
  missing_dates = setdiff(date_names, names(mat_k))
  
  if (length(missing_dates) > 0) {
    mat_k[missing_dates] = NA_real_
  }
  
  mat_k = mat_k %>%
    select(country, all_of(date_names))
  
  Omega[, , k] = as.matrix(mat_k[, date_names])
}

# -----------------------------
# Create Y tensor
# -----------------------------

Y = array(
  NA_real_,
  dim = c(
    length(country_levels),
    length(date_names),
    length(outcomes)
  ),
  dimnames = list(
    country = country_levels,
    date = date_names,
    outcome = outcomes
  )
)

for (k in seq_along(outcomes)) {
  
  outcome_k = outcomes[k]
  
  mat_k = final_df %>%
    select(
      country,
      date,
      value = all_of(outcome_k)
    ) %>%
    mutate(date_chr = as.character(date)) %>%
    select(country, date_chr, value) %>%
    pivot_wider(
      names_from = date_chr,
      values_from = value
    ) %>%
    right_join(
      tibble(country = country_levels),
      by = "country"
    ) %>%
    arrange(match(country, country_levels))
  
  missing_dates = setdiff(date_names, names(mat_k))
  
  if (length(missing_dates) > 0) {
    mat_k[missing_dates] = NA_real_
  }
  
  mat_k = mat_k %>%
    select(country, all_of(date_names))
  
  Y[, , k] = as.matrix(mat_k[, date_names])
}

cat("Omega dimensions:", dim(Omega), "\n")
cat("Y dimensions:", dim(Y), "\n")
cat("NA in Omega:", sum(is.na(Omega)), "\n")
cat("NA in Y:", sum(is.na(Y)), "\n")
cat("Outcome delay in days:", delay_days, "\n")


country_labels = country_levels
country_labels[seq_along(country_labels) %% 1 == 1] = ""
names(country_labels) = country_levels

# -----------------------------
# Plot Omega tensor
#
# Blue = policy OFF
# Red  = policy ON
# -----------------------------

omega_plot_df = final_df %>%
  select(country, date, all_of(policies)) %>%
  pivot_longer(
    cols = all_of(policies),
    names_to = "policy",
    values_to = "policy_on"
  ) %>%
  drop_na(policy_on) %>%
  mutate(
    country = factor(country, levels = country_levels),
    policy = factor(policy, levels = policies)
  )

p_omega = ggplot(
  omega_plot_df,
  aes(x = date, y = country, fill = factor(policy_on))
) +
  geom_tile() +
  facet_wrap(~ policy, nrow = 1) +
  scale_y_discrete(
    breaks = country_levels,
    labels = country_labels,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "0" = "blue",
      "1" = "red"
    ),
    breaks = c("0", "1"),
    labels = c("OFF", "ON"),
    name = "Policy"
  ) +
  labs(
    title = paste0("Omega: policy adoption paths between ", start.date, " to ", time_horizon),
    x = "Date",
    y = "Country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_omega)

# -----------------------------
# Plot Y
# -----------------------------

blue_pal = grDevices::colorRampPalette(c("lightblue", "blue"))(101)
red_pal  = grDevices::colorRampPalette(c("pink", "red"))(101)

outcome_long = final_df %>%
  select(country, date, all_of(outcomes)) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "outcome",
    values_to = "value"
  )

policy_long = final_df %>%
  select(country, date, all_of(policies)) %>%
  pivot_longer(
    cols = all_of(policies),
    names_to = "policy",
    values_to = "policy_on"
  )

y_plot_df = outcome_long %>%
  left_join(policy_outcome_map, by = "outcome") %>%
  left_join(
    policy_long,
    by = c("country", "date", "policy")
  ) %>%
  drop_na(value, policy_on) %>%
  mutate(
    value_capped = pmin(value, outcome_cap),
    value_scaled = value_capped / outcome_cap,
    intensity_index = as.integer(round(value_scaled * 100)) + 1L,
    intensity_index = pmin(101L, pmax(1L, intensity_index)),
    fill_colour = if_else(
      policy_on == 1,
      red_pal[intensity_index],
      blue_pal[intensity_index]
    ),
    country = factor(country, levels = country_levels),
    panel_label = factor(panel_label, levels = policy_outcome_map$panel_label)
  )

p_y = ggplot(
  y_plot_df,
  aes(x = date, y = country, fill = fill_colour)
) +
  geom_tile() +
  facet_wrap(~ panel_label, nrow = 1) +
  scale_y_discrete(
    breaks = country_levels,
    labels = country_labels,
    drop = FALSE
  ) +
  scale_fill_identity() +
  labs(
    title = paste0(
      "Outcome tensor values between ",
      start.date,
      " to ",
      time_horizon,
      " \nusing outcomes at t + ",
      delay_days,
      " days"
    ),
    x = "Policy date",
    y = "Country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

print(p_y)

ggsave(
  filename = paste0(
    "figure/Omega_policy_tensor_until_",
    time_horizon,
    "_delay_",
    delay_days,
    ".png"
  ),
  plot = p_omega,
  width = 18,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = paste0(
    "figure/Y_outcome_tensor_until_",
    time_horizon,
    "_delay_",
    delay_days,
    ".png"
  ),
  plot = p_y,
  width = 18,
  height = 8,
  units = "in",
  dpi = 600
)

tensor_object = list(
  Omega = Omega,
  Y = Y,
  country_levels = country_levels,
  date_names = date_names,
  policies = policies,
  outcomes = outcomes,
  outcome_caps = outcome_caps,
  start_date = start_date,
  start.date = start.date,
  time_horizon = time_horizon,
  delay_days = delay_days,
  policy_outcome_map = policy_outcome_map,
  reliable_iso_codes = reliable_iso_codes,
  selected_countries = selected_countries,
  missing_reliable_countries = missing_reliable_countries
)

saveRDS(
  tensor_object,
  paste0(
    "data/Omega_Y_until_",
    time_horizon,
    "_delay_",
    delay_days,
    ".rds"
  )
)
