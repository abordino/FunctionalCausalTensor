setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CovidOx")

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)

# -----------------------------
# Settings
# -----------------------------

out_dir = "data"
dir.create(out_dir, showWarnings = FALSE)

years = c(2020, 2021, 2022)

start_date = as.Date("2020-01-01")
end_date   = as.Date("2022-12-31")

oxcgrt_urls = paste0(
  "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/",
  "OxCGRT_nat_differentiated_withnotes_", years, ".csv"
)

owid_url = paste0(
  "https://raw.githubusercontent.com/owid/covid-19-data/master/",
  "public/data/owid-covid-data.csv"
)

output_file = file.path(
  out_dir,
  "oxcgrt_owid_merged_daily_2020_2022.csv"
)

# -----------------------------
# Helpers
# -----------------------------

coalesce_existing = function(data, columns) {
  columns = columns[columns %in% names(data)]
  
  if (length(columns) == 0) {
    return(rep(NA_real_, nrow(data)))
  }
  
  dplyr::coalesce(!!!lapply(columns, function(x) data[[x]]))
}

# Proceed backwards row by row.
# Each zero is filled with the latest original nonzero value
# encountered later in time.
fill_zeros_backward_from_original_nonzero = function(x) {
  x_filled = x
  latest_original_nonzero = NA_real_
  
  for (i in seq(length(x), 1)) {
    
    if (is.na(x[i])) {
      next
    }
    
    if (x[i] != 0) {
      latest_original_nonzero = x[i]
      x_filled[i] = x[i]
    } else if (!is.na(latest_original_nonzero)) {
      x_filled[i] = latest_original_nonzero
    }
  }
  
  x_filled
}

# -----------------------------
# Download and prepare OxCGRT
# -----------------------------

oxcgrt_raw = map2_dfr(
  oxcgrt_urls,
  years,
  function(url, year) {
    read_csv(url, show_col_types = FALSE, guess_max = 100000) %>%
      mutate(source_year_file = year)
  }
)

oxcgrt_raw$StringencyIndex = coalesce_existing(
  oxcgrt_raw,
  c(
    "StringencyIndex_SimpleAverage",
    "StringencyIndex_WeightedAverage",
    "StringencyIndex_NonVaccinated",
    "StringencyIndex"
  )
)

oxcgrt_raw$GovernmentResponseIndex = coalesce_existing(
  oxcgrt_raw,
  c(
    "GovernmentResponseIndex_SimpleAverage",
    "GovernmentResponseIndex_WeightedAverage",
    "GovernmentResponseIndex_NonVaccinated",
    "GovernmentResponseIndex"
  )
)

oxcgrt_policy = oxcgrt_raw %>%
  mutate(
    date = ymd(as.character(Date)),
    iso_code = CountryCode
  ) %>%
  filter(date >= start_date, date <= end_date) %>%
  filter(is.na(Jurisdiction) | Jurisdiction == "NAT_TOTAL") %>%
  transmute(
    iso_code,
    date,
    
    C1_school_closing =
      `C1E_School closing`,
    
    C2_workplace_closing =
      `C2E_Workplace closing`,
    
    C6_stay_at_home_requirements =
      `C6E_Stay at home requirements`,
    
    C8_international_travel_controls =
      `C8E_International travel controls`,
    
    H2_testing_policy =
      `H2E_Testing policy`,
    
    H3_contact_tracing =
      `H3E_Contact tracing`,
    
    H6_facial_coverings =
      `H6E_Facial Coverings`,
    
    H7_vaccination_policy =
      `H7E_Vaccination policy`,
    
    StringencyIndex,
    GovernmentResponseIndex
  ) %>%
  distinct(iso_code, date, .keep_all = TRUE)

# -----------------------------
# Download and prepare OWID
# -----------------------------

owid_outcomes = read_csv(
  owid_url,
  show_col_types = FALSE,
  guess_max = 100000
) %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= start_date, date <= end_date) %>%
  transmute(
    iso_code,
    country = location,
    date,
    
    new_deaths_per_million,
    new_cases_per_million,
    new_deaths_smoothed_per_million,
    new_cases_smoothed_per_million,
    people_fully_vaccinated_per_hundred,
    icu_patients_per_million,
    hosp_patients_per_million,
    positive_rate,
    reproduction_rate,    
    new_tests_smoothed_per_thousand,
    weekly_hosp_admissions_per_million
  ) %>%
  distinct(iso_code, date, .keep_all = TRUE)

# -----------------------------
# Merge
# -----------------------------

merged_data = owid_outcomes %>%
  inner_join(
    oxcgrt_policy,
    by = c("iso_code", "date")
  ) %>%
  arrange(iso_code, date) %>%
  group_by(iso_code) %>%
  mutate(
    death_growth_rate =
      log1p(new_deaths_per_million) -
      lag(log1p(new_deaths_per_million)),
    
    case_growth_rate =
      log1p(new_cases_per_million) -
      lag(log1p(new_cases_per_million)),
    
    death_acceleration =
      death_growth_rate - lag(death_growth_rate),
    
    case_acceleration =
      case_growth_rate - lag(case_growth_rate),
    
    death_growth_rate = replace_na(death_growth_rate, 0),
    case_growth_rate = replace_na(case_growth_rate, 0),
    death_acceleration = replace_na(death_acceleration, 0),
    case_acceleration = replace_na(case_acceleration, 0),
    
    death_acceleration =
      fill_zeros_backward_from_original_nonzero(
        death_acceleration
      ),
    
    case_acceleration =
      fill_zeros_backward_from_original_nonzero(
        case_acceleration
      )
  ) %>%
  ungroup() %>%
  mutate(
    week = floor_date(date, "week", week_start = 1),
    month = floor_date(date, "month"),
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5) ~ "Spring",
      month(date) %in% c(6, 7, 8) ~ "Summer",
      month(date) %in% c(9, 10, 11) ~ "Fall"
    ),
    season_year = if_else(month(date) == 12, year(date) + 1L, year(date)),
    season_id = paste0(season_year, "_", season)
  ) %>%
  select(
    iso_code,
    country,
    date,
    week,
    month,
    season_id,
    
    new_deaths_per_million,
    new_cases_per_million,
    new_deaths_smoothed_per_million,
    new_cases_smoothed_per_million,
    death_growth_rate,
    case_growth_rate,
    death_acceleration,
    case_acceleration,
    people_fully_vaccinated_per_hundred,
    icu_patients_per_million,
    hosp_patients_per_million,
    positive_rate,
    reproduction_rate,
    new_tests_smoothed_per_thousand,
    weekly_hosp_admissions_per_million,
    
    C1_school_closing,
    C2_workplace_closing,
    C6_stay_at_home_requirements,
    C8_international_travel_controls,
    H2_testing_policy,
    H3_contact_tracing,
    H6_facial_coverings,
    H7_vaccination_policy,
    StringencyIndex,
    GovernmentResponseIndex
  ) %>%
  arrange(iso_code, date)

# -----------------------------
# Save only final merged CSV
# -----------------------------

write_csv(merged_data, output_file)

cat("Saved:", output_file, "\n")
cat("Rows:", nrow(merged_data), "\n")
cat("Countries:", n_distinct(merged_data$iso_code), "\n")
cat("Date range:", as.character(min(merged_data$date)), "to", as.character(max(merged_data$date)), "\n")