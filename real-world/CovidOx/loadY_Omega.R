setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CovidOx")

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)

# -----------------------------
# User input
# -----------------------------

tensor_file = "data/Omega_Y_until_2020-04-05_delay_28.rds"

# -----------------------------
# Load saved tensors
# -----------------------------

tensor_data = readRDS(tensor_file)

Omega = tensor_data$Omega
Y = tensor_data$Y

country_levels = tensor_data$country_levels
date_names = tensor_data$date_names
policies = tensor_data$policies
outcomes = tensor_data$outcomes
outcome_caps = tensor_data$outcome_caps
policy_outcome_map = tensor_data$policy_outcome_map

start_date = tensor_data$start_date
start.date = tensor_data$start.date
time_horizon = tensor_data$time_horizon
delay_days = tensor_data$delay_days

all_dates = as.Date(date_names)

cat("Loaded Omega dimensions:", dim(Omega), "\n")
cat("Loaded Y dimensions:", dim(Y), "\n")
cat("Countries:", length(country_levels), "\n")
cat("Dates:", length(date_names), "\n")
cat("Policies:", paste(policies, collapse = ", "), "\n")
cat("Outcomes:", paste(outcomes, collapse = ", "), "\n")
cat("Outcome delay in days:", delay_days, "\n")
cat("NA in Omega:", sum(is.na(Omega)), "\n")
cat("NA in Y:", sum(is.na(Y)), "\n")

# -----------------------------
# Country labels
# -----------------------------

country_labels = country_levels
names(country_labels) = country_levels

# -----------------------------
# Omega long dataframe
# -----------------------------

omega_plot_df = as.data.frame.table(
  Omega,
  responseName = "policy_on",
  stringsAsFactors = FALSE
) %>%
  rename(
    country = country,
    date = date,
    policy = policy
  ) %>%
  mutate(
    date = as.Date(date),
    country = factor(country, levels = country_levels),
    policy = factor(policy, levels = policies)
  ) %>%
  drop_na(policy_on)

# -----------------------------
# Plot Omega tensor
#
# Blue = policy OFF
# Red  = policy ON
# Missing values are not plotted.
# -----------------------------

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
    title = paste0(
      "Omega: policy adoption paths between ",
      start.date,
      " to ",
      time_horizon
    ),
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
# Y long dataframe
# -----------------------------

outcome_long = as.data.frame.table(
  Y,
  responseName = "value",
  stringsAsFactors = FALSE
) %>%
  rename(
    country = country,
    date = date,
    outcome = outcome
  ) %>%
  mutate(date = as.Date(date))

policy_long = as.data.frame.table(
  Omega,
  responseName = "policy_on",
  stringsAsFactors = FALSE
) %>%
  rename(
    country = country,
    date = date,
    policy = policy
  ) %>%
  mutate(date = as.Date(date))

# -----------------------------
# Plot Y tensor
# -----------------------------

blue_pal = grDevices::colorRampPalette(c("lightblue", "blue"))(101)
red_pal  = grDevices::colorRampPalette(c("pink", "red"))(101)

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
      "\nusing outcomes at t + ",
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

# -----------------------------
# Save plots
# -----------------------------

dir.create("figure", showWarnings = FALSE)

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
