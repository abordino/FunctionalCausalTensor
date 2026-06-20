setwd("~/Documents/phd/projects/causalMatrix/code/real-world/CovidOx")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(forcats)
library(patchwork)

# -----------------------------
# User input
# -----------------------------

data_file = "data/oxcgrt_owid_merged_daily_2020_2022.csv"

policy_vars = c(
  "C6_stay_at_home_requirements",
  "H3_contact_tracing"
)

# -----------------------------
# Load data
# -----------------------------

df = read_csv(data_file, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

# -----------------------------
# Country labels: show one every 5 countries
# -----------------------------

country_levels = df %>%
  distinct(country) %>%
  arrange(country) %>%
  pull(country)

country_labels = country_levels

country_labels[seq_along(country_labels) %% 5 != 1] = ""

names(country_labels) = country_levels

# -----------------------------
# Prepare data for plotting
# -----------------------------

plot_data = df %>%
  select(
    iso_code,
    country,
    date,
    all_of(policy_vars)
  ) %>%
  pivot_longer(
    cols = all_of(policy_vars),
    names_to = "policy",
    values_to = "policy_value"
  ) %>%
  mutate(
    policy_on = case_when(
      is.na(policy_value) ~ NA_real_,
      policy_value == 0 ~ 0,
      policy_value > 0 ~ 1
    ),
    country = factor(country, levels = country_levels),
    policy = factor(policy, levels = policy_vars)
  )

# -----------------------------
# Plot function
# -----------------------------

plot_policy_matrix = function(policy_name) {
  
  p_data = plot_data %>%
    filter(policy == policy_name)
  
  p = ggplot(
    p_data,
    aes(
      x = date,
      y = country,
      fill = factor(policy_on)
    )
  ) +
    geom_tile() +
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
      na.value = "grey85",
      breaks = c("0", "1"),
      labels = c("OFF", "ON"),
      name = "Policy"
    ) +
    labs(
      title = str_replace_all(policy_name, "_", " "),
      x = "Date",
      y = "Country"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 4),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

# -----------------------------
# Create plots
# -----------------------------

plots = lapply(policy_vars, plot_policy_matrix)
names(plots) = policy_vars

# -----------------------------
# Combine plots in a 2 x 2 grid
# -----------------------------

combined_plot =
  (
      plots[["C6_stay_at_home_requirements"]]
  ) /
  (
    plots[["H3_contact_tracing"]]
  ) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

print(combined_plot)

ggsave(
  filename = "figure/policy_matrices_2x2.pdf",
  plot = combined_plot,
  width = 14,
  height = 10
)

ggsave(
  filename = "figure/policy_matrices_2x2.png",
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 300
)