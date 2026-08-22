library(ggplot2)

setwd("~/Desktop/code/")

data_file = "Real/CovidOx/data/Omega_Y_until_2020-04-05_delay_28.rds"

obj = readRDS(data_file)

Omega = obj$Omega
Y = obj$Y
outcome_caps = obj$outcome_caps

blue_pal = grDevices::colorRampPalette(c("lightblue", "blue"))(101)
red_pal = grDevices::colorRampPalette(c("pink", "red"))(101)

plot_list = vector("list", dim(Y)[3])

for (k in seq_len(dim(Y)[3])) {
  
  outcome = dimnames(Y)$outcome[k]
  cap = unname(outcome_caps[outcome])
  
  y_mat = Y[, , k]
  omega_mat = Omega[, , k]
  
  value_capped = pmin(y_mat, cap)
  value_scaled = value_capped / cap
  
  intensity = round(value_scaled * 100) + 1
  intensity = pmin(101, pmax(1, intensity))
  
  colours = matrix(
    NA_character_,
    nrow = nrow(y_mat),
    ncol = ncol(y_mat)
  )
  
  colours[omega_mat == 0] = blue_pal[intensity[omega_mat == 0]]
  colours[omega_mat == 1] = red_pal[intensity[omega_mat == 1]]
  
  df = expand.grid(
    country = dimnames(Y)$country,
    date = dimnames(Y)$date
  )
  
  df$value = as.vector(y_mat)
  df$colour = as.vector(colours)
  df$outcome = outcome
  
  plot_list[[k]] = df
}

plot_df = do.call(rbind, plot_list)

plot_df$date = as.Date(plot_df$date)

plot_df$country = factor(
  plot_df$country,
  levels = dimnames(Y)$country
)

plot_df$outcome = factor(
  plot_df$outcome,
  levels = dimnames(Y)$outcome
)

plot_df = plot_df[
  !is.na(plot_df$value) &
    !is.na(plot_df$colour),
]

p = ggplot(
  plot_df,
  aes(
    x = date,
    y = country,
    fill = colour
  )
) +
  geom_tile() +
  facet_wrap(
    ~ outcome,
    nrow = 1,
    scales = "free"
  ) +
  scale_fill_identity() +
  labs(
    x = "Policy date",
    y = "Country"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

print(p)

ggsave(
  "Plots/Y_two_layers.png",
  p,
  width = 18,
  height = 8,
  dpi = 600
)
