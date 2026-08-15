setwd("~/Desktop/code")

suppressPackageStartupMessages(library(tidyverse))

source("bilinearTensorAllFunction.R")
source("CYNuclear.R")

options(digits = 10)

# Fix lambda so the numerical check is deterministic
cy_lambda_matlab = function(omega,
                             sigma_squared,
                             monte_carlo_draws = 100L,
                             probability = 0.95,
                             lambda_multiplier = 1.14) {
  list(
    lambda = 1,
    simulated_norms = rep(NA_real_, monte_carlo_draws),
    noise_quantile = NA_real_,
    sigma_squared = sigma_squared
  )
}

# Full counterfactual matrix M0
M0 = outer(c(1, 2, 3, 4), c(1, 2, 3, 4))
rownames(M0) = paste0('row', 1:4)
colnames(M0) = paste0('t', 1:4)

# Rows 1-2 never adopt; rows 3-4 adopt at t=2
omega = rbind(
  c(1, 1, 1, 1),
  c(1, 1, 1, 1),
  c(1, 0, 0, 0),
  c(1, 0, 0, 0)
)

y_obs = M0
y_obs[omega == 0] = 0

fit = cy_complete_staggered_layer(
  y = y_obs,
  omega = omega,
  adoption_time = c(Inf, Inf, 2, 2),
  time = 1:4,
  rank = 1,
  sigma = 1,
  subgroup_size = 1,
  tau = 0.5,
  tolerance = 1e-12,
  max_iterations = 100000L,
  monte_carlo_draws = 1L,
  processing_order = 'matlab',
  verbose = FALSE
)

N_parts = list(c(2, 2))
T_parts = list(c(1, 3))
eta = c(1, 1, 1, -1)

psi0 = function(Mhat) {
  c(
    ATE = pluginPsi1_ATE(Mhat, 1, N_parts, T_parts),
    Trend = pluginPsi1_Trend(Mhat, 1, N_parts, T_parts),
    RowHet = pluginPsi1_RowHet(Mhat, 1, N_parts, T_parts, eta),
    ThirdRow = pluginPsi1_Local(Mhat, 1, N_parts, T_parts, row_index = 3)
  )
}

cat('\nFull M0:\n')
print(M0)
cat('\nMasked input (NA shown for missing):\n')
masked = M0
masked[omega == 0] = NA_real_
print(masked)
cat('\nCY nuclear/regularized completion:\n')
print(fit$completed_nuclear)
cat('\nCY debiased completion:\n')
print(fit$completed_debiased)
cat('\nPlug-in Psi0 functionals:\n')
print(rbind(
  truth = psi0(M0),
  nuclear = psi0(fit$completed_nuclear),
  debiased = psi0(fit$completed_debiased)
))
