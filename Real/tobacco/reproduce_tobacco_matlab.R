setwd("~/Desktop/code")

source("CYNuclear.R")

# Settings from Main_tobacco.m
data_path = "Real/tobacco/data/Data_cigarette_sales (open data).xlsx"
number_iterations = 10L
seed = NULL

if (!is.null(seed)) set.seed(seed)

# Tobacco data: 38 states x 31 years
raw = readxl::read_excel(
  data_path,
  sheet = "for matlab_new",
  range = "A2:AM33",
  col_names = FALSE,
  .name_repair = "minimal"
)

state_names = as.character(unlist(raw[1, 2:39], use.names = FALSE))
years = as.integer(unlist(raw[2:32, 1], use.names = FALSE))

tobacco = t(matrix(
  as.numeric(unlist(raw[2:32, 2:39], use.names = FALSE)),
  nrow = 31,
  ncol = 38
))
rownames(tobacco) = state_names
colnames(tobacco) = years

# Same staggered missing-data pattern as the MATLAB code
omega = rbind(
  matrix(1, 23, 31),
  cbind(matrix(1, 6, 26), matrix(0, 6, 5)),
  cbind(matrix(1, 6, 21), matrix(0, 6, 10)),
  cbind(matrix(1, 3, 16), matrix(0, 3, 15))
)

adoption_time = c(
  rep(Inf, 23),
  rep(27, 6),
  rep(22, 6),
  rep(17, 3)
)

time = 1:31
missing = omega == 0

RMSE_CY = numeric(number_iterations)
RMSE_CY_nuclear = numeric(number_iterations)

for (i in seq_len(number_iterations)) {
  mild_1 = sample.int(6, 3, replace = FALSE) + 20
  mild_2 = setdiff(21:26, mild_1)
  
  moderate_1 = sample.int(6, 3, replace = FALSE) + 26
  moderate_2 = setdiff(27:32, moderate_1)
  
  severe_1 = sample.int(6, 3, replace = FALSE) + 32
  severe_2 = setdiff(33:38, severe_1)
  
  rows = c(
    1:20,
    mild_1, mild_2,
    moderate_1, moderate_2,
    severe_1, severe_2
  )
  
  y = tobacco[rows, , drop = FALSE]
  
  fit = cy_complete_staggered_layer(
    y = y,
    omega = omega,
    adoption_time = adoption_time,
    time = time,
    rank = 2L,
    sigma = NULL,
    subgroup_size = 1L,
    tau = 0.5,
    tolerance = 1e-7,
    sigma_tolerance = 1e-2,
    max_iterations = Inf,
    max_sigma_iterations = Inf,
    monte_carlo_draws = 100L,
    lambda_probability = 0.95,
    lambda_multiplier = 1.14,
    processing_order = "matlab",
    seed = NULL,
    verbose = FALSE
  )
  
  truth = y[missing]
  
  RMSE_CY[i] = sqrt(mean(
    (fit$completed_debiased[missing] - truth)^2
  ))
  
  RMSE_CY_nuclear[i] = sqrt(mean(
    (fit$completed_nuclear[missing] - truth)^2
  ))
  
  cat(
    "Iteration", i,
    "- CY:", signif(RMSE_CY[i], 8),
    "- nuclear:", signif(RMSE_CY_nuclear[i], 8),
    "\n"
  )
}

avgRMSE_CY = mean(RMSE_CY)
avgRMSE_CY_nuclear = mean(RMSE_CY_nuclear)

summary = data.frame(
  estimator = c("CY debiased (M_CFMY)", "CY non-debiased (M_nucl)"),
  average_RMSE = c(avgRMSE_CY, avgRMSE_CY_nuclear),
  sd_RMSE = c(sd(RMSE_CY), sd(RMSE_CY_nuclear)),
  iterations = number_iterations
)

cat("\n")
print(summary, row.names = FALSE)
