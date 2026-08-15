setwd("~/Desktop/code")

suppressPackageStartupMessages(library(yaml))

# ----------------------------------------------------------
# 1. Reproduce synthetic 4-block simulation. 
#    Generate 2x2 panel shown in Figure 5.
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$synthetic4block)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)

  simulation_files = c(
    "Synthetic/4Block/sim_varyK.R",
    "Synthetic/4Block/sim_varyK_SNRZoom.R",
    "Synthetic/4Block/sim_varyKZoom.R",
    "Synthetic/4Block/sim_varyKrobustness.R"
  )

  plot_file = "Synthetic/4Block/combine4plots.R"
  files_to_run = c(simulation_files, plot_file)

  missing_files = files_to_run[!file.exists(files_to_run)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }

  for (file in simulation_files) {
    message("Running ", file)
    source(file, echo = FALSE)
  }

  message("Creating combined plot with ", plot_file)
  source(plot_file, echo = FALSE)

  message("Synthetic/4Block simulations and plotting completed.")
} else {
  message("Synthetic/4Block is disabled in wrapper.yaml")
}

# ----------------------------------------------------------
# 2. Reproduce synthetic data simulation on tau-sensitivity.
#    Generates Figure 6(a).
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$sensitivityTau)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)

  simulation_file = "Synthetic/Staggered/SensitivityTau/SensitivityTau.R"
  plot_file = "Synthetic/Staggered/SensitivityTau/SensitivityTauPlot.R"
  files_to_run = c(simulation_file, plot_file)

  missing_files = files_to_run[!file.exists(files_to_run)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }

  message("Running ", simulation_file)
  source(simulation_file, echo = FALSE)

  message("Creating plots with ", plot_file)
  source(plot_file, echo = FALSE)

  message("Tau-sensitivity simulation and plotting completed.")
} else {
  message("Tau sensitivity is disabled in wrapper.yaml")
}

# ----------------------------------------------------------
# 3. Reproduce runtime-versus-accuracy simulation
#    Generates Figure 6(b,c).
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$runtimeVsAccuracy)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)
  
  simulation_file = "Synthetic/Staggered/runtimeVsAccuracy/runtimeVsAccuracy.R"
  plot_file = "Synthetic/Staggered/runtimeVsAccuracy/runtimeVsAccuracyPlot.R"
  files_to_run = c(simulation_file, plot_file)
  
  missing_files = files_to_run[!file.exists(files_to_run)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  message("Running ", simulation_file)
  source(simulation_file, echo = FALSE)
  
  message("Creating plots with ", plot_file)
  source(plot_file, echo = FALSE)
  
  message("Runtime-versus-accuracy simulation and plotting completed.")
} else {
  message("Runtime versus accuracy is disabled in wrapper.yaml.")
}

# ----------------------------------------------------------
# 4. Reproduce Castle Doctrine simulation results.
#    Generates Figures 7--8 and Table 2. 
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$castle)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)
  
  simulation_file_1 = "Real/CastleDoctrine/CastleSimul.R"
  simulation_file_2 = "Real/CastleDoctrine/CastleSimulMasking.R"
  plot_file = "Real/CastleDoctrine/castlePlot.R"
  files_to_run = c(simulation_file_1, simulation_file_2, plot_file)
  
  missing_files = files_to_run[!file.exists(files_to_run)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  message("Running ", simulation_file_1)
  source(simulation_file_1, echo = FALSE)
  
  message("Running ", simulation_file_2)
  source(simulation_file_2, echo = FALSE)
   
  message("Creating plots with ", plot_file)
  source(plot_file, echo = FALSE)
  
  message("Castle Doctrine simulation and plotting completed.")
} else {
  message("Castle Doctrine is disabled in wrapper.yaml.")
}

# ----------------------------------------------------------
# 5. Reproduce COVID-19 simulation results.
#    Generates Figure 9 and Table 3.
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$covid)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)
  
  simulation_file_1 = "Real/CovidOx/CovidSimul.R"
  simulation_file_2 = "Real/CovidOx/CovidSimulMasking.R"
  plot_file = "Real/CovidOx/CovidPlot.R"
  files_to_run = c(simulation_file_1) #, simulation_file_2, plot_file)
  
  missing_files = files_to_run[!file.exists(files_to_run)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  message("Running ", simulation_file_1)
  source(simulation_file_1, echo = FALSE)
  
  # message("Running ", simulation_file_2)
  # source(simulation_file_2, echo = FALSE)
  
  # message("Creating plots with ", plot_file)
  # source(plot_file, echo = FALSE)
  
  message("Covid-19 simulation and plotting completed.")
} else {
  message("Covid-19 is disabled in wrapper.yaml.")
}


# ----------------------------------------------------------
# 6. Reproduce the synthetic simulations study with artificial 
#    masking included in the Supplement. 
#    Generates Figures S1--S2.
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$masking)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)
  
  simulation_file = "Synthetic/Staggered/maskingError/simulMaskingAVG.R"
  plot_file = "Synthetic/Staggered/maskingError/plotMaskingAVG.R"
  files_to_run = c(simulation_file, plot_file)
  
  missing_files = files_to_run[!file.exists(files_to_run)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  message("Running ", simulation_file)
  source(simulation_file, echo = FALSE)
  
  message("Creating plots with ", plot_file)
  source(plot_file, echo = FALSE)
  
  message("Masking error simulation and plotting completed.")
} else {
  message("Masking error is disabled in wrapper.yaml.")
}


# ----------------------------------------------------------
# 7. (Optional) Reproduce produce results similar to those in the first 
# column of Table A.4 in the supplementary file in
# https://www.tandfonline.com/doi/suppl/10.1080/01621459.2024.2380105?scroll=top
# This is to ensure that the methods in Choi and Yuan are correctly implemented
# in R (rather than MATLAB)
# ----------------------------------------------------------
config = read_yaml("wrapper.yaml")

if (isTRUE(config$simulations$CY)) {
  dir.create("Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Plots", showWarnings = FALSE, recursive = TRUE)
  
  simulation_file_1 = "Synthetic/Staggered/CYzero/CYzero.R"
  simulation_file_2 = "Real/tobacco/reproduce_tobacco_matlab.R"
  files_to_run = c(simulation_file_1, simulation_file_2)
  
  missing_files = files_to_run[!file.exists(files_to_run)]
  
  if (length(missing_files) > 0L) {
    stop(
      "Missing required file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  message("Running ", simulation_file_1)
  source(simulation_file_1, echo = FALSE)
  
  message("Running ", simulation_file_2)
  source(simulation_file_2, echo = FALSE)
  
  message("CY diagnostics done.")
} else {
  message("CY is disabled in wrapper.yaml.")
}