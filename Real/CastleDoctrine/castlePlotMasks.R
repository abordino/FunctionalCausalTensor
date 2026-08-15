# Castle: display saved initial and terminal masks -----------------------------

setwd("~/Desktop/code")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# Configuration ---------------------------------------------------------------

STEP_SIZE = 3L
JUMP = 80L
REP = 1L
RANK_TO_LOAD = 1L
RUN_CY = FALSE

RESULT_FILE = NULL

results_dir = file.path(
  "Results",
  "castle_staggered_masking_rank2",
  paste0("step", STEP_SIZE)
)

cy_file_tag = if (RUN_CY) "" else "_noCY"
rep_file_tag = if (REP == 1L) "" else paste0("_rep", REP)

if (is.null(RESULT_FILE)) {
  RESULT_FILE = file.path(
    results_dir,
    paste0(
      "castle_staggered_masking_results_motor_rank",
      RANK_TO_LOAD,
      "_jump",
      JUMP,
      "_step",
      STEP_SIZE,
      rep_file_tag,
      cy_file_tag,
      ".rds"
    )
  )
}

if (!file.exists(RESULT_FILE)) {
  stop(
    "Saved Castle result file was not found:\n",
    RESULT_FILE
  )
}

# Load saved masks -------------------------------------------------------------

results = readRDS(RESULT_FILE)

if (is.null(results$mask_sequences)) {
  stop(
    "The saved result does not contain `mask_sequences`."
  )
}

if (is.null(results$D_original)) {
  stop(
    "The saved result does not contain `D_original`."
  )
}

mask_sequences = results$mask_sequences
D_original = results$D_original

target_layer = results$metadata$target_layer
target_label = results$metadata$target_label

target_idx = match(
  target_layer,
  dimnames(D_original)[[3]]
)

if (is.na(target_idx)) {
  stop(
    "Could not find the saved target layer in D_original."
  )
}

# M0 is the original target-layer treatment/missing region.
M0 = D_original[, , target_idx] == 1

# Plot helper ------------------------------------------------------------------

make_mask_plot = function(
    M,
    title,
    subtitle
) {
  mask_status = matrix(
    "Observed",
    nrow = nrow(M),
    ncol = ncol(M),
    dimnames = dimnames(M)
  )
  
  mask_status[M0] = "Original missing"
  mask_status[M & !M0] = "Artificially masked"
  
  plot_data = as.data.frame.table(
    mask_status,
    responseName = "mask_status",
    stringsAsFactors = FALSE
  )
  
  names(plot_data)[1:2] = c(
    "state",
    "year"
  )
  
  state_levels = rownames(M0)
  year_levels = colnames(M0)
  
  if (is.null(state_levels)) {
    state_levels = unique(plot_data$state)
  }
  
  if (is.null(year_levels)) {
    year_levels = unique(plot_data$year)
  }
  
  plot_data = as_tibble(plot_data) %>%
    mutate(
      state = factor(
        state,
        levels = rev(state_levels)
      ),
      year = factor(
        year,
        levels = year_levels
      ),
      mask_status = factor(
        mask_status,
        levels = c(
          "Observed",
          "Artificially masked",
          "Original missing"
        )
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = year,
      y = state,
      fill = mask_status
    )
  ) +
    geom_tile(linewidth = 0.05) +
    scale_fill_manual(
      values = c(
        "Observed" = "#f7f7f7",
        "Artificially masked" = "#d62728",
        "Original missing" = "#525252"
      ),
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Year",
      y = "State",
      fill = "Mask status"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 6)
    )
}

# Plot the first saved artificial mask from every replicate --------------------

for (rep_id in seq_along(mask_sequences)) {
  mask_sequence = mask_sequences[[rep_id]]
  
  if (!length(mask_sequence)) {
    warning(
      "Replicate ",
      rep_id,
      " has an empty mask sequence."
    )
    next
  }
  
  initial_mask_object = mask_sequence[[1L]]
  
  initial_plot = make_mask_plot(
    M = initial_mask_object$M,
    title = paste0(
      "Castle initial artificial mask: replicate ",
      rep_id
    ),
    subtitle = paste0(
      target_label,
      "; ",
      initial_mask_object$n_artificial_masked,
      " artificially masked cells; JUMP = ",
      results$metadata$JUMP
    )
  )
  
  print(initial_plot)
}

# Plot the common terminal mask once -------------------------------------------

terminal_masks = lapply(
  mask_sequences,
  function(mask_sequence) {
    mask_sequence[[length(mask_sequence)]]$M
  }
)

terminal_masks_identical = all(
  vapply(
    terminal_masks,
    function(M) {
      identical(M, terminal_masks[[1L]])
    },
    logical(1)
  )
)

if (!terminal_masks_identical) {
  warning(
    paste0(
      "Terminal masks are not identical across replicates. ",
      "The final plot uses replicate 1."
    )
  )
}

terminal_sequence = mask_sequences[[1L]]
terminal_mask_object =
  terminal_sequence[[length(terminal_sequence)]]

terminal_plot = make_mask_plot(
  M = terminal_mask_object$M,
  title = "Castle final artificial mask",
  subtitle = paste0(
    target_label,
    "; ",
    terminal_mask_object$n_artificial_masked,
    " artificially masked cells; common across replicates = ",
    terminal_masks_identical
  )
)

print(terminal_plot)

cat(
  "Loaded: ",
  RESULT_FILE,
  "\nDisplayed ",
  length(mask_sequences),
  " initial masks and one terminal mask.\n",
  sep = ""
)
