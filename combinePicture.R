library(magick)

setwd("~/Desktop/code")

#------------------------------------------------------------

files = c(
  "Plots/synthetic_tucker2_50x11x4_staggered_masking_signalrank2x2/synthetic_tucker2_50x11x4_staggered_masking_summary_layer2_rank1_jump80_step3_rep30.png",
  "Plots/synthetic_tucker2_50x11x4_staggered_masking_signalrank2x2/synthetic_tucker2_50x11x4_staggered_masking_summary_layer2_rank2_jump80_step3_rep30.png",
  "Plots/synthetic_tucker2_50x11x4_staggered_masking_signalrank2x2/synthetic_tucker2_50x11x4_staggered_masking_summary_layer2_rank3_jump80_step3_rep30.png"
)

imgs = image_read(files)

n_imgs = length(imgs)

if (n_imgs == 0) {
  stop("No images provided.")
}

#------------------------------------------------------------

scale_vec = c(1, rep(1, n_imgs - 1))

for (i in seq_len(n_imgs)) {
  
  if (scale_vec[i] != 1) {
    
    info_i = image_info(imgs[i])
    
    new_width_i = round(
      info_i$width * scale_vec[i]
    )
    
    imgs[i] = image_resize(
      imgs[i],
      paste0(new_width_i, "x"),
      filter = "Lanczos"
    )
    
  }
}

#------------------------------------------------------------

img_list = lapply(
  seq_len(n_imgs),
  function(i) {
    image_trim(imgs[i])
  }
)

#------------------------------------------------------------

info_list = lapply(
  img_list,
  image_info
)

cell_width = max(
  vapply(
    info_list,
    function(z) z$width,
    numeric(1)
  )
)

cell_height = max(
  vapply(
    info_list,
    function(z) z$height,
    numeric(1)
  )
)

inner_pad = 4

cell_width = cell_width + 2 * inner_pad

cell_height = cell_height + 2 * inner_pad

#------------------------------------------------------------

pad_to_cell = function(
    img,
    cell_width,
    cell_height
) {
  
  image_extent(
    img,
    geometry = paste0(
      cell_width,
      "x",
      cell_height
    ),
    gravity = "center",
    color = "white"
  )
}

cells = lapply(
  img_list,
  function(img) {
    
    pad_to_cell(
      img,
      cell_width,
      cell_height
    )
    
  }
)

#------------------------------------------------------------

panel_tags = paste0(
  "(",
  letters[seq_len(n_imgs)],
  ")"
)

for (i in seq_len(n_imgs)) {
# for (i in 1) {
  
  tag_file = tempfile(
    fileext = ".png"
  )
  
  png(
    filename = tag_file,
    width = cell_width,
    height = cell_height,
    units = "px",
    bg = "transparent"
  )
  
  par(
    mar = c(0, 0, 0, 0),
    xaxs = "i",
    yaxs = "i"
  )
  
  plot.new()
  
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )
  
  text(
    x = 0.05,
    y = 0.05,
    labels = panel_tags[i],
    adj = c(1, 1),
    cex = 5.8,
    font = 2
  )
  
  dev.off()
  
  tag_layer = image_read(
    tag_file
  )
  
  cells[[i]] = image_composite(
    cells[[i]],
    tag_layer,
    operator = "over",
    gravity = "northwest"
  )
  
  unlink(
    tag_file
  )
  
}

#------------------------------------------------------------

if (n_imgs <= 3) {
  
  nrow = 1
  ncol = n_imgs
  
} else {
  
  ncol = ceiling(
    sqrt(n_imgs)
  )
  
  nrow = ceiling(
    n_imgs / ncol
  )
}

#------------------------------------------------------------

n_total_cells = nrow * ncol

n_missing = n_total_cells - n_imgs

if (n_missing > 0) {
  
  blank_cell = image_blank(
    width = cell_width,
    height = cell_height,
    color = "white"
  )
  
  for (i in seq_len(n_missing)) {
    cells[[length(cells) + 1]] = blank_cell
  }
}

#------------------------------------------------------------

gap = 0

#------------------------------------------------------------

combine_row = function(
    row_cells,
    gap = 0,
    cell_height
) {
  
  if (length(row_cells) == 1) {
    return(row_cells[[1]])
  }
  
  pieces = list()
  
  for (j in seq_along(row_cells)) {
    
    pieces[[length(pieces) + 1]] = row_cells[[j]]
    
    if (
      j < length(row_cells) &&
      gap > 0
    ) {
      
      pieces[[length(pieces) + 1]] = image_blank(
        width = gap,
        height = cell_height,
        color = "white"
      )
    }
    
  }
  
  row_images = do.call(
    c,
    pieces
  )
  
  image_append(
    row_images,
    stack = FALSE
  )
}

#------------------------------------------------------------

rows = vector(
  "list",
  nrow
)

for (r in seq_len(nrow)) {
  
  first = (r - 1) * ncol + 1
  
  last = r * ncol
  
  rows[[r]] = combine_row(
    cells[first:last],
    gap = gap,
    cell_height = cell_height
  )
}

#------------------------------------------------------------

if (nrow == 1) {
  
  combo = rows[[1]]
  
} else {
  
  pieces = list()
  
  row_width = image_info(
    rows[[1]]
  )$width
  
  for (r in seq_len(nrow)) {
    
    pieces[[length(pieces) + 1]] = rows[[r]]
    
    if (
      r < nrow &&
      gap > 0
    ) {
      
      pieces[[length(pieces) + 1]] = image_blank(
        width = row_width,
        height = gap,
        color = "white"
      )
    }
    
  }
  
  all_rows = do.call(
    c,
    pieces
  )
  
  combo = image_append(
    all_rows,
    stack = TRUE
  )
}

#------------------------------------------------------------

#------------------------------------------------------------

outer_margin = 30

combo = image_border(
  combo,
  color = "white",
  geometry = paste0(outer_margin, "x", outer_margin)
)

#------------------------------------------------------------

#------------------------------------------------------------

output_file = "Plots/combined_SyntheticALL.png"

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

image_write(
  combo,
  output_file
)

#------------------------------------------------------------