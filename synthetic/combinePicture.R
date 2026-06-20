library(magick)

setwd("~/Documents/phd/projects/causalMatrix/code/synthetic/4Block")

files = c(
  "figure/compare_mse_vs_K_pooling_oracle_randomXY_SNR1_N1win70-70_T1win60-60.png",
  
  "figure/Zoomcompare_mse_vs_K_pooling_oracle_randomXY_SNR1_N1win70-70_T1win60-60.png",
  
  "figure/compare_mse_vs_K_pooling_oracle_SNR1_N1win30-70_T1win30-60_leadingEigXY_trueRank6_estRank11.png",
  
  "figure/Zoomcompare_logMSE_vs_K_pooling_oracleLocal_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.png"
)

imgs = image_read(files)

# ------------------------------------------------------------
#  Scaling
# ------------------------------------------------------------

scale_vec = c(1.00, 1.00, 1.00, 1.12)

for (i in seq_along(imgs)) {
  if (scale_vec[i] != 1) {
    info_i = image_info(imgs[i])
    new_width_i = round(info_i$width * scale_vec[i])
    imgs[i] = image_resize(imgs[i], paste0(new_width_i, "x"))
  }
}

# ------------------------------------------------------------
# Draw (a)/(b)/(c)/(d)
# ------------------------------------------------------------

label_draw = function(img, lab, x = 15, y = 35) {
  info = image_info(img)
  h = info$height
  cex_val = max(1, h / 400)
  
  img = image_draw(img)
  par(mar = c(0, 0, 0, 0), xpd = NA, family = "HersheySans")
  
  w = strwidth(lab, cex = cex_val)
  ht = strheight(lab, cex = cex_val)
  
  rect(
    x - 6, y - 1.6 * ht,
    x + 6 + 1.1 * w, y + 6,
    col = rgb(1, 1, 1, 0.7),
    border = NA
  )
  
  text(
    x, y,
    labels = lab,
    cex = cex_val,
    col = "black",
    adj = c(0, 1)
  )
  
  dev.off()
  
  img
}

a = label_draw(imgs[1], "(a)")
b = label_draw(imgs[2], "(b)")
c = label_draw(imgs[3], "(c)")
d = label_draw(imgs[4], "(d)")

a = image_trim(a)
b = image_trim(b)
c = image_trim(c)
d = image_trim(d)

global_scale = 3.50

a = image_resize(a, paste0(round(image_info(a)$width * global_scale), "x"))
b = image_resize(b, paste0(round(image_info(b)$width * global_scale), "x"))
c = image_resize(c, paste0(round(image_info(c)$width * global_scale), "x"))
d = image_resize(d, paste0(round(image_info(d)$width * global_scale), "x"))

pad_to_cell = function(img, cell_width, cell_height) {
  image_extent(
    img,
    geometry = paste0(cell_width, "x", cell_height),
    gravity = "center",
    color = "white"
  )
}

info_list = lapply(list(a, b, c, d), image_info)

cell_width = max(sapply(info_list, function(z) z$width))
cell_height = max(sapply(info_list, function(z) z$height))

inner_pad = 4

cell_width = cell_width + 2 * inner_pad
cell_height = cell_height + 2 * inner_pad

a_cell = pad_to_cell(a, cell_width, cell_height)
b_cell = pad_to_cell(b, cell_width, cell_height)
c_cell = pad_to_cell(c, cell_width, cell_height)
d_cell = pad_to_cell(d, cell_width, cell_height)

# ------------------------------------------------------------
# Combine into 2x2
# ------------------------------------------------------------

gap = 0

spacer_v = image_blank(width = gap, height = cell_height, color = "white")
spacer_h = image_blank(width = 2 * cell_width + gap, height = gap, color = "white")

top_row = image_append(image_join(a_cell, spacer_v, b_cell), stack = FALSE)
bottom_row = image_append(image_join(c_cell, spacer_v, d_cell), stack = FALSE)

combo = image_append(image_join(top_row, spacer_h, bottom_row), stack = TRUE)

combo = image_trim(combo)

image_write(combo, "figure/combined_1.png")