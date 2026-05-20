# --------------------------------------
# combine pictures: 2 on top, 1 on bottom
# --------------------------------------

library(magick)

setwd("~/Documents/phd/projects/causalMatrix/code/synthetic/4Block")

files = c(
  "figure/compare_mse_vs_K_pooling_oracle_randomXY_SNR1_N1win70-70_T1win60-60.png",
  "figure/Zoomcompare_mse_vs_K_pooling_oracle_randomXY_SNR1_N1win70-70_T1win60-60.png",
  "figure/Zoomcompare_logMSE_vs_K_pooling_oracleLocal_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.png"
)

imgs = image_read(files)

# Make third picture slightly bigger
scale_third = 1.12

info3 = image_info(imgs[3])
new_width3 = round(info3$width * scale_third)

imgs[3] = image_resize(imgs[3], paste0(new_width3, "x"))

# Draw (a)/(b)/(c)
label_draw = function(img, lab, x = 15, y = 35) {
  info = image_info(img)
  h = info$height
  cex_val = max(1, h / 400)
  
  img = image_draw(img)
  par(mar = c(0, 0, 0, 0), xpd = NA, family = "HersheySans")
  
  # soft white box behind text
  w  = strwidth(lab, cex = cex_val)
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

# ------------------------------------------------------------
# Layout:
#   top row    = (a) and (b)
#   bottom row = (c), centered under the top row
# ------------------------------------------------------------

top_row = image_append(image_join(a, b), stack = FALSE)

top_info = image_info(top_row)
c_info = image_info(c)

# If c is wider than top row, resize c down to top-row width
if (c_info$width > top_info$width) {
  c = image_resize(c, paste0(top_info$width, "x"))
  c_info = image_info(c)
}

# Add horizontal white padding to center c under top row
left_pad = floor((top_info$width - c_info$width) / 2)
right_pad = top_info$width - c_info$width - left_pad

c_centered = image_extent(
  c,
  geometry = paste0(top_info$width, "x", c_info$height),
  gravity = "center",
  color = "white"
)

# Stack top row and centered bottom image
combo = image_append(image_join(top_row, c_centered), stack = TRUE)

image_write(combo, "figure/combined_1.png")