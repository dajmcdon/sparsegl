library(rgl)
rgl_orientation <- readRDS("inst/logo/rgl_view.rds")

f <- function(x, y) {
  pmax(1 - 0.5 * sqrt(2 * (x^2 + y^2)) - 0.5 * (abs(x) + abs(y)), 0)
}

draw_constraint <- function(colour = "cornflowerblue", alpha = .6, bg = "white") {
  persp3d(f, xlim = c(-1, 1), ylim = c(-1, 1), col = colour,
          zlim = c(-1, 1), alpha = alpha, box = FALSE, axes = FALSE,
          xlab = "", ylab = "", zlab = "")
  persp3d(\(x,y) -f(x, y), xlim = c(-1, 1), ylim = c(-1, 1),
          col = colour, add = TRUE, alpha = alpha,
          box = FALSE, axes = FALSE)
  par3d(userMatrix = rgl_orientation, zoom = .575)
  bg3d(color = "white")
  snapshot3d("inst/logo/constraint.png", webshot = FALSE)
  close3d()
}

bg <- "#950f27"
border <- "#c41232"

draw_constraint("#284dab", alpha = .6) # need to crop it
rgl.quit()

library(hexSticker)
sysfonts::font_add_google("Montserrat")
sticker(
  "inst/logo/constraint.png",
  package = "sparsegl", filename = "man/figures/logo.png",
  s_x = 1, s_y = .75, s_width = .6, s_height = .6,
  p_size = 20, p_x = 1, p_family = "Montserrat", p_color = "white",
  h_fill = bg, h_color = border, h_size = 1.5,
)


