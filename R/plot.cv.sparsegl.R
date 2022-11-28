#' Plot cross-validation curves produced from a `cv.sparsegl` object.
#'
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the `lambda` values used.
#'
#' A [ggplot2::ggplot()] plot is produced. Additional user
#' modifications may be added as desired.
#'
#' @param x Fitted `"cv.sparsegl"` object, produced with [cv.sparsegl()].
#' @param log_axis Apply log scaling to the requested axes.
#' @param sign.lambda Either plot against `log(lambda)` (default) or the
#'    reverse if `sign.lambda < 0`.
#' @param ... Not used.
#' @seealso [cv.sparsegl()].
#' @method plot cv.sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' cv_fit <- cv.sparsegl(X, y, groups)
#' plot(cv_fit)
plot.cv.sparsegl <- function(x, log_axis = c("xy", "x", "y", "none"),
                             sign.lambda = 1, ...) {
  rlang::check_dots_empty()
  cvobj <- x
  dat <- data.frame("X" = sign.lambda * cvobj$lambda,
                    "y" = cvobj$cvm,
                    "upper" = cvobj$cvupper,
                    "lower" = cvobj$cvlo)
  log_axis <- match.arg(log_axis)
  sign.lambda <- sign(sign.lambda)
  g <- dat %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$X, y = .data$y)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      color = 'darkgrey') +
    ggplot2::geom_point(color = 'darkblue') +
    ggplot2::xlab("Lambda") +
    ggplot2::ylab(cvobj$name) +
    ggplot2::theme_bw()
  switch(log_axis,
         xy = g + ggplot2::scale_x_log10() + ggplot2::scale_y_log10(),
         x = g + ggplot2::scale_x_log10(),
         y = g + ggplot2::scale_y_log10(),
         g)
}
