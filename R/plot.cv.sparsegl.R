#' Plot cross-validation curves produced from a `cv.sparsegl` object.
#'
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the `lambda` values used. This function is
#' modified based on [glmnet::plot.cv.glmnet()].
#'
#' A plot is produced, a [ggplot2::ggplot()] object. Additional user
#' modifications can be added as desired.
#'
#' @param x Fitted [cv.sparsegl()] object
#' @param sign.lambda Either plot against `log(lambda)` (default) or its
#'   negative if `sign.lambda = -1`.
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
plot.cv.sparsegl <- function(x, sign.lambda = 1, ...) {
    cvobj <- x
    dat <- data.frame("X" = sign.lambda * cvobj$lambda,
                      "y" = cvobj$cvm,
                      "upper" = cvobj$cvupper,
                      "lower" = cvobj$cvlo)
    dat %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$X, y = .data$y)) +
        ggplot2::geom_point(color = 'red') +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
            width = 0.1,
            color = 'darkgrey') +
        ggplot2::xlab("Lambda") +
        ggplot2::ylab(cvobj$name) +
        ggplot2::scale_x_log10(
            labels = scales::trans_format("log10", scales::math_format())) +
        ggplot2::theme_bw()
}
