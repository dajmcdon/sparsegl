#' Get coefficients from a `cv.sparsegl` object.
#'
#' This function gets coefficients from a
#' cross-validated [sparsegl()] model, using the stored `"sparsegl.fit"`
#' object, and the optimal value chosen for `lambda`.
#'
#' @param object Fitted [cv.sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored on the CV `object` (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#' @param ... Not used. Other arguments to [predict()].
#'
#'
#' @return The coefficients at the requested value(s) for `lambda`.
#' @seealso [cv.sparsegl()], and [predict.cv.sparsegl()] methods.
#' @method coef cv.sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' cv_fit <- cv.sparsegl(X, y, groups)
#' coef(cv_fit, s = c(0.02, 0.03))
coef.cv.sparsegl <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
    assertthat::assert_that(
        is.numeric(s) || is.character(s),
        msg = "Invalid form for s.")
    if (is.numeric(s)) lambda <- s
    else {
        s <- match.arg(s)
        lambda <- object[[s]]
    }
    coef(object$sparsegl.fit, s = lambda, ...)
}



#' Make predictions from a `cv.sparsegl` object.
#'
#' This function makes predictions from a cross-validated [cv.sparsegl()] object,
#' using the stored `"sparsegl.fit"` object, and the value chosen for `lambda`.
#'
#'
#' @param object Fitted [cv.sparsegl()] object.
#' @param newx Matrix of new values for `x` at which predictions are to be
#'   made. Must be a matrix. See documentation for [predict.sparsegl()].
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored on the CV `object` (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#' @param ... Not used. Other arguments to [predict()].
#'
#' @return A matrix or vector of predicted values.
#' @seealso [cv.sparsegl()], and [coef.cv.sparsegl()] methods.
#'
#' @method predict cv.sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' cv_fit <- cv.sparsegl(X, y, groups)
#' predict(cv_fit, newx = X[50:60, ], s = "lambda.min")
#'
predict.cv.sparsegl <- function(object, newx,
                                s = c("lambda.1se", "lambda.min"), ...) {
    assertthat::assert_that(
        is.numeric(s) || is.character(s),
        msg = "Invalid form for s.")
    if (is.numeric(s)) lambda <- s
    else {
        s <- match.arg(s)
        lambda <- object[[s]]
    }
    predict(object$sparsegl.fit, newx, s = lambda, ...)
}
