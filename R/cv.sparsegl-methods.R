#' Extract coefficients from a `cv.sparsegl` object.
#'
#' This function etracts coefficients from a
#' cross-validated [sparsegl()] model, using the stored `"sparsegl.fit"`
#' object, and the optimal value chosen for `lambda`.
#'
#' @param object Fitted [cv.sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored in the CV object (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#' @param ... Not used.
#'
#'
#' @return The coefficients at the requested value(s) for `lambda`.
#' @seealso [cv.sparsegl()] and [predict.cv.sparsegl()].
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
  rlang::check_dots_empty()
  if (!(is.numeric(s) || is.character(s)))
    cli::cli_abort("Invalid form for `s`.")
  if (is.numeric(s)) lambda <- s
  else {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  coef(object$sparsegl.fit, s = lambda)
}



#' Make predictions from a `cv.sparsegl` object.
#'
#' This function makes predictions from a cross-validated [cv.sparsegl()] object,
#' using the stored `sparsegl.fit` object, and the value chosen for `lambda`.
#'
#' @inheritParams predict.sparsegl
#' @param object Fitted [cv.sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored in the CV object (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#'
#' @return A matrix or vector of predicted values.
#' @seealso [cv.sparsegl()] and [coef.cv.sparsegl()].
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
predict.cv.sparsegl <- function(
    object, newx,
    s = c("lambda.1se", "lambda.min"),
    type = c("link", "response", "coefficients", "nonzero", "class"), ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (!(is.numeric(s) || is.character(s)))
    cli::cli_abort("Invalid form for `s`.")
  if (is.numeric(s)) lambda <- s
  else {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  predict(object$sparsegl.fit, newx, s = lambda, type = type)
}

#' @method fitted cv.sparsegl
#' @export
fitted.cv.sparsegl <- function(object, ...) {
  cli::cli_abort(c(
    "!" = "Because design matrices are typically large, these are not stored ",
    "!" = "in the estimated `cv.sparsegl` object. Use `predict()` instead, and ",
    "!" = "pass in the original data."))
}


#' @method summary cv.sparsegl
#' @export
summary.cv.sparsegl <- function(object, ...) {
  rlang::check_dots_empty()
  optlams <- c(object$lambda.1se, object$lambda.min)
  optlax <- c(1, match(optlams, object$lambda), length(object$lambda))
  tab <- with(object, data.frame(
    lambda = lambda[optlax],
    index = optlax,
    cvm = cvm[optlax],
    cvsd = cvsd[optlax],
    nnzero = nnzero[optlax],
    active_grps = active_grps[optlax])
  )
  rownames(tab) <- c("Max.", "lambda.1se", "lambda.min", "Min.")
  out <- structure(
    list(
      call = object$call,
      error_measure = object$name,
      table = tab
    ),
    class = "summary.cvsparsegl"
  )
  out

}

#' @method print summary.cvsparsegl
#' @export
print.summary.cvsparsegl <- function(
    x,
    digits = max(3, getOption("digits") - 3), ...
) {

  rlang::check_dots_empty()
  lambda_warning = NULL
  if (x$table$index[2] == 1) lambda_warning = "smallest"
  if (x$table$index[3] == x$table$index[4]) lambda_warning = "largest"
  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)

  cat("Error measure: ", x$error_measure, "\n\n")

  if (!is.null(lambda_warning)) {
    cat("Warning: the CV minimum occurred at the", lambda_warning,
        "lambda in the path.\n\n")
  }

  print(x$tab, digits = digits)
  cat("\n")
}

#' @method print cv.sparsegl
#' @export
print.cv.sparsegl <- function(x, digits = max(3, getOption("digits") - 3),
                              ...) {

  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}

