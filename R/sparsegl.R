#' Regularization paths for sparse group-lasso models
#'
#' Fits regularization paths for sparse group-lasso penalized learning problems at a
#' sequence of regularization parameters `lambda`.
#'
#' Note that the objective function for least squares is
#' \deqn{RSS/(2n) + \lambda penalty}
#' Users can also tweak the penalty by choosing a different penalty factor.
#'
#' For computing speed reason, if models are not converging or running slowly,
#' consider increasing `eps`, decreasing `nlambda`, or increasing
#' `lambda.factor` before increasing `maxit`.
#'
#' @template param_x-template
#' @template param_y-template
#' @template param_group-template
#' @template param_family-template
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.factor The factor for getting the minimal lambda in the
#'   `lambda` sequence, where `min(lambda) = lambda.factor * max(lambda)`.
#'   `max(lambda)` is the smallest value of `lambda` for which all coefficients
#'   are zero. The default depends on the relationship between \eqn{n}
#'   (the number of rows in the matrix of predictors) and \eqn{p}
#'   (the number of predictors). If \eqn{n \geq p}, the
#'   default is `0.0001`.  If \eqn{n < p}, the default is `0.01`.
#'   A very small value of `lambda.factor` will lead to a
#'   saturated fit. This argument has no effect if there is user-defined
#'   `lambda` sequence.
#' @template param_lambda-template
#' @param pf Penalty factor, a vector of the same length as the total number
#'   of groups. Separate penalty weights can be applied to each group of
#' \eqn{\beta}{beta's}s to allow differential shrinkage. Can be 0 for some
#' groups, which implies no shrinkage, and results in that group always being
#' included in the model. Default value for each entry is the square-root of
#' the corresponding size of each group.
#' @param dfmax Limit the maximum number of groups in the model. Default is
#'   no limit.
#' @param pmax Limit the maximum number of groups ever to be nonzero. For
#'   example once a group enters the model, no matter how many times it exits or
#'   re-enters model through the path, it will be counted only once.
#' @param eps Convergence termination tolerance. Defaults value is `1e-8`.
#' @param maxit Maximum number of outer-loop iterations allowed at fixed lambda
#'   value. Default is `3e8`. If models do not converge, consider increasing
#'   `maxit`.
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#' @param asparse The weight to put on the \eqn{\ell_1}-norm in sparse group
#'   lasso. Default is `0.05`.
#' @param standardize Logical flag for variable standardization (scaling) prior
#'   to fitting the model. Default is TRUE.
#' @param lower_bnd Lower bound for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-positive numbers only.
#'   Default value for each entry is `-Inf`.
#' @param upper_bnd Upper for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-negative numbers only.
#'   Default value for each entry is `Inf`.
#'
#'
#' @return An object with S3 class [sparsegl()].
#' * `call` The call that produced this object.
#' * `b0` Intercept sequence of length `length(lambda)`.
#' * `beta` A `p` x `length(lambda)` sparse matrix of coefficients.
#' * `df` The number of features with nonzero coefficients for each value of
#'     `lambda`.
#' * `dim` Dimension of coefficient matrix.
#' * `lambda` The actual sequence of `lambda` values used.
#' * `npasses` Total number of iterations summed over all `lambda` values.
#' * `jerr` Error flag, for warnings and errors, 0 if no error.
#' * `group` A vector of consecutive integers describing the grouping of the
#'     coefficients.
#'
#'
#' @seealso [plot.sparsegl()], [coef.sparsegl()], [predict.sparsegl()]
#'   and [print.sparsegl()] methods.
#'
#' @export
#'
#'
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit <- sparsegl(X, y, group = groups)
#'
#'
sparsegl <- function(
  x, y, group = NULL, family = c("gaussian", "binomial"),
  nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL, pf = sqrt(bs),
  intercept = TRUE, asparse = 0.05, standardize = TRUE,
  lower_bnd = -Inf, upper_bnd = Inf,
  dfmax = as.integer(max(group)) + 1L,
  pmax = min(dfmax * 1.2, as.integer(max(group))),
  eps = 1e-08, maxit = 3e+08) {

  this.call <- match.call()
  family <- match.arg(family)
  if (!is.matrix(x) && !inherits(x, "sparseMatrix"))
    stop("x must be a matrix")

  if (any(is.na(x))) stop("Missing values in x not allowed!")

  y <- drop(y)
  if (!is.null(dim(y))) stop("y must be a vector or 1-column matrix.")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  if (length(y) != nobs) stop("x and y have different numbers of observations.")

  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    assertthat::assert_that(
      length(group) == nvars,
      msg = "group length does not match the number of predictors in x")
  }

  bn <- as.integer(max(group))  # number of groups
  bs <- as.integer(as.numeric(table(group)))  # number of elements in each group

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn)))
    stop("Groups must be consecutively numbered 1, 2, 3, ...")

  assertthat::assert_that(
    asparse <= 1,
    msg = "asparse must be less than or equal to 1, you may want glmnet::glmnet()")

  if (asparse < 0) {
    asparse <- 0
    warning("asparse must be in [0,1], running ordinary group lasso.")
  }


  iy <- cumsum(bs) # last column of x in each group
  ix <- c(0, iy[-bn]) + 1 # first column of x in each group
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)

  #parameter setup
  assertthat::assert_that(
    length(pf) == bn,
    msg = paste("The length of group-lasso penalty factor must be",
                "same as the number of groups"))

  maxit <- as.integer(maxit)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)

  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    assertthat::assert_that(lambda.factor < 1,
                            msg = "lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin = 1 if user define lambda
    flmin <- as.double(1)
    assertthat::assert_that(all(lambda >= 0),
                            msg = "lambdas must be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)

  ### check on upper/lower bounds
  assertthat::assert_that(all(lower_bnd <= 0),
                          msg = "Lower bounds should be non-positive")
  assertthat::assert_that(all(upper_bnd >= 0),
                          msg = "Upper bounds should be non-negative")
  lower_bnd[lower_bnd == -Inf] <- -9.9e30
  upper_bnd[upper_bnd == Inf] <- 9.9e30
  if (length(lower_bnd) < bn) {
    if (length(lower_bnd) == 1) {
      lower_bnd <- rep(lower_bnd, bn)
    } else {
      stop("Lower bounds must be length 1 or length the number of groups")
    }
  } else {
    lower_bnd <- lower_bnd[seq_len(bn)]
  }
  if (length(upper_bnd) < bn) {
    if (length(upper_bnd) == 1) {
      upper_bnd <- rep(upper_bnd, bn)
    } else {
      stop("Upper bounds must be length 1 or length the number of groups")
    }
  } else {
    upper_bnd <- upper_bnd[seq_len(bn)]
  }
  storage.mode(upper_bnd) <- "double"
  storage.mode(lower_bnd) <- "double"
  ### end check on limits

  # call R sub-function
  fit <- switch(
    family,
    gaussian = sgl_ls(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam,
                flmin, ulam,
                eps, maxit, vnames, group, intr, as.double(asparse),
                standardize, lower_bnd, upper_bnd),
    binomial = sgl_logit(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax,
                      nlam, flmin, ulam,
                      eps, maxit, vnames, group, intr, as.double(asparse),
                      standardize, lower_bnd, upper_bnd)
  )

  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  fit$asparse <- asparse
  class(fit) <- c("sparsegl", class(fit))
  fit
}

