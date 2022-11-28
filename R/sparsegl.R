#' Regularization paths for sparse group-lasso models
#'
#' @description
#' Fits regularization paths for sparse group-lasso penalized learning problems at a
#' sequence of regularization parameters `lambda`.
#' Note that the objective function for least squares is
#' \deqn{RSS/(2n) + \lambda penalty}
#' Users can also tweak the penalty by choosing a different penalty factor.
#'
#'
#' @param x Double. A matrix of predictors, of dimension
#'   \eqn{n \times p}{n * p}; each row
#'   is a vector of measurements and each column is a feature. Objects of class
#'   [`Matrix::sparseMatrix`] are supported.
#' @param y Double/Integer/Factor. The response variable.
#'   Quantitative for `family="gaussian"` and for other exponential families.
#'   If `family="binomial"` should be either a factor with two levels or
#'   a vector of integers taking 2 unique values. For a factor, the last level
#'   in alphabetical order is the target class.
#' @param group Integer. A vector of consecutive integers describing the
#'   grouping of the coefficients (see example below).
#' @param family Character or function. Specifies the generalized linear model
#'   to use. Valid options are:
#'   * `"gaussian"` - least squares loss (regression, the default),
#'   * `"binomial"` - logistic loss (classification)
#'
#'   For any other type, a valid [stats::family()] object may be passed. Note
#'   that these will generally be much slower to estimate than the built-in
#'   options passed as strings. So for example, `family = "gaussian"` and
#'   `family = gaussian()` will produce the same results, but the first
#'   will be much faster.
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.factor A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `min(lambda) = lambda.factor * max(lambda)`.
#'   `max(lambda)` is the smallest value of `lambda` for which all coefficients
#'   are zero. The default depends on the relationship between \eqn{n}
#'   (the number of rows in the matrix of predictors) and \eqn{p}
#'   (the number of predictors). If \eqn{n \geq p}, the
#'   default is `0.0001`.  If \eqn{n < p}, the default is `0.01`.
#'   A very small value of `lambda.factor` will lead to a
#'   saturated fit. This argument has no effect if there is user-defined
#'   `lambda` sequence.
#' @param lambda A user supplied `lambda` sequence. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the smallest value
#'   of `lambda` that would give the null model (all coefficient estimates equal
#'   to zero), and `lambda.factor`. Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param pf_group Penalty factor on the groups, a vector of the same
#'   length as the total number of groups. Separate penalty weights can be applied
#'   to each group of \eqn{\beta}{beta's}s to allow differential shrinkage.
#'   Can be 0 for some
#'   groups, which implies no shrinkage, and results in that group always being
#'   included in the model (depending on `pf_sparse`). Default value for each
#'   entry is the square-root of the corresponding size of each group.
#' @param pf_sparse Penalty factor on l1-norm, a vector the same length as the
#'   total number of columns in `x`. Each value corresponds to one predictor
#'   Can be 0 for some predictors, which
#'   implies that predictor will be receive only the group penalty.
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
#' @param asparse The relative weight to put on the \eqn{\ell_1}-norm in
#'   sparse group lasso. Default is `0.05` (resulting in `0.95` on the
#'   \eqn{\ell_2}-norm).
#' @param standardize Logical flag for variable standardization (scaling) prior
#'   to fitting the model. Default is TRUE.
#' @param lower_bnd Lower bound for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-positive numbers only.
#'   Default value for each entry is `-Inf`.
#' @param upper_bnd Upper for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-negative numbers only.
#'   Default value for each entry is `Inf`.
#' @param weights Double vector. Optional observation weights. These can
#'   only be used with a [stats::family()] object.
#' @param offset Double vector. Optional offset (constant predictor without a
#'   corresponding coefficient). These can only be used with a
#'   [stats::family()] object.
#' @param warm List created with [make_irls_warmup()]. These can only be used
#'   with a [stats::family()] object, and is not typically necessary even then.
#' @param trace_it Scalar integer. Larger values print more output during
#'   the irls loop. Typical values are `0` (no printing), `1` (some printing
#'   and a progress bar), and `2` (more detailed printing).
#'   These can only be used with a [stats::family()] object.
#'
#' @return An object with S3 class `"sparsegl"`. Among the list components:
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
#' * `nobs` The number of observations used to estimate the model.
#'
#' If `sparsegl()` was called with a [stats::family()] method, this may also
#' contain information about the deviance and the family used in fitting.
#'
#'
#' @seealso [cv.sparsegl()] and the [`plot()`][plot.sparsegl()],
#'   [`predict()`][predict.sparsegl()], and [`coef()`][coef.sparsegl()]
#'   methods for `"sparsegl"` objects.
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
#' yp <- rpois(n, abs(X %*% beta_star))
#' fit_pois <- sparsegl(X, yp, group = groups, family = poisson())
sparsegl <- function(
  x, y, group = NULL, family = c("gaussian", "binomial"),
  nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL, pf_group = sqrt(bs), pf_sparse = rep(1, nvars),
  intercept = TRUE, asparse = 0.05, standardize = TRUE,
  lower_bnd = -Inf, upper_bnd = Inf,
  weights = NULL, offset = NULL, warm = NULL,
  trace_it = 0,
  dfmax = as.integer(max(group)) + 1L,
  pmax = min(dfmax * 1.2, as.integer(max(group))),
  eps = 1e-08, maxit = 3e+06) {

  this.call <- match.call()
  if (!is.matrix(x) && !inherits(x, "sparseMatrix"))
    abort("`x` must be a matrix.")

  if (any(is.na(x))) abort("Missing values in `x` are not allowed!")

  y <- drop(y)
  if (!is.null(dim(y))) abort("`y` must be a vector or 1-column matrix.")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  if (length(y) != nobs) {
    cli::cli_abort(c("`x` has {nobs} rows while `y` has {length(y)}.",
                     "Both should be the same."))
  }

  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    if (length(group) != nvars) {
      cli::cli_abort(c("The length of group is {length(group)}.",
                       "It must match the number of columns in `x`: {nvars}"))
    }
  }

  bn <- as.integer(max(group))  # number of groups
  bs <- as.integer(as.numeric(table(group)))  # number of elements in each group

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn)))
    abort("Groups must be consecutively numbered 1, 2, 3, ...")

  if (asparse > 1) {
    abort(c("`asparse` must be less than or equal to 1",
            i = "You may want glmnet::glmnet() instead."))}

  if (asparse < 0) {
    asparse <- 0
    rlang::warn("asparse must be in [0,1], running ordinary group lasso.")
  }
  if (any(pf_sparse < 0)) abort("`pf_sparse` must be non-negative.")
  if (any(is.infinite(pf_sparse)))
    abort("`pf_sparse` may not be infinite. Simply remove the column from `x`.")
  if (any(pf_group < 0)) abort("`pf_group` must be non-negative.")
  if (any(is.infinite(pf_group)))
    abort("`pf_group` may not be infinite. Simply remove the group from `x`.")
  if (all(pf_sparse == 0)) {
    if (asparse > 0) {
      abort(c("`pf_sparse` is identically 0 but `asparse` suggests",
            "!" = "some L1 penalty is desired."))
    } else {
      rlang::warn("`pf_sparse` was set to 1 since `asparse == 0`.")
      pf_sparse = rep(1, nvars)
    }
  }

  ## Note: should add checks to see if any columns are completely unpenalized
  ## This is not currently expected.

  iy <- cumsum(bs) # last column of x in each group
  ix <- c(0, iy[-bn]) + 1 # first column of x in each group
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)

  #parameter setup
  if (length(pf_group) != bn) {
    cli::cli_abort(
      "The length of `pf_group` must be the same as the number of groups {bn}."
    )
  }
  if (length(pf_sparse) != nvars) {
    cli::cli_abort(
      "The length of `pf_sparse` must be equal to the number of predictors {nvars}."
    )
  }

  pf_sparse <- pf_sparse / sum(pf_sparse) * nvars
  maxit <- as.integer(maxit)
  pf_group <- as.double(pf_group)
  pf_sparse <- as.double(pf_sparse)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)

  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) abort("`lambda.factor` must be less than 1.")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin = 1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) abort("`lambda` must be non-negative.")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)

  ### check on upper/lower bounds
  if (any(lower_bnd > 0)) abort("`lower_bnd` must be non-positive")
  if (any(upper_bnd < 0)) abort("`upper_bnd` must be non-negative")
  lower_bnd[lower_bnd == -Inf] <- -9.9e30
  upper_bnd[upper_bnd == Inf] <- 9.9e30
  if (length(lower_bnd) < bn) {
    if (length(lower_bnd) == 1) lower_bnd <- rep(lower_bnd, bn)
    else cli::cli_abort("`lower_bnd` must be length 1 or length {bn}.")
  } else {
    lower_bnd <- lower_bnd[seq_len(bn)]
  }
  if (length(upper_bnd) < bn) {
    if (length(upper_bnd) == 1) upper_bnd <- rep(upper_bnd, bn)
    else cli::cli_abort("`upper_bnd` must be length 1 or length {bn}.")
  } else {
    upper_bnd <- upper_bnd[seq_len(bn)]
  }
  storage.mode(upper_bnd) <- "double"
  storage.mode(lower_bnd) <- "double"

  # call R sub-function
  if (is.character(family)) {
    family <- match.arg(family)
    fit <- switch(
      family,
      gaussian = sgl_ls(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        as.double(asparse), standardize, lower_bnd, upper_bnd),
      binomial = sgl_logit(
        bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
        as.double(asparse), standardize, lower_bnd, upper_bnd)
    )
  } else {
    fit <- sgl_irwls(
      bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
      dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
      as.double(asparse), standardize, lower_bnd, upper_bnd, weights,
      offset, family, trace_it, warm
    )
  }

  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  fit$asparse <- asparse
  fit$nobs <- nobs
  fit$pf_group <- pf_group
  fit$pf_sparse <- pf_sparse
  class(fit) <- c(class(fit), "sparsegl")
  fit
}

