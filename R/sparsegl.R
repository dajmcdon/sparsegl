#' Fits the regularization paths for group-lasso penalized learning problems
#'
#' Fits regularization paths for group-lasso penalized learning problems at a
#' sequence of regularization parameters lambda.
#'
#' Note that the objective function for \code{"ls"} least squares is
#' \deqn{RSS/(2*n) + lambda * penalty;} for \code{"hsvm"} Huberized squared
#' hinge loss, \code{"sqsvm"} Squared hinge loss and \code{"logit"} logistic
#' regression, the objective function is \deqn{-loglik/n + lambda * penalty.}
#' Users can also tweak the penalty by choosing different penalty factor.
#'
#' For computing speed reason, if models are not converging or running slow,
#' consider increasing \code{eps}, decreasing \code{nlambda}, or increasing
#' \code{lambda.factor} before increasing \code{maxit}.
#'
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#' is an observation vector.
#' @param y response variable. This argument should be quantitative for
#' regression (least squares), and a two-level factor for classification
#' (logistic model, huberized SVM, squared SVM).
#' @param group a vector of consecutive integers describing the grouping of the
#' coefficients (see example below).
#' @param pen a character string specifying the penalty function to use, valid
#' options are: \itemize{ \item \code{"sparsegl"} for group plus l1 penalty,
#' \item \code{"gglasso"} for group penalty only }Default is \code{"sparsegl"}.
#' @param algorithm a character string specifying the algorithm for sparse
#' group lasso (sparsegl). Valid options are: \itemize{ }
#' @param nlambda the number of \code{lambda} values - default is 100.
#' @param lambda.factor the factor for getting the minimal lambda in
#' \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} *
#' \code{max(lambda)}.  \code{max(lambda)} is the smallest value of
#' \code{lambda} for which all coefficients are zero. The default depends on
#' the relationship between \eqn{n} (the number of rows in the matrix of
#' predictors) and \eqn{p} (the number of predictors). If \eqn{n >= p}, the
#' default is \code{0.001}, close to zero.  If \eqn{n<p}, the default is
#' \code{0.05}.  A very small value of \code{lambda.factor} will lead to a
#' saturated fit. It takes no effect if there is user-defined \code{lambda}
#' sequence.
#' @param lambda a user supplied \code{lambda} sequence. Typically, by leaving
#' this option unspecified users can have the program compute its own
#' \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}.
#' Supplying a value of \code{lambda} overrides this. It is better to supply a
#' decreasing sequence of \code{lambda} values than a single (small) value, if
#' not, the program will sort user-defined \code{lambda} sequence in decreasing
#' order automatically.
#' @param pf penalty factor, a vector in length of bn (bn is the total number
#' of groups). Separate penalty weights can be applied to each group of
#' \eqn{\beta}{beta's}s to allow differential shrinkage. Can be 0 for some
#' groups, which implies no shrinkage, and results in that group always being
#' included in the model. Default value for each entry is the square-root of
#' the corresponding size of each group.
#' @param dfmax limit the maximum number of groups in the model. Useful for
#' very large \code{bs} (group size), if a partial path is desired. Default is
#' \code{bs+1}.
#' @param pmax limit the maximum number of groups ever to be nonzero. For
#' example once a group enters the model, no matter how many times it exits or
#' re-enters model through the path, it will be counted only once. Default is
#' \code{min(dfmax*1.2,bs)}.
#' @param eps convergence termination tolerance. Defaults value is \code{1e-8}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda
#' value. Default is 3e8. If models do not converge, consider increasing
#' \code{maxit}.
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#' @param asparse the weight to put on the ell1 norm in sparse group lasso. Default
#' is 0.05
#' @return An object with S3 class \code{\link{gglasso}}.  \item{call}{the call
#' that produced this object} \item{b0}{intercept sequence of length
#' \code{length(lambda)}} \item{beta}{a \code{p*length(lambda)} matrix of
#' coefficients.} \item{df}{the number of nonzero groups for each value of
#' \code{lambda}.} \item{dim}{dimension of coefficient matrix (ices)}
#' \item{lambda}{the actual sequence of \code{lambda} values used}
#' \item{npasses}{total number of iterations (the most inner loop) summed over
#' all lambda values} \item{jerr}{error flag, for warnings and errors, 0 if no
#' error.} \item{group}{a vector of consecutive integers describing the
#' grouping of the coefficients.}
#' @author Aaron Cohen and Dan Mcdonald\cr Maintainer: Aaron Cohen <cohenaa@indiana.edu>
#' @seealso \code{plot.gglasso}
#' @keywords models regression sparse
#' @examples
#' # load sparsegl library
#' library(sparsegl)
#'
#' # define group index
#' group1 <- rep(1:20,each=5)
#'
#' # fit group lasso penalized least squares
#' # m1 <- sparsegl(x=bardet$x,y=bardet$y,group=group1,pen="gglasso")
#' # define group index
#' group2 <- rep(1:20,each=5)
#'
#'
#' @export
sparsegl <- function(
  x, y, group = NULL,
  nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL, pf = sqrt(bs),
  intercept = TRUE, asparse = 0.05, standardize = TRUE,
  lower_bnd = -Inf, upper_bnd = Inf,
  dfmax = as.integer(max(group)) + 1L,
  pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08) {
  #################################################################################
  #\tDesign matrix setup, error checking
  this.call <- match.call()
  pen = match.arg(pen)
  algorithm = match.arg(algorithm)

  if (!is.matrix(x) && !inherits(x, "sparseMatrix"))
    stop("x has to be a matrix")

  if (any(is.na(x)))
    stop("Missing values in x not allowed!")

  y <- drop(y)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  assertthat::assert_that(length(y) == nobs,
                          msg = "x and y have different number of rows")
  assertthat::assert_that(is.numeric(y),
                          msg = "The response y must be numeric.")
  #################################################################################
  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    assertthat::assert_that(
      length(group) == nvars,
      msg = "group length does not match the number of predictors in x")
  }

  bn <- as.integer(max(group))
  bs <- as.integer(as.numeric(table(group)))

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn)))
    stop("Groups must be consecutively numbered 1,2,3,...")
  #Need to add if(loss=sparsegl)...
  if (pen == "sparsegl" && (asparse>1 || asparse<0)){
    asparse <- 0
    pen <- "gglasso"
    warning("asparse must be in [0,1], running ordinary group lasso.")
  }

  ix <- rep(NA, bn)
  iy <- rep(NA, bn)
  j <- 1
  for (g in seq_len(bn)) {
    ix[g] <- j
    iy[g] <- j + bs[g] - 1
    j <- j + bs[g]
  }
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)
  #################################################################################
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
  #################################################################################
  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    assertthat::assert_that(lambda.factor < 1,
                            msg = "lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin=1 if user define lambda
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


  #################################################################################
  # call R sub-function
  fit <- sgl(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam,
    eps, maxit, vnames, group, intr, as.double(asparse),
    standardize, algorithm, lower_bnd, upper_bnd)
  #################################################################################
  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  class(fit) <- c("sparsegl", class(fit))
  fit
}
