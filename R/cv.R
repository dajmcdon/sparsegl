#' Cross-validation for a `sparsegl` object.
#'
#' Performs k-fold cross-validation for [sparsegl()].
#' This function is largely similar [glmnet::cv.glmnet()].
#'
#' The function runs [sparsegl()] `nfolds + 1` times; the first to
#' get the `lambda` sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard error
#' over the folds are computed.
#'
#' @aliases cv.sparsegl
#' @inheritParams sparsegl
#' @param pred.loss Loss to use for cross-validation error. Valid options are:
#'  * `"default"` the same as deviance (mse for regression and deviance otherwise)
#'  * `"mse"` mean square error
#'  * `"deviance"` the default (mse for Gaussian regression, and negative
#'    log-likelihood otherwise)
#'  * `"mae"` mean absolute error, can apply to any family
#'  * `"misclass"` for classification only, misclassification error.
#' @param nfolds Number of folds - default is 10. Although `nfolds` can be
#'   as large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. Smallest value allowable is `nfolds = 3`.
#' @param foldid An optional vector of values between 1 and `nfolds`
#'   identifying which fold each observation is in. If supplied, `nfolds` can
#'   be missing.
#' @param ... Additional arguments to [sparsegl()].
#'
#' @return An object of class [cv.sparsegl()] is returned, which is a
#'   list with the components describing the cross-validation error.
#'   \item{lambda}{The values of \code{lambda} used in the fits.}
#'   \item{cvm}{The mean cross-validated error - a vector of
#'     length \code{length(lambda)}.}
#'   \item{cvsd}{Estimate of standard error of \code{cvm}.}
#'   \item{cvupper}{Upper curve = \code{cvm + cvsd}.}
#'   \item{cvlower}{Lower curve = \code{cvm - cvsd}.}
#'   \item{name}{A text string indicating type of measure (for plotting
#'     purposes).}
#'   \item{nnzero}{The number of non-zero coefficients for each \code{lambda}}
#'   \item{active_grps}{The number of active groups for each \code{lambda}}
#'   \item{sparsegl.fit}{A fitted [sparsegl()] object for the full data.}
#'   \item{lambda.min}{The optimal value of \code{lambda} that gives
#'     minimum cross validation error \code{cvm}.}
#'   \item{lambda.1se}{The largest value of \code{lambda} such that error
#'     is within 1 standard error of the minimum.}
#'   \item{call}{The function call.}
#'
#'
#' @seealso [sparsegl()], as well as [`plot()`][plot.cv.sparsegl()],
#'   [`predict()`][predict.cv.sparsegl()], and [`coef()`][coef.cv.sparsegl()]
#'   methods for `"cv.sparsegl"` objects.
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
#' cv_fit <- cv.sparsegl(X, y, groups)
#'
cv.sparsegl <- function(
    x, y, group = NULL, family = c("gaussian", "binomial"),
    lambda = NULL,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass"),
    nfolds = 10, foldid = NULL, ...) {

  if (is.character(family)) family <- match.arg(family)
  else validate_family(family)

  pred.loss <- match.arg(pred.loss)
  if (pred.loss == "misclass" && !(is.character(family) && family == "binomial"))
    cli::cli_abort('`pred.loss = "misclass"` only works if `family = "binomial"`.')
  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  sparsegl.object <- sparsegl(x, y, group, lambda = lambda, family = family,
                              ...)
  lambda <- sparsegl.object$lambda
  # predict -> coef
  if (is.null(foldid)) foldid <- sample(rep(seq(nfolds), length = N))
  else nfolds <- max(foldid)
  if (nfolds < 2) {
    cli::cli_abort("`nfolds` must be at least 2; `nfolds = 10` is recommended.")
  }
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    outlist[[i]] <- sparsegl(
      x = x[!test_fold, , drop = FALSE],
      y = y[!test_fold], group = group, lambda = lambda, family = family,
      ...)
  }
  ###What to do depends on the pred.loss and the model fit
  cvstuff <- cverror(sparsegl.object, outlist, lambda, x, y, foldid, pred.loss)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  nz <- predict(sparsegl.object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(group[x])))
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
              cvlo = cvm - cvsd, name = cvname,
              nnzero = nnzero, active_grps = active_grps,
              sparsegl.fit = sparsegl.object,
              call = match.call())
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sparsegl"
  obj
}



cverror <- function(fullfit, outlist, lambda, x, y, foldid, pred.loss) {
  UseMethod("cverror")
}

#' @export
cverror.lsspgl <- function(fullfit, outlist, lambda, x, y, foldid,
                           pred.loss = c("default", "mse", "deviance", "mae")) {
  typenames <- c(default = "Mean squared error", mse = "Mean squared error",
                 deviance = "Mean squared error", mae = "Mean absolute error")
  pred.loss <- match.arg(pred.loss)
  predmat <- matrix(NA, length(y), length(lambda))
  nfolds <- max(foldid)
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- switch(pred.loss, mae = abs(y - predmat), (y - predmat)^2)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(apply(scaled, 2, mean, na.rm = TRUE) / (N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cverror.logitspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass")) {
  typenames <- c(default = "Binomial deviance", mse = "Mean squared error",
                 deviance = "Binomial deviance", mae = "Mean absolute error",
                 misclass = "Missclassification error")
  pred.loss <- match.arg(pred.loss)
  prob_min <- 1e-05
  fmax <- log(1 / prob_min - 1)
  fmin <- -fmax
  y <- as.factor(y)
  y <- as.numeric(y) - 1 # 0 / 1 valued
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  predmat <- pmin(pmax(predmat, fmin), fmax)
  binom_deviance <- function(m) stats::binomial()$dev.resids(y, m, 1)
  cvraw <- switch(
    pred.loss,
    mse = (y - predmat)^2,
    mae = abs(y - predmat),
    misclass = y != ifelse(predmat > 0.5, 1, 0),
    apply(predmat, 2, binom_deviance)
  )
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, NA.RM = TRUE) /
                 (N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cverror.irlsspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae")) {
  typenames <- c(default = "Deviance", mse = "Mean squared error",
                 deviance = "Deviance", mae = "Mean absolute error")
  pred.loss <- match.arg(pred.loss)

  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }

  dev_fun <- function(m) fullfit$family$dev.resids(y, m, 1)
  cvraw <- switch(
    pred.loss,
    mse = (y - predmat)^2,
    mae = abs(y - predmat),
    apply(predmat, 2, dev_fun)
  )
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, NA.RM = TRUE) /
                 (N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
