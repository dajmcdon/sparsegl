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
#'  * `"auc"` for classification only, area under the ROC curve
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
#'   \item{i.min}{The index of \code{lambda.min} in the \code{lambda} sequence.}
#'   \item{i.1se}{The index of \code{lambda.1se} in the \code{lambda} sequence.}
#'   \item{call}{The function call.}
#'
#'
#' @seealso [sparsegl()], as well as [`plot()`][plot.cv.sparsegl()],
#'   [`predict()`][predict.cv.sparsegl()], and [`coef()`][coef.cv.sparsegl()]
#'   methods for `"cv.sparsegl"` objects.
#'
#' @references Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., and
#'   McDonald, D.J. 2024.
#'   \emph{sparsegl: An `R` Package for Estimating Sparse Group Lasso.}
#'   Journal of Statistical Software, Vol. 110(6): 1–23.
#'   \doi{10.18637/jss.v110.i06}.
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
    pred.loss = c("default", "mse", "deviance", "mae", "misclass", "auc"),
    nfolds = 10, foldid = NULL,
    weights = NULL, offset = NULL,
    ...) {
  fam <- validate_family(family)
  is_binomial <- FALSE
  if (fam$check == "char") {
    family <- match.arg(family)
    if (family == "binomial") is_binomial <- TRUE
  }
  if (fam$check == "fam" && fam$family$family == "binomial") is_binomial <- TRUE

  # not allowed for some families
  pred.loss <- match.arg(pred.loss)
  if (pred.loss == "misclass" || pred.loss == "auc") {
    if (!is_binomial) {
      cli_abort(c(
        "When `pred.loss` is {.val {pred.loss}}, `family` must be either:",
        `!` = "{.val {'binomial'}}, or {.fn stats::binomial}."
      ))
    }
  }

  N <- nrow(x)
  if (pred.loss == "auc" && (N / nfolds < 10)) {
    cli_warn(c(
      "For `pred.loss = 'auc'`, at least {.val {10}} observations are needed per fold.",
      `!` = "Here, only {.val {floor(N/nfolds)}} are avaliable.",
      `!` = "The results use `pred.loss = 'deviance'` instead.",
      `!` = "Alternatively, reduce `nfolds`."
    ))
  }
  ### Fit the model once to get dimensions etc of output
  y <- drop(y)
  sparsegl.object <- sparsegl(x, y, group,
    lambda = lambda, family = family,
    weights = weights, offset = offset, ...
  )
  lambda <- sparsegl.object$lambda
  # predict -> coef
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else {
    nfolds <- max(foldid)
  }
  if (nfolds < 2) {
    cli_abort(
      "`nfolds` must be at least {.val {2}}. `nfolds` = {.val {10}} is recommended."
    )
  }
  outlist <- as.list(seq(nfolds))
  ### Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    test_fold <- foldid == i
    outlist[[i]] <- sparsegl(
      x = x[!test_fold, , drop = FALSE],
      y = y[!test_fold], group = group, lambda = lambda, family = family,
      weights = weights[!test_fold], offset = offset[!test_fold], ...
    )
  }
  ### What to do depends on the pred.loss and the model fit
  cvstuff <- cverror(
    sparsegl.object, outlist, lambda, x, y, foldid,
    pred.loss, weights, is_binomial
  )
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  nz <- predict(sparsegl.object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(group[x])))
  out <- list(
    lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
    cvlo = cvm - cvsd, name = cvname,
    nnzero = nnzero, active_grps = active_grps,
    sparsegl.fit = sparsegl.object,
    call = match.call()
  )
  lamin <- getmin(lambda, cvm, cvsd, pred.loss)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sparsegl"
  obj
}

cvall <- function(fullfit, outlist, lambda, x, y, foldid, pred.loss,
                  dev_fun, weights = NULL, is_binomial = FALSE) {
  N <- length(y)
  nlambda <- length(lambda)
  nfolds <- max(foldid)
  predmat <- matrix(NA, N, nlambda)
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq_len(nlami)] <- preds
    nlams[i] <- nlami
  }
  if (is_binomial) {
    prob_min <- 1e-05
    fmax <- log(1 / prob_min - 1)
    fmin <- -fmax
    predmat <- pmin(pmax(predmat, fmin), fmax)
  }

  weights <- weights %||% rep(1, N)
  if (pred.loss == "auc") {
    err <- matrix(NA, nfolds, nlams)
    good <- matrix(0, nfolds, nlams)
    for (i in seq_len(nfolds)) {
      good[i, seq_len(nlams[i])] = 1
      which = foldid == i
      for (j in seq_len(nlams[i])) {
        err[i, j] = auc_mat(y[which], predmat[which, j],  weights[which])
      }
    }
    weights <- tapply(weights, foldid, sum)
    N <- colSums(good)
  } else {
    err <- switch(
      pred.loss,
      mse = (y - predmat)^2,
      mae = abs(y - predmat),
      misclass = y != ifelse(predmat > 0.5, 1, 0),
      deviance = apply(predmat, 2, dev_fun)
    )
    N <- length(y) - apply(is.na(predmat), 2, sum)
  }
  cvm <- apply(err, 2, stats::weighted.mean, na.rm = TRUE, w = weights)
  scaled <- scale(err, cvm, FALSE)^2
  cvsd <- sqrt(apply(
    scaled, 2, stats::weighted.mean, w = weights, na.rm = TRUE
  ) / (N - 1))
  list(cvm = cvm, cvsd = cvsd)
}

cverror <- function(fullfit, outlist, lambda, x, y, foldid, pred.loss, ...) {
  UseMethod("cverror")
}


#' @export
cverror.lsspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae"),
    ...) {
  typenames <- c(
    default = "Mean squared error", mse = "Mean squared error",
    deviance = "Mean squared error", mae = "Mean absolute error"
  )
  pred.loss <- match.arg(pred.loss)
  if (pred.loss != "mae") pred.loss <- "mse"
  dev_fun <- function(m) (y - m)^2
  c(
    cvall(fullfit, outlist, lambda, x, y, foldid, pred.loss, dev_fun,
          is_binomial = FALSE),
    name = typenames[pred.loss]
  )
}

#' @export
cverror.logitspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass", "auc"),
    is_binomial = TRUE,
    ...) {
  typenames <- c(
    default = "Binomial deviance", mse = "Mean squared error",
    deviance = "Binomial deviance", mae = "Mean absolute error",
    misclass = "Missclassification error", auc = "Area under the curve"
  )
  pred.loss <- match.arg(pred.loss)
  if (pred.loss == "default") pred.loss <- "deviance"
  dev_fun <- function(m) stats::binomial()$dev.resids(y, m, 1)
  y <- as.factor(y)
  y <- as.numeric(y) - 1 # 0 / 1 valued
  c(
    cvall(fullfit, outlist, lambda, x, y, foldid, pred.loss, dev_fun,
          is_binomial = TRUE),
    name = typenames[pred.loss]
  )
}

#' @export
cverror.irlsspgl <- function(
    fullfit, outlist, lambda, x, y, foldid,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass", "auc"),
    weights = NULL, is_binomial, ...) {
  typenames <- c(
    default = "Deviance", mse = "Mean squared error",
    deviance = "Deviance", mae = "Mean absolute error",
    misclass = "Missclassification error", auc = "Area under the curve"
  )
  pred.loss <- match.arg(pred.loss)
  if (pred.loss == "default") pred.loss <- "deviance"
  dev_fun <- function(m) fullfit$family$dev.resids(y, m, 1)
  weights <- weights %||% rep(1, length(y))
  c(
    cvall(fullfit, outlist, lambda, x, y, foldid, pred.loss, dev_fun,
          weights = weights, is_binomial = is_binomial),
    name = typenames[pred.loss]
  )
}
