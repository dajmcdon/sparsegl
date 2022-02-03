#' Cross-validation for a `sparsegl` object.
#'
#' Does k-fold cross-validation for [sparsegl()].
#' This function is largely similar [glmnet::cv.glmnet()].
#'
#' The function runs [sparsegl()] `nfolds + 1` times; the first to
#' get the `lambda` sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard deviation
#' over the folds are computed.
#'
#' @aliases cv.sparsegl cv.ls
#' @template param_x-template
#' @template param_y-template
#' @template param_group-template
#' @template param_family-template
#' @template param_lambda-template
#' @param pred.loss Loss to use for cross-validation error. Valid options are:
#'  * `"L2"` for regression, mean square error
#'  * `"L1"` for regression, mean absolute error
#'  * `"binomial"` for classification, binomial deviance loss
#'  * `"misclass"` for classification, misclassification error.
#' @param nfolds Number of folds - default is 10. Although `nfolds` can be
#'   as large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. Smallest value allowable is `nfolds = 3`.
#' @param foldid An optional vector of values between 1 and `nfolds`
#'   identifying which fold each observation is in. If supplied, `nfolds` can
#'   be missing.
#' @param ... Other arguments that can be passed to sparsegl.
#'
#' @return An object of class [cv.sparsegl()] is returned, which is a
#'   list with the ingredients of the cross-validation fit.
#'   \item{lambda}{The values of \code{lambda} used in the fits.}
#'   \item{cvm}{The mean cross-validated error - a vector of
#'     length \code{length(lambda)}.}
#'   \item{cvsd}{Estimate of standard error of \code{cvm}.}
#'   \item{cvupper}{Upper curve = \code{cvm + cvsd}.}
#'   \item{cvlower}{Lower curve = \code{cvm - cvsd}.}
#'   \item{name}{A text string indicating type of measure (for plotting
#'     purposes).}
#'   \item{sparsegl.fit}{A fitted [sparsegl()] object for the full data.}
#'   \item{lambda.min}{The optimal value of \code{lambda} that gives
#'     minimum cross validation error \code{cvm}.}
#'   \item{lambda.1se}{The largest value of \code{lambda} such that error
#'     is within 1 standard error of the minimum.}
#'
#'
#' @seealso [sparsegl()], [plot.cv.sparsegl()],
#' [predict.cv.sparsegl()], and [coef.cv.sparsegl()] methods.
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
cv.sparsegl <- function(x, y, group = NULL, family = c("gaussian", "binomial"),
                        lambda = NULL,
                        pred.loss = c("L2", "L1", "binomial", "misclass"),
                        nfolds = 10, foldid = NULL, ...) {
    family <- match.arg(family)
    pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    sparsegl.object <- sparsegl(x, y, group, lambda = lambda, family = family,
                                ...)
    lambda <- sparsegl.object$lambda
    # predict -> coef
    if (is.null(foldid)) foldid <- sample(rep(seq(nfolds), length = N))
    else nfolds <- max(foldid)
    assertthat::assert_that(
        nfolds > 1, msg = "nfolds must be at least 2; nfolds = 10 recommended")
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
    fun <- paste("cv", class(sparsegl.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
                cvlo = cvm - cvsd, name = cvname,
                sparsegl.fit = sparsegl.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.sparsegl"
    obj
}


cv.ls <- function(outlist, lambda, x, y, foldid,
                        pred.loss = c("L2","L1")) {
    typenames <- c(L2 = "Least-Squares loss", L1 = "Absolute loss")
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
    cvraw <- switch(pred.loss, L2 = (y - predmat)^2, L1 = abs(y - predmat))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    scaled <- scale(cvraw, cvm, FALSE)^2
    cvsd <- sqrt(apply(scaled, 2, mean, na.rm = TRUE) / (N - 1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

cv.logit <- function(outlist, lambda, x, y, foldid,
                        pred.loss = c("binomial", "misclass")) {
    typenames <- c(binomial = "Binomial Deviance Loss",
                   misclass = "Misclassification Error")
    pred.loss <- match.arg(pred.loss)
    prob_min <- 1e-05
    fmax <- log(1/prob_min - 1)
    fmin <- -fmax
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        test_fold <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[test_fold, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- switch(pred.loss, binomial = log(1 + exp(-y * predmat)),
                    misclass = (y != ifelse(predmat > 0, 1, -1)))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, NA.RM = TRUE) /
                     (N - 1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
