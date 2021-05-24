#' Cross-validation for sparsegl
#'
#' Does k-fold cross-validation for sparsegl, produces a plot, and returns a
#' value for \code{lambda}. This function is modified based on the \code{cv}
#' function from the \code{glmnet} package.
#'
#' The function runs \code{\link{sparsegl}} \code{nfolds}+1 times; the first to
#' get the \code{lambda} sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard deviation
#' over the folds are computed.
#'
#' @aliases cv.sparsegl cv.ls
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#' is an observation vector.
#' @param y response variable. This argument should be quantitative for
#' regression (least squares).
#' @param group a vector of consecutive integers describing the grouping of the
#' coefficients (see example below).
#' @param lambda optional user-supplied lambda sequence; default is
#' \code{NULL}, and \code{\link{sparsegl}} chooses its own sequence.
#' @param pred.loss loss to use for cross-validation error. Valid options are:
#' \itemize{ \item \code{"L1"} for regression, mean square
#' error used by least squares regression \code{loss="ls"}, it measure the
#' deviation from the fitted mean to the response.  \item \code{"L2"} for
#' regression, mean absolute error used by least squares regression
#' \code{loss="ls"}, it measure the deviation from the fitted mean to the
#' response.  }. Default is \code{"L2"}.
#' @param nfolds number of folds - default is 5. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}.
#' @param foldid an optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param \dots other arguments that can be passed to sparsegl.
#' @return an object of class \code{\link{cv.sparsegl}} is returned, which is a
#' list with the ingredients of the cross-validation fit.  \item{lambda}{the
#' values of \code{lambda} used in the fits.} \item{cvm}{the mean
#' cross-validated error - a vector of length \code{length(lambda)}.}
#' \item{cvsd}{estimate of standard error of \code{cvm}.} \item{cvupper}{upper
#' curve = \code{cvm+cvsd}.} \item{cvlower}{lower curve = \code{cvm-cvsd}.}
#' \item{name}{a text string indicating type of measure (for plotting
#' purposes).} \item{sparsegl.fit}{a fitted \code{\link{sparsegl}} object for the
#' full data.} \item{lambda.min}{The optimal value of \code{lambda} that gives
#' minimum cross validation error \code{cvm}.} \item{lambda.1se}{The largest
#' value of \code{lambda} such that error is within 1 standard error of the
#' minimum.}
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{sparsegl}}, \code{\link{plot.cv.sparsegl}},
#' \code{\link{predict.cv.sparsegl}}, and \code{\link{coef.cv.sparsegl}} methods.
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression

#' @export
cv.sparsegl <- function(x, y, group, lambda = NULL,
                        pred.loss = c("L2", "L1"),
                        nfolds = 5, foldid, ...) {
    if  (missing(pred.loss) || !match(pred.loss, c("L2", "L1"), FALSE)) {
        warning("Only 'L2' and 'L1'  available for least squares models; 'L2' used")
        pred.loss <- "L2"
    } 
    pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    sparsegl.object <- sparsegl(x, y, group, lambda = lambda, ...)
    lambda <- sparsegl.object$lambda
    # predict -> coef
    if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = N))
    else nfolds <- max(foldid)
    assertthat::assert_that(
        nfolds > 1, msg = "nfolds must be at least 2; nfolds = 10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- sparsegl(x = x[!which, , drop = FALSE], y = y_sub,
                                 group = group, lambda = lambda,  ...)
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


cv.ls <- function(outlist, lambda, x, y, foldid, pred.loss = c("L2","L1")) {
    typenames <- c(L2 = "Least-Squares loss", L1 = "Absolute loss")
    pred.loss <- match.arg(pred.loss)
    predmat <- matrix(NA, length(y), length(lambda))
    nfolds <- max(foldid)
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE])
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, L2 = (y - predmat)^2, L1 = abs(y - predmat))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    scaled <- scale(cvraw, cvm, FALSE)^2
    cvsd <- sqrt(apply(scaled, 2, mean, na.rm = TRUE) / (N - 1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


