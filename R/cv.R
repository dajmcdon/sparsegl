#' Cross-validation for gglasso
#' 
#' Does k-fold cross-validation for gglasso, produces a plot, and returns a
#' value for \code{lambda}. This function is modified based on the \code{cv}
#' function from the \code{glmnet} package.
#' 
#' The function runs \code{\link{gglasso}} \code{nfolds}+1 times; the first to
#' get the \code{lambda} sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard deviation
#' over the folds are computed.
#' 
#' @aliases cv.gglasso cv.ls cv.logit cv.hsvm cv.sqsvm
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#' is an observation vector.
#' @param y response variable. This argument should be quantitative for
#' regression (least squares), and a two-level factor for classification
#' (logistic model, huberized SVM, squared SVM).
#' @param group a vector of consecutive integers describing the grouping of the
#' coefficients (see example below).
#' @param lambda optional user-supplied lambda sequence; default is
#' \code{NULL}, and \code{\link{gglasso}} chooses its own sequence.
#' @param pred.loss loss to use for cross-validation error. Valid options are:
#' \itemize{ \item \code{"loss"} for classification, margin based loss
#' function.  \item \code{"misclass"} for classification, it gives
#' misclassification error.  \item \code{"L1"} for regression, mean square
#' error used by least squares regression \code{loss="ls"}, it measure the
#' deviation from the fitted mean to the response.  \item \code{"L2"} for
#' regression, mean absolute error used by least squares regression
#' \code{loss="ls"}, it measure the deviation from the fitted mean to the
#' response.  } Default is \code{"loss"}.
#' @param nfolds number of folds - default is 5. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}.
#' @param foldid an optional vector of values between 1 and \code{nfold}
#' identifying what fold each observation is in. If supplied, \code{nfold} can
#' be missing.
#' @param delta parameter \eqn{\delta}{delta} only used in huberized SVM for
#' computing log-likelihood on validation set, only available with
#' \code{pred.loss = "loss"}, \code{loss = "hsvm"}.
#' @param \dots other arguments that can be passed to gglasso.
#' @return an object of class \code{\link{cv.gglasso}} is returned, which is a
#' list with the ingredients of the cross-validation fit.  \item{lambda}{the
#' values of \code{lambda} used in the fits.} \item{cvm}{the mean
#' cross-validated error - a vector of length \code{length(lambda)}.}
#' \item{cvsd}{estimate of standard error of \code{cvm}.} \item{cvupper}{upper
#' curve = \code{cvm+cvsd}.} \item{cvlower}{lower curve = \code{cvm-cvsd}.}
#' \item{name}{a text string indicating type of measure (for plotting
#' purposes).} \item{gglasso.fit}{a fitted \code{\link{gglasso}} object for the
#' full data.} \item{lambda.min}{The optimal value of \code{lambda} that gives
#' minimum cross validation error \code{cvm}.} \item{lambda.1se}{The largest
#' value of \code{lambda} such that error is within 1 standard error of the
#' minimum.}
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{gglasso}}, \code{\link{plot.cv.gglasso}},
#' \code{\link{predict.cv.gglasso}}, and \code{\link{coef.cv.gglasso}} methods.
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load data set
#' data(bardet)
#' 
#' # define group index
#' group <- rep(1:20,each=5)
#' 
#' # 5-fold cross validation using group lasso 
#' # penalized logisitic regression
#' cv <- cv.gglasso(x=bardet$x, y=bardet$y, group=group, loss="ls",
#' pred.loss="L2", lambda.factor=0.05, nfolds=5)
#' 
#' @export
cv.gglasso <- function(x, y, group, lambda = NULL, pred.loss = c("misclass", 
    "loss", "L1", "L2"), nfolds = 5, foldid, delta, ...) {
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    if (missing(delta)) 
        delta <- 1
    if (delta < 0) 
        stop("delta must be non-negtive")
    gglasso.object <- gglasso(x, y, group, lambda = lambda, delta = delta, ...)
    lambda <- gglasso.object$lambda
    # predict -> coef
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- gglasso(x = x[!which, , drop = FALSE], y = y_sub, group = group, 
            lambda = lambda, delta = delta, ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(gglasso.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
        cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    obj
}


#' @export
cv.hsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for HHSVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * hubercls(y * predmat, delta), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}



#' @export
cv.logit <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for logistic regression; 'misclass' used")
        pred.loss <- "misclass"
    }
    prob_min <- 1e-05
    fmax <- log(1/prob_min - 1)
    fmin <- -fmax
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- switch(pred.loss, loss = 2 * log(1 + exp(-y * predmat)), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cv.sqsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for squared SVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * ifelse((1 - y * predmat) <= 0, 0, 
        (1 - y * predmat))^2, misclass = (y != ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}

#' @export
cv.ls <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(L2 = "Least-Squared loss", L1 = "Absolute loss")
    if (pred.loss == "default") 
        pred.loss <- "L2"
    if (!match(pred.loss, c("L2", "L1"), FALSE)) {
        warning("Only 'L2' and 'L1'  available for least squares models; 'L2' used")
        pred.loss <- "L2"
    }
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
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
