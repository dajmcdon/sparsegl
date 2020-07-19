######################################################################
## This function is adapted/modified based on the plot.cv function from
## the glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.




#' plot the cross-validation curve produced by cv.gglasso
#' 
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. This function is
#' modified based on the \code{plot.cv} function from the \code{glmnet}
#' package.
#' 
#' A plot is produced.
#' 
#' @param x fitted \code{\link{cv.gglasso}} object
#' @param sign.lambda either plot against \code{log(lambda)} (default) or its
#' negative if \code{sign.lambda=-1}.
#' @param \dots other graphical parameters to plot
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{cv.gglasso}}.
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' 
#' Friedman, J., Hastie, T., and Tibshirani, R. (2010), ``Regularization paths
#' for generalized linear models via coordinate descent,'' \emph{Journal of
#' Statistical Software}, 33, 1.\cr \url{http://www.jstatsoft.org/v33/i01/}
#' @keywords models regression
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load data set
#' data(colon)
#' 
#' # define group index
#' group <- rep(1:20,each=5)
#' 
#' # 5-fold cross validation using group lasso 
#' # penalized logisitic regression
#' cv <- cv.gglasso(x=colon$x, y=colon$y, group=group, loss="logit",
#' pred.loss="misclass", lambda.factor=0.05, nfolds=5)
#' 
#' # make a CV plot
#' plot(cv)
#' @method plot cv.gglasso
#' @export
plot.cv.gglasso <- function(x, sign.lambda = 1, ...) {
    cvobj <- x
    xlab <- "log(Lambda)"
    if (sign.lambda < 0) 
        xlab <- paste("-", xlab, sep = "")
    plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm, ylim = range(cvobj$cvupper, 
        cvobj$cvlo), xlab = xlab, ylab = cvobj$name, type = "n")
    new.args <- list(...)
    if (length(new.args)) 
        plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper, cvobj$cvlo, width = 0.01, 
        col = "darkgrey")
    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
} 
