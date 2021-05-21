######################################################################
## This function is adapted/modified based on the plot.cv function from
## the glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.


utils::globalVariables(c("X", "y", "lower", "upper"))

#' plot the cross-validation curve produced by cv.sparsegl
#' 
#' Plots the cross-validation curve, and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. This function is
#' modified based on the \code{plot.cv} function from the \code{glmnet}
#' package.
#' 
#' A plot is produced.
#' 
#' @param x fitted \code{\link{cv.sparsegl}} object
#' @param sign.lambda either plot against \code{log(lambda)} (default) or its
#' negative if \code{sign.lambda=-1}.
#' @param \dots other graphical parameters to plot
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{cv.sparsegl}}.
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' 
#' Friedman, J., Hastie, T., and Tibshirani, R. (2010), ``Regularization paths
#' for generalized linear models via coordinate descent,'' \emph{Journal of
#' Statistical Software}, 33, 1.\cr \url{http://www.jstatsoft.org/v33/i01/}
#' @keywords models regression
#' @method plot cv.sparsegl
#' @export
plot.cv.sparsegl <- function(x, sign.lambda = 1, ...) {
    cvobj <- x
    dat <- data.frame("X" = sign.lambda * log(cvobj$lambda), 
                      "y" = cvobj$cvm, 
                      "upper" = cvobj$cvupper,
                      "lower" = cvobj$cvlo)
    dat %>% 
        ggplot2::ggplot(ggplot2::aes(x = X, y = y)) +
        ggplot2::geom_point(color = 'red') + 
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.1, color = 'darkgrey') + 
        ggplot2::xlab("Log Lambda") + 
        ggplot2::ylab(cvobj$name)
} 
