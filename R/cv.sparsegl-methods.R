#' Get coefficients or make coefficient predictions from a `cv.sparsegl` object.
#'
#' This function gets coefficients or makes coefficient predictions from a
#' cross-validated \code{sparsegl} model, using the stored \code{"sparsegl.fit"}
#' object, and the optimal value chosen for \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to get
#' coefficients or make coefficient predictions.
#'
#' @param object fitted \code{\link{cv.sparsegl}} object.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the value \code{s="lambda.1se"} stored
#' on the CV \code{object}, it is the largest value of \code{lambda} such that
#' error is within 1 standard error of the minimum. Alternatively
#' \code{s="lambda.min"} can be used, it is the optimal value of \code{lambda}
#' that gives minimum cross validation error \code{cvm}. If \code{s} is
#' numeric, it is taken as the value(s) of \code{lambda} to be used.
#' @param \dots not used. Other arguments to predict.
#' @return The coefficients at the requested values for \code{lambda}.
#' @seealso \code{\link{cv.sparsegl}}, and \code{\link{predict.cv.sparsegl}}
#' methods.
#' @method coef cv.sparsegl
#' @export
coef.cv.sparsegl <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
    assertthat::assert_that(
        is.numeric(s) || is.character(s),
        msg = "Invalid form for s.")
    if (is.numeric(s)) {
        lambda <- s
    } else {
        s <- match.arg(s)
        lambda <- object[[s]]
    }
    coef(object$sparsegl.fit, s = lambda, ...)
}



#' Make predictions from a `cv.sparsegl` object.
#'
#' This function makes predictions from a cross-validated \code{sparsegl} model,
#' using the stored \code{"sparsegl.fit"} object, and the optimal value chosen
#' for \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object fitted \code{\link{cv.sparsegl}} object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix. See documentation for \code{predict.sparsegl}.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the value \code{s="lambda.1se"} stored
#' on the CV object. Alternatively \code{s="lambda.min"} can be used. If
#' \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
#' used.
#' @param \dots not used. Other arguments to predict.
#' @return The returned object depends on the \dots{} argument which is passed
#' on to the \code{\link{predict}} method for \code{\link{sparsegl}} objects.
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{cv.sparsegl}}, and \code{\link{coef.cv.sparsegl}}
#' methods.
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression

#' @method predict cv.sparsegl
#' @export
predict.cv.sparsegl <- function(object, newx,
                                s = c("lambda.1se", "lambda.min"), ...) {
    assertthat::assert_that(
        is.numeric(s) || is.character(s),
        msg = "Invalid form for s.")
    if (is.numeric(s)) {
        lambda <- s
    } else {
        s <- match.arg(s)
        lambda <- object[[s]]
    }
    predict(object$sparsegl.fit, newx, s = lambda, ...)
}
