#' get coefficients or make coefficient predictions from an "sparsegl" object.
#'
#' Computes the coefficients at the requested values for \code{lambda} from a
#' fitted \code{\link{sparsegl}} object.
#'
#' \code{s} is the new vector at which predictions are requested. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{coef}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' \code{lambda} indices.
#'
#' @param object fitted \code{\link{sparsegl}} model object.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the entire sequence used to create the
#' model.
#' @param \dots not used. Other arguments to predict.
#' @return The coefficients at the requested values for \code{lambda}.
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{predict.sparsegl}} method
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @examples
#'
#' # load sparsegl library
#' library(sparsegl)
#'
#' # load data set
#' data(colon)
#'
#' # define group index
#' group <- rep(1:20, each=5)
#'
#' # fit group lasso
#' m1 <- sparsegl(x=colon$x, y=colon$y, group=group)
#'
#' # the coefficients at lambda = 0.01 and 0.02
#' coef(m1, s=c(0.01, 0.02))
#'
#' @export
#' @method coef sparsegl
coef.sparsegl <- function(object, s = NULL, ...) {
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
        nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  return(nbeta)
}




#' Make predictions from a "sparsegl" object.
#'
#' Similar to other predict methods, this functions predicts fitted values and
#' class labels from a fitted [`sparsegl()`] object.
#'
#' `s` is the new vector at which predictions are requested. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{predict}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of predicted values from both left and
#' right \code{lambda} indices.
#'
#' @param object fitted \code{\link{sparsegl}} model object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the entire sequence used to create the
#' model.
#' @param type type of prediction required: \itemize{ \item Type \code{"link"},
#' for regression it returns the fitted response; for classification it gives
#' the linear predictors.  \item Type \code{"class"}, only valid for
#' classification, it produces the predicted class label corresponding to the
#' maximum probability.}
#'
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{\link{coef}} method
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression

#'
#' @method predict sparsegl
#' @export
predict.sparsegl <- function(object, newx, s = NULL, ...) {
  nbeta <- coef(object, s)
  if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
  nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  return(nfit)
}




#' Print a sparsegl object
#'
#' Print the nonzero group counts at each lambda along the sparsegl path.
#'
#' Print the information about the nonzero group counts at each lambda step in
#' the \code{\link{sparsegl}} object. The result is a two-column matrix with
#' columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of
#' the groups that have nonzero within-group coefficients, the \code{Lambda}
#' column is the the corresponding lambda.
#'
#' @param x fitted \code{\link{sparsegl}} object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @return a two-column matrix, the first columns is the number of nonzero
#' group counts and the second column is \code{Lambda}.
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @examples

#' @method print sparsegl
#' @export
print.sparsegl <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}


