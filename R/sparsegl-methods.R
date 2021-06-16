#' Extract model coefficients from a `sparsegl` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' \code{\link{sparsegl}} object.
#'
#' `s` is the new vector at which predictions are requested. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{coef}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' \code{lambda} indices.
#'
#' @param object Fitted \code{\link{sparsegl}} model object.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the entire sequence used to create the
#' model.
#' @param \dots Not used.
#' @seealso \code{\link{sparsegl}}, \code{\link{predict.sparsegl}} and 
#' \code{\link{print.sparsegl}} methods.
#' @return The coefficients at the requested values for \code{lambda}.
#' 
#' @method coef sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' coef(fit1, s = c(0.02, 0.03))
coef.sparsegl <- function(object, s = NULL, ...) {
  b0 <- as.matrix(object$b0)
  # if conflicts happens and throw an error here, remove t() outside as.matrix()
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*% diag(lamlist$frac) +
        nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  return(nbeta)
}




#' Make predictions from a `sparsegl` object.
#'
#' Similar to other predict methods, this functions predicts fitted values and
#' class labels from a fitted [`sparsegl`] object.
#'
#' `s` is the new vector at which predictions are requested. If \code{s}
#' is not in the lambda sequence used for fitting the model, the \code{predict}
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of predicted values from both left and
#' right \code{lambda} indices.
#'
#' @param object Fitted \code{\link{sparsegl}} model object.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the entire sequence used to create the
#' model.
#'
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on type.
#' @seealso \code{\link{sparsegl}}, \code{\link{coef.sparsegl}} and 
#' \code{\link{print.sparsegl}} methods.
#'
#' @method predict sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' predict(fit1, newx = X[10, ], s = fit1$lambda[3:5])
predict.sparsegl <- function(object, newx, s = NULL, ...) {
  nbeta <- coef(object, s)
  if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
  fit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  return(fit)
}




#' Print a `sparsegl` object
#'
#' Prints a few summaries of the fitted \code{\link{sparsegl}} model object.
#'
#'
#' @param x Fitted \code{\link{sparsegl}} object.
#' @param digits Significant digits in printout.
#' @param \dots Additional print arguments.
#' @seealso \code{\link{sparsegl}}, \code{\link{coef.sparsegl}} and 
#' \code{\link{predict.sparsegl}} methods.
#' @method print sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' print(fit1)
print.sparsegl <- function(x, digits = min(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Approx. degrees of freedom: ", round(min(x$df), digits),
      " - ", round(max(x$df), digits), "\n")
  cat("Range of lambda: ", round(max(x$lambda), digits),
      " - ", round(min(x$lambda), digits), "\n")
  nlams <- length(x$lambda)
  cat("Saturated penalty: ",
      round(sp_group_norm(x$beta[,nlams], x$group, x$asparse), digits))
}


