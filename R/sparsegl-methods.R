#' Extract model coefficients from a `sparsegl` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' [sparsegl()] object.
#'
#' `s` is the new vector at which predictions are requested. If `s`
#' is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#' @param object Fitted [sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'  coefficients are required. Default is the entire sequence.
#' @param ... Not used.
#' @seealso [sparsegl()], [predict.sparsegl()] and
#' [print.sparsegl()] methods.
#'
#' @return The coefficients at the requested values for `lambda`.
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
  b0 <- matrix(object$b0, nrow = 1)
  # if conflicts happens and throw an error here, remove t() outside as.matrix()
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    ls <- length(s)
    if (ls == 1) {
      nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(ls, lamlist$frac) +
        nbeta[, lamlist$right, drop = FALSE] %*% Matrix::Diagonal(ls, 1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  return(nbeta)
}




#' Make predictions from a `sparsegl` object.
#'
#' Similar to other predict methods, this function produces fitted values and
#' class labels from a fitted [`sparsegl`] object.
#'
#' `s` is new vector at which predictions are requested. If `s`
#' is not in the lambda sequence used for fitting the model, the
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of predicted values from both left and
#' right `lambda` indices.
#'
#' @param object Fitted [sparsegl()] model object.
#' @param newx Matrix of new values for `x` at which predictions are to be
#'   made. Must be a matrix.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param type Type of prediction required. Type `"link"` gives the linear
#'   predictors for `"binomial"`; for `"gaussian"` models it gives the fitted
#'   values. Type `"response"` gives the fitted probabilities for
#'   `"binomial"`; for `"gaussian"` type
#'   `"response"` is equivalent to type `"link"`. Type
#'   `"coefficients"` computes the coefficients at the requested values for
#'   `s`.  Note that for `"binomial"` models, results are returned only
#'   for the class corresponding to the second level of the factor response.
#'   Type `"class"` applies only to `"binomial"` models, and produces the
#'   class label corresponding to
#'   the maximum probability. Type `"nonzero"` returns a list of the indices
#'   of the nonzero coefficients for each value of \code{s}.
#'
#' @param ... Not used.
#' @return The object returned depends on type.
#' @seealso [sparsegl()], [coef.sparsegl()] and
#' [print.sparsegl()] methods.
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
predict.sparsegl <- function(object, newx, s = NULL,
                             type = c("link","response","coefficients","nonzero","class"),
                             ...) {
  type <- match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'newx'")
  }
  if (class(object)[2] == "ls" && type == "class")
    stop("No class predictions are available for regression.")
  nbeta <- coef(object, s)
  if (type == "coefficients") return(nbeta)
  if (type == "nonzero") return(nonzeroCoef(nbeta[-1, ,drop = FALSE]))
  if (inherits(newx, "sparseMatrix")) newx <- as(newx,"dgCMatrix")
  dx <- dim(newx)
  p <- object$dim[1]
  if (is.null(dx)) newx <- matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p)
    stop(paste0("The number of variables in newx must be ", p))
  fit <- as.matrix(cbind2(1, newx) %*% nbeta)
  if (type == "link") return(fit)
  if (type == "response" && class(object)[2] == "ls") return(fit)
  if (type == "response" && class(object)[2] == "logit")
    return(1 / (1 + exp(-fit)))
  if (type == "class") {
    fit <- ifelse(fit > 0, 2, 1)
    fit <- object$classnames[fit]
    return(fit)
  }
}




#' Print a `sparsegl` object.
#'
#' Prints some summary information about the fitted [sparsegl()] object.
#'
#'
#' @param x Fitted [sparsegl()] object.
#' @param digits Significant digits in printout.
#' @param ... not used
#'
#' @seealso [sparsegl()], [coef.sparsegl()] and
#' [predict.sparsegl()] methods.
#'
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


