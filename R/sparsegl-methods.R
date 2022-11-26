#' Extract model coefficients from a `sparsegl` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' [sparsegl()] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#' @param object Fitted [sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'  coefficients are required. Default is the entire sequence.
#' @param ... Not used.
#' @seealso [sparsegl()] and [predict.sparsegl()].
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
  rlang::check_dots_empty()
  b0 <- matrix(object$b0, nrow = 1)
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
        nbeta[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(ls, 1 - lamlist$frac)
    }
    namess <- names(s) %||% paste0("s", seq_along(s))
    dimnames(nbeta) <- list(vnames, namess)
  }
  return(nbeta)
}




#' Make predictions from a `sparsegl` object.
#'
#' Similar to other predict methods, this function produces fitted values and
#' class labels from a fitted [`sparsegl`] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#'
#' @param object Fitted [sparsegl()] model object.
#' @param newx Matrix of new values for `x` at which predictions are to be
#'   made. Must be a matrix. This argument is mandatory.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param type Type of prediction required. Type `"link"` gives the linear
#'   predictors for `"binomial"`; for `"gaussian"` models it gives the fitted
#'   values. Type `"response"` gives predictions on the scale of the response
#'   (for example, fitted probabilities for `"binomial"`); for `"gaussian"` type
#'   `"response"` is equivalent to type `"link"`. Type
#'   `"coefficients"` computes the coefficients at the requested values for
#'   `s`.
#'   Type `"class"` applies only to `"binomial"` models, and produces the
#'   class label corresponding to
#'   the maximum probability. Type `"nonzero"` returns a list of the indices
#'   of the nonzero coefficients for each value of \code{s}.
#'
#' @param ... Not used.
#' @return The object returned depends on type.
#' @seealso [sparsegl()], [coef.sparsegl()].
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
predict.sparsegl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero", "class"),
    ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      cli::cli_abort(
        "You must supply a value for `newx` when `type` == '{type}'."
      )
  }
  nbeta <- coef(object, s)
  if (type == "coefficients") return(nbeta)
  if (type == "nonzero") return(nonzeroCoef(nbeta[-1, , drop = FALSE]))
  if (inherits(newx, "sparseMatrix")) newx <- as_dgCMatrix(newx)
  dx <- dim(newx)
  p <- object$dim[1]
  if (is.null(dx)) newx <- matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p)
    cli::cli_abort("The number of variables in `newx` must be {p}.")
  fit <- as.matrix(cbind2(1, newx) %*% nbeta)
  fit
}


#' @export
predict.lsspgl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero"),
    ...) {
  type <- match.arg(type)
  NextMethod("predict")
}

#' @export
predict.logitspgl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero", "class"),
    ...) {
  type <- match.arg(type)
  nfit <- NextMethod("predict")
  switch(
    type,
    response = 1 / (1 + exp(-nfit)),
    class = object$classnames[ifelse(nfit > 0, 2, 1)],
    nfit
  )
}

#' @export
predict.irlsspgl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero"),
    ...) {
  type <- match.arg(type)
  nfit <- NextMethod("predict")
  if (type == "response") return(object$family$linkinv(nfit))
  else return(nfit)
}


#' @export
fitted.sparsegl <- function(object, ...) {
  cli::cli_abort(c(
    "!" = "Because design matrices are typically large, these are not stored ",
    "!" = "in the estimated `sparsegl` object. Use `predict()` instead, and ",
    "!" = "pass in the original data."))
}

#' @method summary sparsegl
#' @export
summary.sparsegl <- function(object, ...) {
  rlang::check_dots_empty()
  ns <- length(object$lambda)
  if (ns > 5) {
    xlam <- round(stats::quantile(1:ns))
    names(xlam) <- c("Max.", "3rd Qu.", "Median", "1st Qu.", "Min.")
  } else {
    xlam <- seq_len(ns)
    names(xlam) <- paste0("s", seq_len(ns))
  }
  nz <- predict(object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(object$group[x])))
  tab <- with(object, data.frame(
    lambda = lambda[xlam],
    index = xlam,
    nnzero = nnzero[xlam],
    active_grps = active_grps[xlam])
  )
  rownames(tab) <- names(xlam)
  out <- structure(list(call = object$call, table = tab),
                   class = "summary.sparsegl")
  out
}

#' @method print summary.sparsegl
#' @export
print.summary.sparsegl <- function(
    x, digits = max(3, getOption("digits") - 3), ...) {

  rlang::check_dots_empty()
  lambda_warning <- all(x$table$nnzero == 0)

  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)


  if (lambda_warning) {
    cat("Warning: all regularization parameters resulted in empty models.\n\n")
  }

  cat("Summary of Lambda sequence:\n")
  print(x$tab, digits = digits)
  cat("\n")

}

#' @method print sparsegl
#' @export
print.sparsegl <- function(x, digits = min(3, getOption("digits") - 3), ...) {
  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}


