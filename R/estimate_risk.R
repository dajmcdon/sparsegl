#' Calculate information criteria.
#'
#' This function uses the degrees of freedom to calculate various information
#' criteria. This function uses the "unknown variance" version of the likelihood.
#' Only implemented for Gaussian regression. The constant is ignored (as in
#' [stats::extractAIC()]).
#'
#' @param object fitted object from a call to [sparsegl()].
#' @template param_x-template
#' @template param_y-template
#' @param type one or more of AIC, BIC, or GCV.
#' @param approx_df the `df` component of a [sparsegl()] object is an
#'   approximation (albeit a fairly accurate one) to the actual degrees-of-freedom.
#'   However, the exact value requires inverting a portion of `X'X`. So this
#'   computation may take some time (the default computes the exact df).
#' @seealso [sparsegl()] method.
#' @references Vaiter S, Deledalle C, Peyré G, Fadili J, Dossal C. (2012). \emph{The
#' Degrees of Freedom of the Group Lasso for a General Design}.
#' \url{https://arxiv.org/pdf/1212.6478.pdf}.
#' @return a `data.frame` with as many rows as `object$lambda`. It contains
#'   columns `lambda`, `df`, and the requested risk types.
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
#' estimate_risk(fit1, X, y, type = "AIC")
estimate_risk <- function(object, x, y,
                          type = c("AIC", "BIC", "GCV"),
                          approx_df = FALSE) {
  if (! "ls" %in% class(object)) stop("Only linear regression is supported.")
  type <- match.arg(type, several.ok = TRUE)
  if (is.matrix(y)) {
    stopifnot(dim(y)[2] == 1)
    y <- drop(y)
  }
  preds <- predict(object, x)
  err <- log(colMeans((y - preds)^2))
  n <- length(y)

  if (approx_df) df <- object$df
  else df <- exact_df(object, x)
  out <- data.frame(lambda = object$lambda,
                    df = df,
                    AIC = err + 2 * df / n,
                    BIC = err + log(n) * df / n,
                    GCV = err - 2 * log(1 - df / n)) # actually log(GCV)
  out <- out[c("lambda", "df", type)]
  return(out)
}


exact_df <- function(object, x) {
  # See the correct formula in https://arxiv.org/pdf/1212.6478.pdf
  # Theorem 2
  Iset <- abs(object$beta) > 0
  Imax <- which(apply(Iset, 1, any))
  Iset <- Iset[Imax,]
  group <- object$group[Imax]
  beta <- object$beta[Imax,]
  nlambda <- length(object$lambda)
  xx <- Matrix::crossprod(x[,Imax])
  df <- double(nlambda)
  for (i in seq(nlambda)) {
    Idx <- Iset[, i]
    if (any(Idx)) {
      xx_sub <- xx[Idx, Idx]
      del <- delP(beta[Idx, i], group[Idx])
      df[i] <- sum(solve(xx_sub + object$lambda[i] * del) * xx_sub)
    } else {
      df[i] <- 0
    }
  }
  return(df)
}

delP <- function(beta, group) {
  betas <- split(beta, group)
  mats <- lapply(betas, function(x) {
    p <- length(x)
    bn <- two_norm(x)
    Matrix::diag(1/bn, nrow = p, ncol = p) - outer(x, x) / bn^3

  })
  return(Matrix::bdiag(mats))
}
