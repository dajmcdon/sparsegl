#' Calculate information criteria.
#'
#' This function uses the degrees of freedom to calculate various information
#' criteria. This function uses the "unknown variance" version of the likelihood.
#' Only implemented for Gaussian regression. The constant is ignored (as in
#' [stats::extractAIC()]).
#'
#' @param object fitted object from a call to [sparsegl()].
#' @param x Matrix. The matrix of predictors used to estimate
#'   the `sparsegl` object. May be missing if `approx_df = TRUE`.
#' @param type one or more of AIC, BIC, or GCV.
#' @param approx_df the `df` component of a `sparsegl` object is an
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
#' estimate_risk(fit1, type = "AIC", approx_df = TRUE)
estimate_risk <- function(object, x,
                          type = c("AIC", "BIC", "GCV"),
                          approx_df = FALSE) {
  UseMethod("estimate_risk")
}

#' @export
estimate_risk.default <- function(object, x,
                                  type = c("AIC", "BIC", "GCV"),
                                  approx_df = FALSE) {
  cli::cli_abort("Risk estimation is only available for Gaussian likelihood.")
}

#' @export
estimate_risk.lsspgl <- function(object, x,
                                 type = c("AIC", "BIC", "GCV"),
                                 approx_df = FALSE) {
  type <- match.arg(type, several.ok = TRUE)

  err <- log(object$mse)
  n <- object$nobs

  if (approx_df) df <- object$df + df_correction(object)
  else df <- exact_df(object, x)
  out <- data.frame(
    lambda = object$lambda,
    df = df,
    AIC = err + 2 * df / n,
    BIC = err + log(n) * df / n,
    GCV = err - 2 * log(1 - df / n) # actually log(GCV)
  )
  out <- out[c("lambda", "df", type)]
  return(out)
}

df_correction <- function(obj) {
  # This is based on X being orthogonal. See the proof of corollary 1
  # in https://arxiv.org/pdf/1212.6478.pdf
  group <- obj$group
  beta <- obj$beta
  lambda <- obj$lambda
  pf <- (1 - obj$asparse) * obj$pf_group
  num <- pmax(apply(beta, 2, grouped_zero_norm, gr = group) - 1, 0)
  norms <- apply(beta, 2, grouped_two_norm, gr = group)
  colSums(num / (1 + outer(pf, lambda) / norms))
}



exact_df <- function(object, x) {
  # See the correct formula in https://arxiv.org/pdf/1212.6478.pdf
  # Theorem 2
  if (missing(x))
    stop("Risk estimation with exact df requires the design matrix `x`.")
  Iset <- abs(object$beta) > 0
  Imax <- which(apply(Iset, 1, any))
  Iset <- Iset[Imax,]
  group <- object$group[Imax]
  beta <- object$beta[Imax,]
  pf <- (1 - object$asparse) * object$pf_group[group]
  nlambda <- length(object$lambda)
  xx <- Matrix::crossprod(x[,Imax])
  df <- double(nlambda)
  for (i in seq(nlambda)) {
    Idx <- Iset[, i]
    if (any(Idx)) {
      gr <- group[Idx]
      xx_sub <- as(xx[Idx, Idx], "Matrix")
      del <- delP(beta[Idx, i], gr)
      li <- object$lambda[i] * pf[Idx]
      df[i] <- sum(Matrix::solve(xx_sub + li * del) * xx_sub)
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
