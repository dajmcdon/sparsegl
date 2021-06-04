#' Information criteria
#'
#' @param object fitted object from a call to [`sparsegl()`]
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#'   is a vector of measurement and each column is a feature
#' @param y real-valued response variable.
#' @param type one of AIC, BIC, or GCV
#' @param approx_df the `df` component of a `sparsegl` object is an approximation
#'   (albeit a fairly accurate one) to the actual degrees-of-freedom. However,
#'   the exact value requires inverting a portion of `X'X`. So this may take
#'   some time.
#'
#' @return a vector of the same length as `object$lambda`
#' @export
estimate_risk <- function(object, x, y,
                          type = c("AIC", "BIC", "GCV"),
                          approx_df = FALSE) {
  type <- match.arg(type)
  if (is.matrix(y)) {
    stopifnot(dim(y)[2] == 1)
    y <- drop(y)
  }
  preds <- predict(object, x)
  err <- log(colMeans((y - preds)^2))
  n <- length(y)
  
  if (approx_df) {
    df <- object$df
  } else {
    df <- exact_df(object, x)
  }
  
  pen <- switch(type,
                AIC = 2 * df / n,
                BIC = log(n) * df / n,
                GCV = -2 * log(1 - df / n))
  risk <- err + pen
  if (type == "GCV") risk <- exp(risk)
  return(risk)
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
    Matrix::diag(1/bn, nrow = p, ncol = p) - outer(x / bn, x)
  })
  return(Matrix::bdiag(mats))
}


two_norm <- function(x) sqrt(sum(x^2))
gr_norm <- function(x, gr) sum(as.vector(tapply(x, gr, two_norm)))
sp_group_norm <- function(x, gr, asparse = 0.05) {
  asparse * sum(abs(x)) + (1 - asparse) * gr_norm(x, gr)
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
  xx <- x[,Imax]
  df <- double(nlambda)
  for (i in seq(nlambda)) {
    Idx <- Iset[, i]
    if (any(Idx)) {
      X_I <- xx[, Idx]
      del <- delP(beta[Idx, i], group[Idx])
      df[i] <- sum(diag(X_I %*% solve(t(X_I) %*% X_I + object$lambda[i] * del) %*% t(X_I)))
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
    (Matrix::diag(x / bn, nrow = p, ncol = p) - outer(x, x)) / bn
  })
  return(Matrix::bdiag(mats))
} 
