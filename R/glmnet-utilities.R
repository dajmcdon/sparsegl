######################################################################
## These functions are minor modifications or directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.


#' @importFrom methods as
as_dgCMatrix <- function(x) {
  as(as(x, "sparseMatrix"), "CsparseMatrix")
}


err <- function(n, maxit, pmax) {
  if (n == 0) msg <- ""
  if (n > 0) {
    #fatal error
    if (n < 7777)
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000)
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in sparsegl fortran code -", msg)
  }
  if (n < 0) {
    #non fatal error
    if (n > -10000)
      msg <- paste(
        "Convergence for ", -n, "th lambda value not reached after maxit=",
        maxit, " iterations; solutions for larger lambdas returned",
        sep = "")
    if (n < -10000)
      msg <- paste(
        "Number of nonzero coefficients along the path exceeds pmax=",
        pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned",
        sep = "")
    n <- -1
    msg <- paste("from gglasso fortran code -", msg)
  }
  list(n = n, msg = msg)
}




getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


getoutput <- function(x, group, fit, maxit, pmax, nvars, vnames, eps) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n),
         `1` = stop(errmsg$msg, call. = FALSE),
         `-1` = print(errmsg$msg, call. = FALSE))
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- Matrix::drop0(
      matrix(fit$beta[seq(nvars * nalam)], nvars, nalam,
             dimnames = list(vnames, stepnames)),
      tol = eps^2)
    df <- apply(abs(beta) > 0, 2, sum) ## this is wrong, but fast
  } else {
    beta <- Matrix::Matrix(0, nvars, nalam,
                           dimnames = list(vnames, stepnames))
    df <- rep(0, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}

lambda.interp <- function(lambda, s) {
  ### lambda is the index sequence that is produced by the model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indicies, and a
  #   vector of fractions.
  ### the new values are interpolated between the two using the
  #   fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}


lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

nonzeroCoef <- function(beta) {
  nr <- nrow(beta)
  if (nr == 1) { # degenerate case
    apply(beta, 2, function(x) if (abs(x) > 0) 1 else NULL)
  } else {
    beta <- abs(beta) > 0 # this is sparse
    which <- seq(nr)
    ones <- rep(1, ncol(beta))
    nz <- as.vector((beta %*% ones) > 0)
    which <- which[nz]

    if (length(which) > 0) {
      beta <- as.matrix(beta[which, , drop = FALSE])
      nzel <- function(x, which) if (any(x)) which[x] else NULL
      which <- apply(beta, 2, nzel, which)
      if (!is.list(which)) which <- data.frame(which)
      which
    } else {
      dn <- dimnames(beta)[[2]]
      which <- vector("list", length(dn))
      names(which) <- dn
      which
    }
  }
}
