sgl_irwls <- function(
  bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
  flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
  lower_bnd, upper_bnd, family, weights) {

  # Validate the family object
  varfun <- family$variance
  linkinv <- family$linkinv
  if (!is.function(varfun) || !is.function(linkinv)) {
    rlang::abort("'family' object is not a valid family object. Try `?family`.")
  }
  mu_eta <- family$mu.eta
  valideta <- family$valideta %||% function(eta) TRUE
  validmu <- family$validmu %||% function(mu) TRUE


  is_sparse <- FALSE
  if (!is.numeric(y)) stop("For family = 'gaussian', y must be numeric.")
  if (inherits(x, "sparseMatrix")) {
    is_sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  if (is_sparse) {
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
  }

  if (is.null(weights)) weights <- rep(1, nvars)

  init <- irls_initilizer(x, y, weights, family, intr, pf, pfl1, asparse)
  nulldev <- init$nulldev
  mu <- init$mu
  eta <- family$linkfun(mu)

  # work out lambda values
  nlam <- as.integer(nlam)
  user_lambda <- !((length(ulam) == 1L) & (ulam == 0)) # user provided lambda
  if (!user_lambda) {
    lambda_max <- init$lambda_max
    ulam <- exp(seq(log(lambda_max), log(lambda_max * flmin),
                        length.out = nlam))
    cur_lam <- 9.9e30
  } else {
    cur_lam <- ulam[1]
  }

  # Everything set up, enter the lambda loop
  for (k in 1:nlam) {
    if (k > 1) cur_lam <- ulam[k]

    


  gamma <- calc_gamma(x, ix, iy, bn)

  if (!is_sparse) {
    fit <- dotCall64::.C64(
      "sparse_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "double", "double", "double",
                    "integer", "integer", "integer", "double", "double",
                    "double", "integer", "integer", "integer", "double", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
      pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam),
      activeGroup = integer_dc(pmax), nbeta = integer_dc(nlam),
      alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1), mse = numeric_dc(nlam),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 9), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  } else { # sparse design matrix
    fit <- dotCall64::.C64(
      "spmat_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer",
                    "integer", "integer", "double", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1), mse = numeric_dc(nlam),
      # Read only
      alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 14), rep("rw", 8), rep("w", 9), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  }

  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs
  if (intr) {
    outlist$b0 <- outlist$b0 + ym
  } else {
    outlist$b0 <- rep(0, dim(outlist$beta)[2])
  }
  outlist <- c(outlist,
               list(npasses = fit$npass, jerr = fit$jerr, group = group,
                    mse = fit$mse[seq(fit$nalam)]))
  class(outlist) <- c("ls")
  outlist
}

dev_function <- function(y, mu, weights, family) {
    sum(family$dev.resids(y, mu, weights))
}

irls_initilizer <- function(x, y, weights, family, intercept, pf, pfl1, alpha) {
  nobs <- nrow(x)
  nvars <- ncol(x)

  mu <- family$linkinv(y * 0)
  nulldev <- dev_function(y, mu, weights, family)

  r <- y - mu
  eta <- family$linkfun(mu)
  v <- family$variance(mu)
  mu_eta <- family$mu.eta(eta)
  weights <- weights / sum(weights)
  rv <- r / v * mu_eta * weights
  g <- abs(drop(crossprod(rv, x))) / pfl1
  lambda_max <- max(g)

  list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
}