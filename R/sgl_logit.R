sgl_logit <- function(
  bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1,
  dfmax, pmax, nlam, flmin, ulam, eps,
  maxit, vnames, group, intr, asparse, standardize,
  lower_bnd, upper_bnd) {
  # call Fortran core
  y <- as.factor(y)
  lev <- levels(y)
  ntab <- table(y)
  minclass <- min(ntab)
  if (minclass <= 1)
    stop("Binomial regression: one class has 1 or 0 observations; not allowed")
  if (length(ntab) != 2)
    stop("Binomial regression: more than one class is not supported")
  if (minclass < 8)
    warning(paste0("Binomial regression: one class has fewer than 8\n",
                   "observations; dangerous ground"))
  # TODO, enable prediction with class labels if factor is passed
  y <- 2 * (as.integer(y) - 1) - 1 # convert to -1 / 1

  is.sparse <- FALSE
  if (inherits(x, "sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  if (is.sparse) {
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
  }

  gamma <- 0.25 * calc_gamma(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "log_sparse_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "double", "double", "double",
                    "integer", "integer", "integer", "double", "double",
                    "double", "integer", "integer", "integer", "double",
                    "double", "integer", "integer", "double", "integer",
                    "integer", "double", "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x), y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 8), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  } else {
     fit <- dotCall64::.C64(
      "log_spmat_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer",
                    "integer", "integer", "double", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1),
      # Read only
      alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 14), rep("rw", 8), rep("w", 8), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  }
  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs

  outlist$b0 <- matrix(outlist$b0, nrow = 1)
  outlist <- c(outlist,
               list(npasses = fit$npass, jerr = fit$jerr, group = group,
                    classnames = lev))
  class(outlist) <- c("logit")
  outlist
}
