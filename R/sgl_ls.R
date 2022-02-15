sgl_ls <- function(
  bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps,
  maxit, vnames, group, intr, asparse, standardize,
  lower_bnd, upper_bnd) {
  # call Fortran core
  is.sparse <- FALSE
  if (!is.numeric(y)) stop("For family = 'gaussian', y must be numeric.")
  if (inherits(x,"sparseMatrix")) {
    is.sparse <- TRUE
    x <- methods::as(x,"CsparseMatrix")
    x <- methods::as(x,"dgCMatrix")
  }
  if (intr) {
    ym <- mean(y)
    y <- y - ym
    if (!is.sparse) {
      xm <- colMeans(x)
      x <- sweep(x,2,xm)
    }
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

  gamma <- calc_gamma(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "sparse_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "double", "double",
                    "integer", "integer", "integer", "double", "double",
                    "double", "integer", "integer", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
      # Read / write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
      eps = eps, maxit = maxit,
      # Write only
      nalam = integer_dc(1), beta = numeric_dc(nvars * nlam),
      activeGroup = integer_dc(pmax), nbeta = integer_dc(nlam),
      alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 10), rep("rw", 7), rep("w", 7), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  } else { # sparse design matrix
    fit <- dotCall64::.C64(
      "spmat_four",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "integer", "integer",
                    "integer", "double", "double", "integer", "integer",
                    "integer", "double", "double", "double", "integer",
                    "integer", "integer", "double", "double", "integer",
                    "integer", "double", "integer", "integer", "double",
                    "double", "double"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf,
      # Read write
      dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Write only
      nalam = integer_dc(1), b0 = numeric_dc(nlam),
      beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
      nbeta = integer_dc(nlam), alam = numeric_dc(nlam),
      npass = integer_dc(1), jerr = integer_dc(1),
      # Read only
      alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd,
      INTENT = c(rep("r", 13), rep("rw", 8), rep("w", 8), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  }

  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) outlist$beta <- outlist$beta * xs
  if (intr) {
    if (is.sparse) outlist$b0 <- outlist$b0 + ym
    else outlist$b0 <- ym - xm %*% outlist$beta
  } else {
    if (!is.sparse) outlist$b0 <- rep(0, dim(outlist$beta)[2])
  }
  outlist <- c(outlist,
               list(npasses = fit$npass, jerr = fit$jerr, group = group))
  class(outlist) <- c("ls")
  outlist
}
