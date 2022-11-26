sgl_ls <- function(
  bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
  flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
  lower_bnd, upper_bnd) {
  # call Fortran core
  is.sparse <- FALSE
  if (!is.numeric(y)) rlang::abort("For family = 'gaussian', y must be numeric.")
  if (inherits(x,"sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  ym <- mean(y)
  if (intr) {
    y <- y - ym
    nulldev <- mean(y^2)
  } else {
    nulldev <- mean((y - ym)^2)
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
  outlist$npasses <- fit$npass
  outlist$jerr <- fit$jerr
  outlist$group <- group
  outlist$mse <- fit$mse[seq(fit$nalam)]
  outlist$dev.ratio <- 1 - outlist$mse / nulldev
  outlist$nulldev <- nulldev
  class(outlist) <- c("lsspgl")
  outlist
}
