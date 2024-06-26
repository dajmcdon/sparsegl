#' @importFrom dotCall64 integer_dc vector_dc numeric_dc
sgl_ls <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
    flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd, weights) {
  # call Fortran core
  is.sparse <- FALSE
  if (!is.numeric(y)) cli_abort("For family = 'gaussian', y must be numeric.")
  if (inherits(x, "sparseMatrix")) {
    is.sparse <- TRUE
    x <- as_dgCMatrix(x)
  }

  my <- stats::weighted.mean(y, weights)
  nulldev <- sum((y - my)^2 * weights)
  #if (intr) y <- y - mean(y)

  y <- y * sqrt(weights)
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  x <- sqrt(weights) * x
  gamma <- calc_gamma(x, ix, iy, bn)

  if (!is.sparse) {
    fit <- dotCall64::.C64(
      "sparse_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "double", "double", "double", "double",
        "integer", "integer", "integer", "double", "double",
        "double", "integer", "integer", "integer", "double", "double", "integer",
        "integer", "double", "integer", "integer", "double",
        "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(x),
      y = as.double(y), pf = pf,
      pfl1 = pfl1, w = as.double(sqrt(weights)),
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
      INTENT = c(rep("r", 12), rep("rw", 8), rep("w", 9), rep("r", 3)),
      NAOK = TRUE,
      PACKAGE = "sparsegl"
    )
  } else { # sparse design matrix
    # x <- sqrt(weights) * x
    xidx <- as.integer(x@i + 1)
    xcptr <- as.integer(x@p + 1)
    xval <- as.double(x@x)
    nnz <- as.integer(utils::tail(x@p, 1))
    fit <- dotCall64::.C64(
      "spmat_four",
      SIGNATURE = c(
        "integer", "integer", "integer", "integer", "double",
        "integer", "integer", "double", "integer", "integer",
        "integer", "double", "double", "double", #"double",
        "integer", "integer",
        "integer", "double", "double", "double", "integer",
        "integer", "integer", "double", "double", "integer",
        "integer", "double", "integer", "integer", "double",
        "double", "double", "double"
      ),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(xval), xidx = xidx, xcptr = xcptr,
      nnz = nnz, y = as.double(y), pf = pf, pfl1 = pfl1, #w = as.double(weights),
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
      PACKAGE = "sparsegl"
    )
  }
  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  # print(outlist$b0)
  # print(ym)
  # print(my)
  # print(mean(y))
  # print(weighted.mean(y, weights))
  # print(mean(x %*% outlist$beta))
  # print(weighted.mean(x %*% outlist$beta, weights))
  if (standardize) outlist$beta <- outlist$beta * xs
  if (!intr) {
    outlist$b0 <- rep(0, dim(outlist$beta)[2])
  } #else {
    #outlist$b0 <- outlist$b0 + mean(y)
  #}

  outlist$npasses <- fit$npass
  outlist$jerr <- fit$jerr
  outlist$group <- group
  outlist$mse <- fit$mse[seq(fit$nalam)]
  outlist$dev.ratio <- 1 - outlist$mse / nulldev
  outlist$nulldev <- nulldev
  class(outlist) <- c("lsspgl")
  outlist
}
