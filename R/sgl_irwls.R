# Description of the structure
#
# IRLS is a sequence of loops for unregularized problems, and here, it is
# potentially even more complicated.
#
# Middle loop:
#   irwls_fit - This function looks much like `stats::glm.fit()`. Essentially,
#     if converts a single Exponential family regression problem to a
#     weighted least squares problem. The wls is solved by a call to the FORTRAN
#     code via `spgl_wlsfit()`. However, we want (1) previous solutions at
#     earlier lambda to be warm starts and (2) previous pass through irwls_fit
#     to continue to serve as a warm start.
#
#     This function does all the checking along with step size halving if
#     needed. Inside it repeats irwls until convergence.
#
# Outer loop:
#   sgl_irwls - This function checks and sets up the lambda loop. If the user
#     provided lambda, we use those. Otherwise, we start at a large value
#     and iterate down until we find nonzero groups. Then we construct the
#     sequence. This is the same structure as in the FORTRAN code for
#     the build in families.
#
#     Inside, we create a warm_start object if it doesn't already exist.
#     We also separate into a list of static arguments that doesn't change
#     in any loops ever. Finally, all results are saved.




# Fit a GLM with sparse group regularization for a path of lambda values
#
# Fit a generalized linear model via penalized maximum likelihood for a path of
# lambda values. Can deal with any GLM family.
#
# This organization is based largely off [stats::glm.fit()] with some extras
# to handle the path.
#
# Sometimes the sequence is truncated before \code{nlam} values of lambda
# have been used. This happens when the
# decrease in deviance is marginal (i.e. we are near a saturated fit).
#
#' @importFrom stats gaussian
sgl_irwls <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
    flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd, weights = NULL, offset = NULL, family = gaussian(),
    trace_it = 0, warm = NULL) {

  validate_family(family)

  ## Init the family just to check that it works
  if (is.null(weights)) weights <- rep(1, nobs)
  etastart <- 0
  mustart <- NULL
  start <- NULL
  eval(family$initialize)
  y <- drop(y)
  has_offset <- !is.null(offset)
  if (!has_offset) offset <- as.double(y * 0)

  # standardize x if necessary
  is_sparse <- FALSE
  if (inherits(x, "sparseMatrix")) {
    is_sparse <- TRUE
    x <- as_dgCMatrix(x)
  }
  xs <- NULL
  if (standardize) {
    sx <- sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- x %*% Matrix::Diagonal(x = xs)
  }


  # get null deviance and lambda max, work out lambda values
  # we ALWAYS fit the intercept inside wsgl, so start it at zero
  init <- initializer(x, y, weights, family, intr = intr,
                      has_offset, offset, ulam, flmin)
  # this is supposed to be an upper bound
  # work out lambda values, cur_lambda is lambda_max / 0.99 when appropriate.
  cur_lambda <- init$cur_lambda
  findlambda <- init$findlambda
  no_user_lambda <- init$findlambda
  nulldev <- init$nulldev

  if (trace_it == 1) pb <- utils::txtProgressBar(min = 0, max = nlam, style = 3)

  # preallocate space to store output
  b0 <- double(nlam)
  beta <- matrix(0, nvars, nlam)
  dev.ratio <- rep(NA, length = nlam)
  mnl <- min(nlam, 6L)
  flmin <- max(flmin, 1e-6) # threshold above zero

  static_args <- list(
    nulldev = as.double(nulldev),
    y = as.double(y),
    weights = as.double(weights),
    offset = as.double(offset),
    bn = as.integer(bn),
    bs = as.integer(bs),
    x = x,
    ix = as.integer(ix),
    iy = as.integer(iy),
    xs = as.double(xs),
    nobs = as.integer(nobs),
    nvars = as.integer(nvars),
    pf = as.double(pf),
    pfl1 = as.double(pfl1),
    dfmax = as.integer(dfmax),
    pmax = as.integer(pmax),
    flmin = as.double(flmin),
    eps = as.double(eps),
    maxit = as.integer(maxit),
    vnames = vnames,
    group = group,
    intr = as.integer(intr),
    asparse = asparse,
    lb = as.double(lower_bnd),
    ub = as.double(upper_bnd),
    family = family,
    trace_it = trace_it)

  if (is.null(warm))
    warm <- make_irls_warmup(nobs, nvars, b0 = init$b0, r = init$r)
  if (!inherits(warm, "irlsspgl_warmup")) {
    cli::cli_abort(
      "the `warm` object should be created with `make_irls_warmup()`."
    )
  }

  warm$activeGroup <- integer(pmax)
  warm$activeGroupIndex <- integer(bn)
  warm$sset <- integer(bn)
  if (length(zn <- which(grouped_zero_norm(warm$beta, group) > 0))) {
    warm$activeGroup[seq_along(zn)] <- zn
    warm$activeGroupIndex[zn] <- seq_along(zn)
    warm$sset[zn] <- 1L
  }
  warm <- c(warm, ni = 0L, npass = 0L, me = 0L, findlambda = findlambda,
            eset = list(warm$sset))

  l <- 0L
  while (l < nlam) {
    if (!findlambda) {
      l <- l + 1L
      warm$al0 <- as.double(cur_lambda)
      warm$ulam <- as.double(ulam[l])
      if (trace_it == 2)
        cat("Fitting lambda index", l, ":", ulam[l], fill = TRUE)
    } else {
      # trying to find lambda max, we started too big
      warm$al0 <- as.double(cur_lambda)
      warm$ulam <- as.double(init$lambda_max)
      if (trace_it == 2)
        cat("Trying to find a reasonable starting lambda.", fill = TRUE)
    }
    warm$l <- as.integer(l)

    # browser()
    # here we dispatch to wls
    warm <- irwls_fit(warm, static_args)
    if (trace_it == 1) utils::setTxtProgressBar(pb, l)
    if (warm$jerr != 0) {
      if (l > 1L) {
        cli::cli_warn(c(
          "Convergence for {l}th lambda value not reached after maxit =",
          "{maxit} iterations; solutions for larger lambdas returned."))
        l <- l - 1L
        break
      } else {
        cli::cli_warn(c(
          "Convergence for initial lambda value not reached after maxit =",
          "{maxit} iterations; no solutions available."))
      }
    }

    # browser()
    if (findlambda) { # we searched inside the FORTRAN code, now we found it
      ulam <- double(nlam)
      alf <- flmin^(1/(nlam - 1))
      ulam[1:nlam] <- exp(log(warm$ulam) + -1:(nlam - 2) * log(alf))
      l <- 2L
      findlambda <- FALSE
      dev.ratio[1] <- 0
      b0[1] <- init$b0
    }

    b0[l] <- warm$b0
    beta[, l] <- as.matrix(warm$beta)
    dev.ratio[l] <- warm$dev.ratio

    # early stopping if dev.ratio almost 1 or no improvement
    if (l >= mnl && no_user_lambda) {
      if (dev.ratio[l] > 1 - 1e-4) break
      else if (l > 1) {
        if (family$family == "gaussian") {
          if (dev.ratio[l] - dev.ratio[l - 1] < 1e-5 * dev.ratio[l]) break
        } else if (family$family == "poisson") {
          if (dev.ratio[l] - dev.ratio[l - mnl + 1] < 1e-4 * dev.ratio[l]) break
        } else if (dev.ratio[l] - dev.ratio[l - 1] < 1e-5) break
      }
    }
  }
  if (trace_it == 1) {
    utils::setTxtProgressBar(pb, nlam)
    cat("", fill = TRUE)
  }

  # truncate if we quit early
  if (l < nlam) {
    b0 <- b0[1:l]
    beta <- beta[, 1:l, drop = FALSE]
    dev.ratio <- dev.ratio[1:l]
    ulam <- ulam[1:l]
  }

  # output
  stepnames <- paste0("s", 0:(length(ulam) - 1))
  out <- list()
  out$b0 <- b0
  names(out$b0) <- stepnames
  out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                             dimnames = list(vnames, stepnames))
  if (standardize) out$beta <- out$beta * xs
  out$df <- apply(abs(out$beta) > 0, 2, sum)
  out$dim <- dim(out$beta)
  out$lambda <- ulam
  out$npasses <- warm$npass
  out$jerr <- warm$jerr
  out$group <- group
  out$deviance <- (1 - dev.ratio) * nulldev / nobs
  out$dev.ratio <- dev.ratio
  out$nulldev <- nulldev / nobs
  out$offset <- offset
  out$family <- family
  class(out) <- "irlsspgl"

  return(out)
}

irwls_fit <- function(warm, static) {
  maxit_irls <- 50

  # initialize everything
  variance <- static$family$variance
  linkinv <- static$family$linkinv
  mu.eta <- static$family$mu.eta
  trace_it <- static$trace_it

  fit <- warm
  nulldev <- static$nulldev
  coefold <- fit$beta   # prev value for coefficients
  intold <- fit$b0    # prev value for intercept
  lambda <- fit$ulam
  eta <- get_eta(static$x, static$xs, coefold, intold)
  mu <- linkinv(eta <- eta + static$offset)

  valideta <- static$family$valideta %||% function(eta) TRUE
  validmu <- static$family$validmu %||% function(mu) TRUE


  if (!validmu(mu) || !valideta(eta)) {
    cli::cli_abort(c("Cannot find valid starting values.",
                   "!" = "Please specify some with `make_irls_warmup()`."))
  }

  start <- NULL     # current value for coefficients
  start_int <- NULL # current value for intercept
  obj_val_old <- obj_function(
    static$y, mu, static$group, static$weights, static$family,
    static$pf, static$pfl1, static$asparse,
    coefold, fit$ulam
  )

  if (static$trace_it == 2) {
    cat("Warm Start Objective:", obj_val_old, fill = TRUE)
  }
  conv <- FALSE      # converged?
  iter <- fit$npass

  # IRLS loop
  while (iter < static$maxit && !conv) {
    # some checks for NAs/zeros
    varmu <- variance(mu)
    if (anyNA(varmu)) rlang::abort("NAs in V(mu)")
    if (any(varmu == 0)) rlang::abort("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (anyNA(mu.eta.val)) rlang::abort("NAs in d(mu)/d(eta)")

    # d ell / d beta = X'W (z - mu) / mu.eta.val (theoretically)
    #                = t(wx) %*% r,  (code)
    # can use directly in wls for strong rule / kkt

    # compute working response and weights
    z <- (eta - static$offset) + (static$y - mu) / mu.eta.val
    w <- sqrt((static$weights * mu.eta.val^2) / variance(mu))
    r <- w * (z - eta + static$offset)
    fit$r <- r

    # we updated w, so we need to update gamma
    wx <- Matrix::Diagonal(x = w) %*% static$x
    # because we use eig(xwx) = svd(w^{1/2}x)
    gamma <- calc_gamma(wx, static$ix, static$iy, static$bn)

    # do WLS with warmstart from previous iteration, dispatch to FORTRAN wls
    fit <- spgl_wlsfit(fit, wx, gamma, static)

    if (fit$jerr != 0) return(list(jerr = fit$jerr))

    # update coefficients, eta, mu and obj_val
    start <- fit$beta
    lambda <- fit$ulam
    start_int <- fit$b0
    eta <- get_eta(static$x, static$xs, start, start_int)
    mu <- linkinv(eta <- eta + static$offset)
    obj_val <- obj_function(
      static$y, mu, static$group, static$weights, static$family, static$pf,
      static$pfl1, static$asparse, start, lambda)
    iter <- iter + fit$npass
    fit$npass <- iter
    if (trace_it == 2) cat("Iteration", iter, "Objective:", obj_val, fill = TRUE)

    boundary <- FALSE
    halved <- FALSE  # did we have to halve the step size?
    # if objective function is not finite, keep halving the stepsize until it is finite
    if (!is.finite(obj_val) || obj_val > 9.9e30) {
      rlang::warn("Infinite objective function!")
      if (is.null(coefold) || is.null(intold)) {
        cli::cli_abort(
          c("No valid set of coefficients has been found.",
            "!" = "Please specify some with `make_irls_warmup()`."))
      }
      rlang::warn("step size truncated due to divergence")
      ii <- 1
      while (!is.finite(obj_val) || obj_val > 9.9e30) {
        if (ii > maxit_irls)
          rlang::abort("inner loop 1; cannot correct step size")
        ii <- ii + 1
        start <- (start + coefold) / 2
        start_int <- (start_int + intold) / 2
        eta <- get_eta(static$x, static$xs, start, start_int)
        mu <- linkinv(eta <- eta + static$offset)
        obj_val <- obj_function(
          static$y, mu, static$group, static$weights, static$family, static$pf,
          static$pfl1, static$asparse, start, lambda)
        if (trace_it == 2)
          cat("Iteration", iter, " Halved step 1, Objective:", obj_val, fill = TRUE)
      }
      boundary <- TRUE
      halved <- TRUE
    }
    # if some of the new eta or mu are invalid, keep halving stepsize until valid
    if (!(valideta(eta) && validmu(mu))) {
      rlang::warn("Invalid eta / mu!")
      if (is.null(coefold) || is.null(intold)) {
        cli::cli_abort(c("No valid set of coefficients has been found.",
                         "!" = "Please specify some with `make_irls_warmup()`."
        ))
      }
      rlang::warn("step size truncated: out of bounds")
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > maxit_irls)
          rlang::abort("inner loop 2; cannot correct step size.")
        ii <- ii + 1
        start <- (start + coefold) / 2
        start_int <- (start_int + intold) / 2
        eta <- get_eta(static$x, static$xs, start, start_int)
        mu <- linkinv(eta <- eta + static$offset)
      }
      boundary <- TRUE
      halved <- TRUE
      obj_val <- obj_function(
        static$y, mu, static$group, static$weights, static$family, static$pf,
        static$pfl1, static$asparse, start, lambda)
      if (trace_it == 2) cat("Iteration", iter, " Halved step 2, Objective:",
                             obj_val, fill = TRUE)
    }
    # extra halving step if objective function value actually increased
    if (obj_val > obj_val_old + 1e-7) {
      ii <- 1
      while (obj_val > obj_val_old + 1e-7) {
        if (ii > maxit_irls)
          rlang::abort("inner loop 3; cannot correct step size")
        ii <- ii + 1
        start <- (start + coefold) / 2
        start_int <- (start_int + intold) / 2
        eta <- get_eta(static$x, static$xs, start, start_int)
        mu <- linkinv(eta <- eta + static$offset)
        obj_val <- obj_function(
          static$y, mu, static$group, static$weights, static$family, static$pf,
          static$pfl1, static$asparse, start, lambda)
        if (trace_it == 2) cat("Iteration", iter, " Halved step 3, Objective:",
                               obj_val, fill = TRUE)
      }
      halved <- TRUE
    }

    # if we did any halving, we have to update the coefficients, intercept
    # and weighted residual in the warm_fit object
    if (halved) {
      fit$beta <- start
      fit$b0 <- start_int
      fit$r <- w * (z - eta)
    }

    # test for convergence
    if (abs(obj_val - obj_val_old) / (0.1 + abs(obj_val)) < static$eps) {
      conv <- TRUE
      break
    } else {
      coefold <- start
      intold <- start_int
      obj_val_old <- obj_val
    }
  }
  # end of IRLS loop

  # checks on convergence and fitted values
  if (!conv)
    rlang::warn("sparsgl_irls: algorithm did not converge")
  if (boundary)
    rlang::warn("sparsgl_irls: algorithm stopped at boundary value")

  # some extra warnings, printed only if trace_it == 2
  if (trace_it == 2) {
    tiny <- 10 * .Machine$double.eps
    if (static$family$family == "binomial" &&
        (any(mu > 1 - tiny) || any(mu < tiny))) {
      rlang::warn("sparsgl_irls: fitted probabilities numerically 0 or 1 occurred")
    }
    if (static$family$family == "poisson" && any(mu < static$eps))
      rlang::warn("sparsegl_irls: fitted rates numerically 0 occurred")
  }

  # prepare output object
  # if (save.fit == FALSE) fit$warm_fit <- NULL
  fit$offset <- static$offset
  fit$nulldev <- nulldev
  fit$dev.ratio <- 1 - dev_function(static$y, mu, static$weights,
                                    static$family) / fit$nulldev
  fit$family <- static$family
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val
  fit$npass <- iter

  fit
}


spgl_wlsfit <- function(warm, wx, gamma, static) {

  r <- warm$r
  ulam <- warm$ulam
  activeGroup <- warm$activeGroup
  activeGroupIndex <- warm$activeGroupIndex
  ni <- warm$ni
  npass <- warm$npass
  sset <- warm$sset
  eset <- warm$eset
  b0old <- warm$b0
  betaold <- warm$beta
  al0 <- warm$al0
  findlambda <- warm$findlambda
  l <- warm$l
  me <- warm$me

  bn <- static$bn
  bs <- static$bs
  ix <- static$ix
  iy <- static$iy
  nobs <- static$nobs
  nvars <- static$nvars
  pf <- static$pf
  pfl1 <- static$pfl1
  pmax <- static$pmax
  eps <- static$eps
  maxit <- static$maxit
  intr <- static$intr
  alsparse <- static$asparse
  lb <- static$lb
  ub <- static$ub


  if (inherits(wx, "sparseMatrix")) {
    xidx <- as.integer(wx@i + 1)
    xcptr <- as.integer(wx@p + 1)
    xval <- as.double(wx@x)
    nnz <- as.integer(utils::tail(wx@p, 1))

    wls_fit <- dotCall64::.C64(
      "spmat_wsgl",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "integer", "integer",
                    "integer", "double", "double", "double",
                    "integer", "double", "double", "integer", "integer",
                    "double", "double", "integer", "integer", "integer",
                    "integer", "integer", "double", "double", "double",
                    "integer", "integer", "double", "double", "double",
                    "integer", "integer", "integer"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = xval, xidx = xidx, xcptr = xcptr, nnz = nnz,
      # Read / write
      r = as.double(r),
      # Read only
      pf = pf, pfl1 = pfl1, pmax = pmax,
      # Read write
      ulam = ulam,
      # Read only
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Read / write
      b0 = double(1),
      # Write only
      beta = numeric_dc(nvars),
      # Read / write
      activeGroup = activeGroup, activeGroupIndex = activeGroupIndex, ni = ni,
      npass = npass, jerr = 0L,
      # read only
      alsparse = alsparse, lb = lb, ub = ub,
      # Read / write
      sset = sset, eset = eset,
      # read only
      b0old = b0old, betaold = betaold,
      # read / write
      al0 = al0, findlambda = findlambda, l = l, me = me,
      INTENT = c(rep("r", 11), "rw", rep("r", 3), "rw", rep("r", 3), "rw", "w",
                 rep("rw", 5), rep("r", 3), rep("rw", 2), rep("r", 2),
                 rep("rw", 4)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  } else {
    wls_fit <- dotCall64::.C64(
      "wsgl",
      SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                    "integer", "integer", "double", "double", "double", "double",
                    "integer", "double", "double", "integer", "integer",
                    "double", "double", "integer", "integer", "integer",
                    "integer", "integer", "double", "double", "double",
                    "integer", "integer", "double", "double", "double",
                    "integer", "integer", "integer"),
      # Read only
      bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma, nobs = nobs,
      nvars = nvars, x = as.double(wx),
      # Read / write
      r = as.double(r),
      # Read only
      pf = pf, pfl1 = pfl1, pmax = pmax,
      # Read write
      ulam = ulam,
      # Read only
      eps = eps, maxit = maxit, intr = as.integer(intr),
      # Read / write
      b0 = double(1),
      # Write only
      beta = numeric_dc(nvars),
      # Read / write
      activeGroup = activeGroup, activeGroupIndex = activeGroupIndex, ni = ni,
      npass = npass, jerr = 0L,
      # read only
      alsparse = alsparse, lb = lb, ub = ub,
      # Read / write
      sset = sset, eset = eset,
      # read only
      b0old = b0old, betaold = betaold,
      # read / write
      al0 = al0, findlambda = findlambda, l = l, me = me,
      INTENT = c(rep("r", 8), "rw", rep("r", 3), "rw", rep("r", 3), "rw", "w",
                 rep("rw", 5), rep("r", 3), rep("rw", 2), rep("r", 2),
                 rep("rw", 4)),
      NAOK = TRUE,
      PACKAGE = "sparsegl")
  }

  jerr <- wls_fit$jerr
  # if error code > 0, fatal error occurred: stop immediately
  # if error code < 0, non-fatal error occurred: return error code
  if (jerr > 0) {
    if (jerr > 7777) errmsg <- "Unknown error"
    else errmsg <- "Memory allocation bug; contact pkg maintainer"
    rlang::abort(errmsg, call = rlang::caller_env())
  } else if (jerr < 0) {
    return(list(jerr = jerr))
  }

  # wls_fit
  wls_fit[c("r", "ulam", "b0", "beta", "activeGroup", "activeGroupIndex",
            "ni", "npass", "sset", "eset", "al0", "findlambda", "l", "me",
            "jerr")]
}


initializer <- function(x, y, weights, family, intr, has_offset, offset,
                        ulam, flmin) {
  nobs <- nrow(x)
  nvars <- ncol(x)
  if (intr) {
    suppressWarnings({
      tempfit <- stats::glm(y ~ 1, family = family, weights = weights,
                            offset = offset)})
    b0 <- coef(tempfit)[1]
    mu <- tempfit$fitted.values
  } else {
    mu <- family$linkinv(offset)
    b0 <- as.double(0)
  }
  nulldev <- dev_function(y, mu, weights, family)

  # compute upper bound for lambda max
  r <- y - mu

  no_user_lambda <- flmin < 1
  findlambda <- no_user_lambda
  eta <- family$linkfun(mu)

  if (no_user_lambda) {
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    weights <- weights / sum(weights)
    rv <- r / v * m.e * weights
    g <- abs(drop(crossprod(rv, x)))
    lambda_max <- max(g)
  } else {
    lambda_max <- ulam[1]
  }

  cur_lambda <- max(lambda_max / 0.99, 1e-6)

  list(r = r, mu = mu, findlambda = findlambda,
       lambda_max = lambda_max, nulldev = nulldev,
       cur_lambda = cur_lambda, b0 = b0, eta = eta)
}


obj_function <- function(y, mu, group, weights, family,
                         pf, pfl1, asparse, coefs, lambda) {
  ll <- dev_function(y, mu, weights, family) / 2
  grpen <- lambda * sum(pf * grouped_two_norm(coefs, group)) * (1 - asparse)
  sppen <- lambda * sum(pfl1 * abs(coefs)) * asparse
  ll + grpen + sppen
}


dev_function <- function(y, mu, weights, family) {
  sum(family$dev.resids(y, mu, weights))
}


# Helper function to get etas (linear predictions) on the original scale
get_eta <- function(x, xs, beta, b0) {
  beta <- drop(beta)
  # if (length(xs) > 0) beta <- beta / xs
  drop(x %*% beta + b0)
}



#' Create starting values for iterative reweighted least squares
#'
#' This function may be used to create potentially valid starting
#' values for calling [sparsegl()] with a [stats::family()] object.
#' It is not typically necessary to call this function (as it is used
#' internally to create some), but in some cases, especially with custom
#' generalized linear models, it may improve performance.
#'
#' Occasionally, the irls fitting routine may fail with an admonition to
#' create valid starting values.
#'
#' @param nobs Number of observations in the response (or rows in `x`).
#' @param nvars Number of columns in `x`
#' @param b0 Scalar. Initial value for the intercept.
#' @param beta Vector. Initial values for the coefficients. Must be length
#'   `nvars` (or a scalar).
#' @param r Vector. Initial values for the deviance residuals. Must be length
#'   `nobs` (or a scalar).
#'
#' @return List of class `irlsspgl_warmup`
#' @export
#' @importFrom rlang abort
make_irls_warmup <- function(nobs, nvars, b0 = 0, beta = double(nvars),
                             r = double(nobs)) {

  if (any(length(nobs) != 1, length(nvars) != 1, length(b0) != 1))
    cli::cli_abort("All of `nobs`, `nvars`, and `b0` must be scalars.")
  if ((length(r) > 1 && length(r) != nobs) || length(r) == 0) {
    cli::cli_abort(
      "`r` must have length 1 or `nobs` but has length {length(r)}."
    )
  }
  if ((length(beta) > 1 && length(beta) != nvars) || length(beta) == 0) {
    cli::cli_abort(
      "`beta` must have length 1 or `nvars` but has length {length(beta)}."
    )
  }

  structure(list(b0 = b0, beta = beta, r = r), class = "irlsspgl_warmup")
}
