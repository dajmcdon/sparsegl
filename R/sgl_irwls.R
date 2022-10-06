validate_family <- function(family) {
  if (!is.function(family$variance) || !is.function(family$linkinv))
    rlang::abort(
      "'family' argument seems not to be a valid family object. See `?family`.")
}



#' Fit a GLM with sparse group regularization for a path of lambda values
#'
#' Fit a generalized linear model via penalized maximum likelihood for a path of
#' lambda values. Can deal with any GLM family.
#'
#' This organization is based largely off [stats::glm.fit] with some extras
#' to handle the path.
#'
#' Sometimes the sequence is truncated before \code{nlambda} values of lambda
#' have been used. This happens when \code{glmnet.path} detects that the
#' decrease in deviance is marginal (i.e. we are near a saturated fit).
#'
sgl_irwls <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
    flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd, weights, offset = NULL, family = gaussian(),
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
  init <- initilizer(x, y, weights, family, intr, has_offset, offset, pfl1, ulam)
  cur_lambda <- init$cur_lambda
  # work out lambda values
  findlambda <- init$findlambda
  no_user_lambda <- init$findlambda
  # this is supposed to be an upper bound
  lambda_max <- init$lambda_max
  cur_lambda <- init$cur_lambda
  nulldev <- init$nulldev

  if (trace_it == 1) pb <- utils::txtProgressBar(min = 0, max = nlam, style = 3)

  # preallocate space
  b0 <- double(nlam)
  beta <- matrix(0, nvars, nlam)
  dev.ratio <- rep(NA, length = nlam)
  mnl <- min(nlam, 6L)

  static_args <- list(
    nulldev = nulldev, y = y, weights = weights, offset = offset,
    bn = bn, bs = bs, ix = ix, iy = iy, x = x, xs = xs, nobs = nobs,
    nvars = nvars, pf = pf,
    pfl1 = pfl1, dfmax = dfmax, pmax = pmax, flmin = flmin, epx = eps,
    maxit = maxit, vnames = vnames, group = group, intr = intr,
    asparse = asparse, lower_bnd = lower_bnd, upper_bnd = upper_bnd,
    weights = weights, offset = offset, family = family)

  if (is.null(warm))
    warm <- make_irls_warmup(nobs, nvars, b0 = init$b0, r = init$r)
  if (!inherits(warm, "irwlsspgl_warmup")) {
    rlang::abort(
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
  warm <- c(warm, ni = 0L, npass = 0L, me = 0L)

  l <- 0
  while (l <= nlam) {
    warm_fit$findlambda <- findlambda
    if (!findlambda) {
      l <- l + 1
      warm_fit$al0 <- cur_lambda
      warm_fit$ulam <- ulam[l]
      if (trace_it == 2)
        cat("Fitting lambda index", l, ":", ulam[l], fill = TRUE)
    } else {
      # trying to find lambda max, we started too big
      warm_fit$al0 <- cur_lambda
      warm_fit$ulam <- cur_lambda * 0.99
      if (trace_it == 2)
        cat("Trying to find a reasonable starting lambda.", fill = TRUE)
    }
    warm_fit$k <- k

    # here we dispatch to wls
    fit <- irwls_fit(warm_fit, static_args)
    if (trace.it == 1) utils::setTxtProgressBar(pb, l)
    if (fit$jerr != 0) {
      cli::cli_warn(
        "Convergence for {k}th lambda value not reached after maxit =
          {maxit} iterations; solutions for larger lambdas returned.")
      l <- l - 1
      break
    }

    b0[l] <- fit$a0
    beta[, l] <- as.matrix(fit$beta)
    dev.ratio[l] <- fit$dev.ratio

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

  if (standardize) beta <- beta / xs

  # output
  stepnames <- paste0("s", 0:(length(ulam) - 1))
  out <- list()
  out$b0 <- b0
  names(out$b0) <- stepnames
  out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                             dimnames = list(vnames, stepnames))
  out$lambda <- ulam
  out$dev.ratio <- dev.ratio
  out$nulldev <- start_val$nulldev
  out$npasses <- fit$npasses
  out$jerr <- fit$jerr
  out$offset <- has_offset
  out$family <- family
  out$nobs <- nobs
  class(out) <- "irwlsspgl"

  return(out)
}


make_irls_warmup <- function(nobs, nvars, b0 = 0, beta = double(nvars),
                             r = double(nobs)) {

  stopifnot(is.double(b0), is.double(beta), is.double(r))
  stopifnot(length(b0) == 1, length(beta) == nvars, length(r) == nobs)

  structure(list(b0 = b0, beta = beta, r = r), class = "irwlsspgl_warmup")
}

irwls_fit <- function(warm, static) {

  # initialize everything
  variance <- static$family$variance
  linkinv <- static$family$linkinv
  mu.eta <- static$family$mu.eta

  fit <- warm
  nulldev <- static$nulldev
  coefold <- fit$beta   # prev value for coefficients
  intold <- fit$b0    # prev value for intercept
  lambda <- fit$ulam
  eta <- get_eta(static_args$x, static$xs, coefold, intold)
  mu <- linkinv(eta <- eta + static$offset)

  valideta <- static$family$valideta %||% function(eta) TRUE
  validmu <- static$family$validmu %||% function(mu) TRUE


  if (!validmu(mu) || !valideta(eta)) {
    stop("cannot find valid starting values: please specify some",
         call. = FALSE)
  }

  start <- NULL     # current value for coefficients
  start_int <- NULL # current value for intercept
  obj_val_old <- obj_function(
    static$y, mu, static$group, static$weights, static$family,
    static$pf, static$pfl1, static$asparse,
    coefold, fit$ulam
  )

  if (static$trace.it == 2) {
    cat("Warm Start Objective:", obj_val_old, fill = TRUE)
  }
  conv <- FALSE      # converged?

  # IRLS loop
  for (iter in seq_len(static$maxit)) {
    # some checks for NAs/zeros
    varmu <- variance(mu)
    if (anyNA(varmu)) stop("NAs in V(mu)")
    if (any(varmu == 0)) stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (anyNA(mu.eta.val)) stop("NAs in d(mu)/d(eta)")

    # d ell / d beta = X'W (z - mu) / mu.eta.val (theoretically)
    #                = t(wx) %*% r,  (code)
    # can use directly in wls for strong rule / kkt

    # compute working response and weights
    z <- (eta - static$offset) + (static$y - mu) / mu.eta.val
    w <- sqrt((static$weights * mu.eta.val^2) / variance(mu))
    r <- w * z
    fit$r <- r

    # we updated w, so we need to update gamma
    wx <- Matrix::Diagonal(w) %*% static$x
    # because we use eig(xwx) = svd(w^{1/2}x)
    gamma <- calc_gamma(wx, static$ix, static$iy, static$bn)

    # do WLS with warmstart from previous iteration, dispatch to FORTRAN wls
    # things that need to change:
    #   wx, r, ulam, b0, beta, activeGroup, activeGroupIndex, ni, npass,
    #   jerr, sset, b0old, betaold, al0, findlambda, l, me
    #
    fit <- spgl_wls_fit(fit, wx, gamma, static)
    # fit <- spgl_wls_fit(wx, r, w, lambda, alpha, intercept,
    #                     thresh = thresh, maxit = maxit, penalty.factor = vp,
    #                     exclude = exclude, lower.limits = lower.limits,
    #                     upper.limits = upper.limits, warm = fit,
    #                     from.glmnet.fit = TRUE, save.fit = TRUE)
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
    if (trace_it == 2) cat("Iteration", iter, "Objective:", obj_val, fill = TRUE)

    boundary <- FALSE
    halved <- FALSE  # did we have to halve the step size?
    # if objective function is not finite, keep halving the stepsize until it is finite
    if (!is.finite(obj_val) || obj_val > 9.9e30) {
      warning("Infinite objective function!", call. = FALSE)
      if (is.null(coefold) || is.null(intold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      warning("step size truncated due to divergence", call. = FALSE)
      ii <- 1
      while (!is.finite(obj_val) || obj_val > 9.9e30) {
        if (ii > maxit)
          stop("inner loop 1; cannot correct step size", call. = FALSE)
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
      warning("Invalid eta / mu!", call. = FALSE)
      if (is.null(coefold) || is.null(intold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      warning("step size truncated: out of bounds", call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > maxit)
          stop("inner loop 2; cannot correct step size", call. = FALSE)
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
      if (trace_it == 2) cat("Iteration", iter, " Halved step 2, Objective:", obj_val, fill = TRUE)
    }
    # extra halving step if objective function value actually increased
    if (obj_val > obj_val_old + 1e-7) {
      ii <- 1
      while (obj_val > obj_val_old + 1e-7) {
        if (ii > maxit_irls)
          stop("inner loop 3; cannot correct step size", call. = FALSE)
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
    if (abs(obj_val - obj_val_old) / (0.1 + abs(obj_val)) < eps) {
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
    warning("sparsgl_irls: algorithm did not converge", call. = FALSE)
  if (boundary)
    warning("sparsgl_irls: algorithm stopped at boundary value", call. = FALSE)

  # some extra warnings, printed only if trace_it == 2
  if (trace_it == 2) {
    tiny <- 10 * .Machine$double.eps
    if ((family$family == "binomial") && (any(mu > 1 - tiny) || any(mu < tiny)))
      warning("sparsgl_irls: fitted probabilities numerically 0 or 1 occurred",
              call. = FALSE)
    if ((family$family == "poisson") && (any(mu < eps)))
      warning("sparsegl_irls: fitted rates numerically 0 occurred",
              call. = FALSE)
  }

  # prepare output object
  if (save.fit == FALSE) fit$warm_fit <- NULL
  fit$offset <- has_offset
  fit$nulldev <- static$nulldev
  fit$dev.ratio <- 1 - dev_function(static$y, mu, static$weights,
                                    static$family) / fit$nulldev
  fit$family <- family
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val

  class(fit) <- c("irwls_fit", "sparsegl")
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
  alsparse <- static$alsparse
  lb <- static$lb
  ub <- static$ub


  if (inherits(wx, "sparseMatrix")) {
    rlang::abort("not currently implemented")
  } else {
    # fortran signature wsgl (bn,bs,ix,iy,gam, nobs,nvars,x,r,pf,pfl1,
    #        pmax,ulam,eps,maxit,intr, b0,beta,activeGroup,activeGroupIndex,ni,&
    #        npass,jerr,alsparse,lb,ub, sset,eset,b0old,betaold,al0, findlambda,l,me)
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
      nvars = nvars, x = as.double(wx), r = as.double(r), pf = pf,
      pfl1 = pfl1, pmax = pmax,
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Read / write
      b0 = double(1),
      # Write only
      beta = numeric_dc(nvars),
      # Read / write
      activeGroup = activeGroup, activeGroupIndex = activeGroupIndex, ni = ni,
      npass = npass, jerr = 0L,
      # read only
      alsparse = asparse, lb = lb, ub = ub,
      # Read / write
      sset = sset, eset = eset,
      # read only
      b0old = b0old, betaold = betaold,
      # read / write
      al0 = al0, findlambda = findlambda, l = l, me = me,
      INTENT = c(rep("r", 16), "rw", "w", rep("rw", 5),
                 rep("r", 3), rep("rw", 2), rep("r", 2), rep("rw", 4)),
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

  wls_fit[c("r", "ulam", "b0", "beta", "activeGroup", "activeGroupIndex",
            "ni", "npass", "sset", "eset", "al0", "findlambda", "l", "me")]
}

#' Get null deviance, starting mu and lambda max
#'
#' Return the null deviance, starting mu and lambda max values for
#' initialization. For internal use only.
#'
#' This function is called by \code{glmnet.path} for null deviance, starting mu
#' and lambda max values. It is also called by \code{glmnet.fit} when used
#' without warmstart, but they only use the null deviance and starting mu values.
#'
#' When \code{x} is not sparse, it is expected to already by centered and scaled.
#' When \code{x} is sparse, the function will get its attributes \code{xm} and
#' \code{xs} for its centering and scaling factors.
#'
#' Note that whether \code{x} is centered & scaled or not, the values of \code{mu}
#' and \code{nulldev} don't change. However, the value of \code{lambda_max} does
#' change, and we need \code{xm} and \code{xs} to get the correct value.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param y Quantitative response variable.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' (See \code{\link[stats:family]{family}} for details on family functions.)
#' @param intercept Does the model we are fitting have an intercept term or not?
#' @param is.offset Is the model being fit with an offset or not?
#' @param offset Offset for the model. If \code{is.offset=FALSE}, this should be
#' a zero vector of the same length as \code{y}.
#' @param exclude Indices of variables to be excluded from the model.
#' @param vp Separate penalty factors can be applied to each coefficient.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.

initilizer <- function(x, y, weights, family, intr, has_offset, offset, pfl1,
                       ulam) {
  nobs <- nrow(x)
  nvars <- ncol(x)
  if (intr) {
    if (has_offset) {
      suppressWarnings({
        tempfit <- glm(y ~ 1, family = family, weights = weights,
                       offset = offset)})
      b0 <- coef(tempfit)[1]
      mu <- tempfit$fitted.values
    } else {
      mu <- rep(weighted.mean(y, weights), times = nobs)
      b0 <- mu[1]
    }
  } else {
    mu <- family$linkinv(offset)
    b0 <- 0
  }
  nulldev <- dev_function(y, mu, weights, family)

  # compute upper bound for lambda max
  # can possibly undershoot if any pfl1 < 1e-6, should warn.
  r <- y - mu
  eta <- family$linkfun(mu)
  v <- family$variance(mu)
  m.e <- family$mu.eta(eta)
  weights <- weights / sum(weights)
  rv <- r / v * m.e * weights
  g <- abs(drop(crossprod(rv, x)))
  lambda_max <- max(g / pmax(pfl1, 1e-6))

  no_user_lambda <- (length(ulam) == 1) & (ulam == 0)
  findlambda <- no_user_lambda
  lambda_max <- if (no_user_lambda) lambda_max else ulam[1]
  cur_lambda <- lambda_max / 0.99

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


#' Helper function to get etas (linear predictions) on the original scale
get_eta <- function(x, xs, beta, b0) {
  beta <- drop(beta)
  if (!is.null(xs)) beta <- beta / xs
  drop(x %*% beta + b0)
}
