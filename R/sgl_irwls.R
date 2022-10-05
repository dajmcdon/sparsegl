validate_family <- function(family) {
  if (!is.function(family$variance) || !is.function(family$linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
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
  fit <- NULL
  mnl <- min(nlam, 6L)


  k <- 0
  while (k <= nlam) {
    init$findlambda <- findlambda
    if (!findlambda) {
      k <- k + 1
      init$prev_lambda <- cur_lambda
      init$cur_lambda <- ulam[k]
      if (trace_it == 2)
        cat("Fitting lambda index", k, ":", ulam[k], fill = TRUE)
    } else {
      # trying to find lambda max, we started too big
      init$prev_lambda <- cur_lambda
      init$cur_lambda <- cur_lambda * 0.99
      if (trace_it == 2)
        cat("Trying to find a reasonable starting lambda.", fill = TRUE)
    }
    init$k <- k



    # here we dispatch to wls
    fit <- irwls_fit(
      bn, bs = bs, ix, iy, gamma, nobs, nvars = nvars,
      x, y, pf = pf, pfl1 = pfl1,
      dfmax, pmax, cur_lambda,
      eps, maxit, intr, b0, beta,
      activeGroup = integer_dc(pmax), nbeta = integer_dc(nlam),
      alam = numeric_dc(nlam), npass = integer_dc(1),
      jerr = integer_dc(1), mse = numeric_dc(nlam),
      # read only
      alsparse = asparse, lb = lower_bnd, ub = upper_bnd,
      bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
      flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
      lower_bnd, upper_bnd, weights, offset, family,
      trace_it = 1,

      x, y, weights / sum(weights), cur_lambda, alpha = alpha, offset = offset,
      family = family, intercept = intercept, thresh = thresh,
      maxit = maxit, penalty.factor = vp, exclude = exclude,
      lower.limits = lower.limits, upper.limits = upper.limits,
      warm = fit, from.glmnet.path = TRUE, save.fit = TRUE,
      trace.it = trace.it
    )
    if (trace.it == 1) utils::setTxtProgressBar(pb, k)
    if (fit$jerr != 0) {
      cli::cli_warn(
        "Convergence for {k}th lambda value not reached after maxit =
          {maxit} iterations; solutions for larger lambdas returned.")
      k <- k - 1
      break
    }

    if (k == 0) {}

    b0[k] <- fit$a0
    beta[, k] <- as.matrix(fit$beta)
    dev.ratio[k] <- fit$dev.ratio

    # early stopping if dev.ratio almost 1 or no improvement
    if (k >= mnl && no_user_lambda) {
      if (dev.ratio[k] > 1 - 1e-4) break
      else if (k > 1) {
        if (family$family == "gaussian") {
          if (dev.ratio[k] - dev.ratio[k - 1] < 1e-5 * dev.ratio[k]) break
        } else if (family$family == "poisson") {
          if (dev.ratio[k] - dev.ratio[k - mnl + 1] < 1e-4 * dev.ratio[k]) break
        } else if (dev.ratio[k] - dev.ratio[k - 1] < 1e-5) break
      }
    }
  }
  if (trace_it == 1) {
    utils::setTxtProgressBar(pb, nlam)
    cat("", fill = TRUE)
  }

  # truncate if we quit early
  if (k < nlam) {
    b0 <- b0[1:k]
    beta <- beta[, 1:k, drop = FALSE]
    dev.ratio <- dev.ratio[1:k]
    ulam <- ulam[1:k]
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


irwls_fit <- function(init,
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlam,
    flmin, eps, maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd, weights, offset, family, trace_it, warm = NULL) {

  # initialize everything
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta

  irls_warmup <- function(warm, ...) {
    UseMethod("irls_warmup")
  }
  irls_warmup.default <- function(warm, init, nvars, ...) {
    mu <- init$mu
    beta <- rep(0, nvars)
    b0 <- init$b0
    r <- init$r
    eta <- init$eta
    structure(list(beta = beta, b0 = b0, r = r, eta = eta, mu = mu),
              class = "irls_warm")
  }
  irls_warmup.irls_warm <- function(warm, init, ...) {
    NextMethod()
    beta <- warm$beta
    b0 <- warm$b0
    r <- warm$r
    structure(list(beta = beta, b0 = b0, r = r, eta = eta, mu = mu),
              class = "irls_warm")
  }


  if (is.null(warm)) {
    start_val <- initilizer(x, y, weights, family, intr, has_offset, offset, pfl1)
    nulldev <- start_val$nulldev
    mu <- start_val$mu
    fit <- NULL
    coefold <- rep(0, nvars)   # initial coefs = 0
    eta <- family$linkfun(mu)
    intold <- (eta - offset)[1]
  } else {
    fit <- warm
    nulldev <- fit$nulldev
    coefold <- fit$warm_fit$beta   # prev value for coefficients
    intold <- fit$warm_fit$b0    # prev value for intercept
    eta <- get_eta(x, xs, coefold, intold)
    mu <- linkinv(eta <- eta + offset)
  }

  valideta <- family$valideta %||% function(eta) TRUE
  validmu <- family$validmu %||% function(mu) TRUE


  if (!validmu(mu) || !valideta(eta)) {
    stop("cannot find valid starting values: please specify some",
         call. = FALSE)
  }

  start <- NULL     # current value for coefficients
  start_int <- NULL # current value for intercept
  obj_val_old <- obj_function(y, mu, weights, family, pf, pfl1, asparse, coefold, lambda)

  if (trace.it == 2) {
    cat("Warm Start Objective:", obj_val_old, fill = TRUE)
  }
  conv <- FALSE      # converged?

  # IRLS loop
  for (iter in seq_len(maxit_irls)) {
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
    z <- (eta - offset) + (y - mu) / mu.eta.val
    w <- sqrt((weights * mu.eta.val^2) / variance(mu))
    r <- z * w
    if (!is.null(fit)) fit$warm_fit$r <- r

    # we updated w, so we need to update gamma
    wx <- Matrix::Diagonal(w) %*% x
    gamma <- calc_gamma(wx, ix, iy, bn) # because we use eig(xwx) = svd(w^{1/2}x)

    # do WLS with warmstart from previous iteration, dispatch to FORTRAN wls
    # things that need to change:
    #   wx, r, ulam, b0, beta, activeGroup, activeGroupIndex, ni, npass,
    #   jerr, sset, b0old, betaold, al0, findlambda, l, me
    #
    fit <- spgl_wls_fit(bn, bs, ix, iy, nobs, nvars, wx, r, pf, pfl1, dfmax,
                        pmax, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
                        lower_bnd, upper_bnd, weights, offset, family,
                        cur_lambda, prev_lambda, k, findlambda, trace_it, warm = fit)
    # fit <- spgl_wls_fit(wx, r, w, lambda, alpha, intercept,
    #                     thresh = thresh, maxit = maxit, penalty.factor = vp,
    #                     exclude = exclude, lower.limits = lower.limits,
    #                     upper.limits = upper.limits, warm = fit,
    #                     from.glmnet.fit = TRUE, save.fit = TRUE)
    if (fit$jerr != 0) return(list(jerr = fit$jerr))

    # update coefficients, eta, mu and obj_val
    start <- fit$warm_fit$beta
    start_int <- fit$warm_fit$b0
    eta <- get_eta(x, xs, start, start_int)
    mu <- linkinv(eta <- eta + offset)
    obj_val <- obj_function(y, mu, weights, family, pf, pfl1, asparse, start, lambda)
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
      while (!is.finite(obj_val) || obj_val > control$big) {
        if (ii > maxit_irls)
          stop("inner loop 1; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold) / 2
        start_int <- (start_int + intold) / 2
        eta <- get_eta(x, xs, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, pf, pfl1, asparse, start, lambda)
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
        if (ii > maxit_irls)
          stop("inner loop 2; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold) / 2
        start_int <- (start_int + intold) / 2
        eta <- get_eta(x, xs, start, start_int)
        mu <- linkinv(eta <- eta + offset)
      }
      boundary <- TRUE
      halved <- TRUE
      obj_val <- obj_function(y, mu, weights, family, pf, pfl1, asparse, start, lambda)
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
        eta <- get_eta(x, xs, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, pf, pfl1, asparse, start, lambda)
        if (trace_it == 2) cat("Iteration", iter, " Halved step 3, Objective:",
                               obj_val, fill = TRUE)
      }
      halved <- TRUE
    }

    # if we did any halving, we have to update the coefficients, intercept
    # and weighted residual in the warm_fit object
    if (halved) {
      fit$warm_fit$beta <- start
      fit$warm_fit$b0 <- start_int
      fit$warm_fit$r <- w * (z - eta)
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
  fit$call <- this.call
  fit$offset <- has_offset
  fit$nulldev <- init$nulldev
  fit$dev.ratio <- 1 - dev_function(y, mu, weights, family) / nulldev
  fit$family <- family
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val

  class(fit) <- c("irwls_fit", "sparsegl")
  fit
}

sgl_warmup <- function(object, ...) {

}


spgl_wlsfit <- function(warmobject, ...) {
  UseMethod("spgl_wls_fit")
}

spgl_wlsfit.default <- function(warmobject, ...) {
  rlang::abort("Fitting methods exist only for class `sgl_warmup`.")
}

new_sgl_warmup <- function(
    r, cur_lambda, b0, beta, activeGroup, activeGroupIndex, ni,
    npass, sset, eset, b0old, betaold, prev_lambda, findlambda, me
) {
  stopifnot(is.double(r), is.double(cur_lambda), is.double(b0), is.double(beta),
            is.integer(activeGroup), is.integer(activeGroupIndex),
            is.integer(ni), is.integer(npass), is.integer(sset),
            is.integer(eset), is.double(b0old), is.double(betaold),
            is.double(prev_lambda), is.integer(findlambda), is.integer(me))

  stopifnot(max(length(cur_lambda), length(b0), length(ni), length(npass),
                length(b0old), length(prev_lambda), length(findlambda),
                length(me)) == 1L)

  structure(list(
    r = r, cur_lambda = cur_lambda, b0 = b0, beta = beta,
    activeGroup = activeGroup, activeGroupIndex = activeGroupIndex, ni = ni,
    npass = npass, sset = sset, eset = eset, b0old = b0old,
    betaold = betaold, prev_lambda = prev_lambda, findlambda = findlambda,
    me = me
  ), class = "sgl_warmup")
}

# things that change over iterations
#   gam, r, ulam, b0, beta, activeGroup, activeGroupIndex, ni, npass,
#   jerr, sset, eset, b0old, betaold, al0, findlambda, l, me
# fortran signature wsgl (bn,bs,ix,iy,gam,nobs,nvars,x,r,pf,pfl1,pmax,&
#        ulam,eps,maxit,intr,b0,beta,activeGroup,activeGroupIndex,ni,&
#        npass,jerr,alsparse,lb,ub,sset,eset,b0old,betaold,al0,findlambda,l,me)
spgl_wlsfit.sgl_warmup <- function(
    warmobject, bn, bs, ix, iy, gam, nobs, nvars, wx, pf, pfl1, pmax,
    eps, maxit, intr, asparse, lb, ub, l) {
  rlang::check_dots_empty()


  r <- warmobject$r
  ulam <- warmobject$cur_lambda
  b0 <- warmobject$b0
  beta <- warmobject$beta
  activeGroup <- warmobject$activeGroup
  activeGroupIndex <- warmobject$activeGroupIndex
  ni <- warmobject$ni
  npass <- warmobject$npass
  sset <- warmobject$sset
  eset <- warmobject$eset
  b0old <- warmobject$b0old
  betaold <- warmobject$betaold
  al0 <- warmobject$al0
  findlambda <- warmobject$findlambda
  me <- warmobject$me

  if (inherits(wx, "sparseMatrix")) {

    # wls_fit <- spwls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
    #                      x=x,xm=xm,xs=xs,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
    #                      maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
    #                      nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
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
      pfl1 = pfl1,
      # Read / write
      pmax = pmax,
      # Read only
      ulam = ulam, eps = eps, maxit = maxit, intr = as.integer(intr),
      # Read / write
      b0 = b0,
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
      INTENT = c(rep("r", 11), "rw", rep("r", 4), "rw", "w", rep("rw", 5),
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

  # r, cur_lambda, b0, beta, activeGroup, activeGroupIndex, ni,
  # npass, sset, eset, b0old, betaold, prev_lambda, findlambda, me
  warm_fit <- with(
    wls_fit,
    new_sgl_warmup(
      r, cur_lambda = ulam, b0, beta, activeGroup, activeGroupIndex, ni, npass, sset,
      eset, b0old, betaold, prev_lambda = al0, findlambda, me
    ))

  return(warm_fit)
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
                       offset = offset)
      })
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

  # compute upper bound for lambda max, can possibly undershoot if pfl1 == 0
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

#' Sparse group lasso objective function value
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients (excluding intercept).
#' @param vp Penalty factors for each of the coefficients.
obj_function <- function(y, mu, weights, family,
                         pf, pfl1, asparse, coefs, lambda) {
  ll <- dev_function(y, mu, weights, family) / 2
  grpen <- lambda * sum(pf * grouped_two_norm(coefs)) * (1 - asparse)
  sppen <- lambda * sum(pfl1 * abs(coefs)) * asparse
  ll + grpen + sppen
}

#' Elastic net deviance value
#'
#' Returns the elastic net deviance value.
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
dev_function <- function(y, mu, weights, family) {
  sum(family$dev.resids(y, mu, weights))
}


#' Get predictions from a \code{glmnetfit} fit object
#'
#' Gives fitted values, linear predictors, coefficients and number of non-zero
#' coefficients from a fitted \code{glmnetfit} object.
#'
#' @param object Fitted "glmnetfit" object.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix. This argument is not used for \code{type =
#' c("coefficients","nonzero")}.
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#' required. Default is the entire sequence used to create the model.
#' @param type Type of prediction required. Type "link" gives the linear
#' predictors (eta scale); Type "response" gives the fitted values (mu scale).
#' Type "coefficients" computes the coefficients at the requested values for s.
#' Type "nonzero" returns a list of the indices of the nonzero coefficients for
#' each value of s.
#' @param exact This argument is relevant only when predictions are made at values
#' of \code{s} (lambda) \emph{different} from those used in the fitting of the
#' original model. If \code{exact=FALSE} (default), then the predict function
#' uses linear interpolation to make predictions for values of \code{s} (lambda)
#' that do not coincide with those used in the fitting algorithm. While this is
#' often a good approximation, it can sometimes be a bit coarse. With
#' \code{exact=TRUE}, these different values of \code{s} are merged (and sorted)
#' with \code{object$lambda}, and the model is refit before predictions are made.
#' In this case, it is required to supply the original data x= and y= as additional
#' named arguments to predict() or coef(). The workhorse \code{predict.glmnet()}
#' needs to update the model, and so needs the data used to create it. The same
#' is true of weights, offset, penalty.factor, lower.limits, upper.limits if
#' these were used in the original call. Failure to do so will result in an error.
#' @param newoffset If an offset is used in the fit, then one must be supplied for
#' making predictions (except for type="coefficients" or type="nonzero").
#' @param ... This is the mechanism for passing arguments like \code{x=} when
#' \code{exact=TRUE}; see \code{exact} argument.
#'
#' @return The object returned depends on type.
#'
#' @method predict glmnetfit
#' @export
predict.glmnetfit <- function(object, newx, s = NULL,
                              type = c("link", "response", "coefficients", "nonzero"),
                              exact = FALSE, newoffset, ...) {
  type = match.arg(type)
  nfit <- NextMethod("predict")
  if (type == "response") {
    object$family$linkinv(nfit)
  } else {
    nfit
  }
}

#' Helper function to get etas (linear predictions)
get_eta <- function(x, xs, beta, b0) {
  beta <- drop(beta)
  if (!is.null(xs)) beta <- beta / xs
  drop(x %*% beta + b0)
}

#' Helper function to compute weighted mean and standard deviation
#'
#' Helper function to compute weighted mean and standard deviation.
#' Deals gracefully whether x is sparse matrix or not.
#'
#' @param x Observation matrix.
#' @param weights Optional weight vector.
#'
#' @return A list with components.
#' \item{mean}{vector of weighted means of columns of x}
#' \item{sd}{vector of weighted standard deviations of columns of x}
weighted_mean_sd <- function(x, weights=rep(1,nrow(x))){
  weights <- weights/sum(weights)
  xm <- drop(t(weights)%*%x)
  xv <- drop(t(weights)%*%scale(x,xm,FALSE)^2)
  xv[xv < 10*.Machine$double.eps] <- 0
  list(mean = xm, sd = sqrt(xv))
}
