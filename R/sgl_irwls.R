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
    trace_it = 1) {

  validate_family(family)

  ## Init the family just to check that it works
  if (is.null(weights)) weights <- rep(0, nobs)
  etastart <- 0
  mustart <- NULL
  start <- NULL
  eval(family$initialize)
  y <- drop(y)
  has_offset <- !is.null(offset)
  if (!has_offset) offset <- as.double(y * 0)

  # standardize x if necessary
  is_sparse <- FALSE
  if (inherits(x,"sparseMatrix")) {
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
  # get null deviance and lambda max
  start_val <- initilizer(x, y, weights, family, inter, has_offset, offset, pfl1)

  # work out lambda values
  no_user_lambda <- (length(ulam) == 1) & (ulam == 0) # user didn't provide lambda
  if (no_user_lambda) lambda_max <- start_val$lambda_max # this is supposed to be an upper bound
  if (trace_it == 1) pb <- utils::txtProgressBar(min = 0, max = nlam, style = 3)
  cur_lambda <- lambda_max / 0.99

  b0 <- rep(NA, length = nlam)
  beta <- matrix(0, nvars, nlam)
  dev.ratio <- rep(NA, length = nlam)
  fit <- NULL
  mnl <- min(nlam, 6L)


  k <- 0
  first_pass <- TRUE
  while (k <= nlam) {
    if (k > 1 || !no_user_lambda) {
      k <- k + 1
      prev_lambda <- cur_lambda
      cur_lambda <- ulam[k]
      if (trace_it == 2)
        cat("Fitting lambda index", k, ":", ulam[k], fill = TRUE)
    } else {
      # trying to find lambda max, we started too big
      prev_lambda <- cur_lambda
      cur_lambda <- cur_lambda * 0.99
    }


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


irwls_fit <- function(
    x, y, weights, lambda, alpha = 1.0,
    offset = rep(0, nobs), family = gaussian(),
    intercept = TRUE, thresh = 1e-10, maxit = 100000,
    penalty.factor = rep(1.0, nvars), exclude = c(), lower.limits = -Inf,
    upper.limits = Inf, warm = NULL, search_lambda_max = FALSE,
    save.fit = FALSE, trace.it = 0
) {

  # initialize everything
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta

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
    eta <- get_eta(x, coefold, intold)
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
  for (iter in 1L:maxit_irls) {
    # some checks for NAs/zeros
    varmu <- variance(mu)
    if (anyNA(varmu)) stop("NAs in V(mu)")
    if (any(varmu == 0)) stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (anyNA(mu.eta.val)) stop("NAs in d(mu)/d(eta)")

    # d ell / d beta = (y - mu) / mu.eta.val, can use for strong rule / kkt

    # compute working response and weights
    z <- (eta - offset) + (y - mu) / mu.eta.val
    w <- (weights * mu.eta.val^2) / variance(mu)

    if (!is.null(fit)) fit$warm_fit$r <- w * (z - eta + offset)


    # do WLS with warmstart from previous iteration, dispatch to FORTRAN wls
    fit <- spgl_wls_fit(x, z, w, lambda, alpha, intercept,
                        thresh = thresh, maxit = maxit, penalty.factor = vp,
                        exclude = exclude, lower.limits = lower.limits,
                        upper.limits = upper.limits, warm = fit,
                        from.glmnet.fit = TRUE, save.fit = TRUE)
    if (fit$jerr != 0) return(list(jerr = fit$jerr))

    # update coefficients, eta, mu and obj_val
    start <- fit$warm_fit$beta
    start_int <- fit$warm_fit$b0
    eta <- drop(x %*% start) + start_int
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
        eta <- drop(x %*% start) + start_int
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
        eta <- drop(x %*% start) + start_int
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
        eta <- drop(x %*% start) + start_int
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
  fit$nulldev <- nulldev
  fit$dev.ratio <- 1 - dev_function(y, mu, weights, family) / nulldev
  fit$family <- family
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val

  class(fit) <- c("irwls_fit", "sparsegl")
  fit
}

#' Solve weighted least squares (WLS) problem for a single lambda value
#'
#' Solves the weighted least squares (WLS) problem for a single lambda value. Internal
#' function that users should not call directly.
#'
#' WARNING: Users should not call \code{elnet.fit} directly. Higher-level functions
#' in this package call \code{elnet.fit} as a subroutine. If a warm start object
#' is provided, some of the other arguments in the function may be overriden.
#'
#' \code{elnet.fit} is essentially a wrapper around a C++ subroutine which
#' minimizes
#'
#' \deqn{1/2 \sum w_i (y_i - X_i^T \beta)^2 + \sum \lambda \gamma_j
#' [(1-\alpha)/2 \beta^2+\alpha|\beta|],}
#'
#' over \eqn{\beta}, where \eqn{\gamma_j} is the relative penalty factor on the
#' jth variable. If \code{intercept = TRUE}, then the term in the first sum is
#' \eqn{w_i (y_i - \beta_0 - X_i^T \beta)^2}, and we are minimizing over both
#' \eqn{\beta_0} and \eqn{\beta}.
#'
#' None of the inputs are standardized except for \code{penalty.factor}, which
#' is standardized so that they sum up to \code{nvars}.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed that any standardization needed
#' has already been done.
#' @param y Quantitative response variable.
#' @param weights Observation weights. \code{elnet.fit} does NOT standardize
#' these weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-7}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' (If a warm start object is provided, the number of passes the warm start object
#' performed is included.)
#' @param penalty.factor Separate penalty factors can be applied to each
#' coefficient. This is a number that multiplies \code{lambda} to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. Default is 1 for all variables (and
#' implicitly infinity for variables listed in exclude). Note: the penalty
#' factors are internally rescaled to sum to \code{nvars}.
#' @param exclude Indices of variables to be excluded from the model. Default is
#' none. Equivalent to an infinite penalty factor.
#' @param lower.limits Vector of lower limits for each coefficient; default
#' \code{-Inf}. Each of these must be non-positive. Can be presented as a single
#' value (which will then be replicated), else a vector of length \code{nvars}.
#' @param upper.limits Vector of upper limits for each coefficient; default
#' \code{Inf}. See \code{lower.limits}.
#' @param warm Either a \code{glmnetfit} object or a list (with names \code{beta}
#' and \code{a0} containing coefficients and intercept respectively) which can
#' be used as a warm start. Default is \code{NULL}, indicating no warm start.
#' For internal use only.
#' @param from.glmnet.fit Was \code{elnet.fit()} called from \code{glmnet.fit()}?
#' Default is FALSE.This has implications for computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is FALSE.
#'
#' @return An object with class "glmnetfit" and "glmnet". The list returned has
#' the same keys as that of a \code{glmnet} object, except that it might have an
#' additional \code{warm_fit} key.
#' \item{a0}{Intercept value.}
#' \item{beta}{A \code{nvars x 1} matrix of coefficients, stored in sparse matrix
#' format.}
#' \item{df}{The number of nonzero coefficients.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{Lambda value used.}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{Always FALSE, since offsets do not appear in the WLS problem.
#' Included for compability with glmnet output.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If \code{save.fit=TRUE}, output of C++ routine, used for
#' warm starts. For internal use only.}
#'
elnet.fit <- function(x, y, weights, lambda, alpha = 1.0, intercept = TRUE,
                      thresh = 1e-7, maxit = 100000,
                      penalty.factor = rep(1.0, nvars), exclude = c(),
                      lower.limits = -Inf, upper.limits = Inf, warm = NULL,
                      from.glmnet.fit = FALSE, save.fit = FALSE) {
  this.call <- match.call()
  internal.parms <- glmnet.control()

  # compute null deviance
  ybar <- weighted.mean(y, weights)
  nulldev <- sum(weights * (y - ybar)^2)

  # if class "glmnetfit" warmstart object provided, pull whatever we want out of it
  # else, prepare arguments, then check if coefs provided as warmstart
  # (if only coefs are given as warmstart, we prepare the other arguments
  # as if no warmstart was provided)
  if (!is.null(warm) && "glmnetfit" %in% class(warm)) {
    warm <- warm$warm_fit
    if (class(warm) != "warmfit") stop("Invalid warm start object")

    a <- warm$a
    aint <- warm$aint
    alm0 <- warm$almc
    cl <- warm$cl
    g <- warm$g
    ia <- warm$ia
    iy <- warm$iy
    iz <- warm$iz
    ju <- warm$ju
    m <- warm$m
    mm <- warm$mm
    nino <- warm$nino
    nobs <- warm$no
    nvars <- warm$ni
    nlp <- warm$nlp
    nx <- warm$nx
    r <- warm$r
    rsqc <- warm$rsqc
    xv <- warm$xv
    vp <- warm$vp
  } else {
    nobs <- as.integer(nrow(x))
    nvars <- as.integer(ncol(x))

    # if calling from glmnet.fit(), we do not need to check on exclude
    # and penalty.factor arguments as they have been prepared by glmnet.fit()
    # Also exclude will include variance 0 columns
    if (!from.glmnet.fit) {
      # check and standardize penalty factors (to sum to nvars)
      if(any(penalty.factor == Inf)) {
        exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
        exclude = sort(unique(exclude))
      }
      if(length(exclude) > 0) {
        jd = match(exclude, seq(nvars), 0)
        if(!all(jd > 0)) stop ("Some excluded variables out of range")
        penalty.factor[jd] = 1 # ow can change lambda sequence
      }
      vp = pmax(0, penalty.factor)
      vp = as.double(vp * nvars / sum(vp))
    } else {
      vp <- as.double(penalty.factor)
    }
    # compute ju
    # assume that there are no constant variables
    ju <- rep(1, nvars)
    ju[exclude] <- 0
    ju <- as.integer(ju)

    # compute cl from lower.limits and upper.limits
    lower.limits[lower.limits == -Inf] <- -internal.parms$big
    upper.limits[upper.limits == Inf] <- internal.parms$big
    if (length(lower.limits) < nvars)
      lower.limits = rep(lower.limits, nvars) else
        lower.limits = lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars)
      upper.limits = rep(upper.limits, nvars) else
        upper.limits = upper.limits[seq(nvars)]
    cl <- rbind(lower.limits, upper.limits)
    storage.mode(cl) = "double"

    nx <- as.integer(nvars)

    a <- double(nvars)
    aint <- double(1)
    alm0 <- double(1)
    g <- double(nvars)
    ia <- integer(nx)
    iy <- integer(nvars)
    iz <- integer(1)
    m <- as.integer(1)
    mm <- integer(nvars)
    nino <- integer(1)
    nlp <- integer(1)
    r <- weights * y
    rsqc <- double(1)
    xv <- double(nvars)

    # check if coefs were provided as warmstart: if so, use them
    if (!is.null(warm)) {
      if ("list" %in% class(warm) && "a0" %in% names(warm) &&
          "beta" %in% names(warm)) {
        a <- as.double(warm$beta)
        aint <- as.double(warm$a0)
        mu <- drop(x %*% a + aint)
        r <- weights * (y - mu)
        rsqc <- 1 - sum(weights * (y - mu)^2) / nulldev
      } else {
        stop("Invalid warm start object")
      }
    }
  }

  # for the parameters here, we are overriding the values provided by the
  # warmstart object
  alpha <- as.double(alpha)
  almc <- as.double(lambda)
  intr <- as.integer(intercept)
  jerr <- integer(1)
  maxit <- as.integer(maxit)
  thr <- as.double(thresh)
  v <- as.double(weights)

  a.new <- a
  a.new[0] <- a.new[0] # induce a copy

  # take out components of x and run C++ subroutine
  if (inherits(x, "sparseMatrix")) {
    xm <- as.double(attr(x, "xm"))
    xs <- as.double(attr(x, "xs"))
    wls_fit <- spwls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                         x=x,xm=xm,xs=xs,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                         maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                         nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
  } else {
    wls_fit <- wls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                       x=x,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                       maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                       nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
  }

  # if error code > 0, fatal error occurred: stop immediately
  # if error code < 0, non-fatal error occurred: return error code
  if (wls_fit$jerr > 0) {
    errmsg <- jerr.glmnetfit(wls_fit$jerr, maxit)
    stop(errmsg$msg, call. = FALSE)
  } else if (wls_fit$jerr < 0)
    return(list(jerr = wls_fit$jerr))
  warm_fit <- wls_fit[c("almc", "r", "xv", "ju", "vp", "cl", "nx",
                        "a", "aint", "g", "ia", "iy", "iz", "mm", "nino",
                        "rsqc", "nlp")]
  warm_fit[['m']] <- m
  warm_fit[['no']] <- nobs
  warm_fit[['ni']] <- nvars
  class(warm_fit) <- "warmfit"

  beta <- Matrix::Matrix(wls_fit$a, sparse = TRUE)

  out <- list(a0 = wls_fit$aint, beta = beta, df = sum(abs(beta) > 0),
              dim = dim(beta), lambda = lambda, dev.ratio = wls_fit$rsqc,
              nulldev = nulldev, npasses = wls_fit$nlp, jerr = wls_fit$jerr,
              offset = FALSE, call = this.call, nobs = nobs, warm_fit = warm_fit)
  if (save.fit == FALSE) {
    out$warm_fit <- NULL
  }
  class(out) <- c("glmnetfit", "glmnet")
  out
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
initilizer <- function(x, y, weights, family, intr, has_offset, offset, pfl1) {
  nobs <- nrow(x)
  nvars <- ncol(x)

  if (intercept) {
    if (has_offset) {
      suppressWarnings(
        tempfit <- glm(y ~ 1, family = family, weights = weights,
                       offset = offset)
      )
      mu <- tempfit$fitted.values
    } else {
      mu <- rep(weighted.mean(y, weights), times = nobs)
    }
  } else {
    mu <- family$linkinv(offset)
  }
  nulldev <- dev_function(y, mu, weights, family)

  # compute upper bound for lambda max (assumes some l1 penalty)
  r <- y - mu
  eta <- family$linkfun(mu)
  v <- family$variance(mu)
  m.e <- family$mu.eta(eta)
  weights <- weights / sum(weights)
  rv <- r / v * m.e * weights
  g <- abs(drop(crossprod(rv, x)))
  lambda_max <- max(g / pmax(pfl1, 1e-6))
  list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
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
  if (inherits(x, "sparseMatrix")) {
    beta <- drop(beta) / attr(x, "xs")
    drop(x %*% beta - sum(beta * attr(x, "xm") ) + a0)
  } else {
    drop(x %*% beta + a0)
  }
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
