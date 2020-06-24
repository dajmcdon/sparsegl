#' Fits the regularization paths for group-lasso penalized learning problems
#' 
#' Fits regularization paths for group-lasso penalized learning problems at a
#' sequence of regularization parameters lambda.
#' 
#' Note that the objective function for \code{"ls"} least squares is
#' \deqn{RSS/(2*n) + lambda * penalty;} for \code{"hsvm"} Huberized squared
#' hinge loss, \code{"sqsvm"} Squared hinge loss and \code{"logit"} logistic
#' regression, the objective function is \deqn{-loglik/n + lambda * penalty.}
#' Users can also tweak the penalty by choosing different penalty factor.
#' 
#' For computing speed reason, if models are not converging or running slow,
#' consider increasing \code{eps}, decreasing \code{nlambda}, or increasing
#' \code{lambda.factor} before increasing \code{maxit}.
#' 
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#' is an observation vector.
#' @param y response variable. This argument should be quantitative for
#' regression (least squares), and a two-level factor for classification
#' (logistic model, huberized SVM, squared SVM).
#' @param group a vector of consecutive integers describing the grouping of the
#' coefficients (see example below).
#' @param loss a character string specifying the loss function to use, valid
#' options are: \itemize{ \item \code{"ls"} least squares loss (regression),
#' \item \code{"logit"} logistic loss (classification).  \item \code{"hsvm"}
#' Huberized squared hinge loss (classification), \item \code{"sqsvm"} Squared
#' hinge loss (classification), }Default is \code{"ls"}.
#' @param nlambda the number of \code{lambda} values - default is 100.
#' @param lambda.factor the factor for getting the minimal lambda in
#' \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} *
#' \code{max(lambda)}.  \code{max(lambda)} is the smallest value of
#' \code{lambda} for which all coefficients are zero. The default depends on
#' the relationship between \eqn{n} (the number of rows in the matrix of
#' predictors) and \eqn{p} (the number of predictors). If \eqn{n >= p}, the
#' default is \code{0.001}, close to zero.  If \eqn{n<p}, the default is
#' \code{0.05}.  A very small value of \code{lambda.factor} will lead to a
#' saturated fit. It takes no effect if there is user-defined \code{lambda}
#' sequence.
#' @param lambda a user supplied \code{lambda} sequence. Typically, by leaving
#' this option unspecified users can have the program compute its own
#' \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}.
#' Supplying a value of \code{lambda} overrides this. It is better to supply a
#' decreasing sequence of \code{lambda} values than a single (small) value, if
#' not, the program will sort user-defined \code{lambda} sequence in decreasing
#' order automatically.
#' @param pf penalty factor, a vector in length of bn (bn is the total number
#' of groups). Separate penalty weights can be applied to each group of
#' \eqn{\beta}{beta's}s to allow differential shrinkage. Can be 0 for some
#' groups, which implies no shrinkage, and results in that group always being
#' included in the model. Default value for each entry is the square-root of
#' the corresponding size of each group.
#' @param dfmax limit the maximum number of groups in the model. Useful for
#' very large \code{bs} (group size), if a partial path is desired. Default is
#' \code{bs+1}.
#' @param pmax limit the maximum number of groups ever to be nonzero. For
#' example once a group enters the model, no matter how many times it exits or
#' re-enters model through the path, it will be counted only once. Default is
#' \code{min(dfmax*1.2,bs)}.
#' @param eps convergence termination tolerance. Defaults value is \code{1e-8}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda
#' value. Default is 3e8. If models do not converge, consider increasing
#' \code{maxit}.
#' @param delta the parameter \eqn{\delta}{delta} in \code{"hsvm"} (Huberized
#' squared hinge loss). Default is 1.
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#' @param asparse the weight to put on the ell1 norm in sparse group lasso. Default 
#' is 0.05 
#' @return An object with S3 class \code{\link{gglasso}}.  \item{call}{the call
#' that produced this object} \item{b0}{intercept sequence of length
#' \code{length(lambda)}} \item{beta}{a \code{p*length(lambda)} matrix of
#' coefficients.} \item{df}{the number of nonzero groups for each value of
#' \code{lambda}.} \item{dim}{dimension of coefficient matrix (ices)}
#' \item{lambda}{the actual sequence of \code{lambda} values used}
#' \item{npasses}{total number of iterations (the most inner loop) summed over
#' all lambda values} \item{jerr}{error flag, for warnings and errors, 0 if no
#' error.} \item{group}{a vector of consecutive integers describing the
#' grouping of the coefficients.}
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @seealso \code{plot.gglasso}
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load bardet data set
#' data(bardet)
#' 
#' # define group index
#' group1 <- rep(1:20,each=5)
#' 
#' # fit group lasso penalized least squares
#' m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")
#' 
#' # load colon data set
#' data(colon)
#' 
#' # define group index
#' group2 <- rep(1:20,each=5)
#' 
#' # fit group lasso penalized logistic regression
#' m2 <- gglasso(x=colon$x,y=colon$y,group=group2,loss="logit")
#' 
#' @export
gglasso <- function(x, y, group = NULL, loss = c("ls", "ls_sparse", "logit", "sqsvm", 
    "hsvm","wls"), nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), 
    lambda = NULL, pf = sqrt(bs), weight = NULL, dfmax = as.integer(max(group)) + 
        1, pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08, 
    delta, intercept=TRUE, asparse = 0.05, standardize=TRUE) {
    #################################################################################
    #\tDesign matrix setup, error checking
    this.call <- match.call()
    loss <- match.arg(loss)
    
    if (!is.matrix(x)) 
        stop("x has to be a matrix")
    
    if (any(is.na(x))) 
        stop("Missing values in x not allowed!")
    
    y <- drop(y)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    
    if (length(y) != nobs) 
        stop("x and y have different number of rows")
    
    if (!is.numeric(y)) 
        stop("The response y must be numeric. Factors must be converted to numeric")
    
    c1 <- loss %in% c("logit", "sqsvm", "hsvm")
    c2 <- any(y %in% c(-1, 1) == FALSE)
    if (c1 && c2) 
        stop("Classification method requires the response y to be in {-1,1}")
	
    if (loss=="wls" & !is.matrix(weight)) 
        stop("User must specify weight matrix for (loss='wls')")
    #################################################################################
    #    group setup
    if (is.null(group)) {
        group <- 1:nvars
    } else if (length(group) != nvars) 
        stop("group length does not match the number of predictors in x")
    
    bn <- as.integer(max(group))
    bs <- as.integer(as.numeric(table(group)))
    
    if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) 
        stop("Groups must be consecutively numbered 1,2,3,...")
    
    if (loss=="ls_sparse" && (asparse>1 || asparse<0)){
      asparse = 0
      loss="ls"
      warning("asparse must be in [0,1], running ordinary group lasso.")
    } 
    
    ix <- rep(NA, bn)
    iy <- rep(NA, bn)
    j <- 1
    for (g in 1:bn) {
        ix[g] <- j
        iy[g] <- j + bs[g] - 1
        j <- j + bs[g]
    }
    ix <- as.integer(ix)
    iy <- as.integer(iy)
    group <- as.integer(group)
    #################################################################################
    #parameter setup
    if (missing(delta)) 
        delta <- 1
    if (delta < 0) 
        stop("delta must be non-negtive")
    delta <- as.double(delta)
    if (length(pf) != bn) 
        stop("The size of group-lasso penalty factor must be same as the number of groups")
    maxit <- as.integer(maxit)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    #################################################################################
    #lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
        if (lambda.factor >= 1) 
            stop("lambda.factor should be less than 1")
        flmin <- as.double(lambda.factor)
        ulam <- double(1)
    } else {
        #flmin=1 if user define lambda
        flmin <- as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    intr <- as.integer(intercept)
    #################################################################################
    # call R sub-functions
    fit <- switch(loss, 
	ls = ls(bn, bs, ix, iy, nobs, nvars, x, y, pf, 
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr), 
	ls_sparse = ls_sparse(bn, bs, ix, iy, nobs, nvars, x, y, pf, 
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, 
        asparse, standardize), 
	logit = logit(bn, 
        bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, 
        ulam, eps, maxit, vnames, group, intr), 
	sqsvm = sqsvm(bn, bs, ix, iy, 
        nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, 
        group, intr), 
	hsvm = hsvm(delta, bn, bs, ix, iy, nobs, nvars, x, y, 
        pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr),
	wls = wls(bn, bs, ix, iy, nobs, nvars, x, y, pf, weight, 
        dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr)
		)
    #################################################################################
    # output
    if (is.null(lambda)) 
        fit$lambda <- lamfix(fit$lambda)
    fit$call <- this.call
    class(fit) <- c("gglasso", class(fit))
    fit
} 
