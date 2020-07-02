ls_sparse <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax,
    pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize) {
    #################################################################################
    # call Fortran core
    intercept = 0L
    if(intr){
        ym = mean(y)
        xm = colMeans(x)
        x = sweep(x,2,xm)
        y = y-ym
    }
    if(standardize){
        xs = sqrt(colSums(x^2))
        x = sweep(x,2,xs,"/")
    }
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- RSpectra::svds(x[,ix[g]:iy[g]],1,0,0)$d^2
    gamma <- gamma/nobs
    gamma <- as.double(gamma)
    fit <- .Fortran("ls_f_sparse", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intercept, nalam = integer(1),
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax),
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),
        alsparse = asparse)
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    if(standardize){
        outlist$beta = outlist$beta/xs
    }
    if(intr){
        outlist$b0 = ym - xm %*% outlist$beta
    }
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("ls")
    outlist
}


#' Title
#'
#' @param jxx 
#' @param bn 
#' @param ga 
#' @param pf 
#' @param tlam 
#' @param alsparse 
#'
#' @return
#' @export
#'
strong_rule <- function(jxx, bn, ga, pf, tlam, alsparse) {

    jxx <- as.integer(jxx)
    bn <- as.integer(bn)
    fit <- .Fortran("strong_rule", jxx, bn, ga, pf, tlam, alsparse)
    ###############################################################
    # output
    fit
}

softthresh <- function(vec, thresh, n) {

    fit <- .Fortran("softthresh", vec, thresh, n)
    ##############################################################
    # output
    fit

}

