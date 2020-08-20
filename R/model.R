sgl <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps,
    maxit, vnames, group, intr, asparse, standardize, algorithm) {
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
    fit <- switch(algorithm,
        original = .Fortran(
            "sparse_orig", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = integer(1),
            beta = double(nvars * nlam), idx = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),
            alsparse = asparse),
        threestep = .Fortran(
            "sparse_three", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = integer(1),
            beta = double(nvars * nlam), idx = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),
            alsparse = asparse),
        threestepalt= .Fortran(
            "sparse_three_alt", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = integer(1),
            beta = double(nvars * nlam), idx = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),
            alsparse = asparse)#,
        # fivestep = .Fortran(
        #     "sparse_five", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
        #     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = integer(1),
        #     beta = double(nvars * nlam), idx = integer(pmax),
        #     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),
        #     alsparse = asparse)
        )
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


gglasso <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax,
    pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, standardize) {
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
    fit <- .Fortran("gglasso", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1),
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax),
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
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




