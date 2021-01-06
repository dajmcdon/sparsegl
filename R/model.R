sgl <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps,
    maxit, vnames, group, intr, asparse, standardize, algorithm) {
    # call Fortran core
    is.sparse <- FALSE
    if (inherits(x,"sparseMatrix")) {
        is.sparse <- TRUE
        x <- as(x,"CsparseMatrix")
        x <- as(x,"dgCMatrix")
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
        xs <- 1 / sqrt(Matrix::colSums(x^2))
        x <- x %*% Matrix::Diagonal(x = xs)
    }
    if (is.sparse) {
        xidx <- as.integer(x@i + 1)
        xcptr <- as.integer(x@p + 1)
        xval <- as.double(x@x)
        nnz <- as.integer(tail(x@p, 1))
        algorithm="fourstepspx"
    }
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- RSpectra::svds(x[,ix[g]:iy[g]], 1, 0, 0)$d^2
    gamma <- gamma/nobs
    gamma <- as.double(gamma)

    fit <- switch(algorithm,
        threestep = .Fortran(
            "sparse_three", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = 0L,
            beta = double(nvars * nlam), activeGroup = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = 0L, jerr = 0L,
            alsparse = as.double(asparse)),
        fourstep = .Fortran(
            "sparse_four", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = 0L,
            beta = double(nvars * nlam), activeGroup = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = 0L, jerr = 0L,
            alsparse = as.double(asparse)),
        fourstepspx = .Fortran(
            "spmat_four", bn,bs,ix,iy,gamma,nobs,nvars,xval,xidx,xcptr,nnz,
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit,
            as.integer(intr), nalam=0L, b0=double(nlam),beta=double(nvars*nlam),
            activeGroup=integer(pmax),nbeta=integer(nlam),alam=double(nlam),
            npass=0L,jerr=0L,alsparse=as.double(asparse)),
        stop("Requested algorithm is not implemented.")
        )
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    if (standardize) {
        outlist$beta = outlist$beta*xs
    }
    if (intr) {
        if (is.sparse) {
            outlist$b0 <- outlist$b0 + ym
        } else {
            outlist$b0 <- ym - xm %*% outlist$beta
        }
    }
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("ls")
    outlist
}


gglasso <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax,
    pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, standardize) {
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




