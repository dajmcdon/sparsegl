sgl <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps,
    maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd) {
    # call Fortran core
    is.sparse <- FALSE
    algorithm <- "sgl"
    if (inherits(x,"sparseMatrix")) {
        is.sparse <- TRUE
        x <- methods::as(x,"CsparseMatrix")
        x <- methods::as(x,"dgCMatrix")
        algorithm <- "sp_sgl"
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
        sx <- sqrt(Matrix::colSums(x^2))
        sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
        xs <- 1 / sx
        x <- x %*% Matrix::Diagonal(x = xs)
    }
    if (is.sparse) {
        xidx <- as.integer(x@i + 1)
        xcptr <- as.integer(x@p + 1)
        xval <- as.double(x@x)
        nnz <- as.integer(tail(x@p, 1))
    }

    gamma <- calc_gamma(x, ix, iy, bn)

    fit <- switch(algorithm,
        sgl = .Fortran(
            "sparse_four", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam = 0L,
            beta = double(nvars * nlam), activeGroup = integer(pmax),
            nbeta = integer(nlam), alam = double(nlam), npass = 0L, jerr = 0L,
            alsparse = as.double(asparse), lb = lower_bnd, ub = upper_bnd),
        sp_sgl = .Fortran(
            "spmat_four", bn,bs,ix,iy,gamma,nobs,nvars,xval,xidx,xcptr,nnz,
            as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit,
            as.integer(intr), nalam = 0L, b0 = double(nlam),
            beta = double(nvars * nlam),
            activeGroup = integer(pmax), nbeta = integer(nlam),
            alam = double(nlam),
            npass = 0L, jerr = 0L, alsparse = as.double(asparse),
            lb = lower_bnd, ub = upper_bnd),
        stop("Requested algorithm is not implemented.")
        )
    # output
    outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
    if (standardize) {
        outlist$beta = outlist$beta * xs
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


calc_gamma <- function(x, ix, iy, bn) {
    gamma <- rep(NA, bn)
    for (g in seq_len(bn)) {
        grabcols <- ix[g]:iy[g]
        # when ncol is larger than 2
        if (length(grabcols) > 2) {
            gamma[g] <- RSpectra::svds(x[ ,grabcols], 1, 0, 0)$d^2
        } 
        # when ncol is equal to 2
        else if (length(grabcols) == 2) {
            gamma[g] <- maxeig2(x[,grabcols])
        } 
        # when ncol is 1
        else {
            gamma[g] <- sum(x[,grabcols]^2)
        }
    }
    return(as.double(gamma / nrow(x)))
}

# returns the largest squared singular value of n-by-2 matrix x (the largest
# eigenvalue of 2-by-2 matrix mat)
maxeig2 <- function(x) {
    mat <- crossprod(x)
    tr <- mat[1] + mat[4]
    dt <- mat[1] * mat[4] - mat[2]^2
    return((tr + sqrt(tr^2 - 4 * dt)) / 2)
}


