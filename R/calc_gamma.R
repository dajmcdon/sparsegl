calc_gamma <- function(x, ix, iy, bn) {
    gamma <- rep(NA, bn)
    for (g in seq_len(bn)) {
        grabcols <- ix[g]:iy[g]
        ncols <- length(grabcols)
        if (ncols > 2) {
            gamma[g] <- RSpectra::svds(x[,grabcols], 1, 0, 0)$d^2
        } else {
            if (ncols == 2) gamma[g] <- maxeig2(x[,grabcols])
            else gamma[g] <- sum(x[,grabcols]^2)
        }
    }
    return(as.double(gamma / nrow(x)))
}

maxeig2 <- function(x) {
    # returns the largest squared singular value of n-by-2 matrix x
    # (the largest eigenvalue of 2-by-2 matrix mat)
    mat <- crossprod(x)
    tr <- mat[1] + mat[4]
    dt <- mat[1] * mat[4] - mat[2]^2
    return((tr + sqrt(tr^2 - 4 * dt)) / 2)
}


