calc_gamma <- function(x, ix, iy, bn) {
  gamma <- rep(NA, bn)
  for (g in seq_len(bn)) {
    grabcols <- ix[g]:iy[g]
    ncols <- length(grabcols)
    if (ncols > 2) gamma[g] <- RSpectra::svds(x[,grabcols], 1, 0, 0)$d^2
    else {
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


# Slight improvements -----------------------------------------------------

#
#
# calc_gamma_irwls_dense <- function(x, ix, iy, bn, sqrw) {
#     gamma <- rep(NA, bn)
#     nrows <- nrow(x)
#     for (g in seq_len(bn)) {
#         grabcols <- ix[g]:iy[g]
#         ncols <- length(grabcols)
#         if (ncols > 2) {
#             if (nrows > ncols) {
#                 gamma[g] <- RSpectra::svds(t(x[, grabcols]), 1, 0, 0,
#                     opts = list(scale = sqrw))$d^2
#             } else {
#                 gamma[g] <- RSpectra::svds(x[, grabcols] / sqrw, 1, 0, 0)$d^2
#             }
#         } else {
#             if (ncols == 2) gamma[g] <- maxeig2(x[, grabcols] / sqrw)
#             else gamma[g] <- sum(x[, grabcols]^2 / sqrw^2)
#         }
#     }
#     return(as.double(gamma / nrow(x)))
# }
#
# calc_gamma_irwls_sparse <- function(x, ix, iy, bn, sqrw) {
#     gamma <- rep(NA, bn)
#     nrows <- nrow(x)
#     for (g in seq_len(bn)) {
#         grabcols <- ix[g]:iy[g]
#         ncols <- length(grabcols)
#         if (ncols > 2) {
#             if (nrows < ncols) { # taller is faster for sparse matrices
#                 gamma[g] <- RSpectra::svds(t(x[, grabcols]), 1, 0, 0,
#                     opts = list(scale = sqrw))$d^2
#             } else {
#                 gamma[g] <- RSpectra::svds(x[, grabcols] / sqrw, 1, 0, 0)$d^2
#             }
#         } else {
#             if (ncols == 2) gamma[g] <- maxeig2(x[, grabcols] / sqrw)
#             else gamma[g] <- sum(x[, grabcols]^2 / sqrw^2)
#         }
#     }
#     return(as.double(gamma / nrow(x)))
# }
