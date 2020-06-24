#' gglasso: A package for group-lasso penalized learning problems
#'
#' There are two main functions in the gglasso: \code{\link{gglasso}} and
#' \code{\link{cv.gglasso}}
#'
#' @useDynLib gglasso, .registration = TRUE
#' @importFrom utils head tail packageDescription
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis par plot.default points segments
#' @importFrom methods cbind2 rbind2
#' @importFrom stats approx coef predict
#' @importFrom RSpectra svds
#' @docType package
NULL
