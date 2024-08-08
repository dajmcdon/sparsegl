#' @useDynLib sparsegl, .registration = TRUE
#' @importFrom utils packageDescription
#' @importFrom methods cbind2 rbind2 as
#' @importFrom stats approx coef predict fitted
#' @importFrom rlang .data abort warn
#' @importFrom rlang := %||%
#' @importFrom cli cli_abort cli_warn
#' @importFrom dotCall64 .C64 integer_dc
#' @importFrom RSpectra svds
#' @import Matrix
#' @keywords internal
#' @references Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., and
#'   McDonald, D.J. 2024.
#'   \emph{sparsegl: An `R` Package for Estimating Sparse Group Lasso.}
#'   Journal of Statistical Software, Vol. 110(6): 1–23.
#'   \doi{10.18637/jss.v110.i06}.
"_PACKAGE"
