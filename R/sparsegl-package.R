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
"_PACKAGE"
