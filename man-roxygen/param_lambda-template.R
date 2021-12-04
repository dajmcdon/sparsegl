#' @param lambda A user supplied `lambda` sequence. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the smallest value
#'   of `lambda` that would give the null model (all coefficient estimates equal
#'   to zero), and `lambda.factor`. Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
