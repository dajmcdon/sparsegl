#' @describeIn grouped_sp_norm l0-norm (number of nonzero entries).
#' @export
zero_norm <- function(x) sum(abs(x) > 0)


#' @describeIn grouped_sp_norm l1-norm (Absolute-value norm).
#' @export
one_norm <- function(x) sum(abs(x))

#' @describeIn grouped_sp_norm l2-norm (Euclidean norm).
#' @export
two_norm <- function(x) sqrt(sum(x^2))


#' @describeIn grouped_sp_norm A vector of group-wise l0-norms.
#' @export
grouped_zero_norm <- function(x, gr) as.vector(tapply(x, gr, zero_norm))

#' @describeIn grouped_sp_norm A vector of group-wise l1-norms.
#' @export
grouped_one_norm <- function(x, gr) as.vector(tapply(x, gr, one_norm))

#' @describeIn grouped_sp_norm A vector of group-wise l2-norms.
#' @export
grouped_two_norm <- function(x, gr) as.vector(tapply(x, gr, two_norm))



#' Calculate common norms
#'
#' Calculate different norms of vectors with or without grouping structures.
#'
#' @param x A numeric vector.
#' @param gr An integer (or factor) vector of the same length as x.
#' @param asparse Scalar. The weight to put on the l1 norm when calculating the
#'   group norm.
#'
#' @return A numeric scalar or vector
#' @export
#'
#' @examples
#' x <- c(rep(-1, 5), rep(0, 5), rep(1,5))
#' gr <- c(rep(1,5), rep(2,5), rep(3,5))
#' asparse <- 0.05
#' grouped_sp_norm(x, gr, asparse)
#'
#' @describeIn grouped_sp_norm A vector of length `unique(gr)` consisting of
#'   the `asparse` convex combination of the l1 and l2-norm for each group.
grouped_sp_norm <- function(x, gr, asparse) {
  stopifnot(length(x) == length(gr), asparse >= 0, asparse <= 1)
  asparse * grouped_one_norm(x, gr) + (1 - asparse) * grouped_two_norm(x, gr)
}

#' @describeIn grouped_sp_norm The l1-norm norm of a vector (a scalar).
#' @export
gr_one_norm <- function(x, gr) one_norm(x)

#' @describeIn grouped_sp_norm The sum of the group-wise l2-norms of a vector
#'   (a scalar).
#' @export
gr_two_norm <- function(x, gr) sum(grouped_two_norm(x, gr))

#' @describeIn grouped_sp_norm The sum of the `asparse` convex combination of
#'   group l1 and l2-norms vectors (a scalar).
#' @export
sp_group_norm <- function(x, gr, asparse = 0.05) sum(grouped_sp_norm(x, gr, asparse))
