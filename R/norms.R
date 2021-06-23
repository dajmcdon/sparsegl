#' @describeIn grouped_sp_norm l1-norm (Absolute-value norm).
one_norm <- function(x) sum(abs(x))

#' @describeIn grouped_sp_norm l2-norm (Euclidean norm).
two_norm <- function(x) sqrt(sum(x^2))


#' @describeIn grouped_sp_norm A vector of group l1-norms.
grouped_one_norm <- function(x, gr) as.vector(tapply(x, gr, one_norm))
#' @describeIn grouped_sp_norm A vector of group l1-norms.
grouped_two_norm <- function(x, gr) as.vector(tapply(x, gr, two_norm))



#' Norm
#' 
#' Norm calculation
#'
#' @param x A numeric vector.
#' @param gr A numeric vector of the same length as x.
#' @param asparse The weight to put on the l1 norm when calculating the group 
#' norm.
#'
#' @return A numeric value or vector
#' @export
#'
#' @examples
#' x <- c(rep(-1, 5), rep(0, 5), rep(1,5))
#' gr <- c(rep(1,5), rep(2,5), rep(3,5))
#' asparse <- 0.05
#' grouped_sp_norm(x, gr, asparse)
#' @describeIn grouped_sp_norm A vector of the combinations of group l1 and l2-norms.
grouped_sp_norm <- function(x, gr, asparse) {
  stopifnot(length(x) == length(gr), asparse >= 0, asparse <= 1)
  asparse * grouped_one_norm(x, gr) + (1 - asparse) * grouped_two_norm(x, gr)
}

#' @describeIn grouped_sp_norm The sum of group l1-norms vector.
gr_one_norm <- function(x, gr) one_norm(x)
#' @describeIn grouped_sp_norm The sum of group l2-norm vector.
gr_two_norm <- function(x, gr) sum(grouped_two_norm(x, gr))
#' @describeIn grouped_sp_norm The sum of combinations of group l1 and l2-norms vectors.
sp_group_norm <- function(x, gr, asparse = 0.05) sum(grouped_sp_norm(x, gr, asparse))
