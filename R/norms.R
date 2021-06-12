two_norm <- function(x) sqrt(sum(x^2))
one_norm <- function(x) sum(abs(x))
grouped_one_norm <- function(x, gr) as.vector(tapply(x, gr, one_norm)) # vector
grouped_two_norm <- function(x, gr) as.vector(tapply(x, gr, two_norm)) # vector
grouped_sp_norm <- function(x, gr, asparse) {
  stopifnot(length(x) == length(gr), asparse >= 0, asparse <= 1)
  asparse * grouped_one_norm(x, gr) + (1 - asparse) * grouped_two_norm(x, gr)
}
gr_two_norm <- function(x, gr) sum(grouped_two_norm(x, gr))
gr_one_norm <- function(x, gr) one_norm(x)
sp_group_norm <- function(x, gr, asparse = 0.05) sum(grouped_sp_norm(x, gr, asparse))
