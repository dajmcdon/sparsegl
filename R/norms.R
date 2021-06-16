# Returns L^1-norm (Absolute-value norm).
one_norm <- function(x) sum(abs(x))
# Returns L^2-norm (Euclidean norm).
two_norm <- function(x) sqrt(sum(x^2))

# Returns a vector of group L^1-norms.
grouped_one_norm <- function(x, gr) as.vector(tapply(x, gr, one_norm))
# Returns a vector of group L^2-norms.
grouped_two_norm <- function(x, gr) as.vector(tapply(x, gr, two_norm))
# Returns a vector of combinations of group L^1 and L^2-norms.
grouped_sp_norm <- function(x, gr, asparse) {
  stopifnot(length(x) == length(gr), asparse >= 0, asparse <= 1)
  asparse * grouped_one_norm(x, gr) + (1 - asparse) * grouped_two_norm(x, gr)
}

# Returns the sum of the vector of group L^1-norms.
gr_one_norm <- function(x, gr) one_norm(x)
# Returns the sum of the vector of group L^2-norm.
gr_two_norm <- function(x, gr) sum(grouped_two_norm(x, gr))
# Returns the sum of the vector of combinations of group L^1 and L^2-norms.
sp_group_norm <- function(x, gr, asparse = 0.05) sum(grouped_sp_norm(x, gr, asparse))
