auc_mat <- function(y, prob, weights = rep(1, N)) {
  N <- length(y)
  y <- cbind(1 - y, y)
  w <- as.vector(weights * y) # length 2N
  Y <- rep(c(0, 1), each = N)
  Prob <- c(prob, prob)
  survival::concordance(Y ~ Prob, weights = w)$concordance
}
