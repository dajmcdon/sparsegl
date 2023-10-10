library(gglasso)

test_that("results match gglasso up to tolerance", {
  n <- 100
  p <- c(50, 100, 150)
  ngroups <- 10

  make_problem <- function(nobs, nvars, SNR = 1, seed = 12345) {
    set.seed(seed)
    grp_size <- nvars / ngroups
    gr <- rep(1:ngroups, each = grp_size)
    coefs <- rep(c(1, 0), each = grp_size)
    X <- matrix(rnorm(nobs * nvars), nobs, nvars)
    beta <- rep(coefs, length.out = nvars)
    mu <- drop(X %*% beta)
    noise_sd <- sqrt(sum(beta^2)) / SNR # in expectation
    epsilon <- rnorm(nobs, 0, noise_sd)
    y <- mu + epsilon
    list(y = y, X = X, seed = seed, gr = gr)
  }

  p1 <- make_problem(n, p[1])
  s <- with(
    p1, sparsegl(X, y, gr, asparse = 0, standardize = FALSE, intercept = FALSE)
  )
  g <- with(p1, gglasso(X, y, gr, intercept = FALSE, lambda = s$lambda))
  expect_lt(mean(abs(s$beta - g$beta)), 1e-4)

  g <- with(p1, gglasso(X, y, gr, intercept = FALSE))
  s <- with(
    p1, sparsegl(X, y, gr, asparse = 0, standardize = FALSE,
                 intercept = FALSE, lambda = g$lambda)
  )
  expect_lt(mean(abs(s$beta - g$beta)), 1e-4)
})
