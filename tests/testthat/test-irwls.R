test_that("get_eta works", {
  n <- 100
  p <- 10
  x <- matrix(rnorm(100*10), 100, 10)
  beta <- p:1
  b0 <- 3
  xs <- sqrt(colSums(x^2))
  ge <- get_eta(x, xs, beta, b0)
  expect_length(ge, n)
  ge <- get_eta(x, NULL, beta, b0)
  expect_length(ge, n)
  ge <- get_eta(x, NULL, as.matrix(beta, ncol = 1), b0)
  expect_length(ge, n)
})

test_that("deviance works", {
  n <- 100
  p <- 20
  gr <- rep(1:4, each = 5)
  y <- rbinom(n, 1, rbeta(100, 2, 2))
  mu <- runif(n)
  w <- rep(1, n)
  pf <- runif(4)
  pfl1 <- runif(p)
  asparse <- .5
  lambda <- .2
  coefs <- rnorm(p)
  expect_silent(dev_function(y, runif(n), rep(1, n), binomial()))
  expect_silent(dev_function(y, runif(n), rep(1, n), poisson()))
  expect_silent(dev_function(y, runif(n), rep(1, n), gaussian()))

  expect_silent(obj_function(y, mu, gr, w, binomial(), pf, pfl1, asparse,
                             coefs, lambda))

  expect_silent(obj_function(y, mu, gr, w, gaussian(), pf, pfl1, asparse,
                             coefs, lambda))

})
