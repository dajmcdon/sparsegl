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


test_that("initializer works", {
  n <- 100
  p <- 20
  gr <- rep(1:4, each = 5)
  y <- rbinom(n, 1, rbeta(100, 2, 2))
  x <- matrix(rnorm(n*p), n, p)
  xsp <- x
  xsp[abs(xsp) < 1] <- 0
  xsp <- Matrix::Matrix(xsp, sparse = TRUE)

  weights <- rnorm(n)^2
  weights <- weights / sum(weights)
  offset <- rep(1, n)
  pfl1 <- rep(1, p)
  pfl2 <- rnorm(p)^2
  pfl2 <- pfl2 / sum(pfl2)
  ulam <- 0.5

  usr_lambda <- initilizer(x, y, weights, gaussian(), TRUE, FALSE, NULL,
                           pfl1, ulam)
  expect_false(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)
  usr_lambda_offset <- initilizer(x, y, weights, gaussian(), TRUE, TRUE, offset,
                           pfl1, ulam)
  expect_false(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)
  usr_lambda_offset_noint <- initilizer(x, y, weights, gaussian(), FALSE,
                                        TRUE, offset, pfl1, ulam)
  expect_false(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)


  usr_lambda <- initilizer(x, y, weights, gaussian(), TRUE, FALSE, NULL,
                           pfl1, 0)
  expect_true(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)
  usr_lambda_offset <- initilizer(x, y, weights, gaussian(), TRUE, TRUE, offset,
                                  pfl1, 0)
  expect_true(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)
  usr_lambda_offset_noint <- initilizer(x, y, weights, gaussian(), FALSE,
                                        TRUE, offset, pfl1, 0)
  expect_true(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)

  # try sparse
  usr_lambda <- initilizer(xsp, y, weights, binomial(), TRUE, FALSE, NULL,
                           pfl1, 0)
  expect_true(usr_lambda$findlambda)
  expect_equal(usr_lambda$cur_lambda * 0.99, usr_lambda$lambda_max)
})

test_that("call to wls FORTRAN subfun works", {
  n <- 100
  p <- 20
  gr <- rep(1:4, each = 5)
  y <- rbinom(n, 1, rbeta(100, 2, 2))
  x <- matrix(rnorm(n*p), n, p)


})
