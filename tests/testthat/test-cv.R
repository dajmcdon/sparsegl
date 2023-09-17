test_that("cv navigates the options properly", {
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  groups <- rep(1:(p / 5), each = 5)

  expect_silent(cv.sparsegl(X, y, groups))
  expect_silent(cv.sparsegl(X, y, groups, pred.loss = "mse"))
  expect_silent(cv.sparsegl(X, y, groups, pred.loss = "mae"))
  expect_silent(cv.sparsegl(X, y, groups, pred.loss = "deviance"))
  expect_error(cv.sparsegl(X, y, groups, pred.loss = "misclass"))
  expect_error(cv.sparsegl(X, y, groups, family = gaussian(), pred.loss = "misclass"))

  y <- rbinom(n, 1, 0.5)
  expect_silent(cv.sparsegl(X, y, groups, family = "binomial"))
  expect_silent(cv.sparsegl(X, y, groups, family = "binomial", pred.loss = "mse"))
  expect_silent(cv.sparsegl(X, y, groups, family = "binomial", pred.loss = "mae"))
  expect_silent(cv.sparsegl(X, y, groups, family = "binomial", pred.loss = "deviance"))
  expect_silent(cv.sparsegl(X, y, groups, family = "binomial", pred.loss = "misclass"))
  expect_silent(cv.sparsegl(X, y, groups, family = binomial(), pred.loss = "misclass"))

})
