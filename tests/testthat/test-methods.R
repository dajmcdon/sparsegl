test_that("sparsegl print / summary work", {
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  groups <- rep(1:(p / 5), each = 5)
  fit <- sparsegl(X, y, group = groups)
  cv_fit <- cv.sparsegl(X, y, groups)

  expect_silent(summary(fit))
  expect_error(summary(fit, z = 12))
  expect_output(print(fit))
  expect_error(print(fit, z = 12))

  expect_silent(summary(cv_fit))
  expect_error(summary(cv_fit, z = 12))
  expect_output(print(cv_fit))
  expect_error(print(cv_fit, z = 12))

  fit <- sparsegl(X, y, groups, lambda = c(2,1))
  cv_fit <- cv.sparsegl(X, y, groups, lambda = c(2,1))

  expect_silent(summary(fit))
  expect_output(print(fit))
  expect_silent(summary(cv_fit))
  expect_output(print(cv_fit))

  fit <- sparsegl(X, y, groups, lambda = 1)
  cv_fit <- cv.sparsegl(X, y, groups, lambda = 1)

  expect_silent(summary(fit))
  expect_output(print(fit))
  expect_silent(summary(cv_fit))
  expect_output(print(cv_fit))

})
