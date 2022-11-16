
test_that("no penalty", {

  set.seed(1)
  X <- matrix(c(rnorm(6)), nrow = 3)
  beta <- c(2, 1)
  y <- X %*% beta + rnorm(3, sd = .1)
  res1 <- sparsegl(X, y, lambda = 0)
  res2 <- summary(lm(y~X))

  coef1 <- as.numeric(c(res1$b0[1], res1$beta[, 1]))
  coef2 <- as.numeric(res2$coefficients[, 1])
  expect_equal(coef1, coef2, tolerance = 1e-3)

})

test_that("penalty is super large", {

  X <- matrix(c(rnorm(20)), nrow = 4)
  beta <- seq(5)
  y <- X %*% beta + rnorm(4, sd = .1)
  res <- sparsegl(X, y, lambda = 10)
  coef <- as.numeric(c(res$beta[, 1]))
  expect_equal(coef, double(5))
})

test_that("the number of nonzero coefficient features with penalty is less than
          or equal to without penalty", {

  X <- matrix(c(rnorm(100)), nrow = 10)
  beta <- seq(10)
  y <- X %*% beta + rnorm(10, sd = .1)
  res1 <- sparsegl(X, y)
  res2 <- sparsegl(X, y, asparse = 0)
  res3 <- sparsegl(X, y, asparse = 0.99)
  res4 <- sparsegl(X, y, lambda = 0)

  num1 <- apply(res1$beta, 2, function(x) sum(x != 0))
  num2 <- apply(res2$beta, 2, function(x) sum(x != 0))
  num3 <- apply(res3$beta, 2, function(x) sum(x != 0))
  num4 <- sum(res2$beta != 0)

  expect_equal(as.numeric(num1 <= num4), rep(1, 100))
  expect_equal(as.numeric(num2 <= num4), rep(1, 100))
  expect_equal(as.numeric(num3 <= num4), rep(1, 100))
})

