
test_that("no penalty", {
  
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


