
test_that("test for exact_df function", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  out <- sparsegl(X, y, gr)
  
  expect_equal(exact_df(out, X), as.vector(out$df), tolerance = 0.01)
})

test_that("type of y", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y1 <- X %*% beta + rnorm(n)
  y2 <- drop(y1)
  out1 <- sparsegl(X, y1, gr)
  out2 <- sparsegl(X, y2, gr)
  
  expect_equal(risk_estimate(out1, X, y1, type = "AIC"), risk_estimate(out2, X, y2, type = "AIC"))
})