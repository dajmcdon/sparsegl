
test_that("test for exact_df function", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  out <- sparsegl(X, y, gr)
  
  expect_equal(exact_df(out, X), as.vector(out$df), tolerance = 0.01)
})
