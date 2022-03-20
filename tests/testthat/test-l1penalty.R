test_that("the length of pfl1 is not same as the number of predictors", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  pfl1 <- rep(1, 10)
  expect_error(sparsegl(X, y, group = gr, pfl1 = pfl1))
})

test_that("any entry of pfl1 is negative", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  pfl1 <- c(-1, rep(1, 14))
  expect_error(sparsegl(X, y, group = gr, pfl1 = pfl1))
})


test_that("function will rescale each entry such that the sum will be constant", {
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  pfl1_1 <- rep(1, 15)
  pfl1_2 <- rep(2, 15)
  out1 <- sparsegl(X, y, group = gr, pfl1 = pfl1_1)
  out2 <- sparsegl(X, y, group = gr, pfl1 = pfl1_2)
  
  expect_equal(out1$b0, out2$b0)
  expect_equal(out1$beta, out2$beta)
  expect_equal(out1$lambda, out2$lambda)
  
})