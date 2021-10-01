
test_that("check if function delP() works as expected", {
  
  # beta is a vector of zeros
  beta <- double(15)
  gr <- rep(1:5, each = 3)
  list1 <- list()
  for (i in seq(5)) {
    list1[[i]] <- matrix(rep(NaN, 9), nrow = 3)
  }
  expect_equal(delP(beta, gr), Matrix::bdiag(list1))
  
  # beta is a matrix of zeros
  beta <- matrix(double(15))
  expect_equal(delP(beta, gr), Matrix::bdiag(list1))
  
  # beta length is not equal to gr length
  gr <- rep(1:5, 2)
  expect_warning(delP(beta, gr))
})

test_that("check if function exact_df() works as expected", {
  
  set.seed(1)
  n <- 100
  beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
  gr <- rep(1:5, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  out <- sparsegl(X, y, gr)
  
  # check the length of the exact DOF
  expect_equal(length(exact_df(out, X)), 100)
  
  # check if exact_df produce 0 when all betas are zero at some lambda
  for (i in seq(100)) {
    if (sum(abs(out$beta[, i])) == 0) {
      expect_equal(exact_df(out, X)[i], 0)
    }
  }
})
