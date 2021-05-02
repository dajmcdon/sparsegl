test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("maximum of eigenvalues", {
  x1 <- matrix(c(1,2,2,1), nrow = 2)  # 2 by 2 matrix
  expect_equal(maxeig2(x1), max(eigen(crossprod(t(x1)))$values))
  
  x2 <- matrix(c(3,2,2,3,2,-2), nrow = 2)  # 2 by 3 matrix
  expect_equal(maxeig2(x2), max(eigen(crossprod(t(x2)))$values))
  
  x3 <- matrix(c(1,0,1,1,0,1,0,0,0,1,1,0), nrow = 3) # 3 by 4 matrix
  start <- 1
  nrow <- 3
  ncol <- 4
  expect_equal(calc_gamma(x3, ix = start, iy = ncol, bn = 1), RSpectra::svds(x3, 1, 0, 0)$d^2 / nrow)
})