
test_that("largest squared singular value", {
  x1 <- matrix(c(1,2,2,1), nrow = 2)  # 2-by-2 matrix
  expect_equal(maxeig2(x1), max(eigen(crossprod(x1))$values))


  x2 <- matrix(c(1,0,1,1,0,1,0,0,0,1,1,0), nrow = 3) # 3-by-4 matrix
  expect_equal(calc_gamma(x2, ix = 1, iy = 4, bn = 1), RSpectra::svds(x2, 1, 0, 0)$d^2 / 3)

  res <- calc_gamma(x2, ix = c(1,3), iy = c(2,4), bn = 2)  # 3-by-4 matrix divided into two 3-by-2 matrices
  res1 <- double(2)
  res1[1] <- maxeig2(x2[, c(1,2)])
  res1[2] <- maxeig2(x2[, c(3,4)])
  expect_equal(res, res1/3)

  A <- matrix(sample(1:50, 48, replace = TRUE), nrow = 6, ncol = 8)  # random 6-by-8 matrix
  expect_equal(calc_gamma(A, ix = 1, iy = 8, bn = 1), RSpectra::svds(A, 1, 0, 0)$d^2/6)

  res <- calc_gamma(A, ix = c(1,3,6), iy = c(2, 5, 8), bn = 3) # 6-by-8 matrix divided into 3 matrices
  res1 <- double(3)
  res1[1] <- maxeig2(A[, 1:2])
  res1[2] <- RSpectra::svds(A[, 3:5], 1, 0, 0)$d^2
  res1[3] <- RSpectra::svds(A[, 6:8], 1, 0, 0)$d^2
  expect_equal(res, res1/6)

  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)  # 8-by-10 sparse matrix
  A <- Matrix::sparseMatrix(i, j, x = x)
  expect_equal(calc_gamma(A, ix = 1, iy = 10, bn = 1), RSpectra::svds(A, 1, 0, 0)$d^2 / 8)

  res <- calc_gamma(A, ix = c(1,6), iy = c(5, 10), bn = 2)  # 8-by-10 sparse matrix divided into two 8-by-5 matrices
  res1 <- double(2)
  res1[1] <- RSpectra::svds(A[, 1:5], 1, 0, 0)$d^2
  res1[2] <- RSpectra::svds(A[, 6:10], 1, 0, 0)$d^2
  expect_equal(res, res1/8)
})
