
test_that("input parameter mismatches cause failures", {

  x <- c(1,2)
  y <- c(1,2)
  expect_error(sparsegl(x, y) )

  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  y <- c(1,2)
  expect_error(sparsegl(x, y), NA)

  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  x[1,1] <- NA
  y <- c(1,2)
  expect_error(sparsegl(x, y))

  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- Matrix::sparseMatrix(i, j, x = x)
  y = c(1,2)
  expect_error(sparsegl(A,y))

  x <- matrix(c(1,2,3,4,5,6), nrow =2)
  y <- c("1", "2")
  expect_error(sparsegl(x, y))

  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- Matrix::sparseMatrix(i, j, x = x)
  y <- c(1:8)
  expect_error(sparsegl(A, y), NA)


  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  y <- c(1,2)
  expect_error(sparsegl(x, y, group = c(1,2)))
  expect_error(sparsegl(x, y, group = c(1,2,4)))
  expect_error(sparsegl(x, y, group = c(2,3,4)))
  expect_error(sparsegl(x, y, group = c(1,2,3), asparse = 1.1))
  expect_warning(sparsegl(x, y, group = c(1,2,3), asparse = -1))
  expect_error(sparsegl(x, y, pf_group = c(1,2)))
  expect_error(sparsegl(x, y, lambda.factor = 2))
  expect_error(sparsegl(x, y, lambda = c(-1, 1)))
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5, 2)))
  expect_error(sparsegl(x, y, upper_bnd = c(-1, -5, 2)))
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5)))
  expect_error(sparsegl(x, y, upper_bnd = c(1,5)))
})

