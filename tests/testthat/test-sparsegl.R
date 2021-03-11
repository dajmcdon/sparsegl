library(Matrix)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("x is not a matrix", {
  
  x <- c(1,2)
  y <- c(1,2)
  expect_false(is.matrix(x))
  expect_false(inherits(x, "sparseMatrix"))
  expect_error(sparsegl(x, y),"x has to be a matrix" )
})

test_that("x is not a sparse matrix",{
  
  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  y <- c(1,2)
  expect_true(is.matrix(x))
  expect_false(inherits(x, "sparseMatrix"))
  expect_error(sparsegl(x, y), NA)
})


test_that("some entries in x is missing", {
  
  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  x[1,1] <- NA
  y <- c(1,2)
  expect_true(any(is.na(x)))
  expect_error(sparsegl(x, y), "Missing values in x not allowed!")
})


test_that("x and y have different dimensions", {
  
  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  y = c(1,2)
  expect_error(sparsegl(A,y), "x and y have different number of rows")
})

test_that("y is not numeric", {
  
  x <- matrix(c(1,2,3,4,5,6), nrow =2)
  y <- c("1", "2")
  expect_type(y, "character")
  expect_false(is.numeric(y))
  expect_error(sparsegl(x, y), "The response y must be numeric.")
})


test_that("x is a sparse matrix with no missings", {
  
  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  y <- c(1:8)
  expect_false(is.matrix(A))
  expect_true(inherits(A, "sparseMatrix"))
  expect_error(sparsegl(A, y), NA)
})

test_that("parameter issues", {
  
  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  y <- c(1,2)
  expect_error(sparsegl(x, y, group = c(1,2)), "group length does not match the number of predictors in x")
  expect_error(sparsegl(x, y, group = c(1,2,4)), "Groups must be consecutively numbered 1,2,3,...")
  expect_error(sparsegl(x, y, group = c(2,3,4)), "Groups must be consecutively numbered 1,2,3,...")
  expect_warning(sparsegl(x, y, group = c(1,2,3), asparse = 1.1))
  expect_warning(sparsegl(x, y, group = c(1,2,3), asparse = -1))
  expect_error(sparsegl(x, y, pf = c(1,2)), "The length of group-lasso penalty factor must be same as the number of groups")
  expect_error(sparsegl(x, y, lambda.factor = 2), "lambda.factor should be less than 1")
  expect_error(sparsegl(x, y, lambda = c(-1, 1)),"lambdas must be non-negative" )
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5, 2)), "Lower bounds should be non-positive")
  expect_error(sparsegl(x, y, upper_bnd = c(-1, -5, 2)), "Upper bounds should be non-negative")
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5)), "Lower bounds must be length 1 or length the number of groups")
  expect_error(sparsegl(x, y, upper_bnd = c(1,5)), "Upper bounds must be length 1 or length the number of groups")
})








