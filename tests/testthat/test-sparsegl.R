library(Matrix)
library(glmnet)
library(gglasso)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("x is not a matrix", {
  
  x <- c(1,2)
  y <- c(1,2)
  expect_false(is.matrix(x))
  expect_false(inherits(x, "sparseMatrix"))
  expect_error(sparsegl(x, y) )
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
  expect_error(sparsegl(x, y))
})


test_that("x and y have different dimensions", {
  
  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  y = c(1,2)
  expect_error(sparsegl(A,y))
})

test_that("y is not numeric", {
  
  x <- matrix(c(1,2,3,4,5,6), nrow =2)
  y <- c("1", "2")
  expect_type(y, "character")
  expect_false(is.numeric(y))
  expect_error(sparsegl(x, y))
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
  expect_error(sparsegl(x, y, group = c(1,2)))
  expect_error(sparsegl(x, y, group = c(1,2,4)))
  expect_error(sparsegl(x, y, group = c(2,3,4)))
  expect_warning(sparsegl(x, y, group = c(1,2,3), asparse = 1.1))
  expect_warning(sparsegl(x, y, group = c(1,2,3), asparse = -1))
  expect_error(sparsegl(x, y, pf = c(1,2)))
  expect_error(sparsegl(x, y, lambda.factor = 2))
  expect_error(sparsegl(x, y, lambda = c(-1, 1)))
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5, 2)))
  expect_error(sparsegl(x, y, upper_bnd = c(-1, -5, 2)))
  expect_error(sparsegl(x, y, lower_bnd = c(-1, -5)))
  expect_error(sparsegl(x, y, upper_bnd = c(1,5)))
})

test_that("no penalty", {
  
  X <- matrix(c(1,2,3,4,5,6), nrow = 3)
  y <- c(1.5, 3.5, 5.5)
  res1 <- sparsegl(X, y, lambda = 0)
  res2 <- summary(lm(y~X))
  
  expect_equal(res1$b0[1], res2$coefficients[1,1])
  expect_equal(res1$beta[1], res2$coefficients[2,1])
  expect_equal(res1$beta[2], 0)
  
  # print(res1$b0)
  # print(res1$beta)
  # print(res2$coefficients)
  
})

test_that("lasso penalty", {
  
  X <- matrix(c(1,2,3,4,5,6), nrow = 3)
  y <- c(2,4,6)
  
  # try large lambda
  expect_equal(sparsegl(X, y, lambda = 2, asparse = 1)$beta[1:2],c(0,0))
  
  fit <- glmnet(X, y, alpha = 1)
  lambda <- fit$lambda[c(1:2, (length(fit$lambda) - 1) :length(fit$lambda))]
  # check the first two lambdas and last two lambdas extracted from glmnet
  
  for (i in 1 : length(lambda)){
    fit1 <- sparsegl(X, y, asparse = 1, lambda = lambda[i])
    fit2 <- glmnet(X, y, alpha = 1)
    print(paste("lambda = ", lambda[i]))
    print(fit1$b0)
    print(fit1$beta)
    print(coef(fit2, s = lambda[i]))
  }
})

test_that("group lasso penalty", {
  
  set.seed(20200110)
  n <- 100
  p <- 48
  g <- 6
  X <- matrix(rnorm(n*p), nrow = n)
  for (i in seq(1, p, by = 6))
    X[, i:(i + g - 1)] <- svd(X[, i:(i + g - 1)])$u
  beta <- c(rep(5, 6), rep(-5, 6), rep(0, p - 12))
  y <- X %*% beta + rnorm(n, sd = .1)
  
  group <-  rep(1:16,each=3)
  fit <- gglasso(x = X, y = y, group = group)
  lambda = fit$lambda[1:5]
  
  # try large lambda
  expect_equal(sparsegl(X, y, lambda = 2, asparse = 0, group = group)$beta[1:48], c(rep(0,48)))


  for (i in 1:length(lambda)){
    fit1 <- sparsegl(X, y, asparse = 0, lambda = lambda[i], group = group)
    fit2 <- gglasso(X, y, group = group)
    print(paste("lambda = ", lambda[i]))
    print(fit1$b0)
    print(fit1$beta)
    print(coef(fit2, s = lambda[i]))
  }
})



