
test_that("no penalty", {
  
  skip("skip tests for no penalty")
  X <- matrix(c(rnorm(6)), nrow = 3)
  beta <- c(2, 1)
  y <- X %*% beta + rnorm(3, sd = .1)
  res1 <- sparsegl(X, y, lambda = 0)
  res2 <- summary(lm(y~X))
  
  coef <- c(res1$b0, res1$beta)
  coeff <- as.numeric(res2$coefficients[, 1])
  expect_equal(coef, coeff, tolerance = 1e-3)
  
})

test_that("lasso penalty as asparse = 1", {
  
  skip("skip tests or lasso only")
  X <- matrix(c(1,2,3,4,5,6), nrow = 3)
  y <- c(2,4,6)
  
  # try large lambda
  expect_equal(sparsegl(X, y, lambda = 2, asparse = 1)$beta[1:2],c(0,0))
  
  fit <- glmnet::glmnet(X, y, alpha = 1)
  lambda <- fit$lambda[c(1:2, (length(fit$lambda) - 1) :length(fit$lambda))]
  # check the first two lambdas and last two lambdas extracted from glmnet
  
  fit1 <- sparsegl(X, y, asparse = 1, lambda = lambda)
  fit2 <- glmnet::glmnet(X, y, alpha = 1, lambda = lambda)
  
  coef <- as.numeric(rbind(fit1$b0, fit1$beta))
  coeff <- as.numeric(coef(fit2, s= lambda))
  expect_equal(coef,coeff, tolerance = 1e-1)
  
  # try without intercept
  fit1 <- sparsegl(X, y, asparse = 1, lambda = lambda, intercept = FALSE)
  fit2 <- glmnet::glmnet(X, y, alpha = 1, intercept = FALSE)
  coef <- as.numeric(fit1$beta)
  coeff <- as.numeric(coef(fit2, s = lambda))
  expect_equal(coef, coeff, tolerance = 1e-1)
  
  set.seed(20200110)
  n <- 100
  p <- 48
  g <- 6
  X <- matrix(rnorm(n*p), nrow = n)
  
  beta <- c(rep(5, 6), rep(-5, 6), rep(0, p - 12))
  y <- X %*% beta + rnorm(n, sd = .1)
  data = list(x = X, y = y)
  group <- rep(1:8,each=6)
  
  fit <- glmnet::glmnet(X, y, alpha = 1)
  lambda <- fit$lambda[c(1:2, (length(fit$lambda) - 1) :length(fit$lambda))]
  
  
  fit1 <- sparsegl(X, y, asparse = 1, lambda = 2*lambda)
  fit2 <- glmnet::glmnet(X, y, alpha = 1)
  coef <- as.numeric(rbind(fit1$b0, fit1$beta))
  coeff <- as.numeric(coef(fit2, s= lambda))
  expect_equal(coef,coeff, tolerance = 1e-1)
  

  fit1 <- sparsegl(X, y, asparse = 1, lambda = lambda, intercept = FALSE)
  fit2 <- glmnet::glmnet(X, y, alpha = 1, intercept = FALSE)
  coef <- as.numeric(fit1$beta)
  coeff <- as.numeric(coef(fit2, s = lambda))
  expect_equal(coef, coeff, tolerance = 1e-1)
  
})

test_that("group lasso penalty as asparse = 0", {
  
  skip("skip tests for group lasso only")
  set.seed(20200110)
  n <- 100
  p <- 48
  g <- 6
  X <- matrix(rnorm(n*p), nrow = n)
  beta <- c(rep(5, 6), rep(-5, 6), rep(0, p - 12))
  y <- X %*% beta + rnorm(n, sd = .1)
  
  group <-  rep(1:8,each=6)
  fit <- gglasso(x = X, y = y, group = group)
  lambda = fit$lambda[1:5]
  
  # try large lambda
  expect_equal(sparsegl(X, y, lambda = 2, asparse = 0, group = group)$beta[1:48], c(rep(0,48)))
  
  
  fit1 <- sparsegl(X, y, asparse = 0, lambda = lambda, group = group)
  fit2 <- gglasso(X, y, group = group)
  coef <- as.numeric(rbind(fit1$b0, fit1$beta))
  coeff <- as.numeric(coef(fit2, s= lambda))
  expect_equal(coef, coeff, tolerance = 1e-2)
  
  
  fit1 <- sparsegl(X, y, asparse = 0, lambda = lambda, group = group, intercept = FALSE)
  fit2 <- gglasso(X, y, group = group, intercept = FALSE)
  coef <- as.numeric(fit1$beta)
  coeff <- as.numeric(coef(fit2, s = lambda))
  expect_equal(coef, coeff, tolerance = 1e-1)
  
})


test_that("sgl", {
  
  skip("skip tests for sparse group lasso as 0 < asparse < 1")
  set.seed(20200110)
  n <- 100
  p <- 48
  g <- 6
  X <- matrix(rnorm(n*p), nrow = n)
  for (i in seq(1, p, by = 6))
    X[, i:(i + g - 1)] <- svd(X[, i:(i + g - 1)])$u
  beta <- c(rep(5, 6), rep(-5, 6), rep(0, p - 12))
  y <- X %*% beta + rnorm(n, sd = .1)
  
  data = list(x = X, y = y)
  group <- rep(1:8,each=6)
  
  fit1 <- SGL(data, index = group, type = "linear")
  lambda <- fit1$lambda[1:5]
  
  
  fit1 <- sparsegl(X, y, asparse = 0.95, lambda = lambda, group = group)
  fit2 <- SGL(data, index = group, type = "linear", lambda = lambda)
  coef <- as.numeric(rbind(fit1$b0, fit1$beta))
  coeff <- as.numeric(rbind(fit2$intercept, fit2$beta))
  expect_equal(coef, coeff, tolerance = 1e-1)
  
  
  set.seed(20200110)
  n <- 100
  p <- 48
  g <- 6
  X <- matrix(rnorm(n*p), nrow = n)
  for (i in seq(1, p, by = 6))
    X[, i:(i + g - 1)] <- svd(X[, i:(i + g - 1)])$u
  beta <- c(rep(5, 6), rep(-5, 6), rep(4, 6), rep(-4, 6), rep(3, 6), rep(-3, 6), rep(2, 6), rep(-2, 6))
  y <- X %*% beta + rnorm(n, sd = .1)
  
  data = list(x = X, y = y)
  group <- rep(1:8,each=6)
  
  fit1 = SGL(data, index = group, type = "linear")
  lambda <- fit1$lambda[1:5]
  
  fit1 <- sparsegl(X, y, asparse = 0.95, lambda = lambda, group = group)
  fit2 <- SGL(data, index = group, type = "linear", lambda = lambda)
  coef <- as.numeric(rbind(fit1$b0, fit1$beta))
  coeff <- as.numeric(rbind(fit2$intercept, fit2$beta))
  expect_equal(coef, coeff, tolerance = 1e-1)

  
})




