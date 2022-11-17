library(glmnet)

test_that("tests for linear regression model", {

  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n)
  eps <- rnorm(n)
  beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
  y <- X %*% beta_star + eps
  groups <- rep(1:(p / 5), each = 5)
  fit1 <- sparsegl(X, y, group = groups, eps = 1e-8, maxit = 3e8)

  # expect error if type is not nonzero or coefficients
  expect_error(predict(fit1, type = "response"))
  expect_error(predict(fit1, type = "link"))

  # expect error when type is class
  expect_error(predict(fit1, newx = X[10, ], type = "class"))

  # expect setting type with response and link would give the same results
  expect_equal(predict(fit1, newx = X[10, ], type = "response"),
               predict(fit1, newx = X[10, ], type = "link"))
  expect_equal(predict(fit1, newx = X[10:15, ], type = "response"),
               predict(fit1, newx = X[10:15, ], type = "link"))

  # test for a simplist case: set lambda equal to zero
  set.seed(1)
  X <- matrix(c(rnorm(6)), nrow = 3)
  beta <- c(2, 1)
  y <- X %*% beta + rnorm(3, sd = .1)
  res1 <- sparsegl(X, y, lambda = 0, eps = 1e-8, maxit = 3e8)
  res2 <- summary(lm(y~X))

  expect_equal(as.numeric(predict(res1, type = "coefficients")),
               as.numeric(res2$coefficients[, 1]), tolerance = 1e-3)
  expect_equal(predict(res1, type = "nonzero")[, 1], c(1, 2))
  expect_equal(predict(res1, newx = X, type = "response")[, 1],
               as.numeric(predict(lm(y~X))), tolerance = 1e-3)
  expect_equal(predict(res1, newx = X, type = "link")[, 1],
               as.numeric(predict(lm(y~X))), tolerance = 1e-3)
})

test_that("tests for logistic regression model", {
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n)
  beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
  groups <- rep(1:(p / 5), each = 5)
  pr <- 1 / (1 + exp(-X %*% beta_star))
  y0 <- rbinom(n, 1, pr)
  y1 <- y0
  y1[which(y1 == 0)] <- -1
  fit2 <- sparsegl(X, y1, group = groups, family = "binomial", eps = 1e-8)

  # expect error if type is not nonzero or coefficients
  expect_error(predict(fit2, type = "response"))
  expect_error(predict(fit2, type = "link"))


  set.seed(1)
  X <- matrix(c(rnorm(50)), nrow = 25)
  beta <- c(2, 1)
  y <- X %*% beta
  pr <- 1 / (1 + exp(-X %*% beta))
  y0 <- rbinom(25, 1, pr)
  y1 <- y0
  y1[which(y1 == 0)] <- -1
  res1 <- sparsegl(X, y1, lambda = 0, family = "binomial", eps = 1e-8)
  res2 <- glmnet(X, y1, family = "binomial", lambda = 0)

  # check whether y can be corrected if y is a numeric vector
  res3 <- sparsegl(X, y0, lambda = 0, family = "binomial", eps = 1e-8)
  chars <- y0
  chars[chars == 0] <- 'a'
  chars[chars == 1] <- 'b'
  res4 <- sparsegl(X, chars, lambda = 0, family = "binomial", eps = 1e-8)

  expect_equal(as.numeric(predict(res1, type = "coefficients")[, 1]),
               as.numeric(c(res2$a0, res2$beta[, 1])), tolerance = 1e-3)
  expect_equal(predict(res1, newx = X, type = "link"), predict(res2, newx = X),
               tolerance = 1e-3)
  expect_equal(predict(res1, newx = X, type = "response"),
               1 / (1 + exp(-predict(res2, newx = X))), tolerance = 1e-3)
  expect_equal(predict(res1, newx = X, type = "class"),
               res2$classnames[ifelse(predict(res2, newx = X) > 0, 2, 1)])

  expect_equal(as.numeric(predict(res1, type = "coefficients")[, 1]),
               as.numeric(predict(res3, type = "coefficients")[, 1]))
  expect_equal(as.numeric(predict(res1, type = "coefficients")[, 1]),
               as.numeric(predict(res4, type = "coefficients")[, 1]))
})

n <- 100
p <- 20
X <- matrix(rnorm(n * p), nrow = n)
Xs <- X
Xs[abs(X) < .05] = 0
Xs <- Matrix::Matrix(Xs, sparse = TRUE)
eps <- rnorm(n)
beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
y <- X %*% beta_star + eps
groups <- rep(1:(p / 5), each = 5)
pr <- 1 / (1 + exp(-X %*% beta_star))
ybin <- rbinom(n, 1, pr)


test_that("intercept works in for all Gaussian cases", {
  fit <- sparsegl(X, y, groups)
  expect_silent(predict(fit, X))
  fit <- sparsegl(X, y, groups, intercept = FALSE)
  expect_silent(predict(fit, X))
  fit <- sparsegl(Xs, y, groups)
  expect_silent(predict(fit, Xs))
  fit <- sparsegl(Xs, y, groups, intercept = FALSE)
  expect_silent(predict(fit, Xs))
})

test_that("intercept works in for all Logistic cases", {
  fit <- sparsegl(X, ybin, groups, family = "binomial")
  expect_silent(predict(fit, X))
  fit <- sparsegl(X, ybin, groups, family = "binomial", intercept = FALSE)
  expect_silent(predict(fit, X))
  fit <- sparsegl(Xs, ybin, groups, family = "binomial")
  expect_silent(predict(fit, Xs))
  fit <- sparsegl(Xs, ybin, groups, family = "binomial", intercept = FALSE)
})
