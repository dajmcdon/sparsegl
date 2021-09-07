test_that("compare for logistic regression", {
  skip("test later")
  set.seed(123)
  n <- 100
  p <- 20
  X1 <- matrix(data = rnorm(n/4*p, mean = 2, sd = 5), nrow = n/4, ncol = p)
  X2 <- matrix(data = rnorm(n/4*p, mean = 5, sd = 1), nrow = n/4, ncol = p)
  X3 <- matrix(data = rnorm(n/4*p, mean = 10, sd = 3), nrow = n/4, ncol = p)
  X4 <- matrix(data = rnorm(n/4*p, mean = 15, sd = 2), nrow = n/4, ncol = p)
  X <- rbind(X1, X2, X3, X4)
  # X <- matrix(data = rnorm(n*p, mean = 2, sd = 5), nrow = n, ncol = p)
  beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, 5))
  groups <- rep(1:(p / 5), each = 5)
  
  pr <- 1 / (1 + exp(-X %*% beta_star))
  y0 <- rbinom(n, 1, pr)
  
  y1 <- y0
  y1[which(y1 == 0)] <- -1
  
  # false standardization, false intercept
  fit_sparsegl <- sparsegl(X, y1, group = groups, loss = "logit", standardize = FALSE, asparse = 0, intercept = FALSE)
  fit_gglasso <- gglasso::gglasso(X, y1, group = groups, loss = "logit", intercept = FALSE, lambda = fit_sparsegl$lambda)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # false standardization, true intercept
  fit_sparsegl <- sparsegl(X, y1, group = groups, loss = "logit", standardize = FALSE, asparse = 0, intercept = TRUE)
  fit_gglasso <- gglasso::gglasso(X, y1, group = groups, loss = "logit", intercept = TRUE, lambda = fit_sparsegl$lambda)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  expect_equal(coef1, coef2, tolerance = 1e-2)

  # true standardization, false intercept
  fit_sparsegl <- sparsegl(X, y1, group = groups, loss = "logit", standardize = TRUE, asparse = 0, intercept = FALSE)
  sx <- sqrt(Matrix::colSums(X^2))
  sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
  xs <- 1 / sx
  xt <- matrix(X %*% Matrix::Diagonal(x = xs), nrow = 100, byrow = FALSE)
  fit_gglasso <- gglasso::gglasso(xt, y1, group = groups, loss = "logit", intercept = FALSE, lambda = fit_sparsegl$lambda)
  fit_gglasso$beta <- fit_gglasso$beta * xs
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  expect_equal(coef1, coef2, tolerance = 1e-2)
  
  # true standardization, true intercept
})