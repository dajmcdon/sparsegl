test_that("compare for linear regression model", {
  
  set.seed(123)
  n <- 100
  p <- 20
  X1 <- matrix(data = rnorm(n/4*p, mean = 2, sd = 5), nrow = n/4, ncol = p)
  X2 <- matrix(data = rnorm(n/4*p, mean = 5, sd = 1), nrow = n/4, ncol = p)
  X3 <- matrix(data = rnorm(n/4*p, mean = 10, sd = 3), nrow = n/4, ncol = p)
  X4 <- matrix(data = rnorm(n/4*p, mean = 15, sd = 2), nrow = n/4, ncol = p)
  X <- rbind(X1, X2, X3, X4)
  beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, 5))
  groups <- rep(1:(p / 5), each = 5)
 
  eps <- rnorm(n, mean = 0, sd = 1)
  y <- X %*% beta_star + eps
  
  # no penalty
  fit_sparsegl <- sparsegl(X, y, lambda = 0)
  fit_lm <- lm(y~X)
  coef1 <- as.numeric(c(fit_sparsegl$b0[1], fit_sparsegl$beta[, 1]))
  coef2 <- as.numeric(fit_lm$coefficients)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # false standardization, false intercept
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = FALSE,
                           asparse = 0, intercept = FALSE)
  fit_gglasso <- gglasso::gglasso(X, y, group = groups, loss = "ls", 
                                  lambda = fit_sparsegl$lambda, intercept = FALSE)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # false standardization, true intercept
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = FALSE,
                           asparse = 0, intercept = TRUE)
  fit_gglasso <- gglasso::gglasso(X, y, group = groups, loss = "ls", 
                                  lambda = fit_sparsegl$lambda, intercept = TRUE)
  fit_SGL <- SGL::SGL(list(x = X, y = y), index = groups, type = "linear", 
                      standardize = FALSE, alpha = 0, lambdas = fit_sparsegl$lambda,
                      nlam = 100)
  intr1 <- as.numeric(fit_sparsegl$b0)
  intr2 <- as.numeric(fit_gglasso$b0)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  coef3 <- as.numeric(fit_SGL$beta)
  expect_equal(intr1, intr2, tolerance = 1e-3)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  expect_equal(coef1, coef3, tolerance = 1e-3)
  
    # with both penalties:
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = FALSE,
                           asparse = 0.5, intercept = TRUE)
  fit_SGL <- SGL::SGL(list(x = X, y = y), index = groups, type = "linear", 
                      standardize = FALSE, alpha = 0.5, lambdas = fit_sparsegl$lambda,
                      nlam = 100)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_SGL$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # true standardization, false intercept
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = TRUE,
                           asparse = 0, intercept = FALSE)
    # manually standardization for gglasso
  sx <- sqrt(Matrix::colSums(X^2))
  sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
  xs <- 1 / sx
  xt <- matrix(X %*% Matrix::Diagonal(x = xs), nrow = 100, byrow = FALSE)
  fit_gglasso <- gglasso::gglasso(xt, y, group = groups, loss = "ls", 
                                  lambda = fit_sparsegl$lambda, intercept = FALSE)
  fit_gglasso$beta <- fit_gglasso$beta * xs
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # true standardization, true intercept
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = TRUE,
                           asparse = 0, intercept = TRUE)
    # manually standardization
  ym <- mean(y)
  yt <- y - ym
  xm <- colMeans(X)
  xt <- sweep(X,2,xm)
  sx <- sqrt(Matrix::colSums(xt^2))
  sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
  xs <- 1 / sx
  xt <- matrix(xt %*% Matrix::Diagonal(x = xs), nrow = 100, byrow = FALSE)
  fit_gglasso <- gglasso::gglasso(xt, yt, group = groups, loss = "ls", 
                                  lambda = fit_sparsegl$lambda, intercept = TRUE)
  fit_gglasso$beta <- fit_gglasso$beta * xs
  fit_gglasso$b0 <- ym - xm %*% fit_gglasso$beta
  fit_SGL <- SGL::SGL(list(x = xt, y = yt), index = groups, type = "linear",
                      standardize = FALSE, alpha = 0, lambdas = fit_sparsegl$lambda,
                      nlam = 100)
  fit_SGL$beta <- fit_SGL$beta*xs
  intr1 <- as.numeric(fit_sparsegl$b0)
  intr2 <- as.numeric(fit_gglasso$b0)
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_gglasso$beta)
  coef3 <- as.numeric(fit_SGL$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  expect_equal(coef1, coef3, tolerance = 1e-3)
  
    # with both penalties:
  fit_sparsegl <- sparsegl(X, y, group = groups, loss = "ls", standardize = TRUE,
                           asparse = 0.5, intercept = TRUE)
  fit_SGL <- SGL::SGL(list(x = xt, y = yt), index = groups, type = "linear", 
                      standardize = FALSE, alpha = 0.5, lambdas = fit_sparsegl$lambda,
                      nlam = 100)
  fit_SGL$beta <- fit_SGL$beta * xs
  coef1 <- as.numeric(fit_sparsegl$beta)
  coef2 <- as.numeric(fit_SGL$beta)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  
  # sparse matrix case
  Xs <- Matrix::rsparsematrix(nrow = 100, ncol = 20, density = 0.4)
  ys <- Xs %*% beta_star + eps
  ys <- as.vector(ys)
  fit_sparse <- sparsegl(Xs, ys, group = groups, loss = "ls", standardize = TRUE,
                         intercept = TRUE)
  fit_nonsparse <- sparsegl(as.matrix(Xs), ys, group = groups, loss = "ls", 
                            standardize = TRUE, intercept = TRUE)
  coef1 <- as.numeric(fit_sparse$beta)
  coef2 <- as.numeric(fit_nonsparse$beta)
  intr1 <- as.numeric(fit_sparse$b0)
  intr2 <- as.numeric(fit_nonsparse$b0)
  expect_equal(coef1, coef2, tolerance = 1e-3)
  expect_equal(intr1, intr2, tolerance = 1e-2)
})
