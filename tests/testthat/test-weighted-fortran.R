test_that("weights operate as expected for gaussian, binomial", {
  set.seed(1)
  nobs <- 1000L
  beta_star <- c(5, 5, 5, -5, -5, -5, 1, 0, 1, 0, 0, 0, 0, 2, 0)
  nvars <- length(beta_star)

  group <- rep(1:5, each = 3)
  x <- matrix(rnorm(nobs * length(beta_star)), nobs)
  x[abs(x) < 0.2] <- 0
  xsp <- as(x, "sparseMatrix")

  eps <- rnorm(nobs)
  y <- x %*% beta_star + eps

  pr <- 1 / (1 + exp(-x %*% beta_star))
  ybin <- rbinom(nobs, 1, pr)

  out <- sparsegl(x, y, group = group, lambda = 0, eps = 1e-16)
  expect_equal(as.vector(coef(out)), unname(coef(lm(y ~ x))))
  out <- sparsegl(x, y, group = group, lambda = 0, eps = 1e-16, intercept = FALSE)
  expect_equal(as.vector(coef(out)), c(0, unname(coef(lm(y ~ x - 1)))))

  out <- sparsegl(xsp, y, group = group, lambda = 0, eps = 1e-16)
  expect_equal(as.vector(coef(out)), unname(coef(lm(y ~ x))))
  out <- sparsegl(xsp, y, group = group, lambda = 0, eps = 1e-16, intercept = FALSE)
  expect_equal(as.vector(coef(out)), c(0, unname(coef(lm(y ~ x - 1)))))

  w <- runif(nobs)
  w <- w / sum(w) * nobs
  out <- sparsegl(x, y, group = group, lambda = 0, eps = 1e-16, weights = w)
  expect_equal(as.vector(coef(out)), unname(coef(lm(y ~ x, weights = w))))
  out <- sparsegl(x, y, group = group, lambda = 0, eps = 1e-16, intercept = FALSE, weights = w)
  expect_equal(as.vector(coef(out)), c(0, unname(coef(lm(y ~ x - 1, weights = w)))))

  out <- sparsegl(xsp, y, group = group, lambda = 0, eps = 1e-16, weights = w)
  expect_equal(as.vector(coef(out)), unname(coef(lm(y ~ x, weights = w))))
  out <- sparsegl(xsp, y, group = group, lambda = 0, eps = 1e-16, intercept = FALSE, weights = w)
  expect_equal(as.vector(coef(out)), c(0, unname(coef(lm(y ~ x - 1, weights = w)))))
})
