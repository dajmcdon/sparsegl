set.seed(1)
nobs <- 100L
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


test_that("wls provides the same result as sparsegl, Gaussian family", {
  res1 <- sparsegl(x, y, group, lambda = .1)
  res2 <- sparsegl(x, y, group, family = gaussian(), lambda = .1)

  expect_lt(max(abs(coef(res1) - coef(res2))), 1e-4)

  res1 <- sparsegl(x, y, group, lambda = .025)
  res2 <- sparsegl(x, y, group, family = gaussian(), lambda = .025)

  expect_lt(max(abs(coef(res1) - coef(res2))), 1e-4)

  res1_lam <- sparsegl(x, y, group)
  res2_lam <- sparsegl(x, y, group, family = gaussian())
  nlam <- length(res2_lam$lambda)
  expect_lt(max(abs(res1_lam$lambda[1:nlam] - res2_lam$lambda)), 1e-10)

  expect_lt(mean(abs(coef(res1_lam)[,1:nlam] - coef(res2_lam))), 1e-4)

  ## sparse case
  res1 <- sparsegl(xsp, y, group, lambda = .1)
  res2 <- sparsegl(xsp, y, group, family = gaussian(), lambda = .1)
  expect_lt(max(abs(coef(res1) - coef(res2))), 1e-4)

  res1 <- sparsegl(xsp, y, group, lambda = .025)
  res2 <- sparsegl(xsp, y, group, family = gaussian(), lambda = .025)
  expect_lt(max(abs(coef(res1) - coef(res2))), 1e-4)

  res1_lam <- sparsegl(xsp, y, group)
  res2_lam <- sparsegl(xsp, y, group, family = gaussian())
  nlam <- seq_along(res2_lam$lambda)
  expect_lt(max(abs(res1_lam$lambda[nlam] - res2_lam$lambda)), 1e-10)

  tt <- abs(coef(res1_lam)[,nlam] - coef(res2_lam)) /
    (1 + abs(coef(res1_lam)[,nlam]))
  expect_lt(max(tt), 1e-3)
})

test_that("wls provides the same result as sparsegl, binomial family", {
  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .01)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .01, eps = 1e-10)
  tt <- abs(coef(res1)[-1] - coef(res2)[-1]) / (1 + abs(coef(res1)[-1]))
  expect_true(all(tt < 1e-3))

  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .0025)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .0025, eps = 1e-10)
  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_true(mean(tt) < 3e-3)

  nlam <- 35
  res1_lam <- sparsegl(x, ybin, group, family = "binomial", nlambda = nlam,
                       lambda.factor = 0.01)
  res2_lam <- sparsegl(x, ybin, group, family = binomial(), eps = 1e-10,
                       nlambda = nlam, lambda.factor = 0.01)
  expect_lt(max(abs(res1_lam$lambda - res2_lam$lambda)), 1e-8)
  tt <- abs(coef(res1_lam) - coef(res2_lam)) / (1 + abs(coef(res1_lam)))
  expect_true(all(colMeans(tt) < 0.1))


  ## sparse case
  res1 <- sparsegl(xsp, ybin, group, family = "binomial", lambda = .01)
  res2 <- sparsegl(xsp, ybin, group, family = binomial(), lambda = .01, eps = 1e-9)
  tt <- abs(coef(res1)[-1] - coef(res2)[-1]) / (1 + abs(coef(res1)[-1]))
  expect_true(max(tt) < 1e-3)

  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .0025)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .0025, eps = 1e-10)
  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_true(mean(tt) < 3e-3)

  nlam <- 35
  res1_lam <- sparsegl(x, ybin, group, family = "binomial", nlambda = nlam,
                       lambda.factor = 0.01, eps = 1e-10)
  res2_lam <- sparsegl(x, ybin, group, family = binomial(), nlambda = nlam,
                       eps = 1e-10, lambda.factor = 0.01)
  expect_lt(max(abs(res1_lam$lambda - res2_lam$lambda)), 1e-8)
  tt <- abs(coef(res1_lam) - coef(res2_lam)) /
    (1 + abs(coef(res1_lam)))
  expect_true(all(colMeans(tt) < 0.1))
})

test_that("wls sparse and dense cases give the same results, Gaussian", {

  res1 <- sparsegl(x, y, group, family = gaussian())
  res2 <- sparsegl(xsp, y, group, family = gaussian())

  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_lt(max(tt), 1e-7)

  res1 <- sparsegl(x, ybin, group, family = binomial(), nlambda = 10,
                   lambda.factor = .01)
  res2 <- sparsegl(xsp, ybin, group, family = binomial(), nlambda = 10,
                   lambda.factor = .01)
  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_lt(max(tt), 1e-6)
})
