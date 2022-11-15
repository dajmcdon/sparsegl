set.seed(1)
nobs <- 100L
beta_star <- c(5, 5, 5, -5, -5, -5, 1, 0, 1, 0, 0, 0, 0, 2, 0)
nvars <- length(beta_star)

group <- rep(1:5, each = 3)
x <- matrix(rnorm(nobs * length(beta_star)), nobs)
x[abs(x) < 0.2] <- 0
xsp <- as(x, "sparseMatrix")

y <- x %*% beta_star + rnorm(nobs)
ysp <- xsp %*% beta_star + rnorm(nobs)

pr <- 1 / (1 + exp(-x %*% beta_star))
ybin <- rbinom(nobs, 1, pr)

pr <- as.double(1 / (1 + exp(-xsp %*% beta_star)))
ybinsp <- rbinom(nobs, 1, as.matrix(pr))


test_that("wls provides the same result as sparsegl, Gaussian family", {
  res1 <- sparsegl(x, y, group, lambda = .1)
  res2 <- sparsegl(x, y, group, family = gaussian(), lambda = .1)

  expect_equal(as.numeric(coef(res1)), as.numeric(coef(res2)), tolerance = 1e-4)

  res1 <- sparsegl(x, y, group, lambda = .025)
  res2 <- sparsegl(x, y, group, family = gaussian(), lambda = .025)

  expect_equal(as.numeric(coef(res1)), as.numeric(coef(res2)), tolerance = 1e-4)

  res1_lam <- sparsegl(x, y, group)
  res2_lam <- sparsegl(x, y, group, family = gaussian())
  nlam <- length(res2_lam$lambda)
  expect_equal(res1_lam$lambda[1:nlam], res2_lam$lambda, tolerance = 1e-4)

  expect_equal(as.numeric(coef(res1_lam)[,1:nlam]),
               as.numeric(coef(res2_lam)),
               tolerance = 1e-3)

  ## sparse case
  res1 <- sparsegl(xsp, ysp, group, lambda = .1)
  res2 <- sparsegl(xsp, ysp, group, family = gaussian(), lambda = .1)

  expect_equal(as.numeric(coef(res1)), as.numeric(coef(res2)), tolerance = 1e-4)

  res1 <- sparsegl(xsp, ysp, group, lambda = .025)
  res2 <- sparsegl(xsp, ysp, group, family = gaussian(), lambda = .025)

  expect_equal(as.numeric(coef(res1)), as.numeric(coef(res2)), tolerance = 1e-4)

  nlam <- 1:60
  res1_lam <- sparsegl(xsp, ysp, group)
  res2_lam <- sparsegl(xsp, ysp, group, family = gaussian(),
                       lambda = res1_lam$lambda[nlam])
  tt <- abs(coef(res1_lam)[,nlam] - coef(res2_lam)) /
    (1 + abs(coef(res1_lam)[,nlam]))
  expect_lt(max(tt), 1e-4)
})

test_that("wls provides the same result as sparsegl, binomial family", {
  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .01)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .01, eps = 1e-10)
  tt <- abs(coef(res1)[-1] - coef(res2)[-1]) / (1 + abs(coef(res1)[-1]))
  expect_true(all(tt < 1e-3))

  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .0025)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .0025, eps = 1e-10)
  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_true(mean(tt) < 1e-3)

  nlam <- 1:20
  res1_lam <- sparsegl(x, ybin, group, family = "binomial")
  res2_lam <- sparsegl(x, ybin, group, family = binomial(),
                       lambda = res1_lam$lambda[nlam], eps = 1e-8)
  tt <- abs(coef(res1_lam)[,nlam] - coef(res2_lam)) /
    (1 + abs(coef(res1_lam)[,nlam]))
  expect_true(all(colMeans(tt) < 1e-3))


  ## sparse case
  res1 <- sparsegl(xsp, ybinsp, group, family = "binomial", lambda = .01)
  res2 <- sparsegl(xsp, ybinsp, group, family = binomial(), lambda = .01, eps = 1e-9)
  tt <- abs(coef(res1)[-1] - coef(res2)[-1]) / (1 + abs(coef(res1)[-1]))
  expect_true(mean(tt) < 1e-3)

  res1 <- sparsegl(x, ybin, group, family = "binomial", lambda = .0025)
  res2 <- sparsegl(x, ybin, group, family = binomial(), lambda = .0025, eps = 1e-10)
  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_true(mean(tt) < 1e-3)

  nlam <- 1:20
  res1_lam <- sparsegl(x, ybin, group, family = "binomial")
  res2_lam <- sparsegl(x, ybin, group, family = binomial(),
                       lambda = res1_lam$lambda[nlam], eps = 1e-8)
  tt <- abs(coef(res1_lam)[,nlam] - coef(res2_lam)) /
    (1 + abs(coef(res1_lam)[,nlam]))
  expect_true(all(colMeans(tt) < 1e-3))

})

test_that("wls sparse and dense cases give the same results, Gaussian", {
  xsp1 <- x
  xsp1[abs(xsp1) < 0.2] <- 0
  xsp2 <- as(xsp1, "sparseMatrix")
  y <- rnorm(nobs)

  res1 <- sparsegl(xsp1, y, group, family = gaussian())
  res2 <- sparsegl(xsp2, y, group, family = gaussian())

  tt <- abs(coef(res1) - coef(res2)) / (1 + abs(coef(res1)))
  expect_lt(max(tt), 1e-6)
})
