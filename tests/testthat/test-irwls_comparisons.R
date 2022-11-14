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

pr <- 1 / (1 + exp(-xsp %*% beta_star))
ybinsp <- rbinom(nobs, 1, as.matrix(pr))


test_that("sgl_irwls provides the same result as sparsegl, gaussian family", {
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


})

test_that("sgl_irwls provides the same result as sparsegl, gaussian family", {
    # Dense matrix
    res1 <- sgl_irwls(
        bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1, dfmax, pmax, nlambda,
        flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
        lower_bnd, upper_bnd,
        family = gaussian()
    )

    res2 <- sparsegl(
        x, y, group, "gaussian", nlambda, flambda, NULL, pf, pfl1, intercept,
        asparse, standardize, lower_bnd, upper_bnd, eps, maxit
    )

    expect_equal(
        as.numeric(res1$coefficients),
        as.numeric(res2$coefficients),
        tolerance = 1e-10
    )

    # Sparse matrix
    res1 <- sgl_irwls(
        bn, bs, ix, iy, nobs, nvars, xsp, ysp, pf, pfl1, dfmax, pmax, nlambda,
        flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
        lower_bnd, upper_bnd,
        family = gaussian()
    )

    res2 <- sparsegl(
        xsp, ysp, group, "gaussian", nlambda, flambda, NULL, pf, pfl1,
        intercept, asparse, standardize, lower_bnd, upper_bnd, eps, maxit
    )

    expect_equal(
        as.numeric(res1$coefficients),
        as.numeric(res2$coefficients),
        tolerance = 1e-6
    )
})

test_that("sgl_irwls provides the same result as sparsegl, binomial family", {
    # Dense matrix
    res1 <- sgl_irwls(
        bn, bs, ix, iy, nobs, nvars, x, ybin, pf, pfl1, dfmax, pmax, nlambda,
        flmin, ulam, eps, maxit, vnames, group, intr, asparse, standardize,
        lower_bnd, upper_bnd,
        family = binomial()
    )

    res2 <- sparsegl(
        x, ybin, group, "binomial", nlambda, flambda, NULL, pf, pfl1, intercept,
        asparse, standardize, lower_bnd, upper_bnd, eps, maxit
    )

    expect_equal(
        as.numeric(res1$coefficients),
        as.numeric(res2$coefficients),
        tolerance = 1e-6
    )

    # Sparse matrix
    res1 <- sgl_irwls(
        bn, bs, ix, iy, nobs, nvars, xsp, ybinsp, pf, pfl1, dfmax, pmax,
        nlambda, flmin, ulam, eps, maxit, vnames, group, intr, asparse,
        standardize, lower_bnd, upper_bnd,
        family = binomial()
    )

    res2 <- sparsegl(
        xsp, ybinsp, group, "binomial", nlambda, flambda, NULL, pf, pfl1,
        intercept, asparse, standardize, lower_bnd, upper_bnd, eps, maxit
    )

    expect_equal(
        as.numeric(res1$coefficients),
        as.numeric(res2$coefficients),
        tolerance = 1e-6
    )
})
