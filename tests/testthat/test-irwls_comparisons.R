set.seed(1)
nobs <- 100L
beta_star <- c(5, 5, 5, -5, -5, -5, 1, 0, 1, 0, 0, 0, 0, 2, 0)
nvars <- length(beta_star)

group <- rep(1:5, each = 3)
x <- matrix(rnorm(nobs * length(beta_star)), nobs)
y <- x %*% beta_star + rnorm(nobs)

pr <- 1 / (1 + exp(-x %*% beta_star))
ybin <- rbinom(nobs, 1, pr)

bn <- as.integer(max(group))
bs <- as.integer(as.numeric(table(group)))
iy <- cumsum(bs)
ix <- c(0, iy[-bn]) + 1
ix <- as.integer(ix)
iy <- as.integer(iy)

pf <- as.double(sqrt(bs))
pfl1 <- rep(as.double(pf / sum(pf) * nvars), 3)
dfmax <- as.integer(max(group)) + 1L
pmax <- as.integer(min(dfmax * 1.2, as.integer(max(group))))

nlambda <- 100L
flambda <- ifelse(nobs < nvars, 0.01, 1e-04)
flmin <- as.double(flambda)
ulam <- double(1)
eps <- as.double(1e-08)
maxit <- as.integer(3e+08)

vnames <- colnames(x)
intercept <- as.integer(TRUE)
asparse <- as.double(0.05)
standardize <- TRUE
lower_bnd <- as.double(rep(-9.9e30, bn))
upper_bnd <- as.double(rep(9.9e30, bn))
intr <- as.integer(intercept)

test_that("sgl_irwls provides the same result as sparsegl, gaussian family", {
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
})

test_that("sgl_irwls provides the same result as sparsegl, binomial family", {
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
        tolerance = 1e-10
    )
})
