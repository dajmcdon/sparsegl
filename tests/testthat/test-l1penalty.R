test_that("Test if function throws error as expected or automatically perform
  adjustment to the given data", {
    n <- 100
    beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
    gr <- rep(1:5, each = 3)
    X <- matrix(rnorm(n * length(beta)), n)
    y <- X %*% beta + rnorm(n)

    # length of pfl1 is not the same as the number of predictors.
    pfl1 <- rep(1, 10)
    expect_error(sparsegl(X, y, group = gr, pf_sparse = pfl1))


    # any entry of plf1 is negative
    index <- sample(seq(15), size = 1)
    pfl1[index] <- -1
    expect_error(sparsegl(X, y, group = gr, pf_sparse = pfl1))


    # function will rescale each entry such that the sum will be constant
    pfl1_1 <- rep(1, 15)
    pfl1_2 <- rep(2, 15)
    out1 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_1)
    out2 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_2)
    expect_equal(out1$b0, out2$b0)
    expect_equal(out1$beta, out2$beta)
    expect_equal(out1$lambda, out2$lambda)

    # test if logistic model works fine
    n <- 100
    beta <- c(5,5,5,-5,-5,-5,1,0,1,0,0,0,0,2,0)
    gr <- rep(1:5, each = 3)
    X <- matrix(rnorm(n * length(beta)), n)
    pr <- 1 / (1 + exp(-X %*% beta))
    y0 <- rbinom(n, 1, pr)

    out1 <- sparsegl(X, y0, group = gr, family = "binomial")
    out2 <- sparsegl(X, y0, group = gr, family = "binomial", pf_sparse = rep(2, 15))
    expect_equal(out1$beta, out2$beta)

  })

test_that("function behaviors changed by pfl1", {
  n <- 100
  beta <- c(5, 5 , 5, -5, -5, -5)
  gr <- rep(1:2, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)
  pfl1_1 <- c(rep(0, 3), rep(10, 3))
  pfl1_2 <- c(rep(0, 3), rep(20, 3))
  out1 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_1)
  out2 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_2)

  # out1 and out2 are the same models.
  expect_equal(out1$b0, out2$b0)
  expect_equal(out1$beta, out2$beta)
  expect_equal(out1$lambda, out2$lambda)

  # will ignore group sparse penalty from this point.

  pfl1_2 <- c(rep(1, 5), 10)
  out1 <- sparsegl(X, y, group = gr, asparse = 1)
  out2 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_2, asparse = 1)


  # lambda is larger when the coefficient of predictor6 turns nonzero in out1.
  expect_true(out1$lambda[which(out1$beta[6, ] != 0)[1]] >=
                out2$lambda[which(out2$beta[6, ] != 0)[1]])


  # coefficient of predictor6 in out1 is larger than out1 with same lambdas.

  out2 <- sparsegl(X, y, group = gr, pf_sparse = pfl1_2, lambda = out1$lambda, asparse = 1)
  expect_equal(as.numeric(abs(as.numeric(out1$beta[6, ])) >= abs(as.numeric(out2$beta[6, ]))),
               rep(1, 100))


  # compare in a single model.

  n <- 100
  beta <- rep(5, 6)
  gr <- rep(1:2, each = 3)
  X <- matrix(rnorm(n * length(beta)), n)
  y <- X %*% beta + rnorm(n)

  # assign increasing l1-penalty along the predictors
  out <- sparsegl(X, y, group = gr, pf_sparse = seq(1, 60, by = 10), asparse = 1)
  expect_true(out$lambda[which(out$beta[1, ] != 0)[1]] >=
                out$lambda[which(out$beta[6, ] != 0)[1]])

})
