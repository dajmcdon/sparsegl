
test_that("cv_sparsegl", {
  set.seed(1)
  
  data(bardet)
  bardet
  group <- rep(1:20,each=5)
  cv <- cv.sparsegl(x=bardet$x, y=bardet$y, group=group,pred.loss="L2", lambda.factor=0.05, nfolds=5)
  cv1 <- gglasso::cv.gglasso(x = bardet$x, y = bardet$y, group = group, pred.loss = "L2", lambda.factor=0.05, nfolds=5)
  expect_equal(cv$lambda.min, cv1$lambda.min, info = paste(cv$lambda.min, " vs ", cv1$lambda.min))
  expect_equal(cv$lambda.lse, cv1$lambda.lse)
  
  cv <- cv.sparsegl(x=bardet$x, y=bardet$y, group=group,pred.loss="L1", lambda.factor=0.05, nfolds=5)
  cv1 <- gglasso::cv.gglasso(x = bardet$x, y = bardet$y, group = group, pred.loss = "L1", lambda.factor=0.05, nfolds=5)
  expect_equal(cv$lambda.min, cv1$lambda.min, info = paste(cv$lambda.min, " vs ", cv1$lambda.min))
  expect_equal(cv$lambda.lse, cv1$lambda.lse)

  n <- 100
  p <- 75
  g <- 5
  X <- matrix(rnorm(n*p), nrow = n)
  for (i in seq(1, p, by = 5))
    X[, i:(i + g - 1)] <- svd(X[, i:(i + g - 1)])$u
  beta <- c(rep(5, 5), rep(-5, 5), rep(0, p - 10))
  y <- X %*% beta + rnorm(n, sd = .1)
  data = list(x = X, y = y)
  group <- rep(1:15,each=5)
  
  cv <- cv.sparsegl(x = X, y = y, group = group, pred.loss = "L2",lambda.factor=0.001, nfold=5)
  cv1 <- gglasso::cv.gglasso(x = X, y = y, group = group, pred.loss = "L2", lambda.factor=0.001, nfolds=5)
  expect_equal(cv$lambda.min, cv1$lambda.min, info = paste(cv$lambda.min, " vs ", cv1$lambda.min))
  expect_equal(cv$lambda.lse, cv1$lambda.lse)
  
  cv <- cv.sparsegl(x = X, y = y, group = group, pred.loss = "L1",  lambda.factor=0.001, nfold=5)
  cv1 <- gglasso::cv.gglasso(x = X, y = y, group = group, pred.loss = "L1", lambda.factor=0.001, nfolds=5)
  expect_equal(cv$lambda.min, cv1$lambda.min, info = paste(cv$lambda.min, " vs ", cv1$lambda.min))
  expect_equal(cv$lambda.lse, cv1$lambda.lse)
})

