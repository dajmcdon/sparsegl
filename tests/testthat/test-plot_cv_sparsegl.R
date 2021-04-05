
test_that("plot_cv_sparsegl", {
  set.seed(1)
  
  skip("skip this test")
  data(bardet)
  bardet
  group <- rep(1:20,each=5)
  cv <- cv.sparsegl(x=bardet$x, y=bardet$y, group=group,pred.loss="L2", lambda.factor=0.05, nfolds=5)
  plot(cv)
  cv1 <- gglasso::cv.gglasso(x = bardet$x, y = bardet$y, group = group, pred.loss = "L2", lambda.factor=0.05, nfolds=5)
  plot(cv1)
  
  print(visualTest::isSimilar(file = "/Users/xiaoxuanliang/desktop/plots/sgl_l2.png", fingerprint = visualTest::getFingerprint(file = "/Users/xiaoxuanliang/desktop/plots/gglasso_l2.png"),
            threshold = 0.1))
  
  
  cv <- cv.sparsegl(x=bardet$x, y=bardet$y, group=group,pred.loss="L1", lambda.factor=0.05, nfolds=5)
  plot(cv)
  cv1 <- gglasso::cv.gglasso(x = bardet$x, y = bardet$y, group = group, pred.loss = "L1", lambda.factor=0.05, nfolds=5)
  plot(cv1)
  
  print(visualTest::isSimilar(file = "/Users/xiaoxuanliang/desktop/plots/sgl_l1.png", fingerprint = visualTest::getFingerprint(file = "/Users/xiaoxuanliang/desktop/plots/gglasso_l1.png"),
                              threshold = 0.1))
  
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
  
  cv <- cv.sparsegl(x = X, y = y, group = group, pred.loss = "L2", lambda.factor = 0.001, nfolds = 5)
  plot(cv)
  cv1 <- gglasso::cv.gglasso(x = X, y = y, group=group, pred.loss = "L2", lambda.factor=0.001, nfolds=5)
  plot(cv1)
  
  print(visualTest::isSimilar(file ="/Users/xiaoxuanliang/desktop/plots/sgl2_l2.png", fingerprint = visualTest::getFingerprint(file = "/Users/xiaoxuanliang/desktop/plots/gglasso2_l2.png"),
                        threshold = 0.1))
  
  cv <- cv.sparsegl(x = X, y = y, group = group, pred.loss = "L1", lambda.factor = 0.001, nfolds = 5)
  plot(cv)
  cv1 <- gglasso::cv.gglasso(x = X, y = y, group=group, pred.loss = "L1", lambda.factor=0.001, nfolds=5)
  plot(cv1)
  
  print(visualTest::isSimilar(file = "/Users/xiaoxuanliang/desktop/plots/sgl2_l1.png", fingerprint = visualTest::getFingerprint(file = "/Users/xiaoxuanliang/desktop/plots/gglasso2_l1.png"),
                              threshold = 0.1))
  
})


