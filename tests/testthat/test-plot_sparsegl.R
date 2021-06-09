set.seed(1010)
n <- 100
p <- 200
X <- matrix(data = rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
eps <- rnorm(n, mean = 0, sd = 1)
beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), c(2, -3, 8, 0, 0), rep(0, (p - 20)))
y <- X %*% beta_star + eps
groups <- rep(1:(p / 5), each = 5)
fit1 <- sparsegl(X, y, group = groups)

test_that("no error in plotting for sparsegl object", {
  
  p1 <- try(plot(fit1, grouped = FALSE, x_axis = "log_lambda"))
  p2 <- try(plot(fit1, grouped = FALSE, x_axis = "norm"))
  p3 <- try(plot(fit1, grouped = TRUE, y_axis = "group_norm",
                 x_axis = "log_lambda"))
  p4 <- try(plot(fit1, grouped = TRUE, y_axis = "group_norm", x_axis = "norm"))
  p5 <- try(plot(fit1, grouped = TRUE, y_axis = "coef", x_axis = "log_lambda"))
  p6 <- try(plot(fit1, grouped = TRUE, y_axis = "coef", x_axis = "norm"))
  p <- list(p1, p2, p3, p4, p5, p6)
  
  for (i in seq(6)) {
    expect_false(inherits(p[[i]], "try-error"))
  }
})