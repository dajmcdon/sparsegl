library(Matrix)
library(tibble)
library(SGL)
library(sparsegl)


nrepls <- 5
n <- 500
p <- c(100, 250, 500, 1000, 2500)
b <- c(rep(1, 5), rep(-1, 5), c(-1,-1,0,0,0), c(1,1,0,0,0), rep(1, 5),
       rep(0, 75)) * 10
s <- .1
gr <- rep(1:(max(p) / 5), each = 5)
x <- matrix(rnorm(n*max(p)), nrow = n)

x[as.logical(rbinom(n*max(p), 1, 1 - s))] <- 0
xsp <- Matrix(x)
mu <- sapply(p, function(z) x[ ,1:z] %*% rep(b, length.out = z) / sqrt(z))
# prob <- 1 / (1 + exp(-mu))

signal <- sqrt(colSums(mu^2))
noise_sd <- sqrt(signal)

y <- mu + rnorm(n*length(p), sd = rep(noise_sd, each = n))

res <- tibble(method = "a", time = proc.time()["elapsed"], .rows = 0)
for (i in seq_along(p)) {
  pp <- seq(p[i])
  dat <- list(y = y[ ,i], x = x[ ,pp])
  xxsp <- xsp[ ,pp]
  g <- gr[pp]
  for (j in seq(nrepls)) {
    s1 <- system.time(
      SGL(dat, nlam = 100, alpha = 0.05, index = g, standardize = FALSE)
    )
    res <- add_row(res, method = "SGL", time = s1["elapsed"])
    s2 <- system.time(
      sparsegl(dat$x, dat$y, g, standardize = FALSE)
    )
    res <- add_row(res, method = "sparsegl", time = s2["elapsed"])
    s3 <- system.time(
      sparsegl(xxsp, dat$y, g, standardize = FALSE)
    )
    res <- add_row(res, method = "sparsegl_sp", time = s3["elapsed"])
    print(paste("Done with p = ", p[i], "repl = ", j))
  }
}

res$p <- rep(p, each = nrepls * 3)
saveRDS(res, "jss-sparsegl/large-data/sparsegl-timing.rds")
