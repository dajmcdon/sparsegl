set.seed(1010)
n <- 100
p <- 10
X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)
eps <- rnorm(n, mean = 0, sd=1)
beta_star = c(rep(5,5),rep(0,5))
# beta_star <- c(rep(5,5), c(5,-5,2,0,0), rep(-5,5), c(2,-3,8,0,0), rep(0,(p-20)))
y <- X%*%beta_star + eps
groups <- rep(1:(p/5), each=5)
sgl <- sparsegl(x = X, y = y, group = groups, pen = "sparsegl", algorithm = "fourstep")
group_norm <- function(x,grp) sum(by(x, grp, function(x) sqrt(sum(x^2))))
sparse_grp_norm <- function(x, grp, alp=.05) {
  group_norm(x,grp)*(1-alp) + alp*sum(abs(x))
}

b = apply(sgl$beta, 2, sparse_grp_norm, grp=groups)
matplot(b/max(b), t(sgl$beta),ty='l',lty=1,col=groups)

Xs = X
Xs[abs(Xs) < 1.5] = 0
ys = Xs %*% beta_star + eps
sgl_sx <- sparsegl(x = Xs, y = ys, group = groups, pen = "sparsegl",
                   algorithm = "fourstep")
bs = apply(sgl_sx$beta, 2, sparse_grp_norm, grp=groups)
matplot(bs/max(bs), t(sgl_sx$beta),ty='l',lty=1,col=groups)

library(Matrix)
Xs_sp = Matrix(Xs, sparse=TRUE)

t1 <- sparsegl(Xs_sp, ys, groups, pen="sparsegl")
bt <- apply(t1$beta, 2, sparse_grp_norm, grp=groups)
