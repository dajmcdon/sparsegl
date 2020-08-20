# For testing the lasso, group lasso, sparse group lasso
# In other words, the gglasso DA project

# Set the n and p values
set.seed(1010)
n <- 100
p <- 200

# There are different ways to generate a data matrix, I'm gonna be a dummy

X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)

eps <- rnorm(n, mean = 0, sd=1)

beta_star <- c(rep(5,5), c(5,-5,2,0,0), rep(-5,5), c(2,-3,8,0,0), rep(0,(p-20)))

y <- X%*%beta_star + eps

groups <- rep(1:(p/5), each=5)

# Compare threestep and threestep alt
mysparse1 <- sparsegl(x = X, y = y, group = groups, pen = "sparsegl", algorithm = "threestep")
mysparse2 <- sparsegl(x = X, y = y, group = groups, pen = "sparsegl", algorithm = "threestepalt")
beta1 <- mysparse1$beta
beta2 <- mysparse2$beta
max(beta2-beta1)
# Now we need to try the lassos

# out <- gglasso(X, y, group = groups, loss = 'ls')
out <- sparsegl(x = X, y = y, group = groups, pen = "gglasso", algorithm = "threestep")

out_sp <- gglasso(X,y, group = groups, loss="ls_sparse")

out_sparsegl <- gglasso(X,y, group = groups, loss="sparsegl")

## Dan's updated (better) versions
grp_norm <- function(x, grp, pf = sqrt(table(grp))){
  out = by(x,grp,function(x) sqrt(sum(x^2)))
  sum(out*pf)
}
sparse_grp_norm <- function(x, grp, alpha = .05, pf = sqrt(table(grp))){
  out1 = grp_norm(x,grp,pf)
  (1-alpha)*out1 + alpha*sum(abs(x))
}

# Dan's code
par(mfrow=c(2,2))
plot(out)
plot(out_sp)
plot(out_sparsegl)
b1 = apply(out$beta, 2, sd)
b2 = apply(out_sp$beta, 2, sd)
matplot(b1, t(out$beta), ty='l', lty=1, col = grp)
matplot(b2, t(out_sp$beta), ty='l', lty=1, col = grp)




###########################################################################
# Dan's code, redux (using group_norm instead of sd)
#group_norm <- function(x) sum(by(x, grp, function(x) sqrt(sum(x^2))))
#sp_group_norm <- function(x, alp=.05) group_norm(x)*(1-alp) + alp*sum(abs(x))


# 
# par(mfrow=c(1,2))
# b1 = apply(out$beta, 2, group_norm)
# b2 = apply(out_sp$beta, 2, sp_group_norm)
# matplot(b1, t(out$beta), ty='l', lty=1, col = grp)
# matplot(b2, t(out_sp$beta), ty='l', lty=1, col = grp)
# 

###########################################################################
###############################################################################
###########################################################################
# OK, let's play again with some data, for both the n<p and p<n cases
# Recall Dan said there are 4 cases: (sparsity/no sparcity within groups) X (p < n vs p > n)

set.seed(1010)
n <- 10000
p <- 1000



X <- matrix(data = rnorm(n*p, mean=0, sd=1), nrow = n, ncol = p)

eps <- rnorm(n, mean = 0, sd=1)

############This beta is NOT sparse#####################################
########################################################################
beta_star <- rep(0,p)
for (i in 0:(20-1)) {
  beta_star[(50*i + 1):(50*(i+1))] <- rep((-1)^i*i, 50)
  beta_star[(50*i + 20):(50*i + 29)] <- rep(0,10)
}
#########################################################################



# As per Dan's suggestion, we make a new beta with fewer nonzero coefficients 
# than observations. This is for n = 10,000 and p = 1000...
########################################################################
beta_star2 <- rep(0,p)
for (i in 0:(20-1)) {
  beta_star2[(50*i + 1):(50*(i+1))] <- rep((-1)^i*i, 50)
  beta_star2[(50*i + 10):(50*i + 39)] <- rep(0,30)
}
#########################################################################


grp <- rep(1:20, each=50)

y <- X%*%beta_star2 + eps


# Now we need to try the lassos
out <- gglasso(X,y,group=grp, loss='ls')

# Next we try ls_sparse
out_sp <- gglasso(X,y,group=grp, loss = 'ls_sparse')

#The below code uses grp_norm etc, not group_norm. See above for grp_norm
group_norm <- function(x) sum(by(x, grp, function(x) sqrt(sum(x^2))))
sp_group_norm <- function(x, alp=.05) group_norm(x)*(1-alp) + alp*sum(abs(x))

b1 = apply(out$beta, 2, group_norm, grp=grp)
b2 = apply(out_sp$beta, 2, sparse_group_norm, grp=grp)

par(mfcol=c(2,2))
matplot(b1/max(b1), t(out$beta[beta_star2 ==0,]), ty='l', lty=1, col = grp)
matplot(b2/max(b2), t(out_sp$beta[beta_star2==0,]), ty='l', lty=1, col = grp)
plot(b1/max(b1), apply(out$beta == 0 & beta_star2 == 0, 2, sum)/sum(beta_star2==0), ty='l')
lines(b2/max(b2), apply(out_sp$beta == 0 & beta_star2 == 0, 2, sum)/sum(beta_star2==0), col=2)
par(mfcol=c(1,1))
