# setup for simple binomial case


nobs <- 100L
nvars <- 20L
group <- rep(1:4, each = 5)
y <- as.double(rbinom(nobs, 1, rbeta(100, 2, 2)))
x <- matrix(rnorm(nobs * nvars), nobs, nvars)
family <- binomial()

# objects created by sparsegl()
nlambda = 100L
lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04)
lambda = NULL
pf_sparse = rep(1, nvars)
intercept = as.integer(TRUE)
asparse = as.double(0.05)
standardize = TRUE
dfmax = as.integer(max(group)) + 1L
pmax = as.integer(min(dfmax * 1.2, as.integer(max(group))))
eps = as.double(1e-08)
maxit = as.integer(3e+08)
bn <- as.integer(max(group))
bs <- as.integer(as.numeric(table(group)))
iy <- cumsum(bs) # last column of x in each group
ix <- c(0, iy[-bn]) + 1 # first column of x in each group
ix <- as.integer(ix)
iy <- as.integer(iy)
group <- as.integer(group)
pf = as.double(sqrt(bs))
pfl1 <- as.double(pf / sum(pf) * nvars)
flmin <- as.double(lambda.factor)
ulam <- double(1)
intr <- as.integer(intercept)
lower_bnd <- as.double(rep(-9.9e30, bn))
upper_bnd <- as.double(rep(9.9e30, bn))
# end sparsegl() objects

# enter sgl_irwls()
weights <- rep(1, nobs)
etastart <- 0
mustart <- NULL
start <- NULL
eval(family$initialize)
y <- drop(y)
has_offset <- !is.null(NULL) # no offset for now, it's NULL
if (!has_offset) offset <- as.double(y * 0)
sx <- sqrt(Matrix::colSums(x^2))
sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
xs <- 1 / sx
x <- x %*% Matrix::Diagonal(x = xs)

init <- initializer(x, y, weights, family, intr = FALSE,
                   has_offset, offset, pfl1, ulam)
cur_lambda <- init$cur_lambda
findlambda <- init$findlambda
no_user_lambda <- init$findlambda
nulldev <- init$nulldev
trace_it <- 0L

vnames <- colnames(x)

static <- list(
  nulldev = as.double(nulldev),
  y = as.double(y),
  weights = as.double(weights),
  offset = as.double(offset),
  bn = as.integer(bn),
  bs = as.integer(bs),
  x = x,
  ix = as.integer(ix),
  iy = as.integer(iy),
  xs = as.double(xs),
  nobs = as.integer(nobs),
  nvars = as.integer(nvars),
  pf = as.double(pf),
  pfl1 = as.double(pfl1),
  dfmax = as.integer(dfmax),
  pmax = as.integer(pmax),
  flmin = as.double(flmin),
  eps = as.double(eps),
  maxit = as.integer(maxit),
  vnames = vnames,
  group = group,
  intr = as.integer(intr),
  asparse = asparse,
  lb = as.double(lower_bnd),
  ub = as.double(upper_bnd),
  family = family,
  trace_it = trace_it
)

warm <- make_irls_warmup(nobs, nvars, b0 = init$b0, r = init$r)
warm <- c(warm,
          activeGroup = list(integer(pmax)),
          activeGroupIndex = list(integer(bn)),
          sset = list(integer(bn)),
          ni = 0L, npass = 0L, me = 0L,
          findlambda = findlambda,
          eset = list(integer(bn)))
l <- 0L
warm$al0 <- as.double(cur_lambda)
warm$ulam <- as.double(init$lambda_max)
warm$l <- l

wx <- static$x
gamma <- calc_gamma(wx, static$ix, static$iy, static$bn)
# end sgl_irwls()


test_that("calling sgl_irwlsfit works", {
  skip("skipping all macro irwls tests.")

  s <- spgl_wlsfit(warm, wx, gamma, static) # stuck



})
