library(gglasso)
library(sparsegl)
library(CVXR)
library(purrr)
library(dplyr)

make_problem <- function(nobs, nvars, SNR = 1,
                         coefs = c(rep(1, 5), rep(0, 5)),
                         seed = 12345) {
  set.seed(seed)
  X <- matrix(rnorm(nobs * nvars), nobs, nvars)
  beta <- rep(coefs, length.out = nvars)
  mu <- drop(X %*% beta)
  noise_sd <- sqrt(sum(beta^2)) / SNR # in expectation
  epsilon <- rnorm(nobs, 0, noise_sd)
  y <- mu + epsilon
  list(y = y, mu = mu, epsilon = epsilon, beta = beta, X = X, seed = seed)
}

create_algos <- function(
    X, y, groups, asparse = 0.05, nlambdas = 20, lambdas = NULL
) {

  bn <- max(groups)
  ngroups <- length(groups) / bn
  pfg <- rep(1, bn)
  spgl <- sparsegl(
    X, y, groups, lambda = lambdas,
    nlambda = nlambdas, asparse = asparse,
    standardize = FALSE, pf_group = pfg, intercept = FALSE)
  # eps = .Machine$double.eps)
  lambdas <- spgl$lambda
  ggl <- gglasso(X, y, groups, lambda = lambdas, pf = pfg, intercept = FALSE)
                 # eps = .Machine$double.eps)

  # CVXR implementation, need same size groups
  beta <- Variable(ncol(X))
  loss <- sum((y - X %*% beta)^2) / (2 * nrow(X))
  pen <- function(beta, asparse) {
    l1 <- p_norm(beta, 1)
    l2 <- sum_entries(p_norm(reshape_expr(beta, c(bn, ngroups)), 2, axis = 2))
    asparse * l1 + (1 - asparse) * l2
  }
  cvx_temp <- sapply(lambdas, function(lambda) {
    obj <- loss + lambda * pen(beta, asparse)
    prob <- Problem(Minimize(obj))
    res <- solve(prob)
    res$getValue(beta)
  })

  betahat <- list(
    sparsegl = as.matrix(spgl$beta),
    gglasso = ggl$beta,
    cvxr = cvx_temp
  )
  muhat <- map(betahat, \(b) X %*% b)
  loss <- map(muhat, \(z) colMeans((y - z)^2) / 2) |>
    bind_cols()
  penalty <- map(
    betahat,
    \(b) lambdas * apply(b, 2, sp_group_norm, gr = groups, asparse = asparse)
  ) |>
    bind_cols()
  obj <- loss + penalty
  list(
    obj = obj, loss = loss, penalty = penalty, lambda = lambdas,
    betahat = betahat
  )
}



# design ------------------------------------------------------------------

SNR <- c(.1, 1, 10)
nobs <- 100
ngroups <- 10
nvars <- c(50, 100, 150)
nlambdas <- 20

make_and_run <- function(nvars, SNR, nobs = 100, ngroups = 10, seed = 12345) {
  stopifnot(nvars %% ngroups == 0)
  stopifnot(ngroups %% 2 == 0)
  cat("nvars = ", nvars, ", SNR = ", SNR, "\n")
  grp_size <- nvars / ngroups
  gr <- rep(1:ngroups, each = grp_size)
  prob <- make_problem(
    nobs, nvars, SNR,
    coefs = rep(c(1, 0), each = grp_size)
  )
  out <- with(prob, create_algos(X, y, gr, nlambdas = nlambdas))
  out
}

design <- tidyr::expand_grid(nvars = nvars, SNR = SNR, seed = 12345)

# runs for about 5 minutes
results <- design |>
  rowwise() |>
  mutate(res = list(make_and_run(nvars, SNR, nobs, ngroups, seed = seed)))


saveRDS(results, here::here("jss-sparsegl", "large-data", "accuracy-results.rds"))


# plotting ----------------------------------------------------------------

results <- readRDS(here::here("jss-sparsegl", "large-data", "accuracy-results.rds"))
results <- results |>
  tidyr::hoist("res", objective = list("obj"), lambda = list("lambda")) |>
  tidyr::unnest(c(objective, lambda)) |>
  mutate(`p/n` = nvars / 100) |>
  select(-res, -seed, -nvars)

library(ggplot2)
library(scales)

results |>
  mutate(across(sparsegl:cvxr, ~ .x / sparsegl - 1)) |>
  mutate(sparsegl = NULL) |>
  pivot_longer(gglasso:cvxr) |>
  ggplot(aes(lambda, value, colour = name)) +
  facet_grid(`p/n` ~ SNR, labeller = label_both, scales = "free") +
  geom_line() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_colour_manual(values = c("darkblue", "orange"), name = "") +
  theme(legend.position = "bottom") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    name = "Change in objective relative to sparsegl"
  )
