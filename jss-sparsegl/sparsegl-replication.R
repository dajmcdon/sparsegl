#' **Note on replication**
#' _14 September 2023_
#'
#' This file will reproduce all results, but downloading the large brain data
#' and running the timing simulation at the beginning takes some time/memory.
#'
#'
#' Necessary packages
#'
#' install.packages(c("ggplot2", "tidyr", "dplyr", "tibble", "SGL", "sparsegl"))
#'
#' install.packages(c("CVXR", "purrr", "gglasso"))

## The following are used in the .Rmd for the pdf submission
library("ggplot2")
library("tidyr")
library("dplyr")
theme_set(theme_bw(base_family = "Palatino"))



#' **Timing comparison**
#'
#'
#'
library("Matrix")
library("tibble")
library("SGL")
library("sparsegl")


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
    # print(paste("Done with p = ", p[i], "repl = ", j))
  }
}

res$p <- rep(p, each = nrepls * 3)
if (!dir.exists("large-data")) dir.create("large-data")
# saveRDS(res, "large-data/sparsegl-timing.rds")

# res <- readRDS("large-data/sparsegl-timing.rds")
res <- res |>
  group_by(method, p) |>
  summarise(time = median(time), .groups = "drop")
p <- unique(res$p)
better_labs <- c(
  SGL = "SGL", sparsegl = "sparsegl",
  sparsegl_sp = "sparsegl (sparse X)"
)
ggplot(res, aes(p, time, colour = method)) +
  geom_point() +
  geom_line() +
  scale_y_log10(name = "Median time (seconds)",
                breaks = c(0.01, .1, 1, 10, 100),
                labels = c(0.01, .1, 1, 10, 100)) +
  scale_x_log10(name = "Number of predictors (p)",
                breaks = p) +
  scale_colour_viridis_d(
    begin = .1, end = .9, name = "Method",
    labels = better_labs)


#' **Simple usage example**
#'
#'
## ----install, eval=FALSE, echo=TRUE--------------------------------------------------
## install.packages("sparsegl")
## library("sparsegl")


## ---- data-simulation----------------------------------------------------------------
set.seed(1010)
n <- 100
p <- 200
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- c(
  rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), c(2, -3, 8, 0, 0),
  rep(0, (p - 20))
)
groups <- rep(1:(p / 5), each = 5)
eps <- rnorm(n, mean = 0, sd = 1)
y <- X %*% beta + eps
pr <- 1 / (1 + exp(-X %*% beta))
y0 <- rbinom(n, 1, pr)



## ------------------------------------------------------------------------------------
fit <- sparsegl(X, y, group = groups)


## ----coef-trace
plot(fit, y_axis = "coef", x_axis = "penalty", add_legend = FALSE)
plot(fit, y_axis = "group", x_axis = "lambda", add_legend = FALSE)


## ----
coef(fit, s = c(0.02, 0.03))[c(1, 3, 25, 29), ]
predict(fit, newx = tail(X), s = fit$lambda[2:3])
print(fit)


## ----
cv_fit <- cv.sparsegl(X, y, groups, nfolds = 15)
plot(cv_fit)


## -----
coef(cv_fit, s = "lambda.1se")[c(1, 3, 25, 29), ]
predict(cv_fit, newx = tail(X), s = "lambda.min") |> c()


## -----
fit_logit <- sparsegl(X, y0, groups, family = "binomial")
cv_fit_logit <- cv.sparsegl(
  X, y0, groups, family = "binomial", pred.loss = "misclass"
)
plot(cv_fit_logit, log_axis = "none")


## ----
er <- estimate_risk(fit, X)


## ----re-plot
er <- er |> dplyr::select(-df) |> pivot_longer(-lambda, values_to = "risk")
err <- er |> group_by(name) |> summarise(lambda = lambda[which.min(risk)])
ggplot(er, aes(lambda, risk, color = name)) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2", name = "") +
  geom_vline(
    data = err, aes(xintercept = lambda, color = name),
    linetype = "dashed", show.legend = FALSE
  ) +
  ylab("Estimated risk") +
  xlab("Lambda") +
  scale_x_log10()

#' **Trust in experts**
#'
#'
#'
library("splines")
df <- 10

data("trust_experts", package = "sparsegl")
trust_experts <- trust_experts %>%
  mutate(across(
    where(is.factor),
    ~ `attr<-`(.x, "contrasts", contr.sum(nlevels(.x), FALSE, TRUE))
  ))

x <- Matrix::sparse.model.matrix(
  ~ 0 + region + age + gender + raceethnicity + period +
    bs(cli, df = df) + bs(hh_cmnty_cli, df = df),
  data = trust_experts, drop.unused.levels = TRUE)

gr <- sapply(trust_experts, function(x) ifelse(is.factor(x), nlevels(x), NA))
gr <- rep(seq(ncol(trust_experts) - 1), times = c(gr[!is.na(gr)], df, df))
fit <- cv.sparsegl(x, trust_experts$trust_experts, gr)

cc <- coef(fit, s = "lambda.1se")
reg <- which(substr(rownames(cc), 1, nchar("region")) == "region")
states <- tibble(state = rownames(cc)[reg], coef = cc[reg]) |>
  mutate(state = substring(state, nchar("region") + 1),
         state_name = tolower(covidcast::abbr_to_name(state, TRUE)))
states_map <- map_data("state")
ggplot(states, aes(map_id = state_name)) +
  geom_map(aes(fill = coef), map = states_map) +
  expand_limits(x = states_map$long, y = states_map$lat) +
  scale_fill_gradient2(
    low = "darkorange", high = "darkblue",
    trans = scales::pseudo_log_trans(),
    breaks = c(-10, -5, -2, 0, 2, 5, 10, 20),
  ) +
  coord_map(projection = "albers", parameters = c(30,40)) +
  theme_void() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_blank())



#' **Brain analysis**
#'
#'
#'
download <- function(url, path) {
  file.exists(path) || download.file(url, path, mode = 'wb')
  cat("Done.\n")
}

options(timeout = max(15 * 60)) # 15 minutes, assuming 1MB/s download speed
#' The below files exceed the limits for JSS submission.
#' They are available with a persistent DOI at the following link.
#' https://doi.org/10.6084/m9.figshare.20314917
download("https://figshare.com/ndownloader/files/36288819", "large-data/A.mtx.gz")
download("https://figshare.com/ndownloader/files/36288825", "large-data/Y.rds")
download("https://figshare.com/ndownloader/files/36288822", "large-data/G.rds")
download("https://figshare.com/ndownloader/files/36294165", "large-data/Gpf.rds")
download("https://figshare.com/ndownloader/files/36288828", "large-data/brain-fit.rds")


A <- Matrix::readMM("large-data/A.mtx.gz")
Y <- readRDS("large-data/Y.rds")
G <- readRDS("large-data/G.rds")
gpf <- readRDS("large-data/Gpf.rds")


#' ~ 1 minutes
system.time(fit <- sparsegl(A, Y, group = G, pf_group = gpf, asparse = 0.0))

df <- estimate_risk(fit, A, approx_df = TRUE)
# saveRDS(fit, file = "large-data/brain-fit.rds")
# fit <- readRDS("large-data/brain-fit.rds")

#' This will produce a warning
#'  "is.na() applied to non-(list or vector) of type 'language'"
#' It is due to `annotate()`, annoying, but meaningless
plot(fit, y_axis = "group", x_axis = "penalty", add_legend = FALSE) +
  scale_y_continuous(labels = scales::label_number(scale = 1e5)) +
  coord_cartesian(c(0, 1), c(0, 5e-5), clip = "off") +
  annotate(
    geom = "text", x = -.14, y = 5.2e-5, fontface = "plain",
    label = quote(phantom(0) %*% 10^{-5}), size = 2.5
  )

#' **Accuracy analysis**
#'
#'
#'
library("CVXR")
library("purrr")
library("gglasso")

SNR <- c(.1, 1, 10)
nobs <- 100
ngroups <- 10
nvars <- c(50, 100, 150)
nlambdas <- 20

make_problem <- function(nobs, nvars, SNR = 1, coefs, seed = 12345) {
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
  lambdas <- spgl$lambda
  ggl <- gglasso(X, y, groups, lambda = lambdas, pf = pfg, intercept = FALSE)

  # CVXR implementation, requires same-size groups
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
    CVXR = cvx_temp
  )
  muhat <- map(betahat, \(b) X %*% b)
  loss <- map(muhat, \(z) colMeans((y - z)^2) / 2) |> bind_cols()
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

design <- expand_grid(nvars = nvars, SNR = SNR, seed = 12345)

# runs for about 1 minute
results <- design |>
  rowwise() |>
  mutate(res = list(make_and_run(nvars, SNR, nobs, ngroups, seed = seed)))

# saveRDS(results, "large-data/accuracy-results.rds")

# results <- readRDS("large-data/accuracy-results.rds")
results <- results |>
  hoist("res", objective = list("obj"), lambda = list("lambda")) |>
  unnest(c(objective, lambda)) |>
  mutate(`p/n` = nvars / 100) |>
  select(-res, -seed, -nvars)

lb <- function(labels, sep = " = ") label_both(labels, TRUE, sep)
results <- results |>
  mutate(across(sparsegl:CVXR, ~ .x / sparsegl - 1)) |>
  mutate(sparsegl = NULL) |>
  group_by(`p/n`, SNR) |>
  mutate(
    lambda = log10(lambda),
    lambda = (lambda - min(lambda)) / (max(lambda) - min(lambda))
  )

pivot_longer(results, gglasso:CVXR) |>
  ggplot(aes(lambda, value, colour = name)) +
  facet_grid(`p/n` ~ SNR, labeller = lb, scales = "free_y") +
  geom_line() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_colour_manual(values = c("darkblue", "orange"), name = "") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    name = "Range of lambda, smallest (left) to largest (right)",
    breaks = c(0, .1, .22, .46,  1),
    labels = NULL, minor_breaks = NULL) +
  scale_y_continuous(
    labels = scales::label_percent(),
    name = "Change in objective relative to sparsegl"
  )

sessionInfo()
