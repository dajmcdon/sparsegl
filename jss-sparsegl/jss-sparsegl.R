# _Note on replication_
#
# This file will reproduce all results, but downloading the large brain data
# and running the timing simulation at the beginning takes some time/memory.
# For this reason, we have provided the results of these experiments in
# `large-data.zip`. The reproduction code is commented out, while the
# `large-data.zip` file is expected to be extracted in the same directory as
# this script.


# Necessary packages
# install.packages(c("knitr", "ggplot2", "tibble", "SGL", "dplyr", "covidcast"))
# install.packages(c("scales", "tidyr))
#

## ---- setup, include=FALSE----------------------------------------------------
## Used in the Rmd for the pdf submission
# knitr::opts_chunk$set(echo = TRUE, cache = TRUE, fig.path = "fig/")
# options(prompt = 'R> ', continue = '+ ')
library(ggplot2)
library(tidyr)
ggplot2::theme_set(ggplot2::theme_bw(base_family = "Palatino"))



## ----timing-comparison, echo=FALSE, fig.width=6, fig.height=3, fig.cap="This figure shows the time required to compute sparse group lasso solutions across a number of different problem sizes. In all cases, we use $n=500$ observations and 100 values of the tuning parameter $\\lambda$. The median is taken across 5 replications for each method and problem size. Note that both axes are on the log scale."----
# library(Matrix)
# library(tibble)
# library(SGL)
# library(sparsegl)
#
# nrepls <- 5
# n <- 500
# p <- c(100, 250, 500, 1000, 2500)
# b <- c(rep(1, 5), rep(-1, 5), c(-1,-1,0,0,0), c(1,1,0,0,0), rep(1, 5),
#        rep(0, 75)) * 10
# s <- .1
# gr <- rep(1:(max(p) / 5), each = 5)
# x <- matrix(rnorm(n*max(p)), nrow = n)
# x[as.logical(rbinom(n*max(p), 1, 1 - s))] <- 0
# xsp <- Matrix(x)
# mu <- sapply(p, function(z) x[ ,1:z] %*% rep(b, length.out = z) / sqrt(z))
# # prob <- 1 / (1 + exp(-mu))
#
# signal <- sqrt(colSums(mu^2))
# noise_sd <- sqrt(signal)
#
# y <- mu + rnorm(n*length(p), sd = rep(noise_sd, each = n))
# res <- tibble(method = "a", time = proc.time()["elapsed"], .rows = 0)
# for (i in seq_along(p)) {
#   pp <- seq(p[i])
#   dat <- list(y = y[ ,i], x = x[ ,pp])
#   xxsp <- xsp[ ,pp]
#   g <- gr[pp]
#   for (j in seq(nrepls)) {
#     s1 <- system.time(
#       SGL(dat, nlam = 100, alpha = 0.05, index = g, standardize = FALSE))
#     res <- add_row(res, method = "SGL", time = s1["elapsed"])
#     s2 <- system.time(sparsegl(dat$x, dat$y, g, standardize = FALSE))
#     res <- add_row(res, method = "sparsegl", time = s2["elapsed"])
#     s3 <- system.time(sparsegl(xxsp, dat$y, g, standardize = FALSE))
#     res <- add_row(res, method = "sparsegl_sp", time = s3["elapsed"])
#     print(paste("Done with p = ", p[i], "repl = ", j))
#   }
# }
#
# res$p <- rep(p, each = nrepls * 3)
# saveRDS(res, "large-data/sparsegl-timing.rds")
library(dplyr)
res <- readRDS("large-data/sparsegl-timing.rds")
res <- res %>%
  group_by(method, p) %>%
  summarise(time = median(time), .groups = "drop")
p <- unique(res$p)
better_labs <- c(
  SGL = "SGL", sparsegl = "sparsegl",
  sparsegl_sp = "sparsegl (sparse X)"
)
ggplot(res, aes(p, time, colour = method)) +
  geom_point() +
  geom_line() +
  scale_y_log10(name = "Median time (seconds)") +
  scale_x_log10(name = "Number of predictors (p)",
                breaks = p) +
  scale_colour_viridis_d(
    begin = .1, end = .9, name = "Method",
    labels = better_labs) +
  theme_bw()


## ----install, eval=FALSE, echo=TRUE------------------------------------------------
## install.packages("sparsegl")


## ---- data-simulation--------------------------------------------------------------
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


## ---- eval=TRUE, echo=FALSE--------------------------------------------------------
library(sparsegl)


## ----------------------------------------------------------------------------------
fit <- sparsegl(X, y, group = groups)


## ----coef-trace, fig.width = 4, fig.height = 3, fig.cap="Coefficient trace plots produced with the \\texttt{plot()} method.", fig.show='hold', out.width="45%"----
plot(fit, y_axis = "coef", x_axis = "penalty", add_legend = FALSE)
plot(fit, y_axis = "group", x_axis = "lambda", add_legend = FALSE)


## ---- message=FALSE----------------------------------------------------------------
coef(fit, s = c(0.02, 0.03))[c(1,3,25,29),] # display a few
predict(fit, newx = tail(X), s = fit$lambda[2:3])
print(fit)


## ----cv-plot, fig.width = 4, fig.height = 2, fig.cap = "Cross validation error produced by the \\texttt{plot()} method for a \\texttt{cv.sparsegl} object."----
cv_fit <- cv.sparsegl(X, y, groups, nfolds = 15)
plot(cv_fit)


## ----------------------------------------------------------------------------------
coef(cv_fit, s = "lambda.1se")[c(1,3,25,29),]
predict(cv_fit, newx = tail(X), s = "lambda.min")


## ----logitres, fig.width = 4, fig.height = 2, fig.cap="Cross validation error for logistic regression produced by the \\texttt{plot()} method using misclassification error on the held-out set."----
fit_logit <- sparsegl(X, y0, groups, family = "binomial")
cv_fit_logit <- cv.sparsegl(
  X, y0, groups, family = "binomial", pred.loss = "misclass"
)
plot(cv_fit_logit, log_axis = "none")


## ----risk-estimate, fig.width = 4, fig.height = 2----------------------------------
er <- estimate_risk(fit, X)


## ----re-plot, echo=FALSE, message=FALSE, warning=FALSE, fig.width = 4, fig.height = 2, fig.cap="AIC, BIC, and GCV (solid lines) along with their minima (vertical dashed lines)."----
er <- er %>% dplyr::select(-df) %>% pivot_longer(-lambda, values_to = "risk")
err <- er %>% group_by(name) %>% summarise(lambda = lambda[which.min(risk)])
ggplot(er, aes(lambda, risk, color = name)) +
  geom_line() +
  scale_color_brewer(palette = "Dark2") +
  geom_vline(data = err, aes(xintercept = lambda, color=name),
             linetype="dashed", show.legend = FALSE) +
  theme_bw() +
  ylab("Estimated risk") +
  xlab("Lambda") +
  scale_x_log10() +
  theme(legend.title = element_blank())


## ----trust, echo=FALSE, fig.cap="State-level estimates for the amount of trust in experts about Covid-19. The value displayed represents the change relative the US-wide average.", fig.width=8, fig.height=4, message=FALSE, warning=FALSE, out.width="5in"----
data("trust_experts")

y <- trust_experts[,"y"]
x <- trust_experts[,-1]
gr <- c(rep(1:7, times = c(8, 51, 5, 4, 8, 10, 10)))
fit <- sparsegl(x, y, gr)
er <- estimate_risk(fit, x, approx_df = FALSE)
cc <- coef(fit, s = er$lambda[which.min(er$BIC)])
states <- tibble(state = rownames(cc)[10:60], coef = cc[10:60]) %>%
  mutate(state_name = tolower(covidcast::abbr_to_name(state, TRUE)))
states_map <- map_data("state")
g <- ggplot(states, aes(map_id = state_name)) +
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
g


## ----dwi-model-description, out.width="85%", fig.cap="Pictorial representation illustrating how the streamlines and voxels are converted from a diffusion-weighted image to a linear model. Each voxel is measured on 90 angles, so it occupies 90 rows in the data. When a streamline (column of $\\mathbf{X}$) passes through a voxel, the values within that voxel are given by a physical model based on the direction of passage. Otherwise, if the streamline does not cross the voxel, the respective rows are zero.", echo=FALSE----
# knitr::include_graphics("fig/model-description.pdf")


## ----dwi-fit, echo=FALSE, warning=FALSE, eval=TRUE, fig.cap="The group norm of the 12 groups based on neuroanatomical structure against the magnitude of the penalty. The fit was produced using \\texttt{sparsegl()} to estimate the group lasso ($\\alpha=0$).", fig.height=3----
# Download and load brain data ---------------------------------------------------------
# download <- function(url, path) {
#   file.exists(path) || download.file(url, path, mode = 'wb')
# }
#
# options(timeout = max(15 * 60)) # 15 minutes, assuming 1MB/s download speed
# # https://doi.org/10.6084/m9.figshare.20314917
# download("https://figshare.com/ndownloader/files/36288819", "large-data/A.mtx.gz")
# download("https://figshare.com/ndownloader/files/36288825", "large-data/Y.rds")
# download("https://figshare.com/ndownloader/files/36288822", "large-data/G.rds")
# download("https://figshare.com/ndownloader/files/36294165", "large-data/Gpf.rds")
# download("https://figshare.com/ndownloader/files/36288828", "large-data/brain-fit.rds")
#
# # The below files exceed the limits for JSS submission. We save only the result.
# A <- Matrix::readMM("large-data/A.mtx.gz")
# Y <- readRDS("large-data/Y.rds")
# G <- readRDS("large-data/G.rds")
# gpf <- readRDS("large-data/Gpf.rds")
#
# # Fit and estimate risk ---------------------------------------------------
# system.time(fit <- sparsegl(A, Y, group=G, pf_group=gpf, asparse=0.0)) # ~ 1.5 minutes
#
# df <- estimate_risk(fit, A, approx_df=TRUE)
# write.csv(df, "large-data/fit.csv", row.names=FALSE, quote=FALSE) # unused
# saveRDS(fit, file = "large-data/brain-fit.rds")
fit <- readRDS("large-data/brain-fit.rds")
plot(fit,
     y_axis = "group",
     x_axis = "penalty",
     add_legend = FALSE) +
  scale_y_continuous(labels = scales::label_number(scale = 1e5)) +
  coord_cartesian(c(0,1), c(0,5e-5), clip = "off") +
  annotate("text", x = -.14, y = 5.2e-5, fontface="plain",
           label = quote(phantom(0) %*% 10^{-5}), size = 2.5)

sessionInfo()
