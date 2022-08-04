## ---- setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE, fig.path = "fig/")
options(prompt = 'R> ', continue = '+ ')
ggplot2::theme_set(ggplot2::theme_bw(base_family = "Palatino"))
# Add any packages we want to cite to the list below. Don't edit the bib
# manually, or changes will be overwritten. Add other references to
# `sparsegl.bib`
#
# Note: tags are R-<pkg>
knitr::write_bib(c("devtools", "knitr", "testthat", "usethis",
                   "rticles", "sparsegl", "glmnet",
                   "SGL","gglasso", "msgl", "biglasso",
                   "RSpectra",
                   "ggplot2", "dotCall64"),
                 file = "pkgs.bib")


## ----timing-comparison, echo=FALSE, fig.width=6, fig.height=3, fig.cap="This figure shows the time required to compute sparse group lasso solutions across a number of different problem sizes. In all cases, we use $n=500$ observations and 100 values of the tuning parameter $\\lambda$. The median is taken across 5 replications for each method and problem size. Note that both axes are on the log scale."----
# source("code/timing.R)
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
## library(sparsegl)


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
library(tidyverse)
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
source("code/covid.R")
g


## ----dwi-model-description, out.width="85%", fig.cap="Pictorial representation illustrating how the streamlines and voxels are converted from a diffusion-weighted image to a linear model. Each voxel is measured on 90 angles, so it occupies 90 rows in the data. When a streamline (column of $\\mathbf{X}$) passes through a voxel, the values within that voxel are given by a physical model based on the direction of passage. Otherwise, if the streamline does not cross the voxel, the respective rows are zero.", echo=FALSE----
knitr::include_graphics("fig/model-description.pdf")


## ----dwi-fit, echo=FALSE, warning=FALSE, eval=TRUE, fig.cap="The group norm of the 12 groups based on neuroanatomical structure against the magnitude of the penalty. The fit was produced using \\texttt{sparsegl()} to estimate the group lasso ($\\alpha=0$).", fig.height=3----
# source("code/brains.R") takes 6 GB memory and about 1 minute, once the data is downloaded
fit <- readRDS("large-data/brain-fit.rds")
plot(fit,
     y_axis = "group",
     x_axis = "penalty",
     add_legend = FALSE) +
  scale_y_continuous(labels = scales::label_number(scale = 1e5)) +
  coord_cartesian(c(0,1), c(0,5e-5), clip = "off") +
  annotate("text", x = -.14, y = 5.2e-5, fontface="plain",
           label = quote(phantom(0) %*% 10^{-5}), size = 2.5)

