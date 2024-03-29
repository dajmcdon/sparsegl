% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.sparsegl-methods.R
\name{predict.cv.sparsegl}
\alias{predict.cv.sparsegl}
\title{Make predictions from a \code{cv.sparsegl} object.}
\usage{
\method{predict}{cv.sparsegl}(
  object,
  newx,
  s = c("lambda.1se", "lambda.min"),
  type = c("link", "response", "coefficients", "nonzero", "class"),
  ...
)
}
\arguments{
\item{object}{Fitted \code{\link[=cv.sparsegl]{cv.sparsegl()}} object.}

\item{newx}{Matrix of new values for \code{x} at which predictions are to be
made. Must be a matrix. This argument is mandatory.}

\item{s}{Value(s) of the penalty parameter \code{lambda} at which
coefficients are desired. Default is the single
value \code{s = "lambda.1se"} stored in the CV object (corresponding to
the largest value of \code{lambda} such that CV error estimate is within 1
standard error of the minimum). Alternatively \code{s = "lambda.min"} can be
used (corresponding to the minimum of cross validation error estimate).
If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}

\item{type}{Type of prediction required. Type \code{"link"} gives the linear
predictors for \code{"binomial"}; for \code{"gaussian"} models it gives the fitted
values. Type \code{"response"} gives predictions on the scale of the response
(for example, fitted probabilities for \code{"binomial"}); for \code{"gaussian"} type
\code{"response"} is equivalent to type \code{"link"}. Type
\code{"coefficients"} computes the coefficients at the requested values for
\code{s}.
Type \code{"class"} applies only to \code{"binomial"} models, and produces the
class label corresponding to
the maximum probability. Type \code{"nonzero"} returns a list of the indices
of the nonzero coefficients for each value of \code{s}.}

\item{...}{Not used.}
}
\value{
A matrix or vector of predicted values.
}
\description{
This function makes predictions from a cross-validated \code{\link[=cv.sparsegl]{cv.sparsegl()}} object,
using the stored \code{sparsegl.fit} object, and the value chosen for \code{lambda}.
}
\examples{
n <- 100
p <- 20
X <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n)
beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
y <- X \%*\% beta_star + eps
groups <- rep(1:(p / 5), each = 5)
fit1 <- sparsegl(X, y, group = groups)
cv_fit <- cv.sparsegl(X, y, groups)
predict(cv_fit, newx = X[50:60, ], s = "lambda.min")

}
\seealso{
\code{\link[=cv.sparsegl]{cv.sparsegl()}} and \code{\link[=coef.cv.sparsegl]{coef.cv.sparsegl()}}.
}
