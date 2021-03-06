% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.sparsegl-methods.R
\name{predict.cv.sparsegl}
\alias{predict.cv.sparsegl}
\title{Make predictions from a \code{cv.sparsegl} object.}
\usage{
\method{predict}{cv.sparsegl}(object, newx, s = c("lambda.1se", "lambda.min"), ...)
}
\arguments{
\item{object}{Fitted \code{\link[=cv.sparsegl]{cv.sparsegl()}} object.}

\item{newx}{Matrix of new values for \code{x} at which predictions are to be
made. Must be a matrix. See documentation for \code{\link[=predict.sparsegl]{predict.sparsegl()}}.}

\item{s}{Value(s) of the penalty parameter \code{lambda} at which
predictions are required. Default is the value \code{s = "lambda.1se"} stored
on the CV object. Alternatively \code{s = "lambda.min"} can be used. If
\code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
used.}

\item{\dots}{Not used. Other arguments to predict.}
}
\value{
The returned object depends on the \dots{} argument which is passed
on to the \code{\link[=predict]{predict()}} method for \code{\link[=sparsegl]{sparsegl()}} objects.
}
\description{
This function makes predictions from a cross-validated \code{\link[=sparsegl]{sparsegl()}} model,
using the stored \code{"sparsegl.fit"} object, and the value chosen
for \code{lambda}.
}
\details{
This function makes it easier to use the results of cross-validation to make
a prediction.
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
\code{\link[=cv.sparsegl]{cv.sparsegl()}}, and \code{\link[=coef.cv.sparsegl]{coef.cv.sparsegl()}}
methods.
}
\keyword{models}
\keyword{regression}
