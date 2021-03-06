% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate.R
\name{ciperm.twosample}
\alias{ciperm.twosample}
\title{Confidence interval for the two-sample difference in means}
\usage{
ciperm.twosample(x, y, ...)
}
\arguments{
\item{x}{Observations from sample 1. A numeric vector of length N_1 or a matrix of dim N_1 x K.}

\item{y}{Observations from sample 2. A numeric vector of length N_2 or a matrix of dim N_2 x K.}

\item{...}{Additional arguments to \link{ciperm} or \link{ciperm.multi}}
}
\value{
Vector or array of lower and upper ends of the confidence interval
}
\description{
Calculates the permutation-test based confidence interval for the difference in means.
Works in both for univariate and multivariate cases.
}
\details{
It is the confidence interval for mean(y) - mean(x)
}
\examples{
x <- rnorm(20)
y <- rnorm(15) + 0.3

ciperm.twosample(x,y, M = 2000, level = 0.99) ## 99\% confidence interval, 2000 permutations.

}
\seealso{
\link{ciperm},  \link{ciperm.multi}
}
