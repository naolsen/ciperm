
## Linear model
make.linear.cov.function <- function(x) {
  function(mu) x*mu
}

## Linear regression test statistic
make.linreg.stat <- function(x) {
  function(eps) abs(sum((x - mean(x))*(eps - mean(eps))))
}


#' Confidence interval for the slope parameter in linear regression
#'
#' Calculates the permutation-test based confidence interval for \eqn{\beta_1} in the linear regression model
#' \deqn{y = \beta_0 + \beta_1 x + \epsilon}. Works in both for univariate and multivariate cases (but \code{x} must be univariate).
#'
#' @param x Values of the independent variable. A numeric vector of length N.
#' @param y Values of the dependent variable. A numeric vector of length N or a matrix of dim N x K.
#' @param bounds Bounds forwarded to \code{uniroot} in the search for l and u values.
#' @param ... Additional arguments to \link{ciperm} or \link{ciperm.multi}
#'
#' @return Vector or array of lower and upper ends of the confidence interval
#'
#' @seealso \link{ciperm},  \link{ciperm.multi}
#'
#' @examples
#' x <- runif(20, -1, 1)
#' y <- 3 + 0.8*x + rt(20, df = 4)
#'
#' ciperm.linreg(x, y, bounds = c(-5, 6))
#'
#' ## Multivariate example (requires package \code{mvtnorm})
#' library(mvtnorm)
#' x <- runif(20, -1, 1)
#' sigma <- matrix(0.80, 4, 4)
#' diag(sigma) <- 1
#'
#' y <- 3 + 0.8*x + rmvnorm(20, mean = rep(0,4), sigma = sigma)
#'
#' ciperm.linreg(x, y, bounds = c(-7, 8))
#'
ciperm.linreg <- function(x, y, bounds, ...) {
  ff <- make.linreg.stat(x)
  phif <- make.linear.cov.function(x)

  if (ncol(as.matrix(y)) == 1)
  ciperm(y, phif, statistic = ff, bounds = bounds, ...)

  else ciperm.multi(y, phif, statistic = ff, bounds = bounds, ...)
}