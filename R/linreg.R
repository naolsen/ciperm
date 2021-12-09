
## Linear model
make.linear.cov.function <- function(x) {
  function(mu) x*mu
}

## Linear regression test statistic
make.linreg.stat <- function(x) {
  function(eps) abs(sum((x - mean(x))*(eps - mean(eps))))
}

#' Confidence interval for slope parameter in linear regression
#'
#' Calculates the permutation-test based confidence interval for \eqn{\beta_1} in the linear regression model
#' \deqn{y = \beta_0 + \beta_1 x + \epsilon}
#'
#' @param x Values of the independent variable
#' @param y Values of the dependent variable
#' @param bounds Bounds forwarded to \code{uniroot} in the search for l and u values.
#' @param ... Additional arguments to \link{ciperm}
#'
#' @return Vector of lower and upper ends of the confidence interval
#'
#' @seealso \link{ciperm}
#'
#' @examples
#' x <- runif(20, -1, 1)
#' y <- 3 + 0.8*x + rt(20, df = 4)
#'
#' ciperm.linreg(x, y, bounds = c(-5, 6))
#'
ciperm.linreg <- function(x, y, bounds, ...) {
  ff <- make.linreg.stat(x)
  phif <- make.linear.cov.function(x)

  ciperm(y, phif, statistic = ff, bounds = bounds, ...)
}