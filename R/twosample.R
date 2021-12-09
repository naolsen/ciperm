


## Covariate function for two-sample case
twosample.cov.function <- function(x,y) {
  make.linear.cov.function(c(rep(0, length(x)), rep(1, length(y))))
}


## Test statistic for two-sample case
make.twosample.stat <- function(f) {
  if (!is.factor(f) || length(levels(f)) != 2) stop ("f must be a factor with two levels")

  f <- as.numeric(f)
  g2 <- which(f == 2)
  g1 <- which(f == 1)
  function(x) abs(mean(x[g2]) - mean(x[g1]))
}


#' Confidence interval for the two-sample difference in means
#'
#' Calculates the permutation-test based confidence interval for the difference in means.
#'
#' @param x Observations from sample 1
#' @param y Observations from sample 2
#' @param ... Additional arguments to \link{ciperm}
#'
#' @details It is the confidence interval for mean(y) - mean(x)
#'
#' @return Vector of lower and upper ends of the confidence interval
#'
#' @seealso \link{ciperm}
#'
#' @examples
#' x <- rnorm(20)
#' y <- rnorm(15) + 0.3
#'
#' ciperm.twosample(x,y, M = 2000, level = 0.99) ## 99% confidence interval, 2000 permutations.
#'
ciperm.twosample <- function(x, y, ...) {
  ## Set up covariate function and test statistic
  theta <- c(rep(0, length(x)), rep(1, length(y)))
  f <- factor(theta)
  ff <- make.twosample.stat(f)
  phif <- twosample.cov.function(x,y)

  ciperm(c(x,y), phif, statistic = ff, bounds = c(min(y) - max(x), max(y) - min(x)), ...)
}
