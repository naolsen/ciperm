


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
#' Works in both for univariate and multivariate cases.
#'
#' @param x Observations from sample 1. A numeric vector of length N_1 or a matrix of dim N_1 x K.
#' @param y Observations from sample 2. A numeric vector of length N_2 or a matrix of dim N_2 x K.
#' @param ... Additional arguments to \link{ciperm} or \link{ciperm.multi}
#'
#' @details It is the confidence interval for mean(y) - mean(x)
#'
#' @return Vector or array of lower and upper ends of the confidence interval
#'
#' @seealso \link{ciperm},  \link{ciperm.multi}
#'
#' @examples
#' x <- rnorm(20)
#' y <- rnorm(15) + 0.3
#'
#' ciperm.twosample(x,y, M = 2000, level = 0.99) ## 99% confidence interval, 2000 permutations.
#'
ciperm.twosample <- function(x, y, ...) {
  ## Set up covariate function and test statistic
  theta <- c(rep(0, NROW(x)), rep(1, NROW(y)))

  f <- factor(theta)
  ff <- make.twosample.stat(f)
  phif <- make.linear.cov.function(theta)

  if (NCOL(y) == 1)
    ciperm(c(x,y), phif, statistic = ff, bounds = c(min(y) - max(x), max(y) - min(x)), ...)

  else ciperm.multi(rbind(x,y), phif, statistic = ff, bounds = c(min(y) - max(x), max(y) - min(x)), ...)

}
