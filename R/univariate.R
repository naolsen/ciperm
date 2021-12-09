
#' Permutation scheme for non-parametric confidence intervals
#'
#' The actual permutation scheme. Note that \code{estimate} (if provided), \code{statistic} and \code{phi} must match.
#'
#' @param y Data (vector of length N)
#' @param phi Covariate function. Must return a vector of length N
#' @param M Number of permutations
#' @param statistic Test statistic
#' @param bounds Bounds forwarded to \code{uniroot} in the search for l and u values.
#' @param estimate Parameter estimate. Optional
#' @param infinities If TRUE set l_m to infinity, when uniroot fails to find a root (Not implemented yet).
#' If FALSE, the procedure fails with an error.
#'
#' @return Matrix of dim M x 2
#'
#' @seealso \link{ciperm}
#'
ciperm0 <- function(y, phi, M = 1000, statistic, bounds, estimate, infinities = FALSE) {
  N <- length(y)
  if (length(y) != length(phi(0))) stop("y and covariate lengths differ!")

  lower <- min(bounds)
  upper <- max(bounds)

  ## Missing estimate, optimize and find out
  if (missing(estimate)) {
    estimate <- optimize(function(mu) statistic(y - phi(mu)), interval = bounds)$minimum
    if (estimate < lower + .Machine$double.eps^0.25 || estimate > upper - .Machine$double.eps^0.25 )
      stop("Estimate not in the interior of the interval. Did you specify your function correctly?")
    ## make some check that estimate is sensible..
  }


  ## Permutation scheme
  lus <- matrix(, M, 2)
  for (j in 1:M) {
    s <- sample(N)
    lus[j,] <- c(uniroot(function(mu) statistic((y - phi(mu))[s]) - statistic(y - phi(mu)),
                         lower = lower, upper = estimate)$root,
                 uniroot(function(mu) statistic((y - phi(mu))[s]) - statistic(y - phi(mu)),
                         lower = estimate, upper = upper)$root)
  }
  lus
}

#' Confidence interval using permutation test
#'
#' Constructs confidence intervals at the specified level.
#'
#' @param level Confidence level (i.e. 1- \alpha).
#'
#' @details See \link{ciperm0} for additional parameters.
#'
#' @return Vector of lower and upper ends of the confidence interval
#'
#' @seealso \link{ciperm.twosample} for a user-friendly version for the two-sample case.
#' \link{ciperm.linreg} for a user-friendly version for the linear regression case.
#'
ciperm <-  function(y, phi, M = 1000, statistic, bounds, estimate, level = 0.95) {
  lus <- ciperm0(y, phi, M = M, statistic = statistic, bounds = bounds, estimate = estimate)

  out <- c(quantile(lus[,1], probs = 1-level), quantile(lus[,2], probs = level))
  names(out) <- c("lower", "upper")
  out
}



