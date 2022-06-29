
#' Permutation scheme for non-parametric confidence intervals (multivariate version)
#'
#' The actual permutation scheme, and the multivariate version of \link{ciperm0}.
#'
#' @param y Matrix of dimension N x K with at least two columns.
#' @param estimate Parameter estimate (vector of length K). Optional
#'
#' @details See \link{ciperm0} for additional parameters.
#' \code{phi}, \code{statistic} and \code{bounds} are assumed univariate and similar across coordinates
#'
#' @return array of dim M x K x 2
#' @export
#'
ciperm0.multi <- function(y, phi, M = 1000, statistic, bounds, estimate,  infinities = FALSE) {

  if (!is.matrix(y) || ncol(y) == 1) stop("Data is of incorrect format")

  N <- nrow(y)
  if (N != length(phi(0))) stop("y and covariate lengths differ!")
  K <- ncol(y)

  lower <- min(bounds)
  upper <- max(bounds)

  ## Missing estimate, optimize and find out
  if (missing(estimate)) {
  estimate <- numeric(K)
  for (k in 1:K) {
    estimate[k] <- optimize(function(mu) statistic(y[,k] - phi(mu)), interval = bounds)$minimum
    if (estimate[k] < lower + .Machine$double.eps^0.25 || estimate[k] > upper - .Machine$double.eps^0.25 )
      stop("Estimate not in the interior of the interval. Did you specify your function correctly?")
    ## make some check that estimate is sensible..
    }
  }

  ## Permutation scheme
  lusm <- array(, c(M, K, 2))
  for (j in 1:M) {
    s <- sample(N)
    for (k in 1:K) {
      if (infinities)
        lusm[j,k,] <- c(tryCatch(uniroot(function(mu) statistic((y[,k] - phi(mu))[s]) - statistic(y[,k] - phi(mu)),
                                      lower = lower, upper = estimate[k])$root, error = function(e) -Inf),
                     tryCatch(uniroot(function(mu) statistic((y[,k] - phi(mu))[s]) - statistic(y[,k] - phi(mu)),
                                      lower = estimate[k], upper = upper)$root, error = function(e) Inf))

      else
        lusm[j,k,] <- c(uniroot(function(mu) statistic((y[,k] - phi(mu))[s]) - statistic(y[,k] - phi(mu)),
                              lower = lower, upper = estimate[k])$root,
                      uniroot(function(mu) statistic((y[,k] - phi(mu))[s]) - statistic(y[,k] - phi(mu)),
                              lower = estimate[k], upper = upper)$root)
    }
  }
  lusm
}


#' Multivariate confidence intervals using permutation tests
#'
#' Constructs confidence intervals at the specified level, and optionally calculates the joint confidence level.
#'
#' @param level Coverage level (i.e. 1- \eqn{\alpha}).
#' @param multilevel Calculate the joint coverage level using \link{alpha.multi}?
#' @param adjusted.ci Calculate the adjusted confidence interval using \link{adjusted_ci}?
#'
#' @details See \link{ciperm0.multi} for additional parameters.
#' The outputs from \code{multilevel} and \code{adjusted.ci} are attached as attributes (when the parameters are set to TRUE).
#'
#' @return A 2 x K array of K confidence intervals.
#' @export
#'
#' @seealso \link{ciperm}
#'
ciperm.multi <- function(y, phi, M = 1000, statistic, bounds, estimate, level = 0.95,
                         multilevel = TRUE, adjusted.ci = FALSE, infinities = FALSE) {
  lus <- ciperm0.multi(y, phi, M = M, statistic = statistic, bounds = bounds, estimate = estimate, infinities = infinities)
  dl <- dim(lus)
  K <- dl[2]

  out <- matrix(, 2, K)

  for (j in 1:K) {
    out[,j] <- c(quantile(lus[,j,1], probs = 1-level), quantile(lus[,j,2], probs = level))
  }
  rownames(out) <- c("lower", "upper")
  if (infinities && any(is.infinite(out)))
    warning("Confidence interval contains infinity. Did you make your bounds wide enough?")
  if (multilevel) attr(out, "alpha_multiple") <- alpha.multi(lus[,,1], lus[,,2], 1-level)

  if (adjusted.ci) {
    out2 <- adjusted_ci(lus, level = level, include_unadjusted = FALSE)
    attr(out, 'adjusted.ci') <- out2
  }
  out
}


#' Joint coverage level
#'
#' @param lows matrix of l values
#' @param ups matrix of u values
#' @param alpha Significance level (for unadjusted/single confidence intervals)
#'
#' @return Calculates alpha_multiple from the output of \link{ciperm0.multi}
#' @export
#'
alpha.multi <- function(lows, ups, alpha) {

  stopifnot(dim(lows) == dim(ups))
  K <- ncol(lows)
  M <- nrow(lows)

    c0 <- 0:(2^K-1)
    corners <- matrix(, 2^K, 0)
    for (i in 1:K) corners <- cbind(corners, c0 %% 2^i %/% 2^(i-1))

    dds <- numeric(2^K)
    for (j in 1:2^K) {
      l <- corners[j,]

      dd <- rep(F, M)
      for(u in 1:K) {
        if (l[u] == 0) dd <- dd | lows[,u] < quantile(lows[,u], alpha)
        else dd <- dd | ups[,u] > quantile(ups[,u], 1-alpha)
      }
      dds[j] <-  sum(dd)
    }
    dds
    max(dds)/M
}


#' Adjusted (multivariate) confidence intervals
#'
#' Adjust coverage level such that joint coverage level corresponds
#'
#' @param level
#' @param tol tolerance for reaching \code{level}.
#' @param include_unadjusted Include unadjusted confidence intervals as an attribute?
#'
#' @details This function uses output from \code{ciperm0.multi}.
#' In a simple bisection algorithm, the confidence level is adjusted until the joint coverage level (form \code{alpha.multi}) is \code{code}.
#'
#' @return Adjusted confidence intervals with the adjusted confidence level as an attribute.
#' @export
#'
#' @seealo \link{alpha.multi}
adjusted_ci <- function(lus, level = 0.95, tol = 0.001, include_unadjusted = TRUE) {

  K <- dim(lus)[2]

  out <- matrix(, 2, K)

  for (j in 1:K) {
    out[,j] <- c(quantile(lus[,j,1], probs = 1-level), quantile(lus[,j,2], probs = level))
  }
  rownames(out) <- c("lower", "upper")

  a_target <- 1-level
  a_max <- a_target
  a_min <- 0
  a <- alpha.multi(lus[,,1], lus[,,2], 1-level)
  a_new <- a_max

  while(abs(a -a_target) > tol) {

    if (a_max - a_min < 1e-7) stop("Convergence error!")

    a_new <- (a_min + a_max) / 2
    a <- alpha.multi(lus[,,1], lus[,,2], a_new)

    if (a < a_target) a_min <- a_new
    else a_max <- a_new
  }

  out2 <- matrix(, 2, K)

  for (j in 1:K) {
    out2[,j] <- c(quantile(lus[,j,1], probs = a_new), quantile(lus[,j,2], probs = 1-a_new))
  }
  rownames(out2) <- c("lower", "upper")

  attr(out2, 'alpha_adjusted') <- a_new

  if (include_unadjusted) attr(out2, 'unadjusted.ci') <- out
  out2
}
