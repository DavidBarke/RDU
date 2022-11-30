#' Calculate all possible outcomes and their objective probabilities
#'
#' @param n Number of periods
#' @param mu Mean of logarithmic one-period stock return
#' @param sigma Standard deviation of logarithmic one-period
#'
#' @export
outcomes <- function(n, mu, sigma, p_up = 0.5) {
  u <- exp(mu + sigma)
  d <- exp(mu - sigma)
  n <- n + 1 # there is one outcome more than the number of periods

  u^(n-(1:n)) * d^((1:n)-1)
}

#' Calculate Decision Weights
#'
#' @param x Distinct outcomes sorted in decreasing order
#' @param cdf Non-decreasing, right-continuous, discrete cumulative distributive
#' function
#' @param pwf Probability weighting function
#'
#' @returns Decision weights
#'
#' @export
decision_weights <- function(x, cdf, pwf = identity) {
  stopifnot(!is.unsorted(rev(x)))

  n <- length(x)

  w <- numeric(n)

  for (i in seq_len(n)) {
    if (i < n) {
      w[i] <- pwf(1 - cdf(x[i+1])) - pwf(1 - cdf(x[i]))
    } else {
      w[i] <- 1 - pwf(1 - cdf(x[i]))
    }
  }

  w
}

#' Factory for probability weighting function from Prelec (1998)
#'
#' @param alpha Shape parameter (alpha > 1 -> S shaped, alpha < 1 -> inverse S
#' shaped)
#' @param beta Location of inflexion point parameter
#'
#' @returns Probability weighting function
#'
#' @export
pwf_prelec <- function(alpha = 1, beta = 1) {
  stopifnot(beta > 0 && alpha > 0)

  function(p) {
    exp(-beta*(-log(p))^alpha)
  }
}

#' Factory for Binomial CDF
#'
#' @param size Number of trials (zero or more).
#' @param prob Probability of success on each trial.
#'
#' @returns CDF of binomial distribution
#'
#' @export
cdf_binom <- function(size, prob) {
  function(q) {
    pbinom(q, size, prob)
  }
}

#' Factory for power utility function
#'
#' @param gamma Coefficient of relative risk aversion
#'
#' @returns Power utility function
#'
#' @export
power_utility <- function(gamma = 1) {
  stopifnot(gamma > 0)

  if (gamma == 1) return(log)

  function(x) {
    (x^(1-gamma)-1)/(1-gamma)
  }
}

#' Utility of risky terminal wealth
#'
#' @inheritParams decision_weights
#'
#' @return The utility of risky terminal wealth
terminal_utility <- function(x, cdf, pwf = identity, v = power_utility(1)) {
  w <- decision_weights(
    x = x,
    cdf = cdf,
    pwf = pwf
  )

  sum(w * v(x))
}

