binomial <- function(data) {
  n <- length(data)  # Total number of trials (samples)
  x <- sum(data)     # Total number of successes
  p_hat <- x / n
  return(p_hat)
}
