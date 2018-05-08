## Example in Stan
## https://gist.github.com/mbjoseph/952b807bf5aad4a72a9d865f84d67afa
## lines 74-77:
## // priors for diagonal entries (Leung and Drton 2016)
## for (i in 1:k) {
##   target += (k - i) * log(beta_diag[i]) - .5 * beta_diag[i] ^ 2 / sigma_L;
## }

LD_diag <- function(x, sigma) {
  ## log-likelihood for Leung & Drton prior
  ## x is diag of loadings matrix
  ## sigma is scale parameter
  n <- length(x)
  ## log pdf
  lpdf <- (n - seq(n)) * log(x) - (0.5/sigma*x^2)
  ## log likelihood
  LL <- sum(lpdf)
  return(LL)
}
