## prior for diagonal of factor loadings matrix
## based on eqn 2.2 in Leung and Drton (2016) 
## https://arxiv.org/abs/1409.7672v1

LD_diagv <- function(X, sigma) {
  ## log pdf for Leung & Drton prior
  ## X is diag of loadings matrix; X > 0
  ## sigma is scale parameter
  n <- length(X)
  lpdf <- (n - seq(n)) * log(X) - X^2 / (2 * sigma)
  return(lpdf)
}

LD_diag <- function(xi, i, m, sigma) {
  ## log pdf for Leung & Drton prior
  ## for i in 1:m, xi = X[i] > 0
  ## m is number of factors
  ## sigma is scale parameter
  (m - i) * log(xi) - xi^2 / (2 * sigma)
}