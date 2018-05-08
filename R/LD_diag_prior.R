## prior for diagonal of factor loadings matrix
## based on eqn 2.2 in Leung and Drton (2016) 
## https://arxiv.org/abs/1409.7672v1

LD_diag <- function(x, sigma) {
  ## log pdf for Leung & Drton prior
  ## x is diag of loadings matrix; x > 0
  ## sigma is scale parameter
  n <- length(x)
  ## log pdf
  lpdf <- (n - seq(n)) * log(x) - x^2 / (2 * sigma)
  return(lpdf)
}
