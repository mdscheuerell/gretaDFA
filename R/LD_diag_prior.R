## Example in Stan
## https://gist.github.com/mbjoseph/952b807bf5aad4a72a9d865f84d67afa
## lines 74-77:
## // priors for diagonal entries (Leung and Drton 2016)
## for (i in 1:k) {
##   target += (k - i) * log(beta_diag[i]) - .5 * beta_diag[i] ^ 2 / sigma_L;
## }

nn <- 10
beta <- seq(10)/10
LL <- 0
sigma <- 1

for(i in 1:nn) {
  LL <- LL + (nn-i) * log(beta[i]) - (0.5*beta[i]^2)/sigma
}
