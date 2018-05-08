## example from Stan
# // priors for diagonal entries (Leung and Drton 2016)
# for (i in 1:k) {
#   target += (k - i) * log(beta_diag[i]) - .5 * beta_diag[i] ^ 2 / sigma_L;
# }
