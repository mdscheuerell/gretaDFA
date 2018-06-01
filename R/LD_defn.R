# a temporary hack to test out the distribution
library (greta)
library (R6)
library (tensorflow)
distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node
check_dims <- .internals$utils$checks$check_dims
as.greta_array <- .internals$greta_arrays$as.greta_array
is_scalar <- .internals$utils$misc$is_scalar
fl <- .internals$utils$misc$fl

# function to call 
ld <- function (i, M, sigma, dim = 1) {
  distrib("ld", i, M, sigma, dim)
  }

# Leung & Drton distn for prior of diag(Z)
ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (i, M, sigma, dim) {
      
      # check if (i, M, dim) are in counting set
      # (i)
      if (length(i) > 1 ||
          i <= 0 ||
          !is.finite(i) ||
          i != floor(i)) {
        
        stop ("i must be a scalar positive integer, but was: ",
              capture.output(dput(i)),
              call. = FALSE)
        
      }
      
      # (M)
      if (length(M) > 1 ||
          M <= 0 ||
          !is.finite(M) ||
          M != floor(M)) {
        
        stop ("M must be a scalar positive integer, but was: ",
              capture.output(dput(M)),
              call. = FALSE)
        
      }
      
      # (dim)
      if (length(dim) > 1 ||
          dim <= 0 ||
          !is.finite(dim) ||
          dim != floor(dim)) {
        
        stop ("dim must be a scalar positive integer, but was: ",
              capture.output(dput(dim)),
              call. = FALSE)
        
      }
      
      # check if M > i
      if (M - i < 1) {
        
        stop ("i can be no larger than M - 1",
              call. = FALSE)
        
      }

      # check if sigma is positive real
      # (sigma)
      if (length(sigma) > 1 ||
          sigma <= 0 ||
          !is.finite(sigma)) {
        
        stop ("sigma must be a scalar positive integer, but was: ",
              capture.output(dput(sigma)),
              call. = FALSE)
        
      }
      
      i <- as.greta_array(i)
      M <- as.greta_array(M)
      sigma <- as.greta_array(sigma)
      
      self$bounds <- c(0, Inf)
      super$initialize("ld", dim, truncation = c(0, Inf))
      self$add_parameter(i, "i")
      self$add_parameter(M, "M")
      self$add_parameter(sigma, "sigma")

    },
    
    tf_distrib = function (parameters, dag) {

      i <- parameters$i
      M <- parameters$M
      sigma <- parameters$sigma

      # log pdf(x | i, M, sigma)
      log_prob = function (x) {
        (M - i) * tf$log(x) - x ^ fl(2) / (fl(2) * sigma)
      }

      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)

    },
    
    # no CDF for discrete distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
    
  )
)

M <- 4

zd <- zeros(M-1)
for(i in 1:(M-1)) {
  
  zd[i] <- ld(i, M = M, sigma = 0.5, dim = 1)
  
}
m <- model(zd)
draws <- mcmc(m, chains = 2, n_samples = 2500, warmup = 500, thin = 10)
plot(draws)

coda::effectiveSize(draws)
coda::gelman.diag(draws)

zd <- function(M, sigma) {
  ## initial greta array
  zd <- zeros(M)
  ## elements 1:(M-1) get L&D prior
  for(i in 1:(M-1)) {
    zd[i] <- ld(i, M = M, sigma = sigma, dim = 1)
  }
  ## element M gets truncated normal
  zd[M] <- normal(0, sigma, truncation = c(0, Inf))
  return(zd)
}

m <- model(zd(4, 1))
draws <- mcmc(m, chains = 4, n_samples = 2500, warmup = 500, thin = 10)
plot(draws)

