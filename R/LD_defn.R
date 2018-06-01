library (greta)

# library(future)

# a temporary hack to test out the distribution
library (R6)
library (tensorflow)
distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node
check_dims <- .internals$utils$checks$check_dims
as.greta_array <- .internals$greta_arrays$as.greta_array
is_scalar <- .internals$utils$misc$is_scalar
fl <- .internals$utils$misc$fl

# function to call 
ld <- function (i, M, C0, dim = 1) {
  distrib("ld", i, M, C0, dim)
  }

# Leung & Drton distn for prior of diag(Z)
ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (i, M, C0, dim) {
      
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

      # check if C0 is positive real
      # (C0)
      if (length(C0) > 1 ||
          C0 <= 0 ||
          !is.finite(C0)) {
        
        stop ("C0 must be a scalar positive integer, but was: ",
              capture.output(dput(C0)),
              call. = FALSE)
        
      }
      
      i <- as.greta_array(i)
      M <- as.greta_array(M)
      C0 <- as.greta_array(C0)
      
      self$bounds <- c(0, Inf)
      super$initialize("ld", dim, truncation = c(0, Inf))
      self$add_parameter(i, "i")
      self$add_parameter(M, "M")
      self$add_parameter(C0, "C0")

    },
    
    tf_distrib = function (parameters, dag) {

      i <- parameters$i
      M <- parameters$M
      C0 <- parameters$C0

      log_prob = function (x) {
        (M - i) * tf$log(x) - x ^ fl(2) / (fl(2) * C0)
      }

      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)

    },
    
    # no CDF for discrete distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
    
  )
)


zd <- zeros(M-1)
for(i in 1:(M-1)) {
  
  zd[i] <- ld(i, M = M, C0 = 0.5, dim = 1)
  
}
m <- model(zd)
draws <- mcmc(m, chains = 2, n_samples = 2500, warmup = 500, thin = 10)
plot(draws)

coda::effectiveSize(draws)
coda::gelman.diag(draws)

