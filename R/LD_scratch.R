library (greta)

# a temporary hack to test out the distribution
library (R6)
library (tensorflow)
distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node
check_dims <- .internals$utils$checks$check_dims
as.greta_array <- .internals$greta_arrays$as.greta_array
is_scalar <- .internals$utils$misc$is_scalar
fl <- .internals$utils$misc$fl

ld <- function (sigma, dim = 2)
  distrib("ld", sigma, dim)

ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (sigma, dim) {
      
      sigma <- as.greta_array(sigma)
      
      if (!is_scalar(sigma)) {
        stop ("sigma must be a scalar, but had dimensions: ",
              capture.output(dput(dim(eta))),
              call. = FALSE)
      }
      
      # check dim is a positive scalar integer
      dim_old <- dim
      dim <- as.integer(dim)
      if (length(dim) > 1 || dim <= 0 || !is.finite(dim)) {
        
        stop ("dim must be a scalar positive integer, but was: ",
              capture.output(dput(dim_old)),
              call. = FALSE)
        
      }
      
      self$bounds <- c(0, Inf)
      super$initialize("ld", dim, truncation = c(0, Inf))
      self$add_parameter(sigma, "sigma")
      
    },
    
    tf_distrib = function (parameters, dag) {
      
      sigma <- parameters$sigma
      n <- self$dim[1]
      
      const <- fl(n - seq(n))
      
      log_prob = function (x) {
        const * tf$log(x) - x ^ fl(2) / (fl(2) * sigma)
      }
      
      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)
      
    },
    
    # no CDF for discrete distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
    
  )
)

z <- ld(sigma = 10, dim = 3)
m <- model(z)
draws <- mcmc(m, chains = 2)
plot(draws)

dd <- do.call(rbind, draws)
## approximate CV
apply(dd, 2, sd) / apply(dd, 2, mean)


coda::effectiveSize(draws)
coda::gelman.diag(draws)


