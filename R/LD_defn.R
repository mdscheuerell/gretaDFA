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

ld <- function (i, M, C0, dim = 1) {
  distrib("ld", i, M, C0, dim)
  }

ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (i, M, C0, dim) {
      
      i <- as.greta_array(i)
      M <- as.greta_array(M)
      C0 <- as.greta_array(C0)

      if (!is_scalar(C0)) {
        stop ("C0 must be a scalar, but had dimensions: ",
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

z1 <- ld(i = 1, M = 4, C0 = 1, dim = 1)
z2 <- ld(i = 2, M = 4, C0 = 1, dim = 1)
z3 <- ld(i = 3, M = 4, C0 = 1, dim = 1)

i <- model(z1, z2, z3)

# plan(multicore)

draws <- mcmc(i, chains = 2, n_samples = 1000, warmup = 500, thin = 4)
plot(draws)



zd <- zeros(M-1)
for(i in 1:(M-1)) {
  
  zd[i] <- ld(i, M = M, C0 = 1, dim = 1)
  
}
m <- model(zd)
draws <- mcmc(m, chains = 2, n_samples = 2000, warmup = 500, thin = 4)
plot(draws)

coda::effectiveSize(draws)
coda::gelman.diag(draws)







dd <- do.call(rbind, draws)
## approximate CV
apply(dd, 2, sd) / apply(dd, 2, mean)

apply(dd, 2, quantile)
dd <- apply(dd, 2, sort)



coda::effectiveSize(draws)
coda::gelman.diag(draws)


npdf <- function(x, mu = 0, C0 = 1) {
  pp <- exp(-(x - mu)^2 / (2*C0^2))
  np <- pp/sqrt(2*pi*C0^2)
  return(np)
}

xx <- seq(-500,500)/100
plot(xx, npdf(xx), type="l")




pld <- function(X, C0 = 1) {
  n <- length(X)
  pdf <- X^(n - seq(n)) * exp(- X^2 / (2 * C0))
  return(pdf)
}

xx <- seq(1000)/100
yy <- matrix(xx,length(xx),2)
pp <- apply(yy, 1, pld)

matplot(xx, t(apply(yy, 1, pld)), type="l")


N <- 3
C0 <- 1


const * tf$log(x) - x ^ fl(2) / (fl(2) * C0)
