## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(tensorflow.one_based_extract = TRUE)

## ---- eval=FALSE, message=FALSE------------------------------------------
## devtools::install_github('greta-dev/greta@dev')

## ----load_pkgs_1, message=FALSE------------------------------------------
library(tensorflow)
library(greta)
installed.packages()[c("greta", "tensorflow"), "Version"]
tensorflow::tf_version()

## ----load_pkgs_2, message=FALSE------------------------------------------
library(MASS)
library(viridis)
library(corrplot)

## ----load_greta, message=FALSE-------------------------------------------
op <- greta::.internals$nodes$constructors$op
cumsum_mat <- function(x, a=1L) {
  ## cumsums across matrix rows/cols
  ## rows: a = 1L
  ## cols: a = 0L
  stopifnot(length(dim(x)) == 2)
  op("cumsum_mat",
     x,
     tf_operation = "tf$cumsum",
     operation_args = list(axis = a))
}

## ------------------------------------------------------------------------
zd <- function(M, sigma) {
  ## initial greta array
  zd <- zeros(M-1)
  ## elements 1:(M-1) get L&D prior
  for(i in 1:(M-1)) {
    zd[i] <- ld(i, M = M, sigma = sigma, dim = 1)
  }
  ## element M gets truncated normal
  # zd[M] <- normal(0, sigma, truncation = c(0, Inf))
  return(zd)
}

## ----define_ld_distrib, echo = FALSE-------------------------------------
library (R6)
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

      # check if sigma is scalar
      # (sigma)
      if (length(sigma) > 1) {
        
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

## ----data_inits----------------------------------------------------------
NN <- 21
TT <- 30
MM <- 3

## ----sim_factors, message=FALSE------------------------------------------
set.seed(123)
## MM x TT matrix of innovations
ww <- matrix(rnorm(MM*TT, 0, 1), MM, TT)
ww[,1] <- rnorm(MM, 0, sqrt(5))
## MM x TT matrix of latent trends
xx <- t(apply(ww,1,cumsum))

## ----plot_factors, fig.height=4, fig.width=7, echo=FALSE-----------------
## plot the trends
par(mai=c(0.8,0.8,0.2,0.2))
matplot(t(xx), type="b", lty="solid", cex=0.7,
        xlab="Time", ylab=expression(italic(x)[italic(t)]),
        col=plasma(MM, end=0.8))

## ----create_Z------------------------------------------------------------
## NN x MM loadings matrix
# ZZ <- matrix(NA, NN, MM)
# LL <- seq(-1,1,length.out=MM+1)
# for(i in 1:MM) {
#   zz <- NULL
#   for(j in 1:MM) {
#     zz <- c(zz,runif(NN/MM,LL[j],LL[j+1]))
#   }
#   ZZ[sample(NN,NN),i] <- round(zz,3)
# }

ZZ <- matrix(runif(NN*MM, -1, 1), NN, MM)
diag(ZZ) <- abs(diag(ZZ))
ZZ[upper.tri(ZZ)] <- 0

## ----create_ts-----------------------------------------------------------
## obs errors
ee <- t(mvrnorm(TT, matrix(0,NN,1), diag(0.2,NN,NN)))
## NN x TT matrix of observed data
yy <- ZZ %*% xx + ee

## ----plot_ts, fig.height=4, fig.width=7, echo=FALSE----------------------
par(mai=c(0.8,0.8,0.2,0.2))
matplot(t(yy), type="l", lty="solid",
        xlab="Time", ylab=expression(italic(y)[italic(i)]),
        col=plasma(NN, alpha=0.7, end=0.8))

## ----hist_cor_coefs, fig.height=7, fig.width=7, echo=FALSE---------------
rho <- cor(t(yy))
par(mai=c(0.8,0.8,0.2,0.2))
corrplot(rho, method="ellipse", type="lower", order = "hclust",
         tl.col = "black", tl.srt = 0, tl.cex = 0.6, tl.offset = 0.7,
         cl.cex = 0.8, cl.offset = 0.9, cl.ratio = 0.1)

## ----priors_loadings_hack, eval=FALSE------------------------------------
## ## empty loadings matrix
## ZZ_est <- zeros(NN,MM)
## ## define sigma
## sigma = normal(0, 5, truncation = c(0, Inf))
## ## diagonal
## diag(ZZ_est) = zd(MM, sigma)
## ## sub-diagonal
## idx_s <- lower.tri(ZZ_est)
## ZZ_est[idx_s] = normal(0, 1, dim = sum(idx_s))

## ----priors_loadings-----------------------------------------------------
## empty loadings matrix
ZZ_est <- zeros(NN,MM)
## definea sigma
sigma = normal(0, 5, truncation = c(0, Inf))
## diagonal for elements 1:(M-1)
idx_d <- row(ZZ_est) == col(ZZ_est) & row(ZZ_est) != MM
ZZ_est_raw_d = zd(MM, sigma)
ZZ_est[idx_d] <- ZZ_est_raw_d
## diagonal for element M
ZZ_est_raw_M = normal(0, sigma, dim = 1, truncation = c(0, Inf))
ZZ_est[MM,MM] <- ZZ_est_raw_M
## sub-diagonal
idx_s <- lower.tri(ZZ_est)
ZZ_est_raw_s = normal(0, 1, dim = sum(idx_s))
ZZ_est[idx_s] <- ZZ_est_raw_s

## ----demo_evaluate, cache=TRUE-------------------------------------------
head(calculate(ZZ_est, list(ZZ_est_raw_d = rep(1, sum(idx_d)),
                            ZZ_est_raw_M = matrix(2),
                            ZZ_est_raw_s = rep(3, sum(idx_s)))))

## ----priors_obs_cov------------------------------------------------------
## diagonal of R
RR_est_raw = student(10, 0, 1, truncation=c(0, Inf))

## ----like_factors--------------------------------------------------------
## inital factor
X0 = normal(0, sqrt(10), c(MM, 1))
## factors for t = 2-TT
XX = normal(0, 1, c(MM, TT))
## combine & cumsum over rows for the random walk
xx_est <- cumsum_mat(cbind(X0,XX))[,-1]

## ----like_obs, eval=TRUE, echo=TRUE--------------------------------------
## scale data
yy_z <- t(scale(t(yy)))
## vectorize data
yy_vec <- as_data(matrix(yy_z, 1, NN*TT))
## vectorize mean
Zx_vec <- c(ZZ_est %*% xx_est)
## create expanded cov matrix
RR_star <- zeros(NN*TT, NN*TT)
diag(RR_star) <- rep(RR_est_raw, NN*TT)
## define likelihood
distribution(yy_vec) = multivariate_normal(Zx_vec, RR_star)

## ----greta_model, cache=TRUE---------------------------------------------
mod_fit <- model(xx_est, ZZ_est, RR_est_raw)
mcmc_smpl <- mcmc(mod_fit, n_samples = 1000, thin = 1, warmup = 1000,
                  chains = 1, verbose = FALSE)

## ----traceplots, fig.height=5, fig.width=7-------------------------------
par(mfrow=c(2,2), mai=c(0.8,0.5,0.2,0.2))
for(i in 1:4) {
  coda::traceplot(mcmc_smpl[[1]][,i]) 
}

## ----bayesplot, fig.height=5, fig.width=7, eval=FALSE, echo=FALSE--------
## par(mai=c(0.8,0.5,0.2,0.2))
## bayesplot::mcmc_trace(as.matrix(mcmc_smpl[[1]][,1:4]))

## ----converge_chk--------------------------------------------------------
## number of estimated params/states
n_par <- dim(mcmc_smpl[[1]])[2]
## H&W diag: pass (1) or fail (NA)
halfwidth_test <- coda::heidel.diag(mcmc_smpl[[1]])[,4]
## proportion passing
round(sum(halfwidth_test, na.rm = TRUE) / n_par, 2)

## ----examine_fits--------------------------------------------------------
mod_smry <- summary(mcmc_smpl)
par_means <- mod_smry$statistics[,"Mean"]

## ----get_ZZ_fit----------------------------------------------------------
ZZ_fit <- par_means[grepl("ZZ",rownames(mod_smry$statistics))]
ZZ_fit <- matrix(ZZ_fit, NN, MM, byrow=FALSE)
round(ZZ_fit, 2)
## rotation matrix
HH_inv <- varimax(ZZ_fit)$rotmat
## rotated Z
ZZ_rot <- ZZ_fit %*% HH_inv
round(ZZ_rot, 2)

## ----get_factors---------------------------------------------------------
# fitted factors
xx_fit <- par_means[grepl("xx",rownames(mod_smry$statistics))]
xx_fit <- matrix(xx_fit, MM, TT, byrow=FALSE)

## ----plot_xx_fits, fig.height=7, fig.width=7, echo=FALSE-----------------
par(mfrow=c(2,1), mai=c(0.8,0.8,0.2,0.2))
matplot(t(xx), type="b", lty="solid", cex=0.7,
        xlab="Time", ylab=expression(italic(x)[italic(t)]),
        col=plasma(MM, end=0.8))
matplot(t(xx_fit), type="l", lty="solid", cex=0.7,
        xlab="Time", ylab=expression(italic(x)[italic(t)]),
        col="darkgray")

## ----corrplot_xx, fig.height=6, fig.width=6------------------------------
par(mai=rep(0.1, 4), omi=rep(0.1, 4))
corrplot(cor(t(xx), t(xx_fit)), method="ellipse",
         tl.col = "black", tl.srt = 0, tl.cex = 0.8, tl.offset = 0.7,
         cl.cex = 0.8, cl.offset = 0.9, cl.ratio = 0.2)

## ----get_RR_fit----------------------------------------------------------
RR_fit <- par_means[grepl("RR",rownames(mod_smry$statistics))]
round(RR_fit, 2)

## ----cor_yy--------------------------------------------------------------
## fitted values
yy_fit <- ZZ_rot %*% xx
## corrrelation
cor_yy <- matrix(diag(cor(t(yy_z), t(yy_fit))), NN/10, 10)
## plots
par(mai=rep(0.1, 4), omi=rep(0.1, 4))
corrplot(cor_yy, method="ellipse",
         tl.pos = "n",
         cl.cex = 0.8, cl.offset = 0.9,
         cl.ratio = 0.2, )#cl.lim = c(0, 1))

