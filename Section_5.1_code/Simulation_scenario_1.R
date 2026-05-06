############################################################
##Simulation 1: Regime changes / structural breaks
##
##Data-generating process:
##r_m,t ~ N(0, sigma_m,t^2)
##r_i,t = alpha_i + beta_i * r_m,t + epsilon_i,t
##epsilon_t ~ N(0, Sigma_eps,t)
##
##Regimes (applied to ASSET correlations):
##   t =   1:250   low correlation
##   t = 251:500   crisis jump to high correlation
##   t = 501:750   moderate correlation
##   t = 751:1000  second break with a different pattern
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript simulation1_run.R <seed> [outdir]")
}
seed <- as.integer(args[1])
if (is.na(seed)) stop("seed must be an integer")
outdir <- if (length(args) >= 2) args[2] else getwd()
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

library(mvtnorm)
library(Matrix)

#############################
##Helper functions
#############################

##Creates an N x N equicorrelation matrix with off-diagonal correlations equal to rho.
make_equicorr <- function(N, rho) {
  R <- matrix(rho, nrow = N, ncol = N)
  diag(R) <- 1
  R
}

##Creates a block correlation matrix with stronger within-block correlations 
##and weaker between-block correlations, adjusting to the nearest valid correlation matrix if needed.
make_blockcorr <- function(block_sizes, rho_within = 0.65, rho_between = 0.35) {
  N <- sum(block_sizes)
  R <- matrix(rho_between, nrow = N, ncol = N)
  diag(R) <- 1
  start <- 1
  for (b in seq_along(block_sizes)) {
    idx <- start:(start + block_sizes[b] - 1)
    R[idx, idx] <- rho_within
    diag(R[idx, idx]) <- 1
    start <- max(idx) + 1
  }
  R <- (R + t(R)) / 2
  eigmin <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (eigmin <= 1e-8) {
    R <- as.matrix(nearPD(R, corr = TRUE)$mat)
  }
  R
}

##Converts a 3D array of covariance matrices 
##into a 3D array of correlation matrices, one time slice at a time.
cov_to_corr_array <- function(Sigma_array) {
  N <- dim(Sigma_array)[1]
  T <- dim(Sigma_array)[3]
  R_array <- array(0, dim = c(N, N, T))
  for (t in 1:T) {
    R_array[,,t] <- cov2cor(Sigma_array[,,t])
  }
  R_array
}

##Build a target asset covariance matrix with the desired asset correlation
##structure, then back out the implied idiosyncratic covariance:
##Sigma_eps = Sigma_asset - sigma_m^2 * beta beta'
##
##If needed, inflate all asset standard deviations by a common factor until
##Sigma_eps is positive definite. This preserves the target asset correlations.
make_feasible_asset_cov <- function(R_asset_target,asset_sd,beta,sigma_m2,
                                    tol = 1e-8,inflate_step = 1.02,max_iter = 500) {
  inflate <- 1
  for (iter in 0:max_iter) {
    sd_now <- asset_sd * inflate
    D_now <- diag(sd_now)
    Sigma_asset <- D_now %*% R_asset_target %*% D_now
    Sigma_eps <- Sigma_asset - sigma_m2 * tcrossprod(beta)
    Sigma_eps <- (Sigma_eps + t(Sigma_eps)) / 2
    ee <- tryCatch(
      eigen(Sigma_eps, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    if (all(is.finite(ee)) && min(ee) > tol) {
      return(list(
        Sigma_asset = Sigma_asset,
        Sigma_eps = Sigma_eps,
        asset_sd = sd_now,
        inflate = inflate
      ))
    }
    inflate <- inflate * inflate_step
  }
  stop("Could not construct a positive definite Sigma_eps. Try reducing market_sd_regimes or beta range, or increasing asset_sd_min / asset_sd_max.")
}

#############################
##Main simulator
#############################
simulate_scenario1 <- function(
    T = 1000,
    N = 30,
    seed = seed,
    ##If TRUE, regime 4 uses a block ASSET correlation structure
    ##to reflect a "different pattern".
    ##If FALSE, regime 4 uses equicorrelation with rho = 0.50.
    regime4_block = TRUE,
    ##Regime-specific off-diagonal ASSET correlation levels
    rho_regimes = c(0.10, 0.70, 0.25, 0.50),
    ## Regime-specific market standard deviations (in %)
    market_sd_regimes = c(0.80, 2.20, 1.10, 1.60),
    ## Regime-specific multipliers for ASSET volatilities
    asset_scale_regimes = c(1.00, 1.35, 1.10, 1.25),
    ##Asset-level parameters
    beta_min = 0.7,
    beta_max = 1.4,
    alpha_sd = 0.03,          #small alpha in percentage units
    asset_sd_min = 1.20,
    asset_sd_max = 2.00
) {
  set.seed(seed)
  if (T != 1000) {
    stop("This version assumes T = 1000 so that each regime has length 250.")
  }
  ##Regime labels
  regime <- rep(1:4, each = 250)
  ##Asset-specific CAPM parameters
  alpha <- rnorm(N, mean = 0, sd = alpha_sd)
  beta <- runif(N, min = beta_min, max = beta_max)
  ##Baseline asset standard deviations
  asset_sd_base <- runif(N, min = asset_sd_min, max = asset_sd_max)
  ##block structure for regime 4
  ##With N = 30, this makes three groups of size 10 by default
  if (N %% 3 == 0) {
    block_sizes <- rep(N / 3, 3)
  } else {
    b1 <- floor(N / 3)
    b2 <- floor(N / 3)
    b3 <- N - b1 - b2
    block_sizes <- c(b1, b2, b3)
  }
  ##Precompute regime-specific target ASSET correlation matrices
  R_asset_regimes <- vector("list", length = 4)
  for (g in 1:4) {
    if (g < 4) {
      R_asset_regimes[[g]] <- make_equicorr(N, rho_regimes[g])
    } else {
      if (regime4_block) {
        ##Different pattern in regime 4:
        ##within-block correlation is high,
        ##between-block correlation is more moderate.
        R_asset_regimes[[g]] <- make_blockcorr(
          block_sizes = block_sizes,
          rho_within = 0.65,
          rho_between = 0.35
        )
      } else {
        R_asset_regimes[[g]] <- make_equicorr(N, rho_regimes[4])
      }
    }
  }
  ##Precompute regime-specific covariance objects
  regime_objects <- vector("list", length = 4)
  for (g in 1:4) {
    sigma_m2_g <- market_sd_regimes[g]^2
    asset_sd_g <- asset_sd_base * asset_scale_regimes[g]
    R_asset_g <- R_asset_regimes[[g]]
    tmp <- make_feasible_asset_cov(
      R_asset_target = R_asset_g,
      asset_sd = asset_sd_g,
      beta = beta,
      sigma_m2 = sigma_m2_g
    )
    regime_objects[[g]] <- list(
      R_asset = R_asset_g,
      Sigma_asset = tmp$Sigma_asset,
      Sigma_eps = tmp$Sigma_eps,
      asset_sd = tmp$asset_sd,
      inflate = tmp$inflate
    )
  }
  ##Storage
  market <- numeric(T)
  Y <- matrix(0, nrow = T, ncol = N)
  ##Truth objects
  Sigma_eps_true <- array(0, dim = c(N, N, T))
  Sigma_asset_true <- array(0, dim = c(N, N, T))
  Sigma_joint_true <- array(0, dim = c(N + 1, N + 1, T))
  ##Store regime-specific quantities
  sigma_m2_t <- numeric(T)
  avg_offdiag_corr_asset <- numeric(T)
  for (t in 1:T) {
    g <- regime[t]
    ##Regime-specific market volatility
    sigma_m <- market_sd_regimes[g]
    sigma_m2 <- sigma_m^2
    sigma_m2_t[t] <- sigma_m2
    ##Regime-specific TRUE ASSET covariance / implied idiosyncratic covariance
    R_asset_t <- regime_objects[[g]]$R_asset
    Sigma_asset_t <- regime_objects[[g]]$Sigma_asset
    Sigma_eps_t <- regime_objects[[g]]$Sigma_eps
    ##Simulate market and idiosyncratic shocks
    r_m_t <- rnorm(1, mean = 0, sd = sigma_m)
    eps_t <- as.numeric(rmvnorm(1, mean = rep(0, N), sigma = Sigma_eps_t))
    ##Asset excess returns
    y_t <- alpha + beta * r_m_t + eps_t
    ##True joint covariance of (market, assets)
    ##[ Var(m)              Cov(m, y)' ]
    ##[ Cov(m, y)           Cov(y)      ]
    Sigma_joint_t <- matrix(0, nrow = N + 1, ncol = N + 1)
    Sigma_joint_t[1, 1] <- sigma_m2
    Sigma_joint_t[1, 2:(N + 1)] <- sigma_m2 * beta
    Sigma_joint_t[2:(N + 1), 1] <- sigma_m2 * beta
    Sigma_joint_t[2:(N + 1), 2:(N + 1)] <- Sigma_asset_t
    ##Save everything
    market[t] <- r_m_t
    Y[t, ] <- y_t
    Sigma_eps_true[,,t] <- Sigma_eps_t
    Sigma_asset_true[,,t] <- Sigma_asset_t
    Sigma_joint_true[,,t] <- Sigma_joint_t
    avg_offdiag_corr_asset[t] <- mean(R_asset_t[row(R_asset_t) != col(R_asset_t)])
  }
  ##Build design matrix for your model
  X <- cbind(Intercept = 1, Mkt.RF = market)
  ##Correlation truths
  R_eps_true <- cov_to_corr_array(Sigma_eps_true)
  R_asset_true <- cov_to_corr_array(Sigma_asset_true)
  R_joint_true <- cov_to_corr_array(Sigma_joint_true)
  colnames(Y) <- paste0("Asset", 1:N)
  list(
    X = X,
    Y = Y,
    market = market,
    alpha = alpha,
    beta = beta,
    regime = regime,
    sigma_m2_t = sigma_m2_t,
    regime_asset_inflation = sapply(regime_objects, `[[`, "inflate"),
    avg_offdiag_corr_asset = avg_offdiag_corr_asset,
    Sigma_eps_true = Sigma_eps_true,
    R_eps_true = R_eps_true,
    Sigma_asset_true = Sigma_asset_true,
    R_asset_true = R_asset_true,
    Sigma_joint_true = Sigma_joint_true,
    R_joint_true = R_joint_true
  )
}

#############################
##Generate one dataset
#############################

sim1 <- simulate_scenario1(T = 1000, N = 30,seed = seed,regime4_block = TRUE)
##Mean true pairwise ASSET correlation by time
plot(sim1$avg_offdiag_corr_asset, type = "l", lwd = 2,
     xlab = "Time", ylab = "Mean off-diagonal correlation",
     main = "True regime-switching ASSET correlation")
abline(v = c(250, 500, 750), lty = 2, col = "red")
##Regime-wise average
tapply(sim1$avg_offdiag_corr_asset, sim1$regime, mean)
##Example: true correlation between Asset1 and Asset2 over time
plot(sim1$R_asset_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset2")
abline(v = c(250, 500, 750), lty = 2, col = "red")