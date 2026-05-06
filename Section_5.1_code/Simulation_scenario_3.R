############################################################
##Simulation 3: Changing number of latent factors
##               and switching loadings
##
##Data-generating process:
##r_m,t ~ N(0, sigma_m,t^2)                    # observed market factor
##z_2,t ~ N(0, sigma_sector,t^2)               # latent sector factor
##z_3,t ~ N(0, sigma_block,t^2)                # latent block factor
##z_4,t ~ N(0, sigma_extra_g2,t^2)             # extra latent factor for group 2
##z_5,t ~ N(0, sigma_extra_g3,t^2)             # extra latent factor for group 3
##z_6,t ~ N(0, sigma_extra_all,t^2)            # extra latent factor for all assets
##
##r_i,t = alpha_i
##           + beta_i * r_m,t
##           + gamma_i,t * z_2,t
##           + delta_i,t * z_3,t
##           + eta_i,t   * z_4,t
##           + theta_i,t * z_5,t
##           + kappa_i,t * z_6,t
##           + e_i,t
##
##e_t ~ N(0, D_t^2)  with diagonal covariance
##
##Regimes:
##t = 1:300   1 factor  (market only / calm)
##t = 301:700   6 factors (market + sector + block + 3 extra / severe crisis)
##t = 701:1000  2 factors (market + sector / milder stressed period)
##
##Blocks:
##3 groups of 10 assets each
##
##Interpretation:
##- Market factor loads on all assets
##- Sector factor loads mainly on groups 1 and 2
##- Block factor loads only on group 1
##- Extra factor 4 loads mainly on group 2
##- Extra factor 5 loads mainly on group 3
##- Extra factor 6 loads on all assets
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript simulation3_run.R <seed> [outdir]")
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
##Converts a 3D array of covariance matrices 
##into a 3D array of correlation matrices, one time slice at a time.
cov_to_corr_array <- function(Sigma_array) {
  N <- dim(Sigma_array)[1]
  TT <- dim(Sigma_array)[3]
  R_array <- array(0, dim = c(N, N, TT))
  for (t in 1:TT) {
    R_array[,,t] <- cov2cor(Sigma_array[,,t])
  }
  R_array
}
##Computes summary statistics of a correlation matrix by group structure
group_corr_summary <- function(R, block_sizes) {
  G <- length(block_sizes)
  grp <- rep(1:G, times = block_sizes)
  within_vals <- c()
  for (g in 1:G) {
    idx <- which(grp == g)
    subR <- R[idx, idx]
    within_vals <- c(within_vals, subR[upper.tri(subR)])
  }
  between_means <- matrix(NA, nrow = G, ncol = G)
  for (g1 in 1:(G - 1)) {
    for (g2 in (g1 + 1):G) {
      vals <- R[grp == g1, grp == g2]
      between_means[g1, g2] <- mean(vals)
      between_means[g2, g1] <- between_means[g1, g2]
    }
  }
  c(within = mean(within_vals), between12 = between_means[1, 2],
    between13 = between_means[1, 3],between23 = between_means[2, 3],
    overall = mean(R[upper.tri(R)])
  )
}
##Creates a binary indicator vector for "active" groups
make_group_indicator <- function(block_sizes, active_groups) {
  grp <- rep(seq_along(block_sizes), times = block_sizes)
  as.numeric(grp %in% active_groups)
}
#############################
##Main simulator
#############################
simulate_scenario3 <- function(
    T = 1000,
    N = 30,
    seed = seed,
    ##Three blocks of 10 assets
    block_sizes = c(10, 10, 10),
    ##Regime lengths
    regime_breaks = c(300, 700),
    ##Regime-specific factor standard deviations
    ##Regime 1: market only
    ##Regime 2: market + sector + block + 3 extra factors
    ##Regime 3: market + sector
    market_sd_regimes = c(0.80, 2.00, 1.20),
    sector_sd_regimes = c(0.00, 1.10, 0.65),
    block_sd_regimes = c(0.00, 0.95, 0.00),
    extra_g2_sd_regimes = c(0.00, 0.85, 0.00),
    extra_g3_sd_regimes = c(0.00, 0.80, 0.00),
    extra_all_sd_regimes = c(0.00, 0.75, 0.00),
    ##Regime-specific idiosyncratic volatility multipliers
    idio_scale_regimes = c(1.00, 1.35, 1.10),
    ##Asset-level parameters
    alpha_sd = 0.03,
    beta_min = 0.7,
    beta_max = 1.4,
    idio_sd_min = 0.70,
    idio_sd_max = 1.20,
    ##Sector-factor loading strength
    ##Sector factor mainly affects groups 1 and 2
    sector_loading_mean = 0.75,
    sector_loading_jitter = 0.08,
    ##Block-factor loading strength
    ##Block factor affects only group 1
    block_loading_mean = 0.90,
    block_loading_jitter = 0.08,
    ##Extra factor loading strengths
    ##Extra factor 4 affects group 2
    extra_g2_loading_mean = 0.80,
    extra_g2_loading_jitter = 0.08,
    ##Extra factor 5 affects group 3
    extra_g3_loading_mean = 0.80,
    extra_g3_loading_jitter = 0.08,
    ##Extra factor 6 affects all assets
    extra_all_loading_mean = 0.65,
    extra_all_loading_jitter = 0.08,
    ##Optional weaker sector loadings in regime 3 than regime 2
    sector_loading_scale_regimes = c(0.00, 1.00, 0.70)
) {
  set.seed(seed)
  if (sum(block_sizes) != N) {
    stop("sum(block_sizes) must equal N.")
  }
  if (length(block_sizes) != 3) {
    stop("This version is written for exactly 3 groups.")
  }
  if (T != 1000) {
    stop("This version assumes T = 1000.")
  }
  #############################
  ##Regimes
  #############################
  regime <- c(rep(1, regime_breaks[1]),rep(2, regime_breaks[2] - regime_breaks[1]),
    rep(3, T - regime_breaks[2]))
  regime_name <- c("Calm_1factor", "Severe_6factors", "Stressed_2factors")[regime]
  ##Number of active factors in each regime
  n_factors_regime <- c(1, 6, 2)
  n_factors_t <- n_factors_regime[regime]
  #############################
  ##Asset-specific parameters
  #############################
  alpha <- rnorm(N, mean = 0, sd = alpha_sd)
  ##Market loadings
  beta <- runif(N, min = beta_min, max = beta_max)
  ##Baseline idiosyncratic standard deviations
  idio_sd_base <- runif(N, min = idio_sd_min, max = idio_sd_max)
  ##Group memberships
  group <- rep(1:3, times = block_sizes)
  #############################
  ##Latent-factor loadings
  #############################
  ##Sector factor mainly affects groups 1 and 2
  sector_indicator <- make_group_indicator(block_sizes, active_groups = c(1, 2))
  gamma_base <- sector_indicator * pmax(
    0,
    rnorm(N, mean = sector_loading_mean, sd = sector_loading_jitter)
  )
  ##Block factor affects only group 1
  block_indicator <- make_group_indicator(block_sizes, active_groups = c(1))
  delta_base <- block_indicator * pmax(
    0,
    rnorm(N, mean = block_loading_mean, sd = block_loading_jitter)
  )
  ##Extra factor 4 affects only group 2
  extra_g2_indicator <- make_group_indicator(block_sizes, active_groups = c(2))
  eta_base <- extra_g2_indicator * pmax(
    0,
    rnorm(N, mean = extra_g2_loading_mean, sd = extra_g2_loading_jitter)
  )
  ##Extra factor 5 affects only group 3
  extra_g3_indicator <- make_group_indicator(block_sizes, active_groups = c(3))
  theta_base <- extra_g3_indicator * pmax(
    0,
    rnorm(N, mean = extra_g3_loading_mean, sd = extra_g3_loading_jitter)
  )
  ##Extra factor 6 affects all assets
  extra_all_indicator <- rep(1, N)
  kappa_base <- extra_all_indicator * pmax(
    0,
    rnorm(N, mean = extra_all_loading_mean, sd = extra_all_loading_jitter)
  )
  #############################
  ##Storage
  #############################
  market <- numeric(T)
  sector_factor <- numeric(T)
  block_factor <- numeric(T)
  extra_g2_factor <- numeric(T)
  extra_g3_factor <- numeric(T)
  extra_all_factor <- numeric(T)
  Y <- matrix(0, nrow = T, ncol = N)
  ##Truth objects
  Sigma_eps_true <- array(0, dim = c(N, N, T))
  Sigma_asset_true <- array(0, dim = c(N, N, T))
  Sigma_joint_true <- array(0, dim = c(N + 1, N + 1, T))
  ##Time-varying loadings (truth)
  gamma_t_store <- matrix(0, nrow = T, ncol = N)
  delta_t_store <- matrix(0, nrow = T, ncol = N)
  eta_t_store <- matrix(0, nrow = T, ncol = N)
  theta_t_store <- matrix(0, nrow = T, ncol = N)
  kappa_t_store <- matrix(0, nrow = T, ncol = N)
  ##Summary series
  sigma_m2_t <- numeric(T)
  avg_offdiag_corr_asset <- numeric(T)
  avg_within_corr_asset <- numeric(T)
  avg_between12_corr_asset <- numeric(T)
  avg_between13_corr_asset <- numeric(T)
  avg_between23_corr_asset <- numeric(T)
  #############################
  ##Simulate over time
  #############################
  for (t in 1:T) {
    g <- regime[t]
    ##Regime-specific factor standard deviations
    sigma_m <- market_sd_regimes[g]
    sigma_sector <- sector_sd_regimes[g]
    sigma_block <- block_sd_regimes[g]
    sigma_extra_g2 <- extra_g2_sd_regimes[g]
    sigma_extra_g3 <- extra_g3_sd_regimes[g]
    sigma_extra_all <- extra_all_sd_regimes[g]
    sigma_m2_t[t] <- sigma_m^2
    ##Regime-specific latent-factor loadings
    ##Regime 1: market only -> all latent loadings = 0
    ##Regime 2: market + sector + block + 3 extra factors
    ##Regime 3: market + sector only
    gamma_t <- gamma_base * sector_loading_scale_regimes[g]
    delta_t <- if (g == 2) delta_base else rep(0, N)
    eta_t <- if (g == 2) eta_base   else rep(0, N)
    theta_t <- if (g == 2) theta_base else rep(0, N)
    kappa_t <- if (g == 2) kappa_base else rep(0, N)
    gamma_t_store[t, ] <- gamma_t
    delta_t_store[t, ] <- delta_t
    eta_t_store[t, ] <- eta_t
    theta_t_store[t, ] <- theta_t
    kappa_t_store[t, ] <- kappa_t
    ##Regime-specific idiosyncratic standard deviations
    idio_sd_t <- idio_sd_base * idio_scale_regimes[g]
    D_t <- diag(idio_sd_t^2)
    ##Residual covariance after removing observed market factor:
    ##Sigma_eps_t = Var(gamma z2 + delta z3 + eta z4 + theta z5 + kappa z6 + e)
    Sigma_eps_t <-
      sigma_sector^2    * tcrossprod(gamma_t) +
      sigma_block^2     * tcrossprod(delta_t) +
      sigma_extra_g2^2  * tcrossprod(eta_t) +
      sigma_extra_g3^2  * tcrossprod(theta_t) +
      sigma_extra_all^2 * tcrossprod(kappa_t) +
      D_t
    ##Full asset covariance:
    ##Sigma_asset_t = sigma_m^2 * beta beta' + Sigma_eps_t
    Sigma_asset_t <- sigma_m^2 * tcrossprod(beta) + Sigma_eps_t
    ##Simulate factors and returns
    r_m_t <- rnorm(1, mean = 0, sd = sigma_m)
    z2_t <- rnorm(1, mean = 0, sd = sigma_sector)
    z3_t <- rnorm(1, mean = 0, sd = sigma_block)
    z4_t <- rnorm(1, mean = 0, sd = sigma_extra_g2)
    z5_t <- rnorm(1, mean = 0, sd = sigma_extra_g3)
    z6_t <- rnorm(1, mean = 0, sd = sigma_extra_all)
    e_t <- rnorm(N, mean = 0, sd = idio_sd_t)
    y_t <- alpha +
      beta * r_m_t +
      gamma_t * z2_t +
      delta_t * z3_t +
      eta_t   * z4_t +
      theta_t * z5_t +
      kappa_t * z6_t +
      e_t
    ##Joint covariance of (market, assets)
    Sigma_joint_t <- matrix(0, nrow = N + 1, ncol = N + 1)
    Sigma_joint_t[1, 1] <- sigma_m^2
    Sigma_joint_t[1, 2:(N + 1)] <- sigma_m^2 * beta
    Sigma_joint_t[2:(N + 1), 1] <- sigma_m^2 * beta
    Sigma_joint_t[2:(N + 1), 2:(N + 1)] <- Sigma_asset_t
    ##Save
    market[t] <- r_m_t
    sector_factor[t] <- z2_t
    block_factor[t] <- z3_t
    extra_g2_factor[t] <- z4_t
    extra_g3_factor[t] <- z5_t
    extra_all_factor[t] <- z6_t
    Y[t, ] <- y_t
    Sigma_eps_true[,,t] <- Sigma_eps_t
    Sigma_asset_true[,,t] <- Sigma_asset_t
    Sigma_joint_true[,,t] <- Sigma_joint_t
    ##Correlation summaries based on ASSET correlations
    R_asset_t <- cov2cor(Sigma_asset_t)
    corr_sum <- group_corr_summary(R_asset_t, block_sizes)
    avg_within_corr_asset[t] <- corr_sum["within"]
    avg_between12_corr_asset[t] <- corr_sum["between12"]
    avg_between13_corr_asset[t] <- corr_sum["between13"]
    avg_between23_corr_asset[t] <- corr_sum["between23"]
    avg_offdiag_corr_asset[t] <- corr_sum["overall"]
  }
  #############################
  ##Build outputs
  #############################
  X <- cbind(Intercept = 1, Mkt.RF = market)
  R_eps_true <- cov_to_corr_array(Sigma_eps_true)
  R_asset_true <- cov_to_corr_array(Sigma_asset_true)
  R_joint_true <- cov_to_corr_array(Sigma_joint_true)
  colnames(Y) <- paste0("Asset", 1:N)
  list(
    X = X,
    Y = Y,
    market = market,
    sector_factor = sector_factor,
    block_factor = block_factor,
    extra_g2_factor = extra_g2_factor,
    extra_g3_factor = extra_g3_factor,
    extra_all_factor = extra_all_factor,
    alpha = alpha,
    beta = beta,
    gamma_base = gamma_base,
    delta_base = delta_base,
    eta_base = eta_base,
    theta_base = theta_base,
    kappa_base = kappa_base,
    gamma_t = gamma_t_store,
    delta_t = delta_t_store,
    eta_t = eta_t_store,
    theta_t = theta_t_store,
    kappa_t = kappa_t_store,
    group = group,
    regime = regime,
    regime_name = regime_name,
    n_factors_t = n_factors_t,
    sigma_m2_t = sigma_m2_t,
    avg_offdiag_corr_asset = avg_offdiag_corr_asset,
    avg_within_corr_asset = avg_within_corr_asset,
    avg_between12_corr_asset = avg_between12_corr_asset,
    avg_between13_corr_asset = avg_between13_corr_asset,
    avg_between23_corr_asset = avg_between23_corr_asset,
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
sim3 <- simulate_scenario3(
  T = 1000,
  N = 30,
  seed = seed
)
#############################
##Diagnostics / example plots
#############################
##Regime boundaries
ablines_sim3 <- c(300, 700)
##Number of active factors over time
plot(sim3$n_factors_t, type = "s", lwd = 2,
     xlab = "Time", ylab = "Number of active factors",
     main = "True number of active factors")
abline(v = ablines_sim3, lty = 2, col = "red")
##Mean true ASSET correlation by time
plot(sim3$avg_offdiag_corr_asset, type = "l", lwd = 2,
     xlab = "Time", ylab = "Mean off-diagonal correlation",
     main = "True ASSET correlation under changing factor structure")
abline(v = ablines_sim3, lty = 2, col = "red")
##Group-based asset correlation summaries
plot(sim3$avg_within_corr_asset, type = "l", lwd = 2,
     ylim = range(c(sim3$avg_within_corr_asset,
                    sim3$avg_between12_corr_asset,
                    sim3$avg_between13_corr_asset,
                    sim3$avg_between23_corr_asset)),
     xlab = "Time", ylab = "Correlation",
     main = "True ASSET correlations by regime")
lines(sim3$avg_between12_corr_asset, lwd = 2, lty = 2)
lines(sim3$avg_between13_corr_asset, lwd = 2, lty = 3)
lines(sim3$avg_between23_corr_asset, lwd = 2, lty = 4)
abline(v = ablines_sim3, lty = 2, col = "red")
legend("topleft",
       legend = c("Within-group", "Between G1-G2", "Between G1-G3", "Between G2-G3"),
       lwd = 2, lty = c(1, 2, 3, 4), bty = "n")
##Regime-wise averages
tapply(sim3$avg_offdiag_corr_asset, sim3$regime_name, mean)
tapply(sim3$avg_within_corr_asset, sim3$regime_name, mean)
tapply(sim3$avg_between12_corr_asset, sim3$regime_name, mean)
tapply(sim3$avg_between13_corr_asset, sim3$regime_name, mean)
tapply(sim3$avg_between23_corr_asset, sim3$regime_name, mean)
##Example: true total correlation between two assets in group 1
plot(sim3$R_asset_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset2")
abline(v = ablines_sim3, lty = 2, col = "red")
##Example: group 1 vs group 2
plot(sim3$R_asset_true[1, 12, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset12")
abline(v = ablines_sim3, lty = 2, col = "red")
##Example: group 1 vs group 3
plot(sim3$R_asset_true[1, 25, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset25")
abline(v = ablines_sim3, lty = 2, col = "red")