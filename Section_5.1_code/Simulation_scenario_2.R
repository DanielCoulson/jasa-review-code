############################################################
##Simulation 2: Sparse dependence / crisis clustering
##
##Data-generating process:
##r_m,t ~ N(0, sigma_m,t^2)
##r_i,t = alpha_i + beta_i * r_m,t + epsilon_i,t
##epsilon_t ~ N(0, Sigma_eps,t)
##
##Groups:
##3 groups of 10 assets each
##
##Calm periods (ASSET correlations):
##within-group correlation = 0.20
##between-group correlation = 0.00
##
##Crisis periods (ASSET correlations):
##within-group correlation = 0.75
##group 1 with group 2 correlation = 0.35
##all other between-group correlations = 0.05
##
##Crisis clustering:
##multiple crisis episodes are inserted over time
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Simulation2_run.R <seed> [outdir]")
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
##Construct a block-structured correlation matrix with specified within-
##and between-group correlations, ensuring it is positive definite.
make_sparse_groupcorr <- function(block_sizes,rho_within,rho_between_mat) {
  G <- length(block_sizes)
  N <- sum(block_sizes)
  if (!all(dim(rho_between_mat) == c(G, G))) {
    stop("rho_between_mat must be a G x G matrix where G = length(block_sizes).")
  }
  grp <- rep(1:G, times = block_sizes)
  R <- matrix(0, nrow = N, ncol = N)
  diag(R) <- 1
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        if (grp[i] == grp[j]) {
          R[i, j] <- rho_within
        } else {
          R[i, j] <- rho_between_mat[grp[i], grp[j]]
        }
      }
    }
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
  TT <- dim(Sigma_array)[3]
  R_array <- array(0, dim = c(N, N, TT))
  for (t in 1:TT) {
    R_array[, , t] <- cov2cor(Sigma_array[, , t])
  }
  R_array
}
##Compute summary statistics of correlations by group:
##- Average within-group correlation
##- Pairwise between-group correlations (for first 3 groups)
##- Overall average correlation (upper triangle)
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
  c(within = mean(within_vals),
    between12 = between_means[1, 2],
    between13 = between_means[1, 3],
    between23 = between_means[2, 3],
    overall = mean(R[upper.tri(R)])
  )
}

##Build a target asset covariance matrix with the desired asset correlation
##structure, then back out the implied idiosyncratic covariance:
##Sigma_eps = Sigma_asset - sigma_m^2 * beta beta'
##
##If needed, inflate all asset standard deviations by a common factor until
##Sigma_eps is positive definite. This preserves the target asset correlations.
make_feasible_asset_cov <- function(R_asset_target,asset_sd,beta,sigma_m2,
                                    tol = 1e-8,inflate_step = 1.02, max_iter = 500) {
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
simulate_scenario2 <- function(T = 1000,N = 30,seed = seed,
    ##Three groups of 10 assets
    block_sizes = c(10, 10, 10),
    ##Crisis episodes (start, end)
    crisis_windows = rbind(
      c(251, 350),
      c(601, 700),
      c(851, 925)),
    ##TARGET ASSET correlation structure
    calm_within = 0.20,
    calm_between = 0.00,
    crisis_within = 0.75,
    crisis_between12 = 0.35,
    crisis_between_other = 0.05,
    ##Regime-specific market standard deviations (in %)
    market_sd_regimes = c(0.80, 2.20),   # calm, crisis
    ##Asset-volatility regime multipliers
    asset_scale_regimes = c(1.00, 1.25), # calm, crisis
    ##Asset-level parameters
    beta_min = 0.7,
    beta_max = 1.4,
    alpha_sd = 0.03,
    asset_sd_min = 1.20,
    asset_sd_max = 2.00
) {
  set.seed(seed)
  if (sum(block_sizes) != N) {
    stop("sum(block_sizes) must equal N.")
  }
  if (length(block_sizes) != 3) {
    stop("This version is written for exactly 3 groups.")
  }
  ##Regime labels: 1 = calm, 2 = crisis
  regime <- rep(1, T)
  for (k in 1:nrow(crisis_windows)) {
    s <- crisis_windows[k, 1]
    e <- crisis_windows[k, 2]
    regime[s:e] <- 2
  }
  regime_name <- ifelse(regime == 1, "Calm", "Crisis")
  ##Asset-specific CAPM parameters
  alpha <- rnorm(N, mean = 0, sd = alpha_sd)
  beta <- runif(N, min = beta_min, max = beta_max)
  ##Baseline asset standard deviations
  asset_sd_base <- runif(N, min = asset_sd_min, max = asset_sd_max)
  ##Group memberships
  group <- rep(1:3, times = block_sizes)
  ##Target ASSET correlation matrices
  rho_between_calm <- matrix(c(1,calm_between,calm_between,calm_between,  
                               1,calm_between,calm_between, calm_between, 1), nrow = 3, byrow = TRUE)
  rho_between_crisis <- matrix(c(1, crisis_between12, crisis_between_other, 
                                 crisis_between12, 1, crisis_between_other, crisis_between_other,  
                                 crisis_between_other, 1), nrow = 3, byrow = TRUE)
  R_asset_calm <- make_sparse_groupcorr(
    block_sizes = block_sizes,
    rho_within = calm_within,
    rho_between_mat = rho_between_calm)
  R_asset_crisis <- make_sparse_groupcorr(
    block_sizes = block_sizes,
    rho_within = crisis_within,
    rho_between_mat = rho_between_crisis)
  ##Precompute regime-specific covariance objects
  regime_objects <- vector("list", length = 2)
  for (g in 1:2) {
    sigma_m2_g <- market_sd_regimes[g]^2
    asset_sd_g <- asset_sd_base * asset_scale_regimes[g]
    R_asset_g <- if (g == 1) R_asset_calm else R_asset_crisis
    tmp <- make_feasible_asset_cov(
      R_asset_target = R_asset_g,
      asset_sd = asset_sd_g,
      beta = beta,
      sigma_m2 = sigma_m2_g)
    regime_objects[[g]] <- list(
      R_asset = R_asset_g,
      Sigma_asset = tmp$Sigma_asset,
      Sigma_eps = tmp$Sigma_eps,
      asset_sd = tmp$asset_sd,
      inflate = tmp$inflate)
  }
  ##Storage
  market <- numeric(T)
  Y <- matrix(0, nrow = T, ncol = N)
  ##Truth objects
  Sigma_eps_true <- array(0, dim = c(N, N, T))
  Sigma_asset_true <- array(0, dim = c(N, N, T))
  Sigma_joint_true <- array(0, dim = c(N + 1, N + 1, T))
  ##Summary series: now based on ASSET correlations
  sigma_m2_t <- numeric(T)
  avg_offdiag_corr_asset <- numeric(T)
  avg_within_corr_asset <- numeric(T)
  avg_between12_corr_asset <- numeric(T)
  avg_between13_corr_asset <- numeric(T)
  avg_between23_corr_asset <- numeric(T)
  for (t in 1:T) {
    g <- regime[t]
    sigma_m <- market_sd_regimes[g]
    sigma_m2 <- sigma_m^2
    sigma_m2_t[t] <- sigma_m2
    R_asset_t <- regime_objects[[g]]$R_asset
    Sigma_asset_t <- regime_objects[[g]]$Sigma_asset
    Sigma_eps_t <- regime_objects[[g]]$Sigma_eps
    ##Simulate market and idiosyncratic shocks
    r_m_t <- rnorm(1, mean = 0, sd = sigma_m)
    eps_t <- as.numeric(rmvnorm(1, mean = rep(0, N), sigma = Sigma_eps_t))
    ##Asset excess returns
    y_t <- alpha + beta * r_m_t + eps_t
    ##True joint covariance of (market, assets)
    Sigma_joint_t <- matrix(0, nrow = N + 1, ncol = N + 1)
    Sigma_joint_t[1, 1] <- sigma_m2
    Sigma_joint_t[1, 2:(N + 1)] <- sigma_m2 * beta
    Sigma_joint_t[2:(N + 1), 1] <- sigma_m2 * beta
    Sigma_joint_t[2:(N + 1), 2:(N + 1)] <- Sigma_asset_t
    ##Save everything
    market[t] <- r_m_t
    Y[t, ] <- y_t
    Sigma_eps_true[, , t] <- Sigma_eps_t
    Sigma_asset_true[, , t] <- Sigma_asset_t
    Sigma_joint_true[, , t] <- Sigma_joint_t
    ##Summary of the true ASSET correlation structure
    corr_sum <- group_corr_summary(R_asset_t, block_sizes)
    avg_within_corr_asset[t] <- corr_sum["within"]
    avg_between12_corr_asset[t] <- corr_sum["between12"]
    avg_between13_corr_asset[t] <- corr_sum["between13"]
    avg_between23_corr_asset[t] <- corr_sum["between23"]
    avg_offdiag_corr_asset[t] <- corr_sum["overall"]
  }
  ##Build design matrix
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
    group = group,
    regime = regime,
    regime_name = regime_name,
    crisis_windows = crisis_windows,
    sigma_m2_t = sigma_m2_t,
    regime_asset_inflation = sapply(regime_objects, `[[`, "inflate"),
    ##ASSET-correlation summaries (primary target)
    avg_offdiag_corr_asset = avg_offdiag_corr_asset,
    avg_within_corr_asset = avg_within_corr_asset,
    avg_between12_corr_asset = avg_between12_corr_asset,
    avg_between13_corr_asset = avg_between13_corr_asset,
    avg_between23_corr_asset = avg_between23_corr_asset,
    ##Truth objects
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
sim2 <- simulate_scenario2(
  T = 1000,
  N = 30,
  seed = seed
)
#############################
##Diagnostics / example plots
#############################
##Vertical lines at crisis boundaries
vlines <- sort(unique(as.vector(t(sim2$crisis_windows))))
##Plot true ASSET correlation summaries
plot(sim2$avg_within_corr_asset, type = "l", lwd = 2,
     ylim = c(0, 0.8),
     xlab = "Time", ylab = "Correlation",
     main = "True ASSET correlations: sparse dependence + crisis clustering")
lines(sim2$avg_between12_corr_asset, lwd = 2, lty = 2)
lines(sim2$avg_between13_corr_asset, lwd = 2, lty = 3)
lines(sim2$avg_between23_corr_asset, lwd = 2, lty = 4)
abline(v = vlines, lty = 2, col = "red")
legend("topleft",
       legend = c("Within-group", "Between G1-G2", "Between G1-G3", "Between G2-G3"),
       lwd = 2, lty = c(1, 2, 3, 4), bty = "n")
##Regime-wise averages: should now match the TARGET ASSET correlations
tapply(sim2$avg_within_corr_asset, sim2$regime_name, mean)
tapply(sim2$avg_between12_corr_asset, sim2$regime_name, mean)
tapply(sim2$avg_between13_corr_asset, sim2$regime_name, mean)
tapply(sim2$avg_between23_corr_asset, sim2$regime_name, mean)
##Example pairwise total correlations
##Asset1 and Asset2 are in the same group
plot(sim2$R_asset_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset2 (same group)")
abline(v = vlines, lty = 2, col = "red")
##Asset1 and Asset12 are in groups 1 and 2
plot(sim2$R_asset_true[1, 12, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset12 (G1 vs G2)")
abline(v = vlines, lty = 2, col = "red")
##Asset1 and Asset25 are in groups 1 and 3
plot(sim2$R_asset_true[1, 25, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset25 (G1 vs G3)")
abline(v = vlines, lty = 2, col = "red")