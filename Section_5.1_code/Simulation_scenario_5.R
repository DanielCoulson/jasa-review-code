############################################################
##Simulation 5: Time-varying alpha/beta + stochastic volatility
##               with larger coefficient innovations in crises
##
##Data-generating process:
##
##Market stochastic volatility:
##     h_t = mu_h + phi_h * (h_{t-1} - mu_h) + sigma_h * eta_t
##     r_m,t | h_t ~ N(0, exp(h_t))
##
##Time-varying coefficients:
##     alpha_i,t = alpha_i_bar
##                 + phi_alpha * (alpha_i,t-1 - alpha_i_bar)
##                 + s_alpha(regime_t) * u_alpha_i,t
##
##     beta_i,t  = beta_i_bar
##                 + phi_beta  * (beta_i,t-1  - beta_i_bar)
##                 + s_beta(regime_t) * u_beta_i,t
##
##Asset returns:
##     r_i,t = alpha_i,t + beta_i,t * r_m,t + epsilon_i,t
##     epsilon_t ~ N(0, Sigma_eps,t)
##
##Regimes are still applied to TARGET ASSET correlations:
##   t = 1:250   low correlation
##   t = 251:500   crisis jump to high correlation
##   t = 501:750   moderate correlation
##   t = 751:1000  second break with a different pattern
##
##Main new feature:
##alpha_t and beta_t evolve slowly most of the time, but their
##innovations become larger in crisis regimes.
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript simulation5_run.R <seed> [outdir]")
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
  TT <- dim(Sigma_array)[3]
  R_array <- array(0, dim = c(N, N, TT))
  for (t in 1:TT) {
    R_array[,,t] <- cov2cor(Sigma_array[,,t])
  }
  R_array
}
##Build a target asset covariance matrix with the desired asset correlation
##structure, then back out the implied idiosyncratic covariance:
##Sigma_eps = Sigma_asset - sigma_m^2 * beta_t beta_t'
##
##If needed, inflate all asset standard deviations by a common factor until
##Sigma_eps is positive definite. This preserves the target asset correlations.
make_feasible_asset_cov <- function(R_asset_target,
                                    asset_sd,
                                    beta_t,
                                    sigma_m2,
                                    tol = 1e-8,
                                    inflate_step = 1.02,
                                    max_iter = 500) {
  inflate <- 1
  for (iter in 0:max_iter) {
    sd_now <- asset_sd * inflate
    D_now <- diag(sd_now)
    Sigma_asset <- D_now %*% R_asset_target %*% D_now
    Sigma_eps <- Sigma_asset - sigma_m2 * tcrossprod(beta_t)
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
  stop("Could not construct a positive definite Sigma_eps. Try lowering beta_tv_sd / crisis multipliers / stochastic-volatility parameters, or increasing asset_sd_min / asset_sd_max.")
}
##Stochastic volatility for the market
simulate_stochastic_volatility <- function(
    T,
    sv_mu = log(1.20^2),
    sv_phi = 0.985,
    sv_sigma = 0.12
) {
  h_t <- numeric(T)
  h_t[1] <- rnorm(
    1,
    mean = sv_mu,
    sd = sv_sigma / sqrt(1 - sv_phi^2)
  )
  for (t in 2:T) {
    h_t[t] <- sv_mu +
      sv_phi * (h_t[t - 1] - sv_mu) +
      rnorm(1, mean = 0, sd = sv_sigma)
  }
  sigma_m2_t <- exp(h_t)
  list(
    h_t = h_t,
    sigma_m2_t = sigma_m2_t
  )
}
##Time-varying alpha and beta with regime-dependent innovation size.
##
##Interpretation:
##- phi_alpha and phi_beta close to 1 => slow evolution / persistence
##- alpha_innov_mult[g] and beta_innov_mult[g] control how large the
##   innovations are in regime g
##- regime 2 is the crisis regime by default, so it gets bigger shocks
simulate_time_varying_coefficients <- function(
    T,
    N,
    regime,
    alpha_bar,
    beta_bar,
    phi_alpha = 0.995,
    phi_beta = 0.995,
    alpha_tv_sd = 0.004,
    beta_tv_sd = 0.012,
    alpha_innov_mult = c(1.0, 3.0, 1.2, 1.8),
    beta_innov_mult = c(1.0, 3.0, 1.2, 1.8),
    beta_floor = 0.05,
    beta_cap = 2.50
) {
  alpha_t <- matrix(0, nrow = T, ncol = N)
  beta_t <- matrix(0, nrow = T, ncol = N)
  ##Stationary initialization using base innovation scales
  alpha_init_sd <- alpha_tv_sd / sqrt(1 - phi_alpha^2)
  beta_init_sd <- beta_tv_sd  / sqrt(1 - phi_beta^2)
  alpha_t[1, ] <- alpha_bar + rnorm(N, mean = 0, sd = alpha_init_sd)
  beta_t[1, ] <- beta_bar  + rnorm(N, mean = 0, sd = beta_init_sd)
  beta_t[1, ] <- pmin(beta_cap, pmax(beta_floor, beta_t[1, ]))
  for (t in 2:T) {
    g <- regime[t]
    alpha_sd_t <- alpha_tv_sd * alpha_innov_mult[g]
    beta_sd_t <- beta_tv_sd  * beta_innov_mult[g]
    alpha_t[t, ] <- alpha_bar +
      phi_alpha * (alpha_t[t - 1, ] - alpha_bar) +
      rnorm(N, mean = 0, sd = alpha_sd_t)
    beta_prop <- beta_bar +
      phi_beta * (beta_t[t - 1, ] - beta_bar) +
      rnorm(N, mean = 0, sd = beta_sd_t)
    beta_t[t, ] <- pmin(beta_cap, pmax(beta_floor, beta_prop))
  }
  colnames(alpha_t) <- paste0("Asset", 1:N)
  colnames(beta_t) <- paste0("Asset", 1:N)
  list(
    alpha = alpha_t,
    beta = beta_t
  )
}
#############################
##Main simulator
#############################
simulate_scenario5 <- function(
    T = 1000,
    N = 30,
    seed = seed,
    ##If TRUE, regime 4 uses a block ASSET correlation structure
    regime4_block = TRUE,
    ##Regime-specific off-diagonal ASSET correlation levels
    rho_regimes = c(0.10, 0.70, 0.25, 0.50),
    ##Regime-specific multipliers for ASSET volatilities
    asset_scale_regimes = c(1.00, 1.35, 1.10, 1.25),
    ##Long-run asset-level coefficient means
    beta_min = 0.7,
    beta_max = 1.4,
    alpha_sd = 0.03,
    asset_sd_min = 1.20,
    asset_sd_max = 2.00,
    ##Time-varying coefficient parameters:
    ##persistent by default, small normal-period innovations
    phi_alpha = 0.995,
    phi_beta = 0.995,
    alpha_tv_sd = 0.004,
    beta_tv_sd = 0.012,
    ##Larger coefficient innovations in crisis/stress regimes
    ##regime 1 = calm, regime 2 = crisis, regime 3 = moderate,
    ##regime 4 = second break / stressed
    alpha_innov_mult = c(1.0, 3.0, 1.2, 1.8),
    beta_innov_mult = c(1.0, 3.0, 1.2, 1.8),
    beta_floor = 0.05,
    beta_cap = 2.50,
    ##Stochastic-volatility parameters for the market
    sv_mu = log(1.20^2),
    sv_phi = 0.985,
    sv_sigma = 0.12
) {
  set.seed(seed)
  if (T != 1000) {
    stop("This version assumes T = 1000 so that each regime has length 250.")
  }
  ##Regime labels
  regime <- rep(1:4, each = 250)
  ##Long-run asset-specific means
  alpha_bar <- rnorm(N, mean = 0, sd = alpha_sd)
  beta_bar <- runif(N, min = beta_min, max = beta_max)
  ##Time-varying coefficients with crisis-sensitive innovations
  tv_coef <- simulate_time_varying_coefficients(
    T = T,
    N = N,
    regime = regime,
    alpha_bar = alpha_bar,
    beta_bar = beta_bar,
    phi_alpha = phi_alpha,
    phi_beta = phi_beta,
    alpha_tv_sd = alpha_tv_sd,
    beta_tv_sd = beta_tv_sd,
    alpha_innov_mult = alpha_innov_mult,
    beta_innov_mult = beta_innov_mult,
    beta_floor = beta_floor,
    beta_cap = beta_cap
  )
  alpha_t <- tv_coef$alpha
  beta_t  <- tv_coef$beta
  ##Baseline asset standard deviations
  asset_sd_base <- runif(N, min = asset_sd_min, max = asset_sd_max)
  ##Optional block structure for regime 4
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
  ##Stochastic-volatility path for the market
  sv_obj <- simulate_stochastic_volatility(
    T = T,
    sv_mu = sv_mu,
    sv_phi = sv_phi,
    sv_sigma = sv_sigma
  )
  sigma_m2_t <- sv_obj$sigma_m2_t
  ##Storage
  market <- numeric(T)
  Y <- matrix(0, nrow = T, ncol = N)
  ##Truth objects
  Sigma_eps_true <- array(0, dim = c(N, N, T))
  Sigma_asset_true <- array(0, dim = c(N, N, T))
  Sigma_joint_true <- array(0, dim = c(N + 1, N + 1, T))
  ##Diagnostics / outputs
  avg_offdiag_corr_asset <- numeric(T)
  regime_asset_inflation <- numeric(T)
  for (t in 1:T) {
    g <- regime[t]
    ##Regime-specific TARGET asset correlation structure
    R_asset_t <- R_asset_regimes[[g]]
    ##Regime-specific asset volatility scale
    asset_sd_t_base <- asset_sd_base * asset_scale_regimes[g]
    ##Time-varying coefficients and stochastic market variance
    alpha_now <- alpha_t[t, ]
    beta_now <- beta_t[t, ]
    sigma_m2 <- sigma_m2_t[t]
    sigma_m <- sqrt(sigma_m2)
    ##Build feasible total asset covariance and implied idiosyncratic covariance
    tmp <- make_feasible_asset_cov(
      R_asset_target = R_asset_t,
      asset_sd = asset_sd_t_base,
      beta_t = beta_now,
      sigma_m2 = sigma_m2
    )
    Sigma_asset_t <- tmp$Sigma_asset
    Sigma_eps_t <- tmp$Sigma_eps
    regime_asset_inflation[t] <- tmp$inflate
    ##Simulate market and idiosyncratic shocks
    r_m_t <- rnorm(1, mean = 0, sd = sigma_m)
    eps_t <- as.numeric(rmvnorm(1, mean = rep(0, N), sigma = Sigma_eps_t))
    ##Asset excess returns
    y_t <- alpha_now + beta_now * r_m_t + eps_t
    ##True joint covariance of (market, assets)
    Sigma_joint_t <- matrix(0, nrow = N + 1, ncol = N + 1)
    Sigma_joint_t[1, 1] <- sigma_m2
    Sigma_joint_t[1, 2:(N + 1)] <- sigma_m2 * beta_now
    Sigma_joint_t[2:(N + 1), 1] <- sigma_m2 * beta_now
    Sigma_joint_t[2:(N + 1), 2:(N + 1)] <- Sigma_asset_t
    ##Save everything
    market[t] <- r_m_t
    Y[t, ] <- y_t
    Sigma_eps_true[,,t] <- Sigma_eps_t
    Sigma_asset_true[,,t] <- Sigma_asset_t
    Sigma_joint_true[,,t] <- Sigma_joint_t
    avg_offdiag_corr_asset[t] <- mean(R_asset_t[row(R_asset_t) != col(R_asset_t)])
  }
  ##Build design matrix
  X <- cbind(Intercept = 1, Mkt.RF = market)
  ##Correlation truths
  R_eps_true <- cov_to_corr_array(Sigma_eps_true)
  R_asset_true <- cov_to_corr_array(Sigma_asset_true)
  R_joint_true <- cov_to_corr_array(Sigma_joint_true)
  colnames(Y) <- paste0("Asset", 1:N)
  colnames(alpha_t) <- colnames(Y)
  colnames(beta_t) <- colnames(Y)
  ##Optional attributes with long-run means
  attr(alpha_t, "long_run_mean") <- alpha_bar
  attr(beta_t,  "long_run_mean") <- beta_bar
  list(
    X = X,
    Y = Y,
    market = market,
    alpha = alpha_t,
    beta = beta_t,
    regime = regime,
    sigma_m2_t = sigma_m2_t,
    regime_asset_inflation = regime_asset_inflation,
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
sim5 <- simulate_scenario5(
  T = 1000,
  N = 30,
  seed = seed,
  regime4_block = TRUE
)
##Mean true pairwise ASSET correlation by time
plot(sim5$avg_offdiag_corr_asset, type = "l", lwd = 2,
     xlab = "Time", ylab = "Mean off-diagonal correlation",
     main = "True regime-switching ASSET correlation")
abline(v = c(250, 500, 750), lty = 2, col = "red")
##Regime-wise average
tapply(sim5$avg_offdiag_corr_asset, sim5$regime, mean)
##Example: true total correlation between Asset1 and Asset2 over time
plot(sim5$R_asset_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True total correlation: Asset1 vs Asset2")
abline(v = c(250, 500, 750), lty = 2, col = "red")
##Example: stochastic market variance path
plot(sim5$sigma_m2_t, type = "l", lwd = 2,
     xlab = "Time", ylab = expression(sigma[m,t]^2),
     main = "Stochastic volatility: market variance")
##Example: time-varying beta for Asset 1
plot(sim5$beta[, 1], type = "l", lwd = 2,
     xlab = "Time", ylab = expression(beta[1*t]),
     main = "Time-varying beta for Asset 1")
##Example: time-varying alpha for Asset 1
plot(sim5$alpha[, 1], type = "l", lwd = 2,
     xlab = "Time", ylab = expression(alpha[1*t]),
     main = "Time-varying alpha for Asset 1")