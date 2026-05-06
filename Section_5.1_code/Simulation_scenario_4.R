############################################################
##Simulation 4: True DGP is DCC with heavy-tailed innovations
##
##Joint return vector:
##   w_t = (r_m,t, r_1,t, ..., r_N,t)'
##
##DGP:
##   w_t = mu + D * z_t
##
##where:
##   z_t | F_{t-1} ~ multivariate Student-t(0, R_t, nu)
##   R_t follows a DCC recursion
##
##The first component of w_t is treated as the observed market series,
##so that:
##   X = [1, market]
##   Y = asset excess returns
##
##This scenario is intentionally NOT a CAPM-style decomposition.
##It is a direct DCC-t multivariate return generator.
############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript simulation4_run.R <seed> [outdir]")
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
##Ensures a matrix is positive definite (PD)
##If not PD or numerically unstable, replaces with nearest PD matrix
force_pd_matrix <- function(S, corr = FALSE) {
  S <- (S + t(S)) / 2
  ee <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (any(!is.finite(ee)) || min(ee) <= 1e-8) {
    S <- as.matrix(Matrix::nearPD(S, corr = corr)$mat)
  }
  S
}
##Normalizes a covariance matrix into a correlation matrix
##Ensures positive definiteness before and after normalization
normalize_to_corr <- function(Q) {
  Q <- force_pd_matrix(Q, corr = FALSE)
  d <- sqrt(diag(Q))
  R <- Q / outer(d, d)
  diag(R) <- 1
  R <- force_pd_matrix(R, corr = TRUE)
  diag(R) <- 1
  R
}
#############################
##Main simulator
#############################

simulate_scenario4 <- function(
    T = 1000,
    N = 30,
    seed = seed,
    ##DCC parameters
    dcc_a = 0.03,
    dcc_b = 0.95,
    ##Degrees of freedom for multivariate Student-t innovations
    df_t = 8,
    ##Long-run correlation structure for the JOINT vector (market + assets)
    market_asset_corr = 0.35,
    asset_equicorr = 0.15,
    ##Scale parameters
    market_sd = 1.20,
    asset_sd_min = 0.80,
    asset_sd_max = 1.40,
    ##Means
    market_mean = 0.00,
    alpha_sd = 0.03,
    ##Burn-in for DCC recursion
    burnin = 300
) {
  set.seed(seed)
  if (dcc_a < 0 || dcc_b < 0 || dcc_a + dcc_b >= 1) {
    stop("Need dcc_a >= 0, dcc_b >= 0, and dcc_a + dcc_b < 1.")
  }
  if (df_t <= 2) {
    stop("df_t must be > 2 so that second moments exist.")
  }
  M <- N + 1   # dimension of joint vector = market + assets
  #############################
  ##Joint unconditional correlation target Qbar
  #############################
  ##Asset block: equicorrelation
  Qbar_assets <- make_equicorr(N, asset_equicorr)
  ##Full joint Qbar
  Qbar <- matrix(0, nrow = M, ncol = M)
  Qbar[1, 1] <- 1
  Qbar[1, 2:M] <- market_asset_corr
  Qbar[2:M, 1] <- market_asset_corr
  Qbar[2:M, 2:M] <- Qbar_assets
  Qbar <- force_pd_matrix(Qbar, corr = TRUE)
  diag(Qbar) <- 1
  #############################
  ##Means and scales
  #############################
  mu_assets <- rnorm(N, mean = 0, sd = alpha_sd)
  mu_joint <- c(market_mean, mu_assets)
  asset_sd <- runif(N, min = asset_sd_min, max = asset_sd_max)
  sd_joint <- c(market_sd, asset_sd)
  D_joint <- diag(sd_joint)
  #############################
  ##Storage
  #############################
  total_T <- T + burnin
  joint_returns <- matrix(0, nrow = total_T, ncol = M)
  z_store <- matrix(0, nrow = total_T, ncol = M)
  Q_array <- array(0, dim = c(M, M, total_T))
  R_joint_true_full <- array(0, dim = c(M, M, total_T))
  Sigma_joint_true_full <- array(0, dim = c(M, M, total_T))
  #############################
  ##Initialize DCC
  #############################
  Q_t <- Qbar
  ##Scaling so that rmvt(..., sigma = R_t, df = nu) has covariance R_t
  ##rather than (nu/(nu-2)) * R_t
  t_scale <- sqrt((df_t - 2) / df_t)
  #############################
  ##Simulate DCC-t process
  #############################
  for (t in 1:total_T) {
    ##Current conditional correlation
    R_t <- normalize_to_corr(Q_t)
    ##Current conditional covariance of joint returns
    Sigma_joint_t <- D_joint %*% R_t %*% D_joint
    ##Draw Student-t standardized innovation
    z_t <- as.numeric(mvtnorm::rmvt(1, sigma = R_t, df = df_t)) * t_scale
    ##Joint returns
    w_t <- mu_joint + as.numeric(D_joint %*% z_t)
    ##Save
    joint_returns[t, ] <- w_t
    z_store[t, ] <- z_t
    Q_array[, , t] <- Q_t
    R_joint_true_full[, , t] <- R_t
    Sigma_joint_true_full[, , t] <- Sigma_joint_t
    ##Update DCC recursion
    Q_t <- (1 - dcc_a - dcc_b) * Qbar + dcc_a * tcrossprod(z_t) + dcc_b * Q_t
    Q_t <- (Q_t + t(Q_t)) / 2
  }
  #############################
  ##Drop burn-in
  #############################
  keep <- (burnin + 1):(burnin + T)
  joint_returns <- joint_returns[keep, , drop = FALSE]
  z_store <- z_store[keep, , drop = FALSE]
  Q_true <- Q_array[, , keep, drop = FALSE]
  R_joint_true <- R_joint_true_full[, , keep, drop = FALSE]
  Sigma_joint_true <- Sigma_joint_true_full[, , keep, drop = FALSE]
  #############################
  ##Build outputs
  #############################
  market <- joint_returns[, 1]
  Y <- joint_returns[, 2:M, drop = FALSE]
  X <- cbind(Intercept = 1, Mkt.RF = market)
  ##Asset block truths
  Sigma_asset_true <- Sigma_joint_true[2:M, 2:M, , drop = FALSE]
  R_asset_true <- R_joint_true[2:M, 2:M, , drop = FALSE]
  ##No separate epsilon decomposition is defined in this DGP.
  ##We return NA arrays for compatibility, and you should use:
  ##truth_target = "asset"
  ##or
  ##truth_target = "joint"
  Sigma_eps_true <- array(NA_real_, dim = c(N, N, T))
  R_eps_true <- array(NA_real_, dim = c(N, N, T))
  ##Diagnostics
  avg_offdiag_corr_asset <- sapply(1:T, function(t) {
    R_t <- R_asset_true[, , t]
    mean(R_t[row(R_t) != col(R_t)])
  })
  avg_mkt_asset_corr <- sapply(1:T, function(t) {
    mean(R_joint_true[1, 2:M, t])
  })
  colnames(Y) <- paste0("Asset", 1:N)
  list(
    X = X,
    Y = Y,
    market = market,
    joint_returns = joint_returns,
    standardized_shocks = z_store,
    alpha = mu_assets,
    market_mean = market_mean,
    market_sd = market_sd,
    asset_sd = asset_sd,
    dcc_a = dcc_a,
    dcc_b = dcc_b,
    df_t = df_t,
    Qbar = Qbar,
    Q_true = Q_true,
    avg_offdiag_corr_asset = avg_offdiag_corr_asset,
    avg_mkt_asset_corr = avg_mkt_asset_corr,
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
sim4 <- simulate_scenario4(
  T = 1000,
  N = 30,
  seed = seed,
  dcc_a = 0.03,
  dcc_b = 0.95,
  df_t = 8
)
#############################
##Objects for your model
#############################
design_matrix <- sim4$X
asset_excess_returns <- sim4$Y
#############################
##Diagnostics / example plots
#############################
##Mean asset correlation over time
plot(sim4$avg_offdiag_corr_asset, type = "l", lwd = 2,
     xlab = "Time", ylab = "Mean off-diagonal asset correlation",
     main = "Simulation 4: True DCC-t asset correlation")
##Mean market-asset correlation over time
plot(sim4$avg_mkt_asset_corr, type = "l", lwd = 2,
     xlab = "Time", ylab = "Mean market-asset correlation",
     main = "Simulation 4: True DCC-t market-asset correlation")
##Example asset pair
plot(sim4$R_asset_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True DCC-t asset correlation: Asset1 vs Asset2")
##Example market-asset pair
plot(sim4$R_joint_true[1, 2, ], type = "l", lwd = 2,
     xlab = "Time", ylab = "Correlation",
     main = "True DCC-t correlation: Market vs Asset1")