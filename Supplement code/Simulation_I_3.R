cov_to_corr_array <- function(Sigma_array) {
  N <- dim(Sigma_array)[1]
  TT <- dim(Sigma_array)[3]
  R_array <- array(0,dim = c(N, N, TT))
  for (t in 1:TT) {
    R_array[, , t] <- cov2cor(Sigma_array[, , t])
  }
  R_array
}
simulate_logvol_ar1 <- function(T,mu,phi, sigma) {
  K <- length(mu)
  h <- matrix(0,nrow = T,ncol = K)
  h[1, ] <-rnorm(K,mean = mu,sd = sigma / sqrt(1 - phi^2))
  if (T >= 2) {
    for (t in 2:T) {
      h[t, ] <- mu +phi * (h[t - 1, ] - mu) +rnorm(K,mean = 0,sd = sigma)
    }
  }
  h
}
simulate_kowal_betas <- function(T,N,group,beta_baseline,block_window_1 = c(41L, 80L),
    block_window_2 = c(121L, 160L),block_shift = 2,rw_end = 100L,rw_endpoint_sd = 1) {
  if (length(group) != N) {
    stop("group must have length N.")
  }
  if (length(beta_baseline) != N) {
    stop("beta_baseline must have length N.")
  }
  if (!all(sort(unique(group)) == 1:3)) {
    stop("group must contain groups 1, 2, and 3.")
  }
  if (block_window_1[1] < 1 ||block_window_1[2] > T ||block_window_1[1] > block_window_1[2]) {
    stop("Invalid block_window_1.")
  }
  if (block_window_2[1] < 1 ||block_window_2[2] > T ||block_window_2[1] > block_window_2[2]) {
    stop("Invalid block_window_2.")
  }
  if (block_window_1[2] >= block_window_2[1]) {
    stop("The two block windows must not overlap.")
  }
  if (rw_end < 2 ||rw_end >= T) {
    stop("Need 2 <= rw_end < T.")
  }
  if (block_shift <= 0) {
    stop("block_shift must be positive.")
  }
  if (rw_endpoint_sd <= 0) {
    stop("rw_endpoint_sd must be positive.")
  }
  beta_t <-matrix(beta_baseline,nrow = T,ncol = N,byrow = TRUE)
  group1_assets <- which(group == 1)
  group2_assets <- which(group == 2)
  beta_t[block_window_1[1]:block_window_1[2],group2_assets] <-sweep(
      beta_t[block_window_1[1]:block_window_1[2],group2_assets,drop = FALSE],
      2,rep(block_shift, length(group2_assets)),"+")
  beta_t[block_window_2[1]:block_window_2[2],group2_assets] <-sweep(
      beta_t[block_window_2[1]:block_window_2[2],group2_assets,drop = FALSE],
      2,rep(block_shift, length(group2_assets)),"-")
  group3_assets <- which(group == 3)
  rw_step_sd <- 0.35
  for (i in group3_assets) {
    z <-rnorm(rw_end - 1L,mean = 0,sd = rw_step_sd)
    z <- z - mean(z)
    z <-z *rw_step_sd /sqrt(mean(z^2))
    beta_t[1, i] <-beta_baseline[i]
    beta_t[2:rw_end, i] <-beta_baseline[i] +cumsum(z)
    beta_t[rw_end, i] <-beta_baseline[i]
    beta_t[(rw_end + 1L):T, i] <-beta_baseline[i]
  }
  omega_t <-matrix(0,nrow = T,ncol = N)
  omega_t[2:T, ] <-beta_t[2:T, , drop = FALSE] -beta_t[1:(T - 1), , drop = FALSE]
  beta_activity <-1 * (abs(omega_t) > 0)
  beta_activity[1, ] <- 0
  activity_group <-matrix(0,nrow = T,ncol = 3)
  for (g in 1:3) {
    idx <- which(group == g)
    activity_group[, g] <-as.numeric(rowSums(beta_activity[, idx, drop = FALSE]) > 0)
  }
  beta_pattern <-c(rep("constant",length(group1_assets)),rep("piecewise_constant",
        length(group2_assets)),rep("persistent_variance_regimes",length(group3_assets)))
  block_breaks <-c(block_window_1[1],block_window_1[2] + 1L,block_window_2[1],block_window_2[2] + 1L)
  rw_reset <- NA_integer_
  colnames(beta_t) <-paste0("Asset", seq_len(N))
  colnames(omega_t) <-colnames(beta_t)
  colnames(beta_activity) <-colnames(beta_t)
  colnames(activity_group) <-paste0("Group", 1:3)
  list(beta = beta_t,omega = omega_t,activity = beta_activity,
    activity_group = activity_group,beta_pattern = beta_pattern,
    beta_baseline = beta_baseline,block_window_1 = block_window_1,
    block_window_2 = block_window_2,block_breaks = block_breaks,
    block_shift = block_shift,rw_end = rw_end,rw_reset = rw_reset,
    rw_step_sd = rw_step_sd,rw_endpoint_sd = rw_endpoint_sd)
}
simulate_scenario6 <- function(T = 200,N = 30,seed = seed,block_sizes = c(10, 10, 10),
    block_window_1 = c(41L, 80L),block_window_2 = c(121L, 160L),block_shift = 2,
    rw_end = 100L,rw_endpoint_sd = 1,alpha_sd = 0.02,beta_baseline_min = 2,beta_baseline_max = 2,
    market_sd = 1,target_rsnr = 10,number_of_residual_factors = 3,factor_logvar_phi = 0.97,
    factor_logvar_sigma = 0,factor_sd_means = c(0.40, 0.32, 0.28),idio_sd_min = 0.50,
    idio_sd_max = 0.70,idio_logvar_phi = 0.98,idio_logvar_sigma = 0) {
  set.seed(seed)
  if (T != 200) {
    stop("This version assumes T = 200.")
  }
  if (sum(block_sizes) != N) {
    stop("sum(block_sizes) must equal N.")
  }
  if (length(block_sizes) != 3) {
    stop("This version is written for exactly 3 groups.")
  }
  if (number_of_residual_factors != 3) {
    stop("This version is written for exactly 3 residual latent factors.")
  }
  if (length(factor_sd_means) !=number_of_residual_factors) {
    stop("factor_sd_means must have length number_of_residual_factors.")
  }
  group <-rep(1:3,times = block_sizes)
  alpha <-rnorm(N,mean = 0,sd = alpha_sd)
  beta_baseline <-runif(N, min = beta_baseline_min,max = beta_baseline_max)
  beta_obj <-simulate_kowal_betas(T = T,N = N,group = group,beta_baseline = beta_baseline,
      block_window_1 = block_window_1,block_window_2 = block_window_2,block_shift = block_shift,
      rw_end = rw_end,rw_endpoint_sd = rw_endpoint_sd)
  beta_t <-beta_obj$beta
  beta_omega_t <-beta_obj$omega
  beta_activity <-beta_obj$activity
  activity_group <-beta_obj$activity_group
  Lambda_eps <-matrix(0,nrow = N,ncol = number_of_residual_factors)
  Lambda_eps[, 1] <-rnorm(N,mean = 0.35,sd = 0.04)
  Lambda_eps[group == 1, 2] <-rnorm(sum(group == 1),mean = 0.28,sd = 0.03)
  Lambda_eps[group == 2, 2] <-rnorm(sum(group == 2),mean = -0.25,sd = 0.03)
  Lambda_eps[group == 3, 2] <-rnorm(sum(group == 3),mean = 0.05,sd = 0.02)
  Lambda_eps[group == 1, 3] <-rnorm(sum(group == 1),mean = 0.03,sd = 0.02)
  Lambda_eps[group == 2, 3] <-rnorm(sum(group == 2),mean = 0.20,sd = 0.03)
  Lambda_eps[group == 3, 3] <-rnorm(sum(group == 3),mean = -0.22,sd = 0.03)
  factor_logvar_mu <-log(factor_sd_means^2)
  factor_logvar <-simulate_logvol_ar1(T = T,mu = factor_logvar_mu,phi = factor_logvar_phi,
      sigma = factor_logvar_sigma)
  idio_sd_base <-runif(N,min = idio_sd_min,max = idio_sd_max)
  idio_logvar_mu <-log(idio_sd_base^2)
  idio_logvar <-simulate_logvol_ar1( T = T,mu = idio_logvar_mu,phi = idio_logvar_phi,
      sigma = idio_logvar_sigma)
  market <-rnorm(T,mean = 0,sd = market_sd)
  Y <-matrix(0,nrow = T,ncol = N)
  Sigma_eps_true <-array( 0,dim = c(N, N, T))
  Sigma_asset_true <-array(0,dim = c(N, N, T))
  Sigma_joint_true <-array( 0,dim = c(N + 1, N + 1, T))
  avg_offdiag_corr_asset <-numeric(T)
  sigma_m2_t <-rep(market_sd^2,T)
  beta_dynamic_component <-sweep(beta_t,2,beta_baseline,"-")
  signal_matrix <-sweep(beta_dynamic_component,1,market,"*")
  rsnr_reference_assets <-which(group == 2)
  signal_sd_by_asset <-apply(signal_matrix[,rsnr_reference_assets,drop = FALSE],2,sd)
  reference_signal_sd <-sqrt(mean(signal_sd_by_asset^2))
  base_residual_variance <-exp(factor_logvar) %*%t(Lambda_eps^2) +exp(idio_logvar)
  base_residual_rms_sd <-sqrt(mean(base_residual_variance[,rsnr_reference_assets,drop = FALSE]))
  target_residual_rms_sd <-reference_signal_sd /target_rsnr
  residual_scale_common <-target_residual_rms_sd /base_residual_rms_sd
  residual_scale <-rep(residual_scale_common,N)
  realized_rsnr <-reference_signal_sd /(residual_scale_common *base_residual_rms_sd)
  cat("Reference dynamic-signal SD:",reference_signal_sd,"\n")
  cat("Target residual RMS SD:",target_residual_rms_sd,"\n")
  cat("Residual scale:",residual_scale_common,"\n")
  cat("Realized pooled RSNR:",realized_rsnr,"\n")
  stopifnot(abs(realized_rsnr -target_rsnr) < 1e-10)
  for (t in 1:T) {
    factor_var_t <-exp(factor_logvar[t, ])
    idio_var_t <-exp(idio_logvar[t, ])
    Sigma_eps_t <-Lambda_eps %*%diag(factor_var_t) %*%t(Lambda_eps) +diag(idio_var_t)
    Sigma_eps_t <-tcrossprod(residual_scale) *Sigma_eps_t
    Sigma_eps_t <-(Sigma_eps_t +t(Sigma_eps_t)) / 2
    beta_now <-beta_t[t, ]
    Sigma_asset_t <-market_sd^2 *tcrossprod(beta_now) +Sigma_eps_t
    Sigma_asset_t <-(Sigma_asset_t +t(Sigma_asset_t)) / 2
    r_m_t <- market[t]
    f_t <-rnorm(number_of_residual_factors,mean = 0,sd = sqrt(factor_var_t))
    u_t <-rnorm(N,mean = 0,sd = sqrt(idio_var_t))
    eps_t <-residual_scale *(as.numeric(Lambda_eps %*% f_t) +u_t)
    y_t <-alpha +beta_now * r_m_t +eps_t
    Sigma_joint_t <-matrix(0,nrow = N + 1,ncol = N + 1)
    Sigma_joint_t[1, 1] <-market_sd^2
    Sigma_joint_t[1,2:(N + 1)] <-market_sd^2 *beta_now
    Sigma_joint_t[2:(N + 1),1] <-market_sd^2 *beta_now
    Sigma_joint_t[2:(N + 1),2:(N + 1)] <-Sigma_asset_t
    market[t] <-r_m_t
    Y[t, ] <-y_t
    Sigma_eps_true[, , t] <-Sigma_eps_t
    Sigma_asset_true[, , t] <-Sigma_asset_t
    Sigma_joint_true[, , t] <-Sigma_joint_t
    R_asset_t <-cov2cor(Sigma_asset_t)
    avg_offdiag_corr_asset[t] <-mean(R_asset_t[upper.tri(R_asset_t)])
  }
  X <-cbind(Intercept = 1,Mkt.RF = market)
  R_eps_true <-cov_to_corr_array(Sigma_eps_true)
  R_asset_true <-cov_to_corr_array(Sigma_asset_true)
  R_joint_true <-cov_to_corr_array(Sigma_joint_true)
  colnames(Y) <-paste0("Asset",seq_len(N))
  colnames(beta_t) <-colnames(Y)
  colnames(beta_omega_t) <-colnames(Y)
  colnames(beta_activity) <-colnames(Y)
  rownames(Lambda_eps) <-colnames(Y)
  colnames(Lambda_eps) <-paste0("ResidualFactor",seq_len(number_of_residual_factors))
  list(X = X,Y = Y,market = market,alpha = alpha,beta = beta_t,beta_omega = beta_omega_t,
    beta_baseline = beta_obj$beta_baseline,beta_pattern = beta_obj$beta_pattern,
    beta_activity = beta_activity,activity_group = activity_group,group = group,
    block_window_1 = beta_obj$block_window_1,block_window_2 = beta_obj$block_window_2,
    block_breaks = beta_obj$block_breaks,block_shift = beta_obj$block_shift,
    rw_end = beta_obj$rw_end,rw_reset = beta_obj$rw_reset,
    rw_step_sd = beta_obj$rw_step_sd,rw_endpoint_sd = beta_obj$rw_endpoint_sd,
    sigma_m2_t = sigma_m2_t,target_rsnr =target_rsnr,realized_rsnr =realized_rsnr,
    residual_scale =residual_scale,Lambda_eps = Lambda_eps,factor_logvar = factor_logvar,
    idio_logvar = idio_logvar,avg_offdiag_corr_asset =avg_offdiag_corr_asset,
    Sigma_eps_true =Sigma_eps_true,R_eps_true =R_eps_true,Sigma_asset_true =Sigma_asset_true,
    R_asset_true = R_asset_true,Sigma_joint_true =Sigma_joint_true,R_joint_true =R_joint_true)
}
seed <- 1
sim6 <-simulate_scenario6(T = 200,N = 30,seed = seed)
design_matrix <- sim6$X
asset_excess_returns <- sim6$Y