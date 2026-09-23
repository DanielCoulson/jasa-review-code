#############################
##Fit models
#############################
##Inputs 
design_matrix <- sim$X
asset_excess_returns <- sim$Y

library(Matrix)
library(mvtnorm)
library(HDInterval)
library(factorstochvol)
library(rugarch)
library(rmgarch)
library(loggle)
library(ggplot2)
library(bayesianVARs)

computeDIC <- FALSE

#Ensures a matrix is symmetric positive definite, using nearPD if needed.
force_pd <- function(S, corr = FALSE) {
  S <- (S + t(S)) / 2
  ee <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (any(!is.finite(ee)) || min(ee) <= 1e-8) {
    S <- as.matrix(Matrix::nearPD(S, corr = corr)$mat)
  }
  S
}

#Converts a covariance matrix to a valid correlation matrix, repairing positive definiteness if necessary.
safe_cov2cor <- function(S) {
  S <- force_pd(S, corr = FALSE)
  R <- cov2cor(S)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R <- force_pd(R, corr = TRUE)
  diag(R) <- 1
  R
}
#Computes summary of a single correlation matrix.
score <- function(R) {
  p <- ncol(R)
  R <- (R + t(R)) / 2
  ch <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(ch)) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    ch <- chol(R)
  }
  ld <- 2 * sum(log(diag(ch)))
  1 - exp(ld / p)
}

#Computes row-wise highest density intervals for posterior or bootstrap draws.
row_hdi <- function(draw_mat, credMass = 0.95) {
  out <- t(apply(draw_mat, 1, function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2) return(c(NA_real_, NA_real_))
    HDInterval::hdi(z, credMass = credMass)
  }))
  colnames(out) <- c("lower", "upper")
  out
}

#Computes row-wise percentile intervals for posterior or bootstrap draws.
row_percentile_interval <- function(draw_mat, probs = c(0.025, 0.975)) {
  out <- t(apply(draw_mat, 1, function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2) return(c(NA_real_, NA_real_))
    as.numeric(quantile(z, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
  }))
  colnames(out) <- c("lower", "upper")
  out
}

#Extracts the true time-varying score from the simulated asset, residual, or joint correlation matrices.
extract_true_scores <- function(sim, truth_target = c("asset", "eps", "joint")) {
  truth_target <- match.arg(truth_target)
  R_true <- switch(
    truth_target,
    asset = sim$R_asset_true,
    eps   = sim$R_eps_true,
    joint = sim$R_joint_true
  )
  TT <- dim(R_true)[3]
  sapply(seq_len(TT), function(t) score(R_true[, , t]))
}

#Computes a weighted covariance matrix for observations in X using weights w.
weighted_cov <- function(X, w) {
  mu <- colSums(X * w)
  Xc <- sweep(X, 2, mu, "-")
  crossprod(sqrt(w) * Xc)
}

#Generates moving-block bootstrap indices for a time series of length TT.
moving_block_bootstrap_idx <- function(TT, block_length) {
  n_blocks <- ceiling(TT / block_length)
  starts <- sample.int(TT - block_length + 1, n_blocks, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + block_length - 1)))
  idx[seq_len(TT)]
}

#Fits a scoring model and uses moving-block bootstrap resamples to estimate score intervals.
bootstrap_score_intervals <- function(
    Y,
    fit_fun,
    B = 200,
    block_length = 20,
    seed = seed,
    probs = c(0.025, 0.975)
) {
  if (!is.null(seed)) set.seed(seed)
  
  base_fit <- fit_fun(Y)
  TT <- length(base_fit$mean)
  boot_scores <- matrix(NA_real_, nrow = TT, ncol = B)
  
  for (b in seq_len(B)) {
    idx <- moving_block_bootstrap_idx(nrow(Y), block_length)
    Yb <- Y[idx, , drop = FALSE]
    tmp <- try(fit_fun(Yb), silent = TRUE)
    if (!inherits(tmp, "try-error") && length(tmp$mean) == TT) {
      boot_scores[, b] <- tmp$mean
    }
  }
  
  ints <- row_percentile_interval(boot_scores, probs = probs)
  base_fit$lower <- ints[, "lower"]
  base_fit$upper <- ints[, "upper"]
  base_fit$boot_scores <- boot_scores
  base_fit
}

#Fits a multivariate factor stochastic volatility model and computes time-varying correlation scores.
fit_mfsv_scores <- function(
    Y,
    factors = 3,
    draws = 12000,
    thin = 4,
    burnin = 1500,
    quiet = TRUE,
    interval_k = 2
) {
  fit <- factorstochvol::fsvsample(
    Y,
    factors = factors,
    draws = draws,
    thin = thin,
    burnin = burnin,
    keeptime = "all",
    quiet = quiet
  )
  
  TT <- nrow(Y)
  
  score_mean  <- rep(NA_real_, TT)
  score_lower <- rep(NA_real_, TT)
  score_upper <- rep(NA_real_, TT)
  
  for (t in seq_len(TT)) {
    R_mean <- factorstochvol::runningcormat(fit, i = t, statistic = "mean")
    R_sd   <- factorstochvol::runningcormat(fit, i = t, statistic = "sd")
    
    score_mean[t]  <- score(R_mean)
    score_lower[t] <- score(R_mean - interval_k * R_sd)
    score_upper[t] <- score(R_mean + interval_k * R_sd)
  }
  
  list(
    mean  = score_mean,
    lower = score_lower,
    upper = score_upper,
    fit   = fit
  )
}

#Fits the DSP-MFSV-CAPM and summarizes posterior score draws with means and HDI intervals.
fit_dsp_mfsv_scores <- function(
    Y,
    X,
    number_of_latent_factors = 3,
    nsave = 3000,
    nburn = 1500,
    nskip = 4,
    credMass = 0.95
) {
  fit <- DSP_MFSV(
    y = Y,
    X = X,
    number_of_latent_factors = number_of_latent_factors,
    nsave = nsave,
    nburn = nburn,
    nskip = nskip,
    mcmc_params = c("scores")
  )
  score_draws <- fit$scores
  hdi_mat <- row_hdi(score_draws, credMass = credMass)
  list(
    mean = rowMeans(score_draws, na.rm = TRUE),
    lower = hdi_mat[, "lower"],
    upper = hdi_mat[, "upper"],
    score_draws = score_draws,
    fit = fit
  )
}

#Computes rolling exponentially weighted covariance estimates and converts them to correlation scores.
fit_ew_roll_scores <- function(Y, window = 60, lambda = 0.94) {
  TT <- nrow(Y)
  score_hat <- rep(NA_real_, TT)
  w <- lambda^((window - 1):0)
  w <- w / sum(w)
  
  for (t in window:TT) {
    Xwin <- Y[(t - window + 1):t, , drop = FALSE]
    S <- weighted_cov(Xwin, w)
    R <- safe_cov2cor(S)
    score_hat[t] <- score(R)
  }
  
  list(mean = score_hat, lower = rep(NA_real_, TT), upper = rep(NA_real_, TT))
}

#Computes rolling Ledoit-Wolf shrinkage covariance estimates and converts them to correlation scores.
fit_lw_roll_scores <- function(Y, window = 60) {
  if (!requireNamespace("cvCovEst", quietly = TRUE)) {
    stop("Package 'cvCovEst' is required for the rolling Ledoit-Wolf estimator.")
  }
  
  TT <- nrow(Y)
  score_hat <- rep(NA_real_, TT)
  
  for (t in window:TT) {
    Xwin <- Y[(t - window + 1):t, , drop = FALSE]
    S <- cvCovEst::linearShrinkLWEst(Xwin)
    R <- safe_cov2cor(S)
    score_hat[t] <- score(R)
  }
  
  list(mean = score_hat, lower = rep(NA_real_, TT), upper = rep(NA_real_, TT))
}

#Applies hard thresholding to small off-diagonal correlations and repairs the matrix to be positive definite.
hard_threshold_corr <- function(R, tau = 0.10) {
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  offdiag <- row(R) != col(R)
  R[offdiag & abs(R) < tau] <- 0
  
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

#Computes rolling hard-thresholded correlation matrices and their associated scores.
fit_ht_roll_scores <- function(Y, window = 60, tau = 0.10) {
  TT <- nrow(Y)
  N  <- ncol(Y)
  
  score_hat <- rep(NA_real_, TT)
  R_arr <- array(NA_real_, dim = c(N, N, TT))
  
  for (t in window:TT) {
    Xwin <- Y[(t - window + 1):t, , drop = FALSE]
    S <- cov(Xwin)
    R <- safe_cov2cor(S)
    R_ht <- hard_threshold_corr(R, tau = tau)
    
    R_arr[, , t] <- R_ht
    score_hat[t] <- score(R_ht)
  }
  
  list(
    mean = score_hat,
    lower = rep(NA_real_, TT),
    upper = rep(NA_real_, TT),
    R = R_arr
  )
}

#Fits the hard-thresholded rolling estimator and adds bootstrap-based score intervals.
fit_ht_roll_scores_boot <- function(
    Y,
    window = 60,
    tau = 0.10,
    B = 200,
    block_length = 20,
    seed = seed
) {
  bootstrap_score_intervals(
    Y = Y,
    fit_fun = function(Z) fit_ht_roll_scores(Z, window = window, tau = tau),
    B = B,
    block_length = block_length,
    seed = seed
  )
}

fit_loggle_path <- function(Y_in, h, d, lambda) {
  X_loggle <- t(Y_in)
  fit <- loggle(X_loggle,
    pos = seq_len(nrow(Y_in)),
    h = h,d = d,lambda = lambda,
    fit.type = "pseudo",refit = TRUE,
    num.thread = 1,print.detail = FALSE)
  stopifnot(length(fit$Omega) == nrow(Y_in))
  out <- rep(NA_real_,nrow(Y_in))
  for (tt in seq_len(nrow(Y_in))) {
    Omega <- as.matrix(fit$Omega[[tt]])
    Omega <- (Omega + t(Omega)) / 2
    Sigma <- solve(Omega)
    Sigma <- (Sigma + t(Sigma)) / 2
    R <- cov2cor(Sigma)
    out[tt] <- score(R)
  }
  stopifnot(all(is.finite(out)))
  out
}

fit_loggle_scores <- function(Y, seed = seed){
  TT <- nrow(Y)
  X_loggle <- t(Y)
  pos_cv <- round(seq(0.1,0.9,length.out = 9) *(TT - 1) +1)
  if (!is.null(seed)) {
    set.seed(100000 + seed)
  }
  cv <- loggle.cv(X_loggle,
    pos = pos_cv,h.list = c(0.20, 0.25),
    d.list = c(0, 0.05, 0.15, 1),lambda.list = c(0.20, 0.25),
    cv.fold = 3,fit.type = "pseudo",
    return.select = TRUE,select.type = "all_fixed",
    cv.vote.thres = 0.8,num.thread = 1,
    print.detail = FALSE)
  sel <- cv$cv.select.result
  h.opt <- as.numeric(sel$h.opt[1])
  d.opt <- as.numeric(sel$d.opt[1])
  lambda.opt <- as.numeric(sel$lambda.opt[1])
  loggle_score <- fit_loggle_path(
    Y,h = h.opt,d = d.opt,
    lambda = lambda.opt)
  list(
    mean = loggle_score, 
    lower = rep(NA_real_, TT),
    upper = rep(NA_real_, TT),
    h = h.opt, 
    d = d.opt, 
    lambda = lambda.opt
  )
}
fit_loggle_scores_boot <- function(Y, B=200, block_length = 20, seed = seed,
                                   probs = c(0.025,0.975)){
  TT <- nrow(Y)
  base_fit <- fit_loggle_scores(Y = Y,seed = seed)
  h.opt <- base_fit$h
  d.opt <- base_fit$d
  lambda.opt <- base_fit$lambda
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_indices <- vector("list",B)
  for (b in seq_len(B)) {
    boot_indices[[b]] <- moving_block_bootstrap_idx(
        TT = TT,block_length = block_length)
  }
  boot_scores <- matrix(NA_real_,nrow = TT,ncol = B)
  for (b in seq_len(B)) {
    idx <- boot_indices[[b]]
    Yb <- Y[idx,,drop = FALSE]
    boot_scores[, b] <- fit_loggle_path(Yb,
      h = h.opt,d = d.opt,lambda = lambda.opt)
  }
  stopifnot(ncol(boot_scores) == B,all(is.finite(boot_scores)))
  intervals <- row_percentile_interval(boot_scores,probs = probs)
  base_fit$lower <- intervals[, "lower"]
  base_fit$upper <- intervals[, "upper"]
  base_fit$boot_scores <- boot_scores
  base_fit
}

#Builds univariate GARCH specifications for each asset to be used in a DCC model.
build_uspec <- function(N, margin_dist = c("norm", "std")) {
  margin_dist <- match.arg(margin_dist)
  u_spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = margin_dist
  )
  rugarch::multispec(replicate(N, u_spec, simplify = FALSE))
}

Q_to_R <- function(Q) {
  Q <- (Q + t(Q)) / 2
  d <- sqrt(pmax(diag(Q), 1e-12))
  R <- Q / tcrossprod(d)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  if (is.null(tryCatch(chol(R), error = function(e) NULL))) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    R <- (R + t(R)) / 2
    diag(R) <- 1
  }
  R
}

#Fits a DCC or asymmetric DCC GARCH model and computes scores from the fitted dynamic correlations.
fit_dcc <- function(
    Y,
    model = c("DCC", "aDCC"),
    distribution = c("mvnorm", "mvt")
) {
  
  model <- match.arg(model)
  distribution <- match.arg(distribution)
  N <- ncol(Y)
  
  margin_dist <- if (distribution == "mvt") {
    "std"
  } else {
    "norm"
  }
  
  uspec <- build_uspec(
    N = N,
    margin_dist = margin_dist
  )
  
  spec <- rmgarch::dccspec(
    uspec = uspec,
    dccOrder = c(1, 1),
    model = model,
    distribution = distribution
  )
  
  rmgarch::dccfit(
    spec = spec,
    data = Y,
    fit.control = list(
      eval.se = FALSE
    )
  )
}

############################################################
## Extract parameters from fitted rmgarch model
############################################################

extract_dcc_parameters <- function(
    fit,
    model = c("DCC", "aDCC")
) {
  
  model <- match.arg(model)
  
  N <- ncol(
    fit@model$modeldata$data
  )
  
  mpars <- fit@model$mpars
  
  get_margin_par <- function(name) {
    
    if (!(name %in% rownames(mpars))) {
      stop("Cannot find parameter: ", name)
    }
    
    as.numeric(
      mpars[
        name,
        seq_len(N)
      ]
    )
  }
  
  dcc_pars <- coef(
    fit,
    type = "dcc"
  )
  
  list(
    
    mu = get_margin_par("mu"),
    
    omega = get_margin_par("omega"),
    
    alpha = get_margin_par("alpha1"),
    
    beta = get_margin_par("beta1"),
    
    dcca = unname(
      dcc_pars["dcca1"]
    ),
    
    dccb = unname(
      dcc_pars["dccb1"]
    ),
    
    dccg = if (model == "aDCC") {
      unname(dcc_pars["dccg1"])
    } else {
      0
    },
    
    Qbar = as.matrix(
      fit@mfit$Qbar
    ),
    
    Nbar = if (model == "aDCC") {
      as.matrix(fit@mfit$Nbar)
    } else {
      matrix(0, N, N)
    }
  )
}

extract_bootstrap_innovations <- function(fit) {
  
  Z <- as.matrix(fit@mfit$stdresid)
  R_arr <- rmgarch::rcor(fit)
  TT <- nrow(Z)
  N <- ncol(Z)
  A <- matrix(NA_real_,nrow = TT,ncol = N)
  for (t in seq_len(TT)) {
    R_t <- R_arr[, , t]
    R_t <- (R_t + t(R_t)) / 2
    diag(R_t) <- 1
    ch <- tryCatch(chol(R_t),error = function(e) NULL)
    if (is.null(ch)) {
      R_t <- as.matrix(Matrix::nearPD(R_t,corr = TRUE)$mat)
      ch <- chol(R_t)
    }
    L <- t(ch)
    A[t, ] <- as.numeric(forwardsolve(L,Z[t, ]))
  }
  good <- apply(A,1,function(x) all(is.finite(x)))
  A[good, , drop = FALSE]
}

simulate_dcc_bootstrap <- function(fit,
    innovations,model = c("DCC", "aDCC")) {
  model <- match.arg(model)
  pars <- extract_dcc_parameters(fit,model)
  original_data <- as.matrix(fit@model$modeldata$data)
  TT <- nrow(original_data)
  N <- ncol(original_data)
  ind <- sample.int(nrow(innovations),TT,replace = TRUE)
  A_star <- innovations[ind,,drop = FALSE]
  persistence <-pars$alpha +pars$beta
  if (any(persistence >= 1)) {
    stop("Non-stationary marginal GARCH estimate.")
  }
  sigma2 <-pars$omega /(1 - persistence)
  sigma2 <- pmax(sigma2,1e-12)
  Q <- pars$Qbar
  Q_constant <-(1 - pars$dcca - pars$dccb) *pars$Qbar
  if (model == "aDCC") {
    Q_constant <-Q_constant -pars$dccg *pars$Nbar
  }
  Y_star <- matrix(NA_real_,nrow = TT,ncol = N)
  colnames(Y_star) <-colnames(original_data)
  
  for (t in seq_len(TT)) {
    R <- Q_to_R(Q)
    L <- t(chol(R))
    z <- as.numeric(L %*% A_star[t, ])
    e <- sqrt(sigma2) * z
    Y_star[t, ] <-pars$mu + e
    if (t < TT) {
      sigma2_next <-pars$omega +pars$alpha * e^2 +pars$beta * sigma2
      sigma2_next <- pmax(sigma2_next,1e-12)
      Q_next <-Q_constant +pars$dcca *tcrossprod(z) +pars$dccb *Q
      if (model == "aDCC") {
        z_negative <- ifelse(z < 0,z,0)
        Q_next <-Q_next +pars$dccg *tcrossprod(z_negative)
      }
      sigma2 <- sigma2_next
      Q <- Q_next
    }
  }
  Y_star
}

filter_original_score <- function(Y,
    bootstrap_fit, model = c("DCC", "aDCC")) {
  model <- match.arg(model)
  Y <- as.matrix(Y)
  TT <- nrow(Y)
  pars <- extract_dcc_parameters(bootstrap_fit,model)
  persistence <-pars$alpha +pars$beta
  if (any(persistence >= 1)) {
    stop("Non-stationary bootstrap GARCH estimate.")
  }
  sigma2 <-pars$omega /(1 - persistence)
  sigma2 <- pmax(sigma2,1e-12)
  Q <- pars$Qbar
  Q_constant <-(1 - pars$dcca - pars$dccb) *pars$Qbar
  if (model == "aDCC") {
    Q_constant <-Q_constant -pars$dccg *pars$Nbar
  }
  score_path <- numeric(TT)
  for (t in seq_len(TT)) {
    R <- Q_to_R(Q)
    score_path[t] <-score(R)
    e <-Y[t, ] -pars$mu
    z <-e /sqrt(sigma2)
    if (t < TT) {
      sigma2_next <-pars$omega +pars$alpha * e^2 +pars$beta * sigma2
      sigma2_next <- pmax(sigma2_next,1e-12)
      Q_next <-Q_constant +pars$dcca *tcrossprod(z) +pars$dccb *Q
      if (model == "aDCC") {
        z_negative <- ifelse(z < 0,z,0)
        Q_next <-Q_next +pars$dccg *tcrossprod(z_negative)
      }
      sigma2 <- sigma2_next
      Q <- Q_next
    }
  }
  score_path
}

fit_dcc_with_uncertainty <- function(Y,
    model = c("DCC", "aDCC"),distribution = c("mvnorm", "mvt"),
    B = 200,seed = NULL,probs = c(0.025, 0.975)) {
  model <- match.arg(model)
  distribution <- match.arg(distribution)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  fit <- fit_dcc( Y,
    model = model,distribution = distribution)
  R_hat <- rmgarch::rcor(fit)
  TT <- dim(R_hat)[3]
  score_hat <- vapply(
    seq_len(TT),
    function(t) score(R_hat[, , t]),
    numeric(1)
  )
  innovations <-extract_bootstrap_innovations(fit)
  boot_scores <- matrix(NA_real_,nrow = TT,ncol = B)
  for (b in seq_len(B)) {
    if (b == 1 ||b %% 20 == 0 ||b == B) {
      message(model," / ",distribution,": bootstrap ",b,"/",B)
    }
    Y_star <- try(simulate_dcc_bootstrap(fit = fit,innovations = innovations,
        model = model),silent = TRUE )
    if (inherits(Y_star, "try-error")) {
      next
    }
    fit_star <- try(fit_dcc(Y_star,model = model,distribution = distribution),
      silent = TRUE)
    if (inherits(fit_star, "try-error")) {
      next
    }
    score_star <- try(filter_original_score(Y = Y,
        bootstrap_fit = fit_star,model = model),
      silent = TRUE)
    if ( !inherits(score_star, "try-error") &&
      length(score_star) == TT &&
      all(is.finite(score_star))
    ) {
      boot_scores[, b] <-
        score_star
    }
  }
  successful <-colSums(is.finite(boot_scores)) == TT
  boot_scores <-boot_scores[,successful,drop = FALSE]
  if (ncol(boot_scores) < 2) {
    stop("Fewer than two bootstrap replications succeeded.")
  }
  lower <- apply(boot_scores,1,quantile,probs = probs[1],names = FALSE,type = 8)
  upper <- apply(boot_scores,1,quantile,probs = probs[2],names = FALSE,type = 8)
  list(
    mean = score_hat,
    lower = lower,
    upper = upper, 
    R = R_hat, 
    fit = fit, 
    boot_scores = boot_scores,
    n_boot_success = ncol(boot_scores)
  )
}

#Fits a Bayesian VAR with Cholesky stochastic volatility and computes time-varying correlation scores.
fit_chol_sv_scores <- function(
    Y,
    lags = 1L,
    draws = 1500L,
    burnin = 1500L,
    thin = 1L,
    cholesky_U_prior = "HS",
    center_data = TRUE,
    quiet = TRUE,
    interval_k = 2 
) {
  if (!requireNamespace("bayesianVARs", quietly = TRUE)) {
    stop("Package 'bayesianVARs' is required for the Cholesky SV benchmark.")
  }
  
  Y_in <- if (center_data) scale(Y, center = TRUE, scale = FALSE) else Y
  Y_in <- as.matrix(Y_in)
  
  TT_full <- nrow(Y_in)
  N <- ncol(Y_in)
  
  prior_sigma <- bayesianVARs::specify_prior_sigma(
    data = Y_in,
    type = "cholesky",
    cholesky_U_prior = cholesky_U_prior,
    cholesky_heteroscedastic = TRUE,
    quiet = TRUE
  )
  
  fit <- bayesianVARs::bvar(
    data = Y_in,
    lags = lags,
    draws = draws,
    burnin = burnin,
    thin = thin,
    prior_intercept = FALSE,
    prior_sigma = prior_sigma,
    sv_keep = "all",
    quiet = quiet
  )
  
  TT_fit <- nrow(fit$logvar)
  
  #Repairs a candidate correlation matrix by symmetrizing it, bounding entries, and enforcing positive definiteness.
  fix_corr_candidate <- function(R) {
    R <- (R + t(R)) / 2
    R[R >  1] <-  1
    R[R < -1] <- -1
    diag(R) <- 1
    R <- force_pd(R, corr = TRUE)
    diag(R) <- 1
    R
  }
  
  # Converts a covariance matrix to a correlation matrix with lightweight numerical stabilization.
  safe_cov2cor_light <- function(S) {
    S <- (S + t(S)) / 2
    ee <- tryCatch(
      eigen(S, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    
    if (any(!is.finite(ee))) {
      S <- S + diag(1e-6, nrow(S))
    } else if (min(ee) <= 1e-10) {
      S <- S + diag(abs(min(ee)) + 1e-6, nrow(S))
    }
    
    R <- cov2cor(S)
    R <- (R + t(R)) / 2
    diag(R) <- 1
    R
  }
  
  #Extracts posterior covariance draws at one time point and summarizes them as mean/lower/upper correlation matrices.
  get_corr_summary_at_t <- function(fit, t) {
    Sigma_draws_t <- stats::vcov(fit, t = t)
    dd <- dim(Sigma_draws_t)
    
    if (length(dd) == 4L) {
      ndraws <- dd[4]
      R_draws <- array(NA_real_, dim = c(N, N, ndraws))
      for (s in seq_len(ndraws)) {
        R_draws[, , s] <- safe_cov2cor_light(Sigma_draws_t[1, , , s, drop = TRUE])
      }
    } else if (length(dd) == 3L) {
      ndraws <- dd[3]
      R_draws <- array(NA_real_, dim = c(N, N, ndraws))
      for (s in seq_len(ndraws)) {
        R_draws[, , s] <- safe_cov2cor_light(Sigma_draws_t[, , s, drop = TRUE])
      }
    } else {
      stop("Unexpected dimensions returned by stats::vcov(fit, t = t).")
    }
    
    R_mean <- apply(R_draws, c(1, 2), mean, na.rm = TRUE)
    R_sd   <- apply(R_draws, c(1, 2), sd,   na.rm = TRUE)
    
    R_mean  <- fix_corr_candidate(R_mean)
    R_lower <- fix_corr_candidate(R_mean - interval_k * R_sd)
    R_upper <- fix_corr_candidate(R_mean + interval_k * R_sd)
    
    list(
      mean_mat  = R_mean,
      lower_mat = R_lower,
      upper_mat = R_upper
    )
  }
  
  score_mean  <- rep(NA_real_, TT_full)
  score_lower <- rep(NA_real_, TT_full)
  score_upper <- rep(NA_real_, TT_full)
  
  idx <- (lags + 1):TT_full
  
  for (t in seq_len(TT_fit)) {
    tmp <- get_corr_summary_at_t(fit, t)
    
    score_mean[idx[t]]  <- score(tmp$mean_mat)
    score_lower[idx[t]] <- score(tmp$lower_mat)
    score_upper[idx[t]] <- score(tmp$upper_mat)
  }
  
  list(
    mean = score_mean,
    lower = score_lower,
    upper = score_upper,
    fit = fit
  )
}

#Fits the exponentially weighted rolling estimator and adds bootstrap-based score intervals.
fit_ew_roll_scores_boot <- function(
    Y,
    window = 60,
    lambda = 0.94,
    B = 200,
    block_length = 20,
    seed = seed
) {
  bootstrap_score_intervals(
    Y = Y,
    fit_fun = function(Z) fit_ew_roll_scores(Z, window = window, lambda = lambda),
    B = B,
    block_length = block_length,
    seed = seed
  )
}

#Fits the rolling Ledoit-Wolf estimator and adds bootstrap-based score intervals.
fit_lw_roll_scores_boot <- function(
    Y,
    window = 60,
    B = 200,
    block_length = 20,
    seed = seed
) {
  bootstrap_score_intervals(
    Y = Y,
    fit_fun = function(Z) fit_lw_roll_scores(Z, window = window),
    B = B,
    block_length = block_length,
    seed = seed
  )
}


#Runs all benchmark and Bayesian score models, stores true scores, and returns fitted results by method.
run_all_score_models <- function(
    sim,
    truth_target = c("asset", "eps", "joint"),
    mfsv_factors = 3,
    mfsv_draws = 12000,
    mfsv_thin = 4,
    mfsv_burnin = 1500,
    dsp_factors = 3,
    dsp_nsave = 3000,
    dsp_nburn = 1500,
    dsp_nskip = 4,
    roll_window = 60,
    ew_lambda = 0.94,
    B_nonbayes = 200,
    block_length = 20,
    bootstrap_seed = seed,
    credMass = 0.95, 
    ht_tau = 0.10, 
    cholsv_lags = 1,
    cholsv_draws = 12000,
    cholsv_burnin = 1500,
    cholsv_thin = 4,
    cholsv_u_prior = "HS"
) {
  truth_target <- match.arg(truth_target)
  
  Y <- sim$Y
  X <- sim$X
  TT <- nrow(Y)
  
  out <- list()
  out$truth <- data.frame(
    t = seq_len(TT),
    true_score = extract_true_scores(sim, truth_target = truth_target)
  )
  
  out$methods <- list()
  
  message("Fitting MFSV...")
  out$methods$mfsv <- fit_mfsv_scores(
    Y = Y,
    factors = mfsv_factors,
    draws = mfsv_draws,
    thin = mfsv_thin,
    burnin = mfsv_burnin
  )
  
  message("Fitting Cholesky SV...")
  out$methods$chol_sv <- fit_chol_sv_scores(
    Y = Y,
    lags = cholsv_lags,
    draws = cholsv_draws,
    burnin = cholsv_burnin,
    thin = cholsv_thin,
    cholesky_U_prior = cholsv_u_prior
  )
  
  message("Fitting DSP_MFSV...")
  out$methods$dsp_mfsv <- fit_dsp_mfsv_scores(
    Y = Y,
    X = X,
    number_of_latent_factors = dsp_factors,
    nsave = dsp_nsave,
    nburn = dsp_nburn,
    nskip = dsp_nskip,
    credMass = credMass
  )
  message("Fitting exponentially weighted rolling estimator with bootstrap...")
  out$methods$ew_roll <- fit_ew_roll_scores_boot(
    Y = Y,
    window = roll_window,
    lambda = ew_lambda,
    B = B_nonbayes,
    block_length = block_length,
    seed = bootstrap_seed
  )
  
  message("Fitting rolling Ledoit-Wolf with bootstrap...")
  out$methods$roll_lw <- fit_lw_roll_scores_boot(
    Y = Y,
    window = roll_window,
    B = B_nonbayes,
    block_length = block_length,
    seed = bootstrap_seed
  )
  
  message("Fitting hard-thresholded rolling estimator with bootstrap...")
  out$methods$ht_roll <- fit_ht_roll_scores_boot(
    Y = Y,
    window = roll_window,
    tau = ht_tau,
    B = B_nonbayes,
    block_length = block_length,
    seed = bootstrap_seed
  )
  message("Fitting LOGGLE with bootstrap...")
  out$methods$loggle <- fit_loggle_scores_boot(
    Y = Y,
    B = B_nonbayes,
    block_length = block_length,
    seed = bootstrap_seed
  )
  message("Fitting Gaussian DCC with Fresoli-Ruiz bootstrap uncertainty...")
  out$methods$dcc_norm <- fit_dcc_with_uncertainty(
    Y = Y,
    model = "DCC",
    distribution = "mvnorm",
    B = B_nonbayes,
    seed = bootstrap_seed,
    probs = c((1 - credMass) / 2, 1 - (1 - credMass) / 2)
  )
  
  message("Fitting t-DCC with Fresoli-Ruiz bootstrap uncertainty...")
  out$methods$dcc_t <- fit_dcc_with_uncertainty(
    Y = Y,
    model = "DCC",
    distribution = "mvt",
    B = B_nonbayes,
    seed = bootstrap_seed + 1L,
    probs = c((1 - credMass) / 2, 1 - (1 - credMass) / 2)
  )
  
  message("Fitting asymmetric DCC with Fresoli-Ruiz bootstrap uncertainty...")
  out$methods$adcc <- fit_dcc_with_uncertainty(
    Y = Y,
    model = "aDCC",
    distribution = "mvnorm",
    B = B_nonbayes,
    seed = bootstrap_seed + 2L,
    probs = c((1 - credMass) / 2, 1 - (1 - credMass) / 2)
  )
  message("Done.")
  out
}
results_sim <- run_all_score_models(
  sim = sim,
  truth_target = "asset",
  mfsv_factors = 3,
  mfsv_draws = 12000,
  mfsv_thin = 4,
  mfsv_burnin = 1500,
  dsp_factors = 3,
  dsp_nsave = 3000,
  dsp_nburn = 1500,
  dsp_nskip = 4,
  roll_window = 60,
  ew_lambda = 0.94,
  B_nonbayes = 200,
  block_length = 20,
  bootstrap_seed = seed,
  credMass = 0.95, 
  ht_tau = 0.10, 
  cholsv_lags = 1,
  cholsv_draws = 12000,
  cholsv_burnin = 1500,
  cholsv_thin = 4,
  cholsv_u_prior = "HS"
)

#Converts the model results object into a long-format data frame for plotting or comparison.
as_score_df <- function(res) {
  TT <- nrow(res$truth)
  
  truth_df <- data.frame(
    t = seq_len(TT),
    method = "truth",
    mean = res$truth$true_score,
    lower = NA_real_,
    upper = NA_real_
  )
  
  method_dfs <- lapply(names(res$methods), function(nm) {
    obj <- res$methods[[nm]]
    data.frame(
      t = seq_len(length(obj$mean)),
      method = nm,
      mean = obj$mean,
      lower = obj$lower,
      upper = obj$upper
    )
  })
  
  do.call(rbind, c(list(truth_df), method_dfs))
}

score_df <- as_score_df(results_sim)

true_score_ts <- results_sim$truth$true_score

mfsv_mean  <- results_sim$methods$mfsv$mean
mfsv_lower <- results_sim$methods$mfsv$lower
mfsv_upper <- results_sim$methods$mfsv$upper

dsp_mean   <- results_sim$methods$dsp_mfsv$mean
dsp_lower  <- results_sim$methods$dsp_mfsv$lower
dsp_upper  <- results_sim$methods$dsp_mfsv$upper

ew_mean    <- results_sim$methods$ew_roll$mean
ew_lower   <- results_sim$methods$ew_roll$lower
ew_upper   <- results_sim$methods$ew_roll$upper

lw_mean    <- results_sim$methods$roll_lw$mean
lw_lower   <- results_sim$methods$roll_lw$lower
lw_upper   <- results_sim$methods$roll_lw$upper

dcc_mean   <- results_sim$methods$dcc_norm$mean
dcc_lower  <- results_sim$methods$dcc_norm$lower
dcc_upper  <- results_sim$methods$dcc_norm$upper

dcct_mean  <- results_sim$methods$dcc_t$mean
dcct_lower <- results_sim$methods$dcc_t$lower
dcct_upper <- results_sim$methods$dcc_t$upper

adcc_mean  <- results_sim$methods$adcc$mean
adcc_lower <- results_sim$methods$adcc$lower
adcc_upper <- results_sim$methods$adcc$upper

ht_mean   <- results_sim$methods$ht_roll$mean
ht_lower  <- results_sim$methods$ht_roll$lower
ht_upper  <- results_sim$methods$ht_roll$upper

loggle_mean  <- results_sim$methods$loggle$mean
loggle_lower <- results_sim$methods$loggle$lower
loggle_upper <- results_sim$methods$loggle$upper

chol_mean   <- results_sim$methods$chol_sv$mean
chol_lower  <- results_sim$methods$chol_sv$lower
chol_upper  <- results_sim$methods$chol_sv$upper

final_output <- list(
  true_score_ts = true_score_ts,
  mfsv = list(mean = mfsv_mean, lower = mfsv_lower, upper = mfsv_upper),
  dsp_mfsv = list(mean = dsp_mean, lower = dsp_lower, upper = dsp_upper),
  ew_roll = list(mean = ew_mean, lower = ew_lower, upper = ew_upper),
  roll_lw = list(mean = lw_mean, lower = lw_lower, upper = lw_upper),
  ht_roll = list(mean = ht_mean, lower = ht_lower, upper = ht_upper),
  dcc_norm = list(mean = dcc_mean, lower = dcc_lower, upper = dcc_upper),
  dcc_t = list(mean = dcct_mean, lower = dcct_lower, upper = dcct_upper),
  adcc = list(mean = adcc_mean, lower = adcc_lower, upper = adcc_upper),
  chol_sv = list(mean = chol_mean, lower = chol_lower, upper = chol_upper),
  loggle = list(mean = loggle_mean, lower = loggle_lower, upper = loggle_upper),
  score_df = score_df,
  raw_results = results_sim
)