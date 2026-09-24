##################################################################################
#### Code for 5-factor example######
##################################################################################
####import packages####
library(stochvol) #https://cran.r-project.org/web/packages/stochvol/index.html
library(mvtnorm) #https://cran.r-project.org/web/packages/mvtnorm/index.html
library(Matrix) #https://cran.r-project.org/web/packages/Matrix/index.html
library(dsp) #https://github.com/drkowal/dsp
library(factorstochvol) #https://cran.r-project.org/web/packages/factorstochvol/index.html
library(HDInterval) #https://cran.r-project.org/web/packages/HDInterval/index.html
library(quantmod) #https://cran.r-project.org/web/packages/quantmod/index.html
library(ggplot2) #https://cran.r-project.org/web/packages/ggplot2/index.html
library(tidyr) #https://cran.r-project.org/web/packages/tidyr/index.html
library(dplyr) #https://cran.r-project.org/web/packages/dplyr/index.html 
library(spam) #https://cran.r-project.org/web/packages/spam/index.html
############################################################
## Cluster paths
############################################################

args <- commandArgs(trailingOnly = TRUE)

project_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")
checkpoint_dir <- file.path(project_dir, "checkpoints")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:", project_dir, "\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Checkpoint directory:", checkpoint_dir, "\n")

set.seed(1)

if (!requireNamespace("Rcpp", quietly = TRUE))
  stop("Package 'Rcpp' is required.")

if (!requireNamespace("RcppArmadillo", quietly = TRUE))
  stop("Package 'RcppArmadillo' is required.")

rcpp_cache_dir <- file.path(project_dir, ".rcpp_cache")
dir.create(rcpp_cache_dir, recursive = TRUE, showWarnings = FALSE)

Rcpp::sourceCpp(
  code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List btf_selected_cov_cpp(
    const arma::mat& X,
    const arma::vec& obs_sigma_t2,
    const arma::mat& evol_sigma_t2) {
  const arma::uword T = X.n_rows;
  const arma::uword p = X.n_cols;
  if (T == 0 || p == 0)
    Rcpp::stop("X must have positive dimensions.");
  if (obs_sigma_t2.n_elem != T)
    Rcpp::stop("obs_sigma_t2 must have length nrow(X).");
  if (evol_sigma_t2.n_rows != T ||evol_sigma_t2.n_cols != p)
    Rcpp::stop("evol_sigma_t2 must be nrow(X) x ncol(X).");
  if (!X.is_finite() ||!obs_sigma_t2.is_finite() ||!evol_sigma_t2.is_finite())
    Rcpp::stop("All inputs must be finite.");
  for (arma::uword t = 0; t < T; ++t) {
    if (obs_sigma_t2(t) <= 0.0)
      Rcpp::stop("obs_sigma_t2 must be strictly positive.");
    for (arma::uword j = 0; j < p; ++j) {
      if (evol_sigma_t2(t,j) <= 0.0)
        Rcpp::stop("evol_sigma_t2 must be strictly positive.");
    }
  }
  arma::cube G(p, p, T, arma::fill::zeros);
  arma::cube P(p, p, T, arma::fill::zeros);
  for (arma::uword t = 0; t < T; ++t) {
    arma::vec xt = X.row(t).t();
    arma::mat F = (xt * xt.t()) / obs_sigma_t2(t);
    for (arma::uword j = 0; j < p; ++j) {
      F(j,j) += 1.0 / evol_sigma_t2(t,j);
      if (t + 1 < T)
        F(j,j) += 1.0 / evol_sigma_t2(t + 1,j);
    }
    if (t > 0) {
      arma::vec c(p);
      for (arma::uword j = 0; j < p; ++j)
        c(j) =-1.0 / evol_sigma_t2(t,j);
      F -= G.slice(t - 1) %(c * c.t());
    }
    F = 0.5 * (F + F.t());
    arma::mat invF;
    bool ok = arma::inv_sympd(invF, F);
    if (!ok)
      Rcpp::stop("Selected-inverse recursion failed at time %d.",static_cast<int>(t + 1));
    G.slice(t) = 0.5 * (invF + invF.t());
  }
  P.slice(T - 1) = G.slice(T - 1);
  if (T > 1) {
    for (int tt = static_cast<int>(T) - 2;tt >= 0;--tt) {
      arma::uword t =static_cast<arma::uword>(tt);
      arma::vec c(p);
      for (arma::uword j = 0; j < p; ++j)
        c(j) =-1.0 /evol_sigma_t2(t + 1,j);
      arma::mat GC = G.slice(t);
      for (arma::uword j = 0; j < p; ++j)
        GC.col(j) *= c(j);
      arma::mat Pt =G.slice(t) +GC *P.slice(t + 1) *GC.t();
      P.slice(t) =0.5 * (Pt + Pt.t());
    }
  }
  arma::mat post_var(T,p,arma::fill::zeros);
  for (arma::uword t = 0; t < T; ++t)
    for (arma::uword j = 0; j < p; ++j)
      post_var(t,j) = P(j,j,t);
  return Rcpp::List::create(Rcpp::Named("cov") = P,Rcpp::Named("var") = post_var);
}
',
cacheDir = rcpp_cache_dir,rebuild = FALSE,showOutput = FALSE,verbose = FALSE,echo = FALSE)
#### Define functions for use later ####
sampleBTF_reg = function (y, X, obs_sigma_t2, evol_sigma_t2, XtX, D = 1, chol0 = NULL, return_mean = FALSE, return_var = FALSE) 
{
  #### function to generate a single sample from the posterior distribution of the regression parameters####
  if ((D < 0) || (D != round(D))) 
    stop("D must be a positive integer")
  if (any(is.na(y))) 
    stop("y cannot contain NAs")
  T = nrow(X)
  p = ncol(X)
  if (D == 1) {
    t_evol_prec_lag_mat = matrix(0, nrow = p, ncol = T)
    t_evol_prec_lag_mat[, 1:(T - 1)] = t(1/evol_sigma_t2[-1,])
    Q_diag = matrix(t(1/evol_sigma_t2) + t_evol_prec_lag_mat)
    Q_off = matrix(-t_evol_prec_lag_mat)[-(T * p)]
    Qevol = bandSparse(T * p, k = c(0, p), diagonals = list(Q_diag, Q_off), symmetric = TRUE)
  }
  else {
    if (D == 2) {
      t_evol_prec_lag2 = t(1/evol_sigma_t2[-(1:2), ])
      Q_diag = t(1/evol_sigma_t2)
      Q_diag[, 2:(T - 1)] = Q_diag[, 2:(T - 1)] + 4*t_evol_prec_lag2
      Q_diag[, 1:(T - 2)] = Q_diag[, 1:(T - 2)] + t_evol_prec_lag2
      Q_diag = matrix(Q_diag)
      Q_off_1 = matrix(0, nrow = p, ncol = T)
      Q_off_1[, 1] = -2/evol_sigma_t2[3, ]
      Q_off_1[, 2:(T - 1)] = Q_off_1[, 2:(T - 1)] + -2*t_evol_prec_lag2
      Q_off_1[, 2:(T - 2)] = Q_off_1[, 2:(T - 2)] + -2*t_evol_prec_lag2[, -1]
      Q_off_1 = matrix(Q_off_1)
      Q_off_2 = matrix(0, nrow = p, ncol = T)
      Q_off_2[, 1:(T - 2)] = t_evol_prec_lag2
      Q_off_2 = matrix(Q_off_2)
      Qevol = bandSparse(T * p, k = c(0, p, 2 * p), diagonals = list(Q_diag, Q_off_1, Q_off_2), symmetric = TRUE)
    }
    else stop("sampleBTF_reg() requires D=1 or D=2")
  }
  Qobs = 1/rep(obs_sigma_t2, each = p) * XtX
  Qpost = Qobs + Qevol
  linht = matrix(t(X * as.numeric(y/obs_sigma_t2)))
  post_mean = NULL 
  post_var_mat = NULL
  post_cov_blocks = NULL
  if (!is.null(chol0)) {
    QHt_Matrix = as.spam.dgCMatrix(as(Qpost, "dgCMatrix"))
    if (return_mean){
      post_mean_vec = spam::solve.spam(QHt_Matrix, linht, Rstruct = chol0)
      post_mean     = matrix(post_mean_vec, nrow = T, byrow = TRUE)
    }
    beta = matrix(rmvnorm.canonical(n = 1, b = linht, Q = QHt_Matrix, 
                                    Rstruct = chol0), nrow = T, byrow = TRUE)
  }
  else {
    chQht_Matrix = Matrix::chol(Qpost)
    rhs = Matrix::solve(Matrix::t(chQht_Matrix), linht)
    if (return_mean){
      post_mean_vec = Matrix::solve (chQht_Matrix, rhs)
      post_mean = matrix(post_mean_vec, nrow = T, byrow = TRUE)
    }
    beta = matrix(Matrix::solve(chQht_Matrix,rhs + rnorm(T*p)), nrow = T, byrow = TRUE)
  }
  if (return_var){
    if (D != 1)
      stop("The contemporaneous posterior covariance calculation requires D = 1.")
    selected_cov <- btf_selected_cov_cpp(X = X,obs_sigma_t2 = as.numeric(obs_sigma_t2),
      evol_sigma_t2 = evol_sigma_t2)
    post_cov_blocks <- aperm(selected_cov$cov,c(3, 1, 2))
    post_var_mat <- selected_cov$var
  }
  list(beta = beta,post_mean = post_mean,post_var = post_var_mat,
    post_cov = post_cov_blocks
  )
}
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

DSP_MFSV = function (y, X = NULL, evol_error = "DHS", D = 1,number_of_latent_factors, nsave = 1000, nburn = 1000, nskip = 4, mcmc_params = list("mu", "yhat", "beta", 
                                                                                                                                                "evol_sigma_t2", "obs_sigma_t2", "dhs_phi","dhs_mean","obs_eps_cov_matrix", "scores"), checkpoint_file = NULL, checkpoint_every = 50, resume = TRUE) 
  #####function that fits the DSP-MFSV model#####
### y = input dependent time series, X = design matrix series, evol_error specifies the prior distribution, D = order of differencing ### 
### number_of_latent_factor = number of latent factors used in the MFSV model for the observation error covariance ###
### nsave = number of saved MCMC samples, nburn = burn-in period of the MCMC sampler, nskip = level of thinning for the MCMC sampler###
### mcmc_params = list of parameters to save, and the function returns MCMC samples of these quantities of interest ###
{
  computeDIC = FALSE
  evol_error = toupper(evol_error)
  if (nrow(X) != nrow(y))
    stop("X and y must have the same number of time points.")
  
  if (ncol(X) != 6)
    stop("The FF5 design matrix must contain an intercept plus five factors.")
  checkpoint_state = NULL
  
  if (!is.null(checkpoint_file) &&
      resume &&
      file.exists(checkpoint_file)) {
    
    checkpoint_state = readRDS(checkpoint_file)
    
    cat(
      "Resuming from checkpoint:",
      "iteration", checkpoint_state$nsi,
      "- saved draw", checkpoint_state$isave, "\n"
    )
  }
  
  factor_cov_file =
    if (!is.null(checkpoint_file))
      paste0(checkpoint_file, "_factor_cov.rds")
  else
    NULL
  risk_factors = scale(X[,2:6, drop = FALSE], center = TRUE,scale = FALSE)
  if (!is.null(checkpoint_state)) {
    if (is.null(factor_cov_file) ||
        !file.exists(factor_cov_file))
      stop("Checkpoint exists but factor covariance file is missing.")
    factor_cov = readRDS(factor_cov_file)
  } else {
    factor_model = fsvsample(y = risk_factors,factors = 3,
      draws = nsave*(nskip+1),burnin = nburn,thin = nskip+1,
      keeptime = "all",zeromean = TRUE,quiet = TRUE)
    factor_cov = factorstochvol::covmat(factor_model,timepoints = "all")
    rm(factor_model)
    gc()
    if (!is.null(factor_cov_file))
      saveRDS(factor_cov,factor_cov_file,compress = FALSE
      )
  }
  expected_factor_cov_dim =c(5, 5, nsave, nrow(y))
  if (!identical(as.integer(dim(factor_cov)),as.integer(expected_factor_cov_dim))) {
    stop(paste("Unexpected factor_cov dimensions:",
        paste(dim(factor_cov), collapse = " x "),"- expected",
        paste(expected_factor_cov_dim, collapse = " x ")))
  }
  T = nrow(y)
  t01 = seq(0, 1, length.out = T)
  XtX = build_XtX(X)
  p = ncol(X)
  N = ncol(y)
  save_scores = ("scores" %in% mcmc_params) && evol_error == "DHS"
  sigma_e = apply(y, 2, sd, na.rm = TRUE)
  sigma_et = matrix(nrow = T, ncol = ncol(y))
  sigma_et <- matrix(rep(sigma_e, each = T), nrow = T)
  chol0 = initCholReg.spam(obs_sigma_t2 = abs(rnorm(T)), evol_sigma_t2 = matrix(abs(rnorm(T *p)), nrow = T), XtX = XtX, D = D)
  beta = array(dim = c(ncol(y),T,p))
  varbeta = array(dim = c(ncol(y),T,p))
  meanbeta = array(dim = c(ncol(y),T,p))
  covbeta = array(dim = c(ncol(y),T,p,p))
  for(i in 1:N){
    res = sampleBTF_reg(y[,i], X, obs_sigma_t2 = sigma_et[,i], evol_sigma_t2 = matrix(0.01 * sigma_et[,i], nrow = T, ncol = p), 
      XtX = XtX, D = D, chol0 = chol0, return_mean = FALSE,return_var =FALSE)
    beta[i, , ] <- res$beta
  }
  mu = matrix(nrow = T, ncol = ncol(y))
  for(i in 1:ncol(y)){
    mu[,i] = rowSums(X*beta[i,,])
  }
  omega = array(dim = c(ncol(y),T-1,p))
  for(i in 1:ncol(y)){
    omega[i,,] = diff(beta[i,,],differences = D)
  }
  beta0 = matrix(nrow = ncol(y), ncol = ncol(X))
  for(i in 1:ncol(y)){
    beta0[i,] = matrix(beta[i,,][1:D,], nrow = D)
  }
  evolParams = list()
  for(i in 1:ncol(y)){
    evolParams[[i]] = initEvolParams(omega[i,,], evol_error = evol_error)
  }
  evolParams0 = list()
  for(i in 1:ncol(y)){
    evolParams0[[i]] = initEvol0(beta0[i,], commonSD = FALSE)
  }
  something = y-mu
  svParams = fsvsample(something,factors = number_of_latent_factors, thin = 1, burnin = 0 , quiet = TRUE, draws = 1, keeptime = "all")
  stuff = factorstochvol::covmat(svParams, timepoints = "all")[,,1,]
  sigma_et <- t(apply(stuff, 3, diag))
  mcmc_output = vector("list", length(mcmc_params))
  names(mcmc_output) = mcmc_params
  if (!is.na(match("mu", mcmc_params)) || computeDIC) 
    post_mu = array(NA, c(nsave,ncol(y), T))
  if (!is.na(match("yhat", mcmc_params))) 
    post_yhat = array(NA, c(nsave,ncol(y), T))
  if (!is.na(match("beta", mcmc_params))) 
    post_beta = array(NA, c(nsave,ncol(y), T, p))
  if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) 
    post_obs_sigma_t2 = array(NA, c(nsave,ncol(y), T))
  if (!is.na(match("evol_sigma_t2", mcmc_params))) 
    post_evol_sigma_t2 = array(NA, c(nsave,ncol(y), T, p))
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == 
      "DHS") 
    post_dhs_phi = array(NA, c(nsave,ncol(y), p))
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == 
      "DHS") 
    post_dhs_mean = array(NA, c(nsave,ncol(y), p))
  if(!is.na(match("obs_eps_cov_matrix",mcmc_params))&&evol_error == "DHS")
    post_obs_eps_cov_matrix = array(0,c(ncol(y), ncol(y),T))
  if(!is.na(match("scores",mcmc_params))&&evol_error == "DHS")
    post_scores = matrix(NA,nrow = T, ncol = nsave)
  nstot = nburn + (nskip + 1) * (nsave)
  skipcount = 0
  isave = 0
  start_nsi = 1L
  if (!is.null(checkpoint_state)) {
    beta = checkpoint_state$beta
    mu = checkpoint_state$mu
    omega = checkpoint_state$omega
    beta0 = checkpoint_state$beta0
    evolParams = checkpoint_state$evolParams
    evolParams0 = checkpoint_state$evolParams0
    svParams = checkpoint_state$svParams
    sigma_et = checkpoint_state$sigma_et
    post_scores = checkpoint_state$post_scores
    isave = checkpoint_state$isave
    skipcount = checkpoint_state$skipcount
    start_nsi = checkpoint_state$nsi + 1L
    assign(".Random.seed",checkpoint_state$random_seed,envir = .GlobalEnv)
  }
  iterations_to_run <- if (start_nsi <= nstot) {
    seq.int(start_nsi, nstot)
  } else {
    integer(0)
  }
  for(nsi in iterations_to_run){
    will_save_this_iter = (nsi > nburn) && (skipcount == nskip)
    need_score_moments = save_scores && will_save_this_iter
    if (nsi %% 100 == 0) {
      cat("Iteration", nsi, "of", nstot, "\n")
    }
    for(i in 1:N){
      res = sampleBTF_reg(
        y[,i], X, 
        obs_sigma_t2 = sigma_et[,i], 
        evol_sigma_t2 = rbind(
          matrix(evolParams0[[i]]$sigma_w0^2,nrow = D),                                                                                     
          evolParams[[i]]$sigma_wt^2), 
        XtX = XtX, D = D, chol0 = chol0,
        return_mean = need_score_moments,
        return_var = need_score_moments )
      beta[i, , ] <- res$beta
      if (need_score_moments){
        varbeta[i, , ] <- res$post_var
        meanbeta[i,,] <- res$post_mean
        covbeta[i,,,] = res$post_cov
      }
    }
    for(i in 1:ncol(y)){
      mu[,i] = rowSums(X*beta[i,,])
    }
    for(i in 1:ncol(y)){
      omega[i,,] = diff(beta[i,,],differences = D)
    }
    for(i in 1:ncol(y)){
      beta0[i,] = matrix(beta[i,,][1:D,], nrow = D)
    }
    for(i in 1:ncol(y)){
      evolParams0[[i]] = sampleEvol0(beta0[i,], evolParams0[[i]], A = 1, commonSD = FALSE)
    }
    for(i in 1:ncol(y)){
      evolParams[[i]] = sampleEvolParams(omega[i,,], evolParams[[i]],1/sqrt(T * p), evol_error )
    }
    something = y-mu
    old_svParams = svParams
    attempt = 1
    max_attempts = 20
    repeat {
      new_svParams = tryCatch(
        fsvsample(something,factors = number_of_latent_factors,
          draws = 1,thin = 1,burnin = 0,keeptime = "all",quiet = TRUE,
          startfac = old_svParams$fac[,,1],
          startpara = old_svParams$para[,,1],
          startlogvar = old_svParams$logvar[,,1],
          startlogvar0 = old_svParams$logvar0[,1],
          startfacload = old_svParams$facload[,,1]
        ),
        error = function(e) NULL
      )
      if (!is.null(new_svParams)) {
        svParams = new_svParams
        break
      }
      attempt = attempt + 1
      if (attempt > max_attempts)
        stop("Residual MFSV failed after 20 retry attempts.")
    }
    stuff = factorstochvol::covmat(svParams, timepoints = "all")[,,1,]
    N = ncol(y)
    idx = cbind(rep(seq_len(N), T),rep(seq_len(N), T),rep(seq_len(T), each = N) )
    sigma_et <- matrix(stuff[idx], nrow = T, ncol = N, byrow = TRUE)
    for(j in 1:ncol(y)){
      for(i in 1:T){
        sigma_et[i,j] = stuff[j,j,i]
      }
    }
    if (need_score_moments){
      scores = numeric(T)
      factor_draw = isave + 1L
      for(i in 1:T){
        mean_B_t = meanbeta[,i,2:6]
        H_t = factor_cov[,,factor_draw,i]
        systematic_cov =mean_B_t %*%H_t %*%t(mean_B_t)
        for(j in 1:N){
          V_beta = covbeta[j,i,2:6,2:6]
          systematic_cov[j,j] = systematic_cov[j,j] +sum(H_t * t(V_beta))
        }
        scores[i] = score(cov2cor(diag(covbeta[,i,1,1]) +systematic_cov +stuff[,,i]))
      }
    }
    if (nsi > nburn) {
      skipcount = skipcount + 1
      if (skipcount > nskip) {
        isave = isave + 1
        if (!is.na(match("mu", mcmc_params)) || computeDIC) 
          for(i in 1:ncol(y)){
            post_mu[isave,i,] = mu[,i]
          }
        if (!is.na(match("yhat", mcmc_params))) 
          for(i in 1:ncol(y)){
            post_yhat[isave,i,] = mu[,i] + sqrt(sigma_et[,i])*rnorm(T)
          }
        if (!is.na(match("beta", mcmc_params))) 
          for(i in 1:ncol(y)){
            post_beta[isave,i,,] = beta[i,,]
          }
        if (!is.na(match("obs_sigma_t2", mcmc_params)) || 
            computeDIC) 
          for(i in 1:ncol(y)){
            post_obs_sigma_t2[isave,i,] = sigma_et[,i]
          }
        if (!is.na(match("evol_sigma_t2", mcmc_params))) {
          for(i in 1:ncol(y)){
            post_evol_sigma_t2[isave,i,,] = rbind(matrix(evolParams0[[i]]$sigma_w0^2,nrow = D), evolParams[[i]]$sigma_wt^2)
          }
        }
        if (!is.na(match("dhs_phi", mcmc_params)) && 
            evol_error == "DHS") 
          for( i in 1:ncol(y)){
            post_dhs_phi[isave,i,] = evolParams[[i]]$dhs_phi
          }
        if (!is.na(match("dhs_mean", mcmc_params)) && 
            evol_error == "DHS") 
          for(i in 1:ncol(y)){
            post_dhs_mean[isave,i,] = evolParams[[i]]$dhs_mean
          }
        if(!is.na(match("obs_eps_cov_matrix",mcmc_params))&& evol_error == "DHS")
          post_obs_eps_cov_matrix <- post_obs_eps_cov_matrix + stuff
        if(save_scores)
          post_scores[,isave] = scores
        skipcount = 0
        if (!is.null(checkpoint_file) &&
            isave %% checkpoint_every == 0) {
          checkpoint = list(nsi = nsi,isave = isave,skipcount = skipcount,
            beta = beta,mu = mu,omega = omega,beta0 = beta0,
            evolParams = evolParams,evolParams0 = evolParams0,
            svParams = svParams,sigma_et = sigma_et,
            post_scores = post_scores,random_seed = .Random.seed
          )
          saveRDS(checkpoint,checkpoint_file,compress = FALSE)
          cat("Checkpoint saved:","iteration", nsi,"- saved draw", isave, "\n")
        }
      }
    }
  }
  if (!is.na(match("mu", mcmc_params))) 
    mcmc_output$mu = post_mu
  if (!is.na(match("yhat", mcmc_params))) 
    mcmc_output$yhat = post_yhat
  if (!is.na(match("beta", mcmc_params))) 
    mcmc_output$beta = post_beta
  if (!is.na(match("obs_sigma_t2", mcmc_params))) 
    mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if (!is.na(match("evol_sigma_t2", mcmc_params))) 
    mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == 
      "DHS") 
    mcmc_output$dhs_phi = post_dhs_phi
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == 
      "DHS") 
    mcmc_output$dhs_mean = post_dhs_mean
  if(!is.na(match("obs_eps_cov_matrix",mcmc_params))&evol_error == "DHS")
    mcmc_output$obs_eps_cov_matrix = post_obs_eps_cov_matrix
  if(!is.na(match("scores",mcmc_params))&evol_error == "DHS")
    mcmc_output$scores = post_scores
  return(mcmc_output)
}
scale_to_range <- function(x, a, b) {
  #####function that scales a time series to lie in the range from a to b#####
  a + ((x - min(x)) * (b - a)) / (max(x) - min(x))
}
############################################################
## Prepare data
############################################################

tech_stocks <- c(
  "AAPL", "MSFT", "NVDA", "AMZN", "GOOG", "META", "TSM", "AVGO",
  "TSLA", "ORCL", "ASML", "NFLX", "SAP", "AMD", "CRM", "ADBE",
  "CSCO", "IBM", "QCOM", "INTU", "TXN", "NOW", "BABA", "PANW",
  "ADI", "MU", "LRCX", "INFY", "DELL", "INTC"
)

diversify_stocks <- c(
  "AAPL", "GOOGL", "AMZN", "KO", "ABBV",
  "PGR", "CAT", "XOM", "LIN", "AMT",
  "BOX", "UVV", "SIG", "COTY", "GKOS",
  "RDN", "ALK", "MGY", "AA", "NHI",
  "MITK", "GTN", "CAL", "PPLI", "KIDS",
  "HMN", "WNC", "NBR", "SCL", "DEA"
)

############################################################
## Fama-French data
############################################################
ff_file <- file.path(
  data_dir,
  "F-F_Research_Data_5_Factors_2x3_daily.csv"
)
if (!file.exists(ff_file)) {
  stop("Fama-French data file not found: ", ff_file)
}
daily_data <- read.csv(ff_file)
daily_data <- daily_data[13973:15229, ]
FF5 <- as.matrix(daily_data[, c("Mkt.RF", "SMB", "HML", "RMW", "CMA")])
design_matrix <- cbind(1, FF5)
FF5_MFSV <- scale(FF5,center = TRUE,scale = FALSE)
############################################################
## Load Yahoo stock/VIX data prepared locally
############################################################

market_file <- file.path(data_dir,"market_data_2019_2023.rds")
if (!file.exists(market_file)) {
  stop("Market-data RDS file not found: ", market_file)
}
market_bundle <- readRDS(market_file)
daily_tech_data <- market_bundle$daily_tech_data
daily_diversify_data <- market_bundle$daily_diversify_data
VIX_data_1 <- market_bundle$VIX_data_1
stopifnot(
  length(daily_tech_data) == length(tech_stocks),
  length(daily_diversify_data) == length(diversify_stocks))

############################################################
## Construct excess returns
############################################################
daily_tech_excess_returns <- matrix( nrow = nrow(daily_tech_data[[1]]) - 1,ncol = length(tech_stocks))
for (i in seq_along(tech_stocks)) {
  daily_tech_excess_returns[, i] <-100 * diff(log(daily_tech_data[[i]][, 6]))[-1] -daily_data$RF
}
daily_diversify_excess_returns <- matrix(nrow = nrow(daily_diversify_data[[1]]) - 1,ncol = length(diversify_stocks))
for (i in seq_along(diversify_stocks)) {
  daily_diversify_excess_returns[, i] <-100 * diff(log(daily_diversify_data[[i]][, 6]))[-1] -daily_data$RF
}
############################################################
## Fit DSP-MFSV models
############################################################
dspmfsv_daily_tech <- DSP_MFSV(y = daily_tech_excess_returns,X = design_matrix,
  nsave = 3000,nburn = 1500,number_of_latent_factors = 3,
  mcmc_params = c("scores"),checkpoint_file = file.path(checkpoint_dir,"checkpoint_tech.rds"),
  checkpoint_every = 50)
saveRDS(dspmfsv_daily_tech,file.path(output_dir, "dspmfsv_daily_tech.rds"),compress = FALSE
)
cat("Technology model completed and saved.\n")
dspmfsv_daily_diversifys <- DSP_MFSV(y = daily_diversify_excess_returns,
  X = design_matrix,nsave = 3000,nburn = 1500,number_of_latent_factors = 3,
  mcmc_params = c("scores"),checkpoint_file = file.path(checkpoint_dir,"checkpoint_diversified.rds"),
  checkpoint_every = 50)
saveRDS(dspmfsv_daily_diversifys,file.path(output_dir, "dspmfsv_daily_diversified.rds"),compress = FALSE)
cat("Diversified model completed and saved.\n")