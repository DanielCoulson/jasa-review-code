#########################################
########DSP-MFSV-CAPM model code#########
#########################################
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
  if (!is.null(chol0)) {
    QHt_Matrix = as.spam.dgCMatrix(as(Qpost, "dgCMatrix"))
    if (return_mean){
      post_mean_vec = spam::solve.spam(QHt_Matrix, linht, Rstruct = chol0)
      post_mean     = matrix(post_mean_vec, nrow = T, byrow = TRUE)
    }
    beta = matrix(rmvnorm.canonical(n = 1, b = linht, Q = QHt_Matrix, 
                                    Rstruct = chol0), nrow = T, byrow = TRUE)
    if (return_var){
      chQht_Matrix = Matrix::chol(Qpost)
    }
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
    invChol_I <- Matrix::solve(chQht_Matrix, Diagonal(n = nrow(Qpost)))
    post_var <- rowSums(invChol_I^2)
    post_var_mat <- matrix(post_var, nrow = T, ncol = p, byrow = TRUE)
  }
  list(beta = beta,post_mean = post_mean ,post_var = post_var_mat)
}

score <- function(R) {
##function to compute the scalar summary of a single correlation matrix
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

DSP_MFSV = function (y, X = NULL, D = 1,number_of_latent_factors, nsave = 1000, nburn = 1000, nskip = 4, mcmc_params = list("mu", "yhat", "beta", 
                                                                                                                                                "evol_sigma_t2", "obs_sigma_t2", "dhs_phi","dhs_mean","obs_eps_cov_matrix", "scores")) 
#####function that fits the DSP-MFSV CAPM model#####
### y = input time series, X = design matrix series, D = order of differencing ### 
### number_of_latent_factor = number of latent factors used in the MFSV model for the observation error covariance ###
### nsave = number of saved MCMC samples, nburn = burnin period of the MCMC sampler, nskip = level of thinning for the MCMC sampler###
### mcmc_params = list of parameters to save, and the function returns MCMC samples of these quantities of interest ###
{
  model_variance = svsample(y = X[,2], draws = (nsave*(nskip+1)+nburn), burnin = nburn, thin = 1, quiet = TRUE) 
  variances = (exp(model_variance$latent[[1]]/2)^2)
  T = nrow(y)
  t01 = seq(0, 1, length.out = T)
  XtX = build_XtX(X)
  p = ncol(X)
  N = ncol(y)
  save_scores = ("scores" %in% mcmc_params) 
  sigma_e = apply(y, 2, sd, na.rm = TRUE)
  sigma_et = matrix(nrow = T, ncol = ncol(y))
  sigma_et <- matrix(rep(sigma_e, each = T), nrow = T)
  chol0 = initCholReg.spam(obs_sigma_t2 = abs(rnorm(T)), evol_sigma_t2 = matrix(abs(rnorm(T *p)), nrow = T), XtX = XtX, D = D)
  beta = array(dim = c(ncol(y),T,p))
  varbeta = array(dim = c(ncol(y),T,p))
  meanbeta = array(dim = c(ncol(y),T,p))
  for(i in 1:N){
    res = sampleBTF_reg(
      y[,i], X, 
      obs_sigma_t2 = sigma_et[,i], 
      evol_sigma_t2 = matrix(0.01 * sigma_et[,i], nrow = T, ncol = p), 
      XtX = XtX, D = D, 
      chol0 = chol0, 
      return_mean = FALSE, 
      return_var =FALSE)
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
    evolParams[[i]] = initEvolParams(omega[i,,], evol_error = "DHS")
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
  if (!is.na(match("dhs_phi", mcmc_params)) ) 
    post_dhs_phi = array(NA, c(nsave,ncol(y), p))
  if (!is.na(match("dhs_mean", mcmc_params)) ) 
    post_dhs_mean = array(NA, c(nsave,ncol(y), p))
  if(!is.na(match("obs_eps_cov_matrix",mcmc_params)))
    post_obs_eps_cov_matrix = array(0,c(ncol(y), ncol(y),T))
  if(!is.na(match("scores",mcmc_params)))
    post_scores = matrix(NA,nrow = T, ncol = nsave)
  nstot = nburn + (nskip + 1) * (nsave)
  skipcount = 0
  isave = 0
  for(nsi in 1:nstot){
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
      evolParams[[i]] = sampleEvolParams(omega[i,,], evolParams[[i]],1/sqrt(T * p), "DHS" )
    }
    something = y-mu
    svParams = fsvsample(something,factors = number_of_latent_factors, draws = 1, thin = 1, burnin = 0,keeptime = "all", quiet = TRUE,startfac = svParams$fac[,,1],
                         startpara = svParams$para[,,1],startlogvar = svParams$logvar[,,1],startlogvar0 = svParams$logvar0[,1],startfacload = svParams$facload[,,1])
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
      for(i in 1:T){
        scores[i] =   score(
          cov2cor(
            diag(varbeta[,i,1])
            + variances[nsi,i] * (diag(varbeta[,i,2]) + tcrossprod(meanbeta[,i,2])) 
            +stuff[,,i]
          )
        )
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
            post_yhat[isave,i,] = mu[,i] + sigma_et[,i]*rnorm(T)
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
        if (!is.na(match("dhs_phi", mcmc_params))) 
          for( i in 1:ncol(y)){
            post_dhs_phi[isave,i,] = evolParams[[i]]$dhs_phi
          }
        if (!is.na(match("dhs_mean", mcmc_params))) 
          for(i in 1:ncol(y)){
            post_dhs_mean[isave,i,] = evolParams[[i]]$dhs_mean
          }
        if(!is.na(match("obs_eps_cov_matrix",mcmc_params)))
          post_obs_eps_cov_matrix <- post_obs_eps_cov_matrix + stuff
        if(save_scores)
          post_scores[,isave] = scores
        skipcount = 0
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
  if (!is.na(match("dhs_phi", mcmc_params))) 
    mcmc_output$dhs_phi = post_dhs_phi
  if (!is.na(match("dhs_mean", mcmc_params))) 
    mcmc_output$dhs_mean = post_dhs_mean
  if(!is.na(match("obs_eps_cov_matrix",mcmc_params)))
    mcmc_output$obs_eps_cov_matrix = post_obs_eps_cov_matrix
  if(!is.na(match("scores",mcmc_params)))
    mcmc_output$scores = post_scores
  return(mcmc_output)
}
