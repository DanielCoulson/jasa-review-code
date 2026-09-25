############################################################
## 1. load all 2020 large-cross-section results
############################################################

downloads <- path.expand("~/Downloads")
input_file <- file.path(downloads,"large_cross_section_cluster_input_2020.rds")
results_dir <- file.path(downloads,"large_cross_section_2020_results")
input <- readRDS(input_file)
dsp_rep <- readRDS(file.path(results_dir,"dsp_representatives.rds"))
residual_fits <- lapply(1:39,function(g) {
    readRDS(file.path(results_dir,sprintf("residual_block_%02d.rds",g)))
  }
)
A <- input$A
dates <- as.Date(input$dates)

############################################################
## 2. Function to reconstruct the full 1000 x 1000
##    covariance matrix at any date
############################################################

get_covariance_1000 <- function(tt) {
  Sigma_rep <-dsp_rep$covariance_mean[, , tt]
  Sigma_resid <-matrix(0,nrow = 1000,ncol = 1000)
  for (g in 1:39) {
    idx <- residual_fits[[g]]$indices
    Sigma_resid[idx,idx] <-residual_fits[[g]]$covariance_mean[, , tt]
  }
  Sigma <-A %*%Sigma_rep %*%t(A) +Sigma_resid
  (Sigma + t(Sigma)) / 2
}

############################################################
## 3. Log determinant helper
############################################################

logdet_pd <- function(M) {
  M <- (M + t(M)) / 2
  ch <- chol(M)
  2 * sum(log(diag(ch)))
}
############################################################
## 4. Exact determinant-score time series
############################################################

determinant_score <- numeric(length(dates))
for (tt in seq_along(dates)) {
  if (tt %% 50 == 0) {
    cat("Processing",tt,"of",length(dates),"\n")
  }
  Sigma_rep <-dsp_rep$covariance_mean[, , tt]
  AS <-A %*% Sigma_rep
  variances <-rowSums(AS * A)
  logdet_resid <- 0
  for (g in 1:39) {
    idx <-residual_fits[[g]]$indices
    Sg <-residual_fits[[g]]$covariance_mean[, , tt]
    variances[idx] <-variances[idx] +diag(Sg)
    logdet_resid <-logdet_resid +logdet_pd(Sg)
  }
  logdet_cov <- logdet_pd(Sigma_rep) +logdet_resid
  logdet_cor <-logdet_cov -sum(log(variances))
  determinant_score[tt] <-1 -exp(logdet_cor / 1000)
}

############################################################
## 5. Results
############################################################

score_1000 <- data.frame(Date = dates,Determinant_Score =determinant_score)
