############################################################
## 1. Locate and load all 2020 large-cross-section results
############################################################

downloads <- path.expand("~/Downloads")

input_file <- file.path(
  downloads,
  "large_cross_section_cluster_input_2020.rds"
)

results_dir <- file.path(
  downloads,
  "large_cross_section_2020_results"
)

stopifnot(
  file.exists(input_file),
  dir.exists(results_dir)
)

input <- readRDS(
  input_file
)

dsp_rep <- readRDS(
  file.path(
    results_dir,
    "dsp_representatives.rds"
  )
)

residual_fits <- lapply(
  1:39,
  function(g) {
    readRDS(
      file.path(
        results_dir,
        sprintf(
          "residual_block_%02d.rds",
          g
        )
      )
    )
  }
)

A <- input$A
dates <- as.Date(input$dates)

stopifnot(
  all(dim(A) == c(1000, 39)),
  all(dim(dsp_rep$covariance_mean) == c(39, 39, 1257)),
  length(residual_fits) == 39,
  length(dates) == 1257
)

cat("ALL 2020 RESULTS LOADED\n")


############################################################
## 2. Function to reconstruct the full 1000 x 1000
##    covariance matrix at any date
############################################################

get_covariance_1000 <- function(tt) {
  
  Sigma_rep <-
    dsp_rep$covariance_mean[, , tt]
  
  Sigma_resid <-
    matrix(
      0,
      nrow = 1000,
      ncol = 1000
    )
  
  for (g in 1:39) {
    
    idx <- residual_fits[[g]]$indices
    
    Sigma_resid[
      idx,
      idx
    ] <-
      residual_fits[[g]]$
      covariance_mean[, , tt]
  }
  
  Sigma <-
    A %*%
    Sigma_rep %*%
    t(A) +
    Sigma_resid
  
  (Sigma + t(Sigma)) / 2
}


############################################################
## Example: first 1000 x 1000 covariance matrix
############################################################

Sigma_1 <- get_covariance_1000(1)

dim(Sigma_1)
# should be 1000 1000


############################################################
## 3. Log determinant helper
############################################################

logdet_pd <- function(M) {
  
  M <- (M + t(M)) / 2
  
  ch <- chol(M)
  
  2 * sum(
    log(
      diag(ch)
    )
  )
}


############################################################
## 4. Exact determinant-score time series
##
## det(R_t) = det(Sigma_t) /
##            product(diag(Sigma_t))
##
## With the representative/residual decomposition:
##
## det(Sigma_t) =
## det(Sigma_rep,t) *
## product_g det(Sigma_resid,g,t)
############################################################

determinant_score <- numeric(
  length(dates)
)

for (tt in seq_along(dates)) {
  
  if (tt %% 50 == 0) {
    cat(
      "Processing",
      tt,
      "of",
      length(dates),
      "\n"
    )
  }
  
  Sigma_rep <-
    dsp_rep$covariance_mean[, , tt]
  
  ##########################################################
  ## Representative contribution to variances
  ##########################################################
  
  AS <-
    A %*% Sigma_rep
  
  variances <-
    rowSums(
      AS * A
    )
  
  ##########################################################
  ## Residual contribution + residual log determinants
  ##########################################################
  
  logdet_resid <- 0
  
  for (g in 1:39) {
    
    idx <-
      residual_fits[[g]]$indices
    
    Sg <-
      residual_fits[[g]]$
      covariance_mean[, , tt]
    
    variances[idx] <-
      variances[idx] +
      diag(Sg)
    
    logdet_resid <-
      logdet_resid +
      logdet_pd(Sg)
  }
  
  ##########################################################
  ## Full covariance and correlation log determinants
  ##########################################################
  
  logdet_cov <-
    logdet_pd(Sigma_rep) +
    logdet_resid
  
  logdet_cor <-
    logdet_cov -
    sum(log(variances))
  
  determinant_score[tt] <-
    1 -
    exp(
      logdet_cor / 1000
    )
}


############################################################
## 5. Results
############################################################

score_1000 <- data.frame(
  Date = dates,
  Determinant_Score =
    determinant_score
)

head(score_1000)
tail(score_1000)

summary(
  score_1000$Determinant_Score
)


############################################################
## 6. Plot
############################################################

library(ggplot2)

p_score_1000 <- ggplot(
  score_1000,
  aes(
    x = Date,
    y = Determinant_Score
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  labs(
    x = "Date",
    y = "Determinant score"
  ) +
  theme_minimal(
    base_size = 14
  )

p_score_1000