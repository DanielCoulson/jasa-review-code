############################################################
## Large cross-section cluster fitting worker
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    "Usage: Rscript large_cross_section_fit.R ",
    "<base_dir> <task_id>"
  )
}

base_dir <- normalizePath(
  args[1],
  mustWork = TRUE
)

task_id <- as.integer(args[2])

if (
  !is.finite(task_id) ||
  task_id < 1L ||
  task_id > 51L
) {
  stop("task_id must be an integer from 1 to 51.")
}

input <- readRDS(
  file.path(
    base_dir,
    "large_cross_section_cluster_input.rds"
  )
)

output_dir <- file.path(
  base_dir,
  "output"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(1000L + task_id)


############################################################
## TASK 1:
## DSP-MFSV-CAPM on the 50 representative stocks
############################################################

if (task_id == 1L) {
  
  source(
    file.path(
      base_dir,
      "DSP-MFSV-CAPM_large.R"
    )
  )
  
  Y_rep <-
    input$representative_excess_returns
  
  X_rep <- cbind(
    Intercept = 1,
    Market = input$market_excess_return
  )
  
  stopifnot(
    all(dim(Y_rep) == c(1006, 50)),
    all(dim(X_rep) == c(1006, 2)),
    all(is.finite(Y_rep)),
    all(is.finite(X_rep))
  )
  
  cat(
    "Fitting 50-stock DSP-MFSV-CAPM\n"
  )
  
  fit <- DSP_MFSV(
    y = Y_rep,
    X = X_rep,
    D = 1,
    number_of_latent_factors = 3,
    nsave = 3000,
    nburn = 1500,
    nskip = 4,
    mcmc_params = c(
      "covariance_mean"
    )
  )
  
  stopifnot(
    all(
      dim(fit$covariance_mean) ==
        c(50, 50, 1006)
    ),
    all(is.finite(fit$covariance_mean))
  )
  
  saveRDS(
    list(
      type = "DSP_representatives",
      covariance_mean =
        fit$covariance_mean
    ),
    file.path(
      output_dir,
      "dsp_representatives.rds"
    )
  )
  
  cat(
    "DSP representative fit complete\n"
  )
  
  
  ############################################################
  ## TASKS 2--51:
  ## Ordinary MFSV on residual blocks 1--50
  ############################################################
  
} else {
  
  library(factorstochvol)
  
  g <- task_id - 1L
  
  Y_block <-
    input$residual_blocks[[g]]
  
  n_g <- ncol(Y_block)
  
  stopifnot(
    nrow(Y_block) == 1006,
    n_g >= 2L,
    all(is.finite(Y_block))
  )
  
  # At most 3 latent factors, and strictly fewer
  # factors than series in very small blocks.
  factors_g <- min(
    3L,
    n_g - 1L
  )
  
  cat(
    "Fitting residual block",
    g,
    "with",
    n_g,
    "stocks and",
    factors_g,
    "factors\n"
  )
  
  fit <- factorstochvol::fsvsample(
    Y_block,
    factors = factors_g,
    draws = 12000,
    thin = 4,
    burnin = 1500,
    keeptime = "all",
    quiet = TRUE
  )
  
  T_obs <- nrow(Y_block)
  
  covariance_mean <-
    array(
      NA_real_,
      dim = c(
        n_g,
        n_g,
        T_obs
      )
    )
  
  for (tt in seq_len(T_obs)) {
    
    covariance_mean[, , tt] <-
      factorstochvol::runningcovmat(
        fit,
        i = tt,
        statistic = "mean"
      )
  }
  
  stopifnot(
    all(is.finite(covariance_mean))
  )
  
  saveRDS(
    list(
      type = "MFSV_residual_block",
      block = g,
      indices =
        input$residual_block_indices[[g]],
      factors = factors_g,
      covariance_mean =
        covariance_mean
    ),
    file.path(
      output_dir,
      sprintf(
        "residual_block_%02d.rds",
        g
      )
    )
  )
  
  cat(
    "Residual block",
    g,
    "complete\n"
  )
}