############################################################
## Large cross-section cluster fitting worker
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    "Usage: Rscript large_cross_section_fit_2020.R ",
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
  task_id > 40L
) {
  stop("task_id must be an integer from 1 to 40.")
}

input <- readRDS(
  file.path(
    base_dir,
    "large_cross_section_cluster_input_2020.rds"
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
## DSP-MFSV-CAPM on the 39 representative stocks
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
    all(dim(Y_rep) == c(1257, 39)),
    all(dim(X_rep) == c(1257, 2)),
    all(is.finite(Y_rep)),
    all(is.finite(X_rep))
  )
  
  cat(
    "Fitting 39-stock DSP-MFSV-CAPM\n"
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
        c(39, 39, 1257)
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
  ## TASKS 2--40:
  ## Ordinary MFSV on residual blocks 1--39
  ############################################################
  
} else {
  
  library(factorstochvol)
  
  g <- task_id - 1L
  
  Y_block <-
    input$residual_blocks[[g]]
  
  n_g <- ncol(Y_block)
  
  stopifnot(
    nrow(Y_block) == 1257,
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

setwd("~/Downloads")

# Check worker syntax
parse(file = "large_cross_section_fit_2020.R")
cat("SYNTAX OK\n")

# Check cluster input
x <- readRDS(
  "large_cross_section_cluster_input_2020.rds"
)

stopifnot(
  x$K == 39L,
  length(x$dates) == 1257L,
  nrow(x$universe) == 1000L,
  all(dim(x$representative_excess_returns) == c(1257, 39)),
  all(dim(x$A) == c(1000, 39)),
  all(dim(x$residuals) == c(1257, 1000)),
  length(x$residual_blocks) == 39L,
  sum(lengths(x$residual_block_indices)) == 961L,
  all(
    vapply(
      x$residual_blocks,
      nrow,
      integer(1)
    ) == 1257L
  ),
  all(
    vapply(
      x$residual_blocks,
      function(z) all(is.finite(z)),
      logical(1)
    )
  )
)

cat("ALL CHECKS PASSED\n")