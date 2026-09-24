##################################################################################
#######################Code for 2008 example#####################################
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
set.seed(1)
scale_to_range <- function(x, a, b) {
#####function that scales a time series to lie in the range from a to b#####
  a + ((x - min(x)) * (b - a)) / (max(x) - min(x))
}
#### prepare data ####

setwd("~/Downloads")

input_1000 <- readRDS("large_cross_section_1000_input.rds")
ff_dates <- input_1000$dates
universe_1000 <-input_1000$universe
daily_1000_excess_returns <-input_1000$excess_returns
daily_data <-read.csv("F-F_Research_Data_Factors_daily.csv")
daily_data <-daily_data[21122:22127, ]
ff_dates_check <- as.Date(as.character(daily_data$X),format = "%Y%m%d")
stopifnot(identical(as.Date(ff_dates),ff_dates_check),
  all(dim(daily_1000_excess_returns) ==c(1006, 1000)),
  nrow(universe_1000) == 1000)
############################################################
## STEP 2: Static CAPM for all 1000 stocks
############################################################

design_matrix <- cbind(Intercept = 1,Market = daily_data$Mkt.RF)
capm_coef <- qr.solve(design_matrix,daily_1000_excess_returns)

############################################################
## STEP 3: T x 1000 residual matrix U
############################################################

U <- daily_1000_excess_returns -design_matrix %*% capm_coef
stopifnot(nrow(U) == 1006,ncol(U) == 1000,all(is.finite(U)))

############################################################
## STEP 4: 1000 x 1000 residual correlation matrix
############################################################

R_U <- cor(U)
stopifnot(all(dim(R_U) == c(1000, 1000)))

############################################################
## STEP 5: correlation distance matrix
############################################################

R_U <- pmax(pmin(R_U, 1),-1)
D_U <- sqrt(2 * (1 - R_U))
diag(D_U) <- 0
D_U <- as.dist(D_U)

############################################################
## STEP 6: k-medoids for K = 30,...,50
############################################################

library(cluster)
K_grid <- 30:50
pam_fits <- vector("list",length(K_grid))
silhouette_width <- numeric(length(K_grid))
for (j in seq_along(K_grid)) {
  cat("Fitting K =",K_grid[j],"\n")
  pam_fits[[j]] <- pam(D_U,k = K_grid[j],diss = TRUE)
  silhouette_width[j] <-pam_fits[[j]]$silinfo$avg.width
}
K_results <- data.frame(K = K_grid,Average_Silhouette = silhouette_width)
print(K_results)
plot(K_results$K,K_results$Average_Silhouette,type = "b",xlab = "K",ylab = "Average silhouette width")

############################################################
## STEP 7: Freeze and extract medoids
############################################################

best_idx <- which.max(K_results$Average_Silhouette)
K_star <- K_results$K[best_idx]
cat("Selected K:", K_star, "\n")
cat("Average silhouette width:",K_results$Average_Silhouette[best_idx],"\n")
pam_star <- pam_fits[[best_idx]]
medoid_idx <- pam_star$id.med
stopifnot(length(medoid_idx) == K_star,length(unique(medoid_idx)) == K_star)
representative_permnos <-universe_1000$PERMNO[medoid_idx]
representative_tickers <-universe_1000$Ticker[medoid_idx]
representatives <- data.frame(
  Cluster = seq_len(K_star),
  Index = medoid_idx,
  PERMNO = representative_permnos,
  Ticker = representative_tickers,
  MthCap = universe_1000$MthCap[medoid_idx]
)
print(representatives)

############################################################
## Cluster membership for all 1000 stocks
############################################################

cluster_membership <- pam_star$clustering
stopifnot(length(cluster_membership) == 1000)
cluster_sizes <- table(cluster_membership)
print(cluster_sizes)
cat("\nMinimum cluster size:",min(cluster_sizes),"\nMaximum cluster size:",
  max(cluster_sizes),"\nMedian cluster size:",median(cluster_sizes),"\n")

############################################################
## Representative-stock excess returns
############################################################

representative_excess_returns <-daily_1000_excess_returns[,medoid_idx,drop = FALSE]
colnames(representative_excess_returns) <-as.character(representative_permnos)
stopifnot(all(dim(representative_excess_returns) == c(1006, K_star)))
dim(representative_excess_returns)
saveRDS(list(K = K_star,pam = pam_star,medoid_idx = medoid_idx,
    representatives = representatives,cluster_membership = cluster_membership,
    representative_excess_returns =representative_excess_returns),
  "large_cross_section_kmedoids_selection.rds")

############################################################
## STEP 8: Relate all 1000 stocks to the K_star representatives
############################################################

medoid_cluster <- as.integer(cluster_membership[medoid_idx])
stopifnot(length(unique(medoid_cluster)) == K_star,all(sort(medoid_cluster) == seq_len(K_star)))
medoid_idx_by_cluster <- medoid_idx[match(seq_len(K_star), medoid_cluster)]
stopifnot(all(cluster_membership[medoid_idx_by_cluster] ==seq_len(K_star)))
representative_excess_returns <-daily_1000_excess_returns[,medoid_idx_by_cluster,drop = FALSE]
representative_permnos <-universe_1000$PERMNO[medoid_idx_by_cluster]
colnames(representative_excess_returns) <-as.character(representative_permnos)

Z_rep <- cbind(Intercept = 1,representative_excess_returns)
ols_coef <- qr.solve(Z_rep,daily_1000_excess_returns)
alpha <- as.numeric(ols_coef[1, ])
A <- t(ols_coef[-1,,drop = FALSE])
stopifnot(all(dim(A) == c(1000, K_star)))
rownames(A) <-as.character(universe_1000$PERMNO)
colnames(A) <-as.character(representative_permnos)

alpha[medoid_idx_by_cluster] <- 0
A[medoid_idx_by_cluster,] <- diag(K_star)

fitted_from_reps <-representative_excess_returns %*% t(A)
fitted_from_reps <- sweep(fitted_from_reps,2,alpha,"+")
E <- daily_1000_excess_returns -fitted_from_reps
E[, medoid_idx_by_cluster] <- 0
colnames(E) <-as.character(universe_1000$PERMNO)
stopifnot(all(dim(E) == c(1006, 1000)),all(E[, medoid_idx_by_cluster] == 0))

normal_eq_max <- max(abs(crossprod(Z_rep,E)))
cat("Maximum absolute OLS normal-equation residual:",normal_eq_max,"\n")

############################################################
## STEP 9: Form the K_star residual blocks
############################################################
residual_block_indices <- lapply(seq_len(K_star),function(g) {
    members <- which(cluster_membership == g)
    setdiff(members,medoid_idx_by_cluster[g])
  })
block_sizes <- lengths(residual_block_indices)
cat("Number of residual blocks:",length(block_sizes),"\n")
cat("Minimum residual-block size:",min(block_sizes),"\n")
cat("Maximum residual-block size:",max(block_sizes),"\n")
cat("Median residual-block size:",median(block_sizes),"\n")
print(block_sizes)
residual_blocks <- lapply(residual_block_indices,function(idx) {
    E[,idx,drop = FALSE]
  })
stopifnot( all(vapply(residual_blocks,nrow,integer(1)) == 1006))

representatives_by_cluster <- data.frame(Cluster = seq_len(K_star),
  Index = medoid_idx_by_cluster,PERMNO = universe_1000$PERMNO[medoid_idx_by_cluster],
  Ticker = universe_1000$Ticker[medoid_idx_by_cluster],
  MthCap = universe_1000$MthCap[medoid_idx_by_cluster])
saveRDS(list(dates = ff_dates,market_excess_return =daily_data$Mkt.RF,
    universe =universe_1000,K =K_star,representatives =representatives_by_cluster,
    medoid_idx =medoid_idx_by_cluster,cluster_membership =cluster_membership,
    representative_excess_returns =representative_excess_returns,
    alpha =alpha,A =A,residuals =E,residual_block_indices =residual_block_indices,
    residual_blocks =residual_blocks),"large_cross_section_cluster_input.rds")
cat("\nCreated large_cross_section_cluster_input.rds\n")