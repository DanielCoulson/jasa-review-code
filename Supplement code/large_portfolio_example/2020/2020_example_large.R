setwd("~/Downloads")
library(dplyr)

############################################################
## 2020 large cross-section: select final 1000 stocks
############################################################

crsp_dec <- read.csv("crsp_dec2018.csv")
crsp_daily <- read.csv("crsp_daily_candidates_2020.csv")

daily_data <- read.csv(
  "F-F_Research_Data_Factors_daily.csv"
)

# Exact sample used in original 2020 example
daily_data <- daily_data[24393:25649, ]

ff_dates <- as.Date(
  as.character(daily_data$X),
  format = "%Y%m%d"
)

stopifnot(length(ff_dates) == 1257)


############################################################
## Prepare CRSP daily data
############################################################

crsp_daily$DlyCalDt <- as.Date(
  crsp_daily$DlyCalDt
)

crsp_daily$DlyRet <- suppressWarnings(
  as.numeric(crsp_daily$DlyRet)
)

crsp_daily$PERMNO <- as.numeric(
  crsp_daily$PERMNO
)


############################################################
## Check/collapse duplicate PERMNO-date observations
############################################################

dup_check <- crsp_daily %>%
  group_by(PERMNO, DlyCalDt) %>%
  filter(n() > 1) %>%
  summarise(
    N_Distinct_Returns =
      n_distinct(DlyRet[is.finite(DlyRet)]),
    .groups = "drop"
  )

conflicts <- dup_check %>%
  filter(N_Distinct_Returns > 1)

cat(
  "Conflicting duplicate observations:",
  nrow(conflicts),
  "\n"
)

stopifnot(nrow(conflicts) == 0)

crsp_daily <- crsp_daily %>%
  group_by(PERMNO, DlyCalDt) %>%
  summarise(
    DlyRet = {
      x <- DlyRet[is.finite(DlyRet)]
      
      if (length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },
    .groups = "drop"
  )


############################################################
## Require complete coverage of all 1257 FF dates
############################################################

candidate_ids <- unique(
  crsp_daily$PERMNO
)

complete_check <- sapply(
  candidate_ids,
  function(id) {
    
    z <- crsp_daily[
      crsp_daily$PERMNO == id,
      c("DlyCalDt", "DlyRet")
    ]
    
    idx <- match(
      ff_dates,
      z$DlyCalDt
    )
    
    if (anyNA(idx)) {
      return(FALSE)
    }
    
    r <- z$DlyRet[idx]
    
    all(
      is.finite(r) &
        r > -1
    )
  }
)

complete_permnos <- candidate_ids[
  complete_check
]

cat(
  "Stocks with complete 2019-2023 coverage:",
  length(complete_permnos),
  "\n"
)


############################################################
## Reapply December 2018 eligibility filters
############################################################

crsp_dec$MthCap <- as.numeric(
  crsp_dec$MthCap
)

crsp_dec$PERMNO <- as.numeric(
  crsp_dec$PERMNO
)

eligible_dec <- subset(
  crsp_dec,
  ShareType == "NS" &
    SecurityType == "EQTY" &
    SecuritySubType == "COM" &
    USIncFlg == "Y" &
    IssuerType %in% c("ACOR", "CORP") &
    PrimaryExch %in% c("N", "A", "Q") &
    ConditionalType == "RW" &
    TradingStatusFlg == "A" &
    is.finite(MthCap) &
    MthCap > 0
)

universe_1000 <- eligible_dec[
  eligible_dec$PERMNO %in% complete_permnos,
]

universe_1000 <- universe_1000[
  order(
    universe_1000$MthCap,
    decreasing = TRUE
  ),
]

universe_1000 <- head(
  universe_1000,
  1000
)

cat(
  "Final portfolio size:",
  nrow(universe_1000),
  "\n"
)

stopifnot(
  nrow(universe_1000) == 1000
)

print(
  universe_1000[
    1:20,
    c("PERMNO", "Ticker", "MthCap")
  ]
)


############################################################
## Construct 1257 x 1000 excess-return matrix
############################################################

daily_1000_returns <- sapply(
  universe_1000$PERMNO,
  function(id) {
    
    z <- crsp_daily[
      crsp_daily$PERMNO == id,
      c("DlyCalDt", "DlyRet")
    ]
    
    idx <- match(
      ff_dates,
      z$DlyCalDt
    )
    
    stopifnot(!anyNA(idx))
    
    100 * log1p(
      z$DlyRet[idx]
    )
  }
)

colnames(daily_1000_returns) <-
  as.character(universe_1000$PERMNO)

daily_1000_excess_returns <- sweep(
  daily_1000_returns,
  1,
  daily_data$RF,
  "-"
)

stopifnot(
  all(
    dim(daily_1000_excess_returns) ==
      c(1257, 1000)
  ),
  all(is.finite(daily_1000_excess_returns))
)

print(
  dim(daily_1000_excess_returns)
)

saveRDS(
  list(
    dates = ff_dates,
    universe = universe_1000,
    excess_returns =
      daily_1000_excess_returns
  ),
  "large_cross_section_1000_input_2020.rds"
)

cat(
  "\n2020 LARGE-CROSS-SECTION INPUT CREATED\n"
)

##############################


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
#### prepare data ####

setwd("~/Downloads")

input_1000 <- readRDS(
  "large_cross_section_1000_input_2020.rds"
)

ff_dates <- input_1000$dates

universe_1000 <-
  input_1000$universe

daily_1000_excess_returns <-
  input_1000$excess_returns

daily_data <-
  read.csv(
    "F-F_Research_Data_Factors_daily.csv"
  )

daily_data <-
  daily_data[24393:25649, ]

ff_dates_check <- as.Date(
  as.character(daily_data$X),
  format = "%Y%m%d"
)

stopifnot(
  identical(
    as.Date(ff_dates),
    ff_dates_check
  ),
  all(
    dim(daily_1000_excess_returns) ==
      c(1257, 1000)
  ),
  nrow(universe_1000) == 1000
)







############################################################
## STEP 2: Static CAPM for all 1000 stocks
############################################################

design_matrix <- cbind(
  Intercept = 1,
  Market = daily_data$Mkt.RF
)

capm_coef <- qr.solve(
  design_matrix,
  daily_1000_excess_returns
)


############################################################
## STEP 3: T x 1000 residual matrix U
############################################################

U <- daily_1000_excess_returns -
  design_matrix %*% capm_coef

stopifnot(
  nrow(U) == 1257,
  ncol(U) == 1000,
  all(is.finite(U))
)


############################################################
## STEP 4: 1000 x 1000 residual correlation matrix
############################################################

R_U <- cor(U)

stopifnot(
  all(dim(R_U) == c(1000, 1000))
)


############################################################
## STEP 5: correlation distance matrix
############################################################

R_U <- pmax(
  pmin(R_U, 1),
  -1
)

D_U <- sqrt(
  2 * (1 - R_U)
)

diag(D_U) <- 0

D_U <- as.dist(D_U)


############################################################
## STEP 6: k-medoids for K = 30,...,50
############################################################

library(cluster)

K_grid <- 30:50

pam_fits <- vector(
  "list",
  length(K_grid)
)

silhouette_width <- numeric(
  length(K_grid)
)

for (j in seq_along(K_grid)) {
  
  cat(
    "Fitting K =",
    K_grid[j],
    "\n"
  )
  
  pam_fits[[j]] <- pam(
    D_U,
    k = K_grid[j],
    diss = TRUE
  )
  
  silhouette_width[j] <-
    pam_fits[[j]]$silinfo$avg.width
}

K_results <- data.frame(
  K = K_grid,
  Average_Silhouette = silhouette_width
)

print(K_results)

plot(
  K_results$K,
  K_results$Average_Silhouette,
  type = "b",
  xlab = "K",
  ylab = "Average silhouette width"
)

############################################################
## STEP 7: Freeze selected K and extract medoids/clusters
############################################################

best_idx <- which.max(
  K_results$Average_Silhouette
)

K_star <- K_results$K[best_idx]

cat(
  "Selected K:",
  K_star,
  "\n"
)

cat(
  "Average silhouette width:",
  K_results$Average_Silhouette[best_idx],
  "\n"
)

# Freeze the value selected by the 2020 data
stopifnot(K_star == 39L)

pam_star <- pam_fits[[best_idx]]

# Indices of the medoid stocks within the 1000-stock universe
medoid_idx <- pam_star$id.med

stopifnot(
  length(medoid_idx) == K_star,
  length(unique(medoid_idx)) == K_star
)

# PERMNOs of representative stocks
representative_permnos <-
  universe_1000$PERMNO[
    medoid_idx
  ]

# Corresponding December-2018 tickers
representative_tickers <-
  universe_1000$Ticker[
    medoid_idx
  ]

representatives <- data.frame(
  Cluster = seq_len(K_star),
  Index = medoid_idx,
  PERMNO = representative_permnos,
  Ticker = representative_tickers,
  MthCap =
    universe_1000$MthCap[
      medoid_idx
    ]
)

print(representatives)


############################################################
## Cluster membership for all 1000 stocks
############################################################

cluster_membership <-
  pam_star$clustering

stopifnot(
  length(cluster_membership) == 1000
)

cluster_sizes <-
  table(cluster_membership)

print(cluster_sizes)

cat(
  "\nMinimum cluster size:",
  min(cluster_sizes),
  "\nMaximum cluster size:",
  max(cluster_sizes),
  "\nMedian cluster size:",
  median(cluster_sizes),
  "\n"
)


############################################################
## Representative-stock excess returns
############################################################

representative_excess_returns <-
  daily_1000_excess_returns[
    ,
    medoid_idx,
    drop = FALSE
  ]

colnames(
  representative_excess_returns
) <- as.character(
  representative_permnos
)

stopifnot(
  all(
    dim(
      representative_excess_returns
    ) == c(1257, K_star)
  )
)

print(
  dim(
    representative_excess_returns
  )
)

saveRDS(
  list(
    K = K_star,
    pam = pam_star,
    medoid_idx = medoid_idx,
    representatives = representatives,
    cluster_membership =
      cluster_membership,
    representative_excess_returns =
      representative_excess_returns
  ),
  "large_cross_section_kmedoids_selection_2020.rds"
)


############################################################
## STEP 8: Relate all 1000 stocks to representatives
############################################################

# Make sure medoids are ordered by cluster label
medoid_cluster <- as.integer(
  cluster_membership[
    medoid_idx
  ]
)

stopifnot(
  length(
    unique(medoid_cluster)
  ) == K_star,
  all(
    sort(medoid_cluster) ==
      seq_len(K_star)
  )
)

medoid_idx_by_cluster <-
  medoid_idx[
    match(
      seq_len(K_star),
      medoid_cluster
    )
  ]

stopifnot(
  all(
    cluster_membership[
      medoid_idx_by_cluster
    ] == seq_len(K_star)
  )
)


############################################################
## Representative returns ordered by cluster
############################################################

representative_excess_returns <-
  daily_1000_excess_returns[
    ,
    medoid_idx_by_cluster,
    drop = FALSE
  ]

representative_permnos <-
  universe_1000$PERMNO[
    medoid_idx_by_cluster
  ]

colnames(
  representative_excess_returns
) <- as.character(
  representative_permnos
)


############################################################
## Multivariate OLS
##
## r_t = alpha + A r_rep,t + e_t
############################################################

Z_rep <- cbind(
  Intercept = 1,
  representative_excess_returns
)

ols_coef <- qr.solve(
  Z_rep,
  daily_1000_excess_returns
)

# alpha: length 1000
alpha <- as.numeric(
  ols_coef[1, ]
)

# A: 1000 x K_star
A <- t(
  ols_coef[
    -1,
    ,
    drop = FALSE
  ]
)

stopifnot(
  all(
    dim(A) ==
      c(1000, K_star)
  )
)

rownames(A) <-
  as.character(
    universe_1000$PERMNO
  )

colnames(A) <-
  as.character(
    representative_permnos
  )


############################################################
## Enforce exact identity for representative stocks
############################################################

alpha[
  medoid_idx_by_cluster
] <- 0

A[
  medoid_idx_by_cluster,
] <- diag(K_star)


############################################################
## Construct residual matrix e_t
############################################################

fitted_from_reps <-
  representative_excess_returns %*%
  t(A)

fitted_from_reps <- sweep(
  fitted_from_reps,
  2,
  alpha,
  "+"
)

E <-
  daily_1000_excess_returns -
  fitted_from_reps

# Representatives have exactly zero residual
E[
  ,
  medoid_idx_by_cluster
] <- 0

colnames(E) <-
  as.character(
    universe_1000$PERMNO
  )

stopifnot(
  all(
    dim(E) ==
      c(1257, 1000)
  ),
  all(
    E[
      ,
      medoid_idx_by_cluster
    ] == 0
  )
)


############################################################
## Diagnostic: OLS orthogonality
############################################################

normal_eq_max <- max(
  abs(
    crossprod(
      Z_rep,
      E
    )
  )
)

cat(
  "Maximum absolute OLS normal-equation residual:",
  normal_eq_max,
  "\n"
)


############################################################
## STEP 9: Form residual blocks
############################################################

residual_block_indices <- lapply(
  seq_len(K_star),
  function(g) {
    
    members <- which(
      cluster_membership == g
    )
    
    # Remove that cluster's medoid
    setdiff(
      members,
      medoid_idx_by_cluster[g]
    )
  }
)

block_sizes <- lengths(
  residual_block_indices
)

cat(
  "Number of residual blocks:",
  length(block_sizes),
  "\n"
)

cat(
  "Minimum residual-block size:",
  min(block_sizes),
  "\n"
)

cat(
  "Maximum residual-block size:",
  max(block_sizes),
  "\n"
)

cat(
  "Median residual-block size:",
  median(block_sizes),
  "\n"
)

print(block_sizes)

stopifnot(
  length(
    residual_block_indices
  ) == K_star,
  sum(block_sizes) ==
    1000 - K_star
)


############################################################
## Actual residual matrices
############################################################

residual_blocks <- lapply(
  residual_block_indices,
  function(idx) {
    
    E[
      ,
      idx,
      drop = FALSE
    ]
  }
)

stopifnot(
  all(
    vapply(
      residual_blocks,
      nrow,
      integer(1)
    ) == 1257
  )
)


############################################################
## Save cluster input
############################################################

representatives_by_cluster <-
  data.frame(
    Cluster =
      seq_len(K_star),
    
    Index =
      medoid_idx_by_cluster,
    
    PERMNO =
      universe_1000$PERMNO[
        medoid_idx_by_cluster
      ],
    
    Ticker =
      universe_1000$Ticker[
        medoid_idx_by_cluster
      ],
    
    MthCap =
      universe_1000$MthCap[
        medoid_idx_by_cluster
      ]
  )

saveRDS(
  list(
    dates = ff_dates,
    
    market_excess_return =
      daily_data$Mkt.RF,
    
    universe =
      universe_1000,
    
    K =
      K_star,
    
    representatives =
      representatives_by_cluster,
    
    medoid_idx =
      medoid_idx_by_cluster,
    
    cluster_membership =
      cluster_membership,
    
    representative_excess_returns =
      representative_excess_returns,
    
    alpha =
      alpha,
    
    A =
      A,
    
    residuals =
      E,
    
    residual_block_indices =
      residual_block_indices,
    
    residual_blocks =
      residual_blocks
  ),
  "large_cross_section_cluster_input_2020.rds"
)

cat(
  "\nCreated large_cross_section_cluster_input_2020.rds\n"
)