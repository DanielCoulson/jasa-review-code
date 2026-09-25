###########################################
############Backtest experiment##########
###########################################
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
computeDIC <- FALSE
scale_to_range <- function(x, a, b) {
  #####function that scales a time series to lie in the range from a to b#####
  a + ((x - min(x)) * (b - a)) / (max(x) - min(x))
}
#### prepare data ####
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 2) {
  as.integer(args[2])
} else {
  NA_integer_
}
project_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  getwd()
}
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")
dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)
source(file.path(project_dir,"DSP-MFSV-CAPM.R"))
source(file.path(project_dir,"simulation_model_fitting_code.R"))

example <- "2008"   # choose "2008" or "2020"
if (example == "2008") {
  diversify_stocks <- c(
    "MSFT", "ADSK", "PLXS",
    "JPM", "FITB", "UMBF",
    "JNJ", "DGX", "CHE",
    "PG", "SJM", "CALM",
    "HD", "DRI", "BKE",
    "UNP", "RHI", "AIR",
    "XOM", "HP", "SM",
    "APD", "RPM", "KWR",
    "DUK", "CMS", "AWR",
    "VZ", "CCOI", "SHEN"
  )
  data_start <- as.Date("2006-01-02")
  data_end   <- as.Date("2011-12-31")
  oos_start  <- as.Date("2008-01-01")
  oos_end    <- as.Date("2011-12-31")
  download_start <- "2005-12-29"
  download_end   <- "2012-01-03"
} else if (example == "2020") {
  diversify_stocks <- c(
    "AAPL", "GOOGL", "AMZN", "KO", "ABBV",
    "PGR", "CAT", "XOM", "LIN", "AMT",
    "BOX", "UVV", "SIG", "COTY", "GKOS",
    "RDN", "ALK", "MGY", "AA", "NHI",
    "MITK", "GTN", "CAL", "PPLI", "KIDS",
    "HMN", "WNC", "NBR", "SCL", "DEA"
  )
  data_start <- as.Date("2018-01-02")
  data_end   <- as.Date("2023-12-31")
  oos_start  <- as.Date("2020-01-01")
  oos_end    <- as.Date("2023-12-31")
  download_start <- "2017-12-28"
  download_end   <- "2024-01-02"
}

estimation_window <- 500L
daily_data_all <-read.csv(file.path(data_dir,"F-F_Research_Data_Factors_daily.csv"))
ff_dates_all <-as.Date(as.character(daily_data_all$X),format = "%Y%m%d")
keep_ff <-!is.na(ff_dates_all) &ff_dates_all >= data_start &ff_dates_all <= data_end
daily_data <-daily_data_all[keep_ff,,drop = FALSE]
ff_dates <-ff_dates_all[keep_ff]
ord <- order(ff_dates)
daily_data <-daily_data[ord,,drop = FALSE]
ff_dates <-ff_dates[ord]

daily_diversify_data <- lapply(diversify_stocks,function(symbol) {
    cat("Downloading", symbol, "...\n")
    getSymbols(symbol,src = "yahoo",from = download_start,to = download_end,
      auto.assign = FALSE)
  })

names(daily_diversify_data) <-diversify_stocks
get_stock_returns <- function(stock_data, symbol) {
  prices <-Ad(stock_data)
  log_ret <-na.omit(100 *diff(log(prices)))
  simple_ret <-na.omit(prices /lag(prices) - 1)
  log_dates <-as.Date(index(log_ret))
  simple_dates <-as.Date(index(simple_ret))
  idx_log <-match(ff_dates,log_dates)
  idx_simple <-match(ff_dates,simple_dates)
  if (anyNA(idx_log) ||anyNA(idx_simple)) {
    bad_dates <-ff_dates[is.na(idx_log) |is.na(idx_simple)]
    stop( paste("Missing observations for",symbol,"- first missing date:",
                bad_dates[1]) )
  }
  list(log_return =as.numeric(log_ret[idx_log]),simple_return =as.numeric(simple_ret[idx_simple]))
}
return_list <-lapply(seq_along(diversify_stocks),
                     function(i) {
                       get_stock_returns(daily_diversify_data[[i]],diversify_stocks[i])
                     }
)
log_return_matrix <-do.call(cbind,lapply(return_list,function(x) x$log_return))
simple_return_matrix <-do.call(cbind,lapply(return_list,function(x) x$simple_return))
colnames(log_return_matrix) <-diversify_stocks
colnames(simple_return_matrix) <-diversify_stocks
daily_diversify_excess_returns <-sweep(log_return_matrix,1,daily_data$RF,"-")
design_matrix <-cbind(Intercept = 1,Market = daily_data$Mkt.RF)
stopifnot(nrow(daily_diversify_excess_returns) ==length(ff_dates),
          nrow(simple_return_matrix) ==length(ff_dates),
  nrow(design_matrix) ==length(ff_dates),ncol(daily_diversify_excess_returns) ==30)
month_id <-format(ff_dates,"%Y-%m")
oos_months <-unique(month_id[ff_dates >= oos_start &ff_dates <= oos_end])
stopifnot(length(oos_months) == 48L)
rebalance_table <-do.call(rbind,lapply(seq_along(oos_months),function(k) {
        this_month <-oos_months[k]
        hold_idx <-which(month_id ==this_month)
        hold_start <-min(hold_idx)
        hold_end <-max(hold_idx)
        train_end <-hold_start -1L
        train_start <-train_end -estimation_window +1L
        if (train_start <1L) {
          stop(paste("Not enough observations for",this_month))
        }
        data.frame(Task = k,Month = this_month,Train_Start_Index =train_start,
          Train_End_Index =train_end,Hold_Start_Index =hold_start,Hold_End_Index =hold_end,
          Train_Start =ff_dates[train_start],Train_End =ff_dates[train_end],
          Hold_Start =ff_dates[hold_start],Hold_End =ff_dates[hold_end],
           N_Train =train_end -train_start +1L,
          N_Hold =length(hold_idx))}))
stopifnot(all(rebalance_table$N_Train ==500L))
print(rebalance_table)
if (is.na(task_id) ||task_id < 1L ||task_id > 48L) {
  stop("Supply OOS month number 1,...,48 as the second command-line argument.")
}

info <-rebalance_table[task_id,,drop = FALSE]
train_idx <-info$Train_Start_Index:info$Train_End_Index
hold_idx <-info$Hold_Start_Index:info$Hold_End_Index
Y_train <-daily_diversify_excess_returns[train_idx,,drop = FALSE]
X_train <-design_matrix[train_idx,,drop = FALSE]
stopifnot(nrow(Y_train) ==500L)
cat("\n========================================\n","OOS MONTH: ",info$Month,
  "\nTRAINING SAMPLE: ", as.character(info$Train_Start)," to ",as.character(info$Train_End),
  "\nHOLDING SAMPLE: ",as.character(info$Hold_Start),
  " to ", as.character(info$Hold_End),"\n========================================\n",
  sep = "")
oos_dir <-file.path(output_dir,"oos_monthly")
dir.create(oos_dir,recursive = TRUE,showWarnings = FALSE)
month_file <-file.path(oos_dir,sprintf("oos_%02d_%s.rds",task_id,gsub("-",
        "",info$Month)))
if (file.exists(month_file)) {
  cat("Loading existing checkpoint:\n",month_file,"\n")
  month_result <-readRDS(month_file)
} else {
  month_result <-list(task_id =task_id,month =info$Month,train_start =info$Train_Start,
      train_end =info$Train_End,hold_dates =ff_dates[hold_idx],
      holding_returns =simple_return_matrix[hold_idx,,drop = FALSE],
      holding_rf =daily_data$RF[hold_idx] /100,
      models =list())
}
