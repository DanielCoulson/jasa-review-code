###########################################
############Code for 2020 example##########
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

scale_to_range <- function(x, a, b) {
#####function that scales a time series to lie in the range from a to b#####
  a + ((x - min(x)) * (b - a)) / (max(x) - min(x))
}

#### prepare data ####
setwd("~/Desktop/Graduate Work /Research /Project 1/Code")
tech_stocks =  c("AAPL", "MSFT", "NVDA", "AMZN", "GOOG", "META", "TSM", "AVGO", "TSLA", "ORCL","ASML", "NFLX", "SAP", "AMD", "CRM", "ADBE", "CSCO","IBM" , "QCOM","INTU", "TXN","NOW" , "BABA", "PANW", "ADI", "MU", "LRCX", "INFY","DELL", "INTC")
daily_data = read.csv("F-F_Research_Data_Factors_daily.csv")
daily_data = daily_data[24393:25649,]
design_matrix = matrix(1,nrow = 1257, ncol = 2)
design_matrix[,2] = daily_data$Mkt.RF
daily_tech_data = lapply(tech_stocks, function(symbol){
  stock_data <- getSymbols(symbol, src = "yahoo", from = "2019-01-02", to = "2023-12-31", auto.assign = FALSE)
  return(stock_data)
})
daily_tech_excess_returns = matrix(nrow = (nrow(daily_tech_data[[1]])-1), ncol = length(tech_stocks))
for(i in 1:length(tech_stocks)){
  daily_tech_excess_returns[,i] = 100*diff( log( daily_tech_data[[i]][,6] ) )[-1] - daily_data[,5]
}
diversify_stocks = c("AAPL", "GOOGL", "AMZN", "KO", "ABBV",
                     "PGR", "CAT", "XOM", "LIN", "AMT", 
                     "BOX", "UVV", "SIG", "COTY", "GKOS",
                     "RDN", "ALK", "MGY", "AA", "NHI", 
                     "MITK", "GTN", "CAL", "IAC", "KIDS",
                     "HMN", "WNC", "NBR", "SCL", "DEA"
)
daily_diversify_data = lapply(diversify_stocks, function(symbol){
  stock_data <- getSymbols(symbol, src = "yahoo", from = "2019-01-02", to = "2023-12-31", auto.assign = FALSE)
  return(stock_data)
})
daily_diversify_excess_returns = matrix(nrow = (nrow(daily_diversify_data[[1]])-1), ncol = length(diversify_stocks))
for(i in 1:length(diversify_stocks)){
  daily_diversify_excess_returns[,i] = 100*diff( log( daily_diversify_data[[i]][,6] ) )[-1] - daily_data[,5]
}

#### Fit DSP-MFSV models ####
dspmfsv_daily_tech = DSP_MFSV(y = daily_tech_excess_returns, X = design_matrix, nsave = 3000,nburn = 1500, number_of_latent_factors = 3)
dspmfsv_daily_diversifys = DSP_MFSV(y = daily_diversify_excess_returns, X = design_matrix, nsave = 3000,nburn = 1500, number_of_latent_factors = 3)

#### Compare the correlations in each portfolio ####
df <- data.frame(date = as.Date(as.character(daily_data$X), format = "%Y%m%d"),
                 Diversified = colMeans(t(dspmfsv_daily_diversifys$scores)),Technology   = colMeans(t(dspmfsv_daily_tech$scores))
)
df_long <- pivot_longer(df, cols = c("Diversified", "Technology"),
                        names_to = "Series",values_to = "Mean")
ggplot(df_long, aes(x = date, y = Mean, color = Series)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Diversified" = "#000000",  
                                "Technology"   = "#D55E00")) +  
  labs(x = "Day", y = "Score", color = "Series") +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18)
  )

#### Compared diversified score with VIX ####
VIX_data_1 = getSymbols("^VIX",src = "yahoo", from = "2019-01-03", to = "2023-12-29", auto.assign = FALSE)
VIX_data_1 = na.omit(VIX_data_1)
common_dates = intersect(index(100*diff( log( daily_diversify_data[[1]][,6] ) )[-1]), index(VIX_data_1))
VIX_data_1_actual = VIX_data_1[as.Date(common_dates)]
VIX_data_1_actual[,6] = scale_to_range(VIX_data_1_actual[,6],min(colMeans(t(dspmfsv_daily_diversifys$scores))),max(colMeans(t(dspmfsv_daily_diversifys$scores))))
common_dates <- as.Date(common_dates) 
df <- data.frame(
  date = common_dates,
  mean = colMeans(t(dspmfsv_daily_diversifys$scores))[which(index(100 * diff(log(daily_diversify_data[[1]][, 6]))[-1]) %in% common_dates)],
  lower = hdi(t(dspmfsv_daily_diversifys$scores))[1, which(index(100 * diff(log(daily_diversify_data[[1]][, 6]))[-1]) %in% common_dates)],
  upper = hdi(t(dspmfsv_daily_diversifys$scores))[2, which(index(100 * diff(log(daily_diversify_data[[1]][, 6]))[-1]) %in% common_dates)],
  VIX.Adjusted = VIX_data_1_actual[, 6]
)
ggplot(df, aes(x = date)) +
  geom_line(aes(y = mean), color = "black", size = 1) +  
  geom_line(aes(y = VIX.Adjusted), color = "red", size = 1) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = "blue", alpha = 0.3) +  
  labs(x = "Day", y = "Score") +
  theme_minimal(base_size = 18) +  
  theme(
    axis.title = element_text(size = 18),      
    axis.text = element_text(size = 18)  )