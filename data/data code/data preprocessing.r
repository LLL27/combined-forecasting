library(dplyr)
library(lubridate)
library(seasonal)
library(seasonalview)
library(TSA)

rm(list=ls())

getwd()
setwd("/Users/liumeishao/paper")
source("code/data process funcs.R")
monthly <- read.csv("raw data/all_monthly_1031.csv")
Names_CN <- colnames(monthly)

Name_ENG <- c("time", "Exp_YoY", "Imp_YoY", "Exp", "Imp",
              "Trade", "MEC", "CFCI", "CMI", 
              "RECI", "PMI_man", "PMI", 
              "CPI_YoY", "CPI_MoM", "CPI_food_YoY",
              "CPI_food_MoM", "CPI_health_YoY", "CPI_health_MoM", "CPI_trans_YoY", "CPI_trans_MoM", 
              "PPI_YoY", "PPI_MoM", "PPI_pro_YoY", "PPI_pro_MoM", "PPI_ex_YoY", "PPI_ex_MoM", "PPI_raw_YoY", "PPI_raw_MoM", 
              "PPI_fac_YoY", "PPI_fac_MoM", "PPI_liv_YoY", "PPI_liv_MoM", "PPI_food_YoY", "PPI_food_MoM", "PPI_clo_YoY", "PPI_clo_MoM", 
              "PPI_gen_YoY", "PPI_gen_MoM", "PPI_dur_YoY", "PPI_dur_MoM", "IPI_YoY", "IPI_MoM", 
              "Commodity", "Consumption", "NEERI", "SHPE", "SZPE", 
              "Deposit", "Loan", "M0_YoY", "M1_YoY", "M2_YoY", "M0", "M1", "M2")
cbind(Names_CN, Name_ENG)
colnames(monthly) <- Name_ENG
# month starts from 2000-01
monthly <- data.frame(monthly[c(1:308), ])
for(i in 2:ncol(monthly)){
  monthly[, i] <- as.numeric(monthly[, i])
}

#Obtain index value from YoY and MoM values


monthly <- monthly %>% mutate(CPI = get_index(Time = time, YoY = CPI_YoY,
                                              MoM = CPI_MoM, index_month = "2011-01"),
                              CPI_food = get_index(Time = time, YoY = CPI_food_YoY,
                                                   MoM = CPI_food_MoM, index_month = "2011-01"),
                              CPI_health = get_index(Time = time, YoY = CPI_health_YoY,
                                                     MoM = CPI_health_MoM, index_month = "2011-01"),
                              CPI_trans = get_index(Time = time, YoY = CPI_trans_YoY,
                                                    MoM = CPI_trans_MoM, index_month = "2011-01"),
                              PPI = get_index(Time = time, YoY = PPI_YoY,
                                              MoM = PPI_MoM, index_month = "2011-01"),
                              PPI_pro = get_index(Time = time, YoY = PPI_pro_YoY,
                                                  MoM = PPI_pro_MoM, index_month = "2011-01"),
                              PPI_ex = get_index(Time = time, YoY = PPI_ex_YoY,
                                                 MoM = PPI_ex_MoM, index_month = "2011-01"),
                              PPI_raw = get_index(Time = time, YoY = PPI_raw_YoY,
                                                  MoM = PPI_raw_MoM, index_month = "2011-01"),
                              PPI_fac = get_index(Time = time, YoY = PPI_fac_YoY,
                                                  MoM = PPI_fac_MoM, index_month = "2011-01"),
                              PPI_liv = get_index(Time = time, YoY = PPI_liv_YoY,
                                                  MoM = PPI_liv_MoM, index_month = "2011-01"),
                              PPI_food = get_index(Time = time, YoY = PPI_food_YoY,
                                                   MoM = PPI_food_MoM, index_month = "2011-01"),
                              PPI_clo = get_index(Time = time, YoY = PPI_clo_YoY,
                                                  MoM = PPI_clo_MoM, index_month = "2011-01"),
                              PPI_gen = get_index(Time = time, YoY = PPI_gen_YoY,
                                                  MoM = PPI_gen_MoM, index_month = "2011-01"),
                              PPI_dur = get_index(Time = time, YoY = PPI_dur_YoY,
                                                  MoM = PPI_dur_MoM, index_month = "2011-01"),
                              IPI = get_index(Time = time, YoY = IPI_YoY,
                                              MoM = IPI_MoM, index_month = "2011-01")) %>%
  select(-c(CPI_YoY, CPI_MoM, CPI_food_YoY, CPI_food_MoM, CPI_health_YoY, CPI_health_MoM, CPI_trans_YoY, CPI_trans_MoM, 
            PPI_YoY, PPI_MoM, PPI_pro_YoY, PPI_pro_MoM, PPI_ex_YoY, PPI_ex_MoM, PPI_raw_YoY, PPI_raw_MoM, 
            PPI_fac_YoY, PPI_fac_MoM, PPI_liv_YoY, PPI_liv_MoM, PPI_food_YoY, PPI_food_MoM, PPI_clo_YoY, PPI_clo_MoM, 
            PPI_gen_YoY, PPI_gen_MoM, PPI_dur_YoY, PPI_dur_MoM, IPI_YoY, IPI_MoM))

#Impute RECI using ARIMA； Commodity, Consumption in Jan and Feb using Seasonal ARIMA (S-ARIMA)


#RECI
RECI_log <- ts(log(monthly$RECI), start = c(2000, 1), frequency = 12)
plot(RECI_log)
na_ind <- c(1:length(RECI_log))[is.na(RECI_log)]
acf(diff(RECI_log[1:(na_ind[1] - 1)]), lag = 48)
pacf(diff(RECI_log[1:(na_ind[1] - 1)]), lag = 48)
## After 1st order differece , the acf is significant nonzero at lag = 1 and the pacf is significant nonzero at lag = 1,so ARIMA(1,1,1) is suitable
RECI_log_imputed <- RECI_log
for(i in na_ind){
  mdl <- arima(RECI_log_imputed[1:(i-1)], order = c(1, 1, 1))
  mdl_pred <- predict(mdl, n.ahead = 1)
  RECI_log_imputed[i] <- mdl_pred$pred[1]
}
plot(RECI_log_imputed, col = "red")
lines(RECI_log)
legend("bottomright", c("original", "imputed"), 
       lty = c(1,1), col = c("black", "red"))
monthly$RECI_imputed <- exp(RECI_log_imputed)


#Commodity
Commodity_log <- ts(log(monthly$Commodity), start = c(2000, 1), frequency = 12)
plot(Commodity_log)
na_ind <- c(1:length(Commodity_log))[is.na(Commodity_log)]
acf(diff(diff(Commodity_log[121:(na_ind[121] - 1)]), lag = 12), lag = 48)
pacf(diff(diff(Commodity_log[121:(na_ind[121] - 1)]), lag = 12), lag = 48)
## After 1st order differece and 12nd order difference, the acf is significant nonzero at lag = 1, 11,12,13,34,35 so S-ARIMA(0,1,1)x(0,1,1)_{12} is suitable2
Commodity_log_imputed <- Commodity_log
for(i in na_ind[121:148]){
  mdl <- arima(Commodity_log_imputed[1:(i-1)], order = c(0, 1, 0),
               seasonal = c(0, 1, 0))
  mdl_pred <- predict(mdl, n.ahead = 1)
  Commodity_log_imputed[i] <- mdl_pred$pred[1]
}
plot(Commodity_log_imputed, col = "red")
lines(Commodity_log)
legend("bottomright", c("original", "imputed"), 
       lty = c(1,1), col = c("black", "red"))
monthly$Commodity_imputed <- exp(Commodity_log_imputed)


#Consumption
consumption_log <- ts(log(monthly$Consumption), start = c(2000, 1), frequency = 12)
plot(consumption_log)
na_ind <- c(1:length(consumption_log))[is.na(consumption_log)]
acf(diff(diff(consumption_log[1:(na_ind[1] - 1)]), lag = 12), lag = 48)
pacf(diff(diff(consumption_log[1:(na_ind[1] - 1)]), lag = 12), lag = 48)
## After 1st order differece and 12nd order difference, the acf is significant nonzero at lag = 1, 11,12,13,34,35 so S-ARIMA(0,1,1)x(0,1,1)_{12} is suitable2
consumption_log_imputed <- consumption_log
for(i in na_ind){
  mdl <- arima(consumption_log_imputed[1:(i-1)], order = c(0, 1, 1),
               seasonal = c(0, 1, 1))
  mdl_pred <- predict(mdl, n.ahead = 1)
  consumption_log_imputed[i] <- mdl_pred$pred[1]
}
plot(consumption_log_imputed, col = "red")
lines(consumption_log)
legend("bottomright", c("original", "imputed"), 
       lty = c(1,1), col = c("black", "red"))
monthly$Consumption_imputed <- exp(consumption_log_imputed)

head(monthly)
tail(monthly)
write.csv(monthly, "processed data/all_monthly_preprocessing.csv", row.names = FALSE)