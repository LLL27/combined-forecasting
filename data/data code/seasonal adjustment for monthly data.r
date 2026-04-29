library(dplyr)
library(tsoutliers)
library(TSA)
library(seasonal)
library(seasonalview)
library(zoo)
library(lubridate)

rm(list=ls())
getwd()
setwd("/Users/liumeishao/paper")
all_monthly <- read.csv("processed data/all_monthly_preprocessing.csv", header = T)
colnames(all_monthly)

#We separate them into flow-type data, stock-type data and data which does not require seasonal adjustment:

#flow: Loan, Deposit, Import, Export, Trade, CPI, CPI_food, CPI_health, CPI_trans, PPI, PPI_pro, PPI_ex, PPI_raw, 
#      PPI_fac, PPI_liv, PPI_food, PPI_clo, PPI_gen, PPI_dur, IPI, 
#     Commodity_imputed, Consumption_imputed
#stock: M0, M1, M2, NEERI, SHPE, SZPE
#May not need seasonal adjustment:  MEC, CFCI, CMI, PMI_man, PMI, RECI_imputed

# 1 Seasonal Adjustment for CPI
# CPI starts from 2000-01. From the following plots, CPI is non-stationary and has obvious seasonality. 
# We do seasonal adjustment for it.

source("code/seas with SF.R", local = knitr::knit_global())
cpi <- ts(all_monthly$CPI, start = c(2000, 01), frequency = 12)
plot_acf_diff(vara = cpi, name = "CPI")
par(mfrow = c(1, 1))

# cpi_seas <- seas_flow(vara = cpi)
# saveRDS(cpi_seas, "output/seas output/cpi_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
cpi_seas <- readRDS("output/seas output/cpi_seas.RDS")
cpi_seas$SF_ind
summary(cpi_seas$x13res)
head(cpi_seas$x13res$data)
plot(cpi_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(cpi_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

cpi_final <- cpi_seas$x13res$data[, 1]
cpi_season <- cpi_seas$x13res$data[, 2]
cpi_trend <- cpi_seas$x13res$data[, 4]


source("code/SARIMAX funcs.R", local = knitr::knit_global())
cpi_sarima <- arima(cpi, order = c(3, 1, 1), 
                    seasonal = list(order = c(0, 1, 1), period = 12))
# x13out <- seas(cpi,arima.model = "(3 1 1)(0 1 1)",
#                outlier = NULL, transform.function = "auto",
#                x11 = '', regression.aictest = NULL)
# summary(x13out)
# x13out$data
# plot(x13out, outliers = TRUE, trend = FALSE,
#      main = "CPI: Index Data", transform = c("none"))
# legend("topleft", c("original","seasonal adjusted"),
#        lty=c(1, 1), col = c(1, 2), lw = 2)
# plot(x13out, outliers = TRUE, trend = FALSE,
#      main = "CPI: Percentage Change of Year", transform = c("PCY"))
# legend("topright", c("original", "seasonal adjusted"),
#        lty = c(1, 1), col = c(1, 2),lw = 2)
# h = ts(x, frequency = 12, start = c(start.year, 01))
# import_seas = seas(import, 
#                    xreg = h, 
#                    transform.function = "auto",
#                    outlier = NULL,
#                    regression.aictest = NULL, 
#                    regression.usertype = "holiday",
#                    forecast.maxlead = 0, 
#                    x11 = "",
#                    regression.variables = c("AO1995.Jan"))  # 手动指定异常值

cpi_sarima
cpi_resid <- residuals(cpi_sarima)

infile <- paste("processed data/flow SF regressor/", cpi_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x <- x[1:length(cpi), ]
h = data.frame(x)

cpi_sarimaSF <- arimax(cpi, order = c(3, 1, 1), xreg = h,
                       seasonal = list(order = c(0, 1, 1), period = 12))
cpi_residSF <- residuals(cpi_sarimaSF)
cpi_parsSF <- coefs2poly(cpi_sarimaSF)

plot(cpi_resid, type = "l")
lines(cpi_residSF, col = 2)
abline(0, 0, lty = 6)

cpi_o1SF <- locate.outliers(resid = cpi_residSF, pars = cpi_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(cpi, start_year = 2000)[cpi_o1SF$ind, ], cpi_o1SF)



# 2 Seasonal Adjustment for CPI Food Category

cpi_food <- ts(all_monthly$CPI_food, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(cpi_food), name = "log(CPI_food)")
par(mfrow = c(1, 1))

# cpi_food_seas <- seas_flow(vara = cpi_food)
# saveRDS(cpi_food_seas, "output/seas output/cpi_food_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
cpi_food_seas <- readRDS("output/seas output/cpi_food_seas.RDS")
summary(cpi_food_seas$x13res)
cpi_food_seas$SF_ind
head(cpi_food_seas$x13res$data)

plot(cpi_food_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Food Category: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(cpi_food_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Food Category: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
cpi_food_final <- cpi_food_seas$x13res$data[, 1]
cpi_food_season <- cpi_food_seas$x13res$data[, 2]
cpi_food_trend <- cpi_food_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
cpi_food_sarima <- arima(log(cpi_food), order = c(3,1,1),
                         seasonal = list(order = c(0, 1, 1), period = 12))
cpi_food_sarima

cpi_food_resid <- residuals(cpi_food_sarima)

start.year <- year(as.Date(time(cpi_food)))[1]
infile <- paste("processed data/flow SF regressor/", cpi_food_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(cpi_food)), ] 
h = data.frame(x)


cpi_food_sarimaSF <- arimax(log(cpi_food), order = c(3, 1, 1), xreg = h,
                            seasonal = list(order = c(0, 1, 1), period = 12))
cpi_food_residSF <- residuals(cpi_food_sarimaSF)
cpi_food_parsSF <- coefs2poly(cpi_food_sarimaSF)

plot(cpi_food_resid, type = "l")
lines(cpi_food_residSF, col = 2)
abline(0, 0, lty = 6)

cpi_food_o1SF <- locate.outliers(resid = cpi_food_residSF, pars = cpi_food_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(cpi_food, start_year = 2000)[cpi_food_o1SF$ind, ], cpi_food_o1SF)


# 3 Seasonal Adjustment for CPI Health Category

cpi_health <- ts(all_monthly$CPI_health, start = c(2000, 01), frequency = 12)
plot_acf_diff(cpi_health, name = "CPI_health")
par(mfrow = c(1, 1))

# cpi_health_seas <- seas_flow(vara = cpi_health)
# saveRDS(cpi_health_seas, "output/seas output/cpi_health_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
cpi_health_seas <- readRDS("output/seas output/cpi_health_seas.RDS")
summary(cpi_health_seas$x13res)
cpi_health_seas$SF_ind
head(cpi_health_seas$x13res$data)

plot(cpi_health_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Health Category: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(cpi_health_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Health Category: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
cpi_health_final <- cpi_health_seas$x13res$data[, 1]
cpi_health_season <- cpi_health_seas$x13res$data[, 2]
cpi_health_trend <- cpi_health_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
cpi_health_sarima <- arima(cpi_health, order = c(1, 1, 0),
                         seasonal = list(order = c(0, 1, 1), period = 12))
cpi_health_sarima

cpi_health_resid <- residuals(cpi_health_sarima)

start.year <- year(as.Date(time(cpi_health)))[1]
infile <- paste("processed data/flow SF regressor/", cpi_health_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(cpi_health)), ] 
h = data.frame(x)


cpi_health_sarimaSF <- arimax(cpi_health, order = c(1, 1, 0), xreg = h,
                            seasonal = list(order = c(0, 1, 1), period = 12))
cpi_health_residSF <- residuals(cpi_health_sarimaSF)
cpi_health_parsSF <- coefs2poly(cpi_health_sarimaSF)

plot(cpi_health_resid, type = "l")
lines(cpi_health_residSF, col = 2)
abline(0, 0, lty = 6)

cpi_health_o1SF <- locate.outliers(resid = cpi_health_residSF, pars = cpi_health_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(cpi_health, start_year = 2000)[cpi_health_o1SF$ind, ], cpi_health_o1SF)


# 4 Seasonal Adjustment for CPI Trans Category

cpi_trans <- ts(all_monthly$CPI_trans, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(cpi_trans), name = "log(CPI_trans)")
par(mfrow = c(1, 1))

# cpi_trans_seas <- seas_flow(vara = cpi_trans)
# saveRDS(cpi_trans_seas, "output/seas output/cpi_trans_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
cpi_trans_seas <- readRDS("output/seas output/cpi_trans_seas.RDS")
summary(cpi_trans_seas$x13res)
cpi_trans_seas$SF_ind
head(cpi_trans_seas$x13res$data)

plot(cpi_trans_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Trans Category: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(cpi_trans_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "CPI Trans Category: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
cpi_trans_final <- cpi_trans_seas$x13res$data[, 1]
cpi_trans_season <- cpi_trans_seas$x13res$data[, 2]
cpi_trans_trend <- cpi_trans_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
cpi_trans_sarima <- arima(log(cpi_trans), order = c(0, 1, 1),
                           seasonal = list(order = c(0, 1, 1), period = 12))
cpi_trans_sarima

cpi_trans_resid <- residuals(cpi_trans_sarima)

start.year <- year(as.Date(time(cpi_trans)))[1]
infile <- paste("processed data/flow SF regressor/", cpi_trans_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(cpi_trans)), ] 
h = data.frame(x)


cpi_trans_sarimaSF <- arimax(log(cpi_trans), order = c(0, 1, 1), xreg = h,
                              seasonal = list(order = c(0, 1, 1), period = 12))
cpi_trans_residSF <- residuals(cpi_trans_sarimaSF)
cpi_trans_parsSF <- coefs2poly(cpi_trans_sarimaSF)

plot(cpi_trans_resid, type = "l")
lines(cpi_trans_residSF, col = 2)
abline(0, 0, lty = 6)

cpi_trans_o1SF <- locate.outliers(resid = cpi_trans_residSF, pars = cpi_trans_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(cpi_trans, start_year = 2000)[cpi_trans_o1SF$ind, ], cpi_trans_o1SF)



# 5 Seasonal Adjustment for PPI
PPI <- ts(all_monthly$PPI, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI), name = "log(PPI)")
par(mfrow = c(1, 1))

# PPI_seas <- seas_flow(vara = PPI)
# saveRDS(PPI_seas, "output/seas output/PPI_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_seas <- readRDS("output/seas output/PPI_seas.RDS")

summary(PPI_seas$x13res)
PPI_seas$SF_ind
head(PPI_seas$x13res$data)

plot(PPI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

# The seasonal adjustment output says “Series should not be a candidate for seasonal 
# adjustment because the spectrum of the prior adjusted series (Table B1) has no visually 
# significant seasonal peaks”. In addition, the plots do not show obvious difference 
# between original data and seasonal adjusted data.

PPI_final <- PPI_seas$x13res$data[, 1]
PPI_season <- PPI_seas$x13res$data[, 2]
PPI_trend <- PPI_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
PPI_sarima <- arima(log(PPI), order = c(3, 1, 1))
PPI_sarima

PPI_resid <- residuals(PPI_sarima)

start.year <- year(as.Date(time(PPI)))[1]
infile <- paste("processed data/flow SF regressor/", PPI_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(PPI)), ] 
h = data.frame(x)


PPI_sarimaSF <- arimax(log(PPI), order = c(3, 1, 1), xreg = h)
PPI_residSF <- residuals(PPI_sarimaSF)
PPI_parsSF <- coefs2poly(PPI_sarimaSF)

plot(PPI_resid, type = "l")
lines(PPI_residSF, col = 2)
abline(0, 0, lty = 6)

# After looking at the model fitting performance before and after adjusting SF effects, 
# no obvious difference in the residual plots. So there is nothing needed to be adjusted for this series.

# 6 Seasonal Adjustment for PPI_pro
PPI_pro <- ts(all_monthly$PPI_pro, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI_pro), name = "log(PPI_pro)")
par(mfrow = c(1, 1))


# PPI_pro_seas <- seas_flow(vara = PPI_pro)
# saveRDS(PPI_pro_seas, "output/seas output/PPI_pro_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_pro_seas <- readRDS("output/seas output/PPI_pro_seas.RDS")

summary(PPI_pro_seas$x13res)
PPI_pro_seas$SF_ind
head(PPI_pro_seas$x13res$data)

plot(PPI_pro_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_pro: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_pro_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_pro: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

# The seasonal adjustment output says “Series should not be a candidate for seasonal 
# adjustment because the spectrum of the prior adjusted series (Table B1) has no visually 
# significant seasonal peaks”. In addition, the plots do not show obvious difference 
# between original data and seasonal adjusted data.

PPI_pro_final <- PPI_pro_seas$x13res$data[, 1]
PPI_pro_season <- PPI_pro_seas$x13res$data[, 2]
PPI_pro_trend <- PPI_pro_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
PPI_pro_sarima <- arima(log(PPI_pro), order = c(3, 1, 1))
PPI_pro_sarima

PPI_pro_resid <- residuals(PPI_pro_sarima)

start.year <- year(as.Date(time(PPI_pro)))[1]
infile <- paste("processed data/flow SF regressor/", PPI_pro_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(PPI_pro)), ] 
h = data.frame(x)


PPI_pro_sarimaSF <- arimax(log(PPI_pro), order = c(3, 1, 1), xreg = h)
PPI_pro_residSF <- residuals(PPI_pro_sarimaSF)
PPI_pro_parsSF <- coefs2poly(PPI_pro_sarimaSF)

plot(PPI_pro_resid, type = "l")
lines(PPI_pro_residSF, col = 2)
abline(0, 0, lty = 6)

# 7 Seasonal Adjustment for PPI_ex
PPI_ex <- ts(all_monthly$PPI_ex, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI_ex), name = "log(PPI_ex)")
par(mfrow = c(1, 1))

# PPI_ex_seas <- seas_flow(vara = PPI_ex)
# saveRDS(PPI_ex_seas, "output/seas output/PPI_ex_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_ex_seas <- readRDS("output/seas output/PPI_ex_seas.RDS")

summary(PPI_ex_seas$x13res)
PPI_ex_seas$SF_ind
head(PPI_ex_seas$x13res$data)

plot(PPI_ex_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_ex: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_ex_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_ex: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

# The seasonal adjustment output says “Series should not be a candidate for seasonal 
# adjustment because the spectrum of the prior adjusted series (Table B1) has no visually 
# significant seasonal peaks”. In addition, the plots do not show obvious difference 
# between original data and seasonal adjusted data. Also, SF_ind = 1 means no SF effect.



# 8 Seasonal Adjustment for PPI_raw
PPI_raw <- ts(all_monthly$PPI_raw, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI_raw), name = "log(PPI_raw)")
par(mfrow = c(1, 1))


# PPI_raw_seas <- seas_flow(vara = PPI_raw)
# saveRDS(PPI_raw_seas, "output/seas output/PPI_raw_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_raw_seas <- readRDS("output/seas output/PPI_raw_seas.RDS")

summary(PPI_raw_seas$x13res)
PPI_raw_seas$SF_ind
head(PPI_raw_seas$x13res$data)

plot(PPI_raw_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_raw: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_raw_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_raw: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
# The seasonal adjustment output says “Series should not be a candidate for seasonal 
# adjustment because the spectrum of the prior adjusted series (Table B1) has no visually 
# significant seasonal peaks”. In addition, the plots do not show obvious difference 
# between original data and seasonal adjusted data.

PPI_raw_final <- PPI_raw_seas$x13res$data[, 1]
PPI_raw_season <- PPI_raw_seas$x13res$data[, 2]
PPI_raw_trend <- PPI_raw_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
PPI_raw_sarima <- arima(log(PPI_raw), order = c(2, 1, 0))
PPI_raw_sarima

PPI_raw_resid <- residuals(PPI_raw_sarima)

start.year <- year(as.Date(time(PPI_raw)))[1]
infile <- paste("processed data/flow SF regressor/", PPI_raw_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(PPI_raw)), ] 
h = data.frame(x)


PPI_raw_sarimaSF <- arimax(log(PPI_raw), order = c(2, 1, 0), xreg = h)
PPI_raw_residSF <- residuals(PPI_raw_sarimaSF)
PPI_raw_parsSF <- coefs2poly(PPI_raw_sarimaSF)

plot(PPI_raw_resid, type = "l")
lines(PPI_raw_residSF, col = 2)
abline(0, 0, lty = 6)

PPI_raw_o1SF <- locate.outliers(resid = PPI_raw_residSF, pars = PPI_raw_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(PPI_raw, start_year = 2000)[PPI_raw_o1SF$ind, ], PPI_raw_o1SF)

# 9 Seasonal Adjustment for PPI_fac
PPI_fac <- ts(all_monthly$PPI_fac, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI_fac), name = "log(PPI_fac)")
par(mfrow = c(1, 1))


# PPI_fac_seas <- seas_flow(vara = PPI_fac)
# saveRDS(PPI_fac_seas, "output/seas output/PPI_fac_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_fac_seas <- readRDS("output/seas output/PPI_fac_seas.RDS")

summary(PPI_fac_seas$x13res)
PPI_fac_seas$SF_ind
head(PPI_fac_seas$x13res$data)

plot(PPI_fac_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_fac_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

#The seasonal adjustment shows PPI_fac does not need adjustment and no SF effects.

# 10 Seasonal Adjustment for PPI_liv
PPI_liv <- ts(all_monthly$PPI_liv, start = c(2000, 01), frequency = 12)
plot_acf_diff(PPI_liv, name = "PPI_liv")
par(mfrow = c(1, 1))


# PPI_liv_seas <- seas_flow(vara = PPI_liv)
# saveRDS(PPI_liv_seas, "output/seas output/PPI_liv_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_liv_seas <- readRDS("output/seas output/PPI_liv_seas.RDS")

summary(PPI_liv_seas$x13res)
PPI_liv_seas$SF_ind
head(PPI_liv_seas$x13res$data)

plot(PPI_liv_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_liv: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_liv_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_liv: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

#The seasonal adjustment shows PPI_liv does not need adjustment and no SF effects.

# 11 Seasonal Adjustment for PPI_food
PPI_food <- ts(all_monthly$PPI_food, start = c(2000, 01), frequency = 12)
plot_acf_diff(PPI_food, name = "PPI_food")
par(mfrow = c(1, 1))


# PPI_food_seas <- seas_flow(vara = PPI_food)
# saveRDS(PPI_food_seas, "output/seas output/PPI_food_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_food_seas <- readRDS("output/seas output/PPI_food_seas.RDS")

summary(PPI_food_seas$x13res)
PPI_food_seas$SF_ind
head(PPI_food_seas$x13res$data)

plot(PPI_food_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_food: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_food_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_food: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

#The seasonal adjustment shows PPI_food does not need adjustment and no SF effects.

# 12 Seasonal Adjustment for PPI_clo
PPI_clo <- ts(all_monthly$PPI_clo, start = c(2000, 01), frequency = 12)
plot_acf_diff(PPI_clo, name = "PPI_clo")
par(mfrow = c(1, 1))


# PPI_clo_seas <- seas_flow(vara = PPI_clo)
# saveRDS(PPI_clo_seas, "output/seas output/PPI_clo_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_clo_seas <- readRDS("output/seas output/PPI_clo_seas.RDS")

summary(PPI_clo_seas$x13res)
PPI_clo_seas$SF_ind
head(PPI_clo_seas$x13res$data)

plot(PPI_clo_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_clo: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_clo_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_clo: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

PPI_clo_final <- PPI_clo_seas$x13res$data[, 1]
PPI_clo_season <- PPI_clo_seas$x13res$data[, 2]
PPI_clo_trend <- PPI_clo_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
PPI_clo_sarima <- arima(log(PPI_clo), order = c(3, 1, 1),
                        seasonal = list(order = c(1, 0, 1), period = 12))
PPI_clo_sarima

PPI_clo_resid <- residuals(PPI_clo_sarima)

start.year <- year(as.Date(time(PPI_clo)))[1]
infile <- paste("processed data/flow SF regressor/", PPI_clo_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(PPI_clo)), ] 
h = data.frame(x)


PPI_clo_sarimaSF <- arimax(log(PPI_clo), order = c(3, 1, 1), xreg = h,
                           seasonal = list(order = c(1, 0, 1), period = 12))
PPI_clo_residSF <- residuals(PPI_clo_sarimaSF)
PPI_clo_parsSF <- coefs2poly(PPI_clo_sarimaSF)

plot(PPI_clo_resid, type = "l")
lines(PPI_clo_residSF, col = 2)
abline(0, 0, lty = 6)

PPI_clo_o1SF <- locate.outliers(resid = PPI_clo_residSF, pars = PPI_clo_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(PPI_clo, start_year = 2000)[PPI_clo_o1SF$ind, ], PPI_clo_o1SF)

#The seasonal adjustment shows PPI_clo does not need adjustment and no SF effects.


# 13 Seasonal Adjustment for PPI_gen
PPI_gen <- ts(all_monthly$PPI_gen, start = c(2000, 01), frequency = 12)
plot_acf_diff(PPI_gen, name = "PPI_gen")
par(mfrow = c(1, 1))


# PPI_gen_seas <- seas_flow(vara = PPI_gen)
# saveRDS(PPI_gen_seas, "output/seas output/PPI_gen_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_gen_seas <- readRDS("output/seas output/PPI_gen_seas.RDS")

summary(PPI_gen_seas$x13res)
PPI_gen_seas$SF_ind
head(PPI_gen_seas$x13res$data)

plot(PPI_gen_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_gen: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_gen_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_gen: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

PPI_gen_final <- PPI_gen_seas$x13res$data[, 1]
PPI_gen_season <- PPI_gen_seas$x13res$data[, 2]
PPI_gen_trend <- PPI_gen_seas$x13res$data[, 4]

source("code/SARIMAX funcs.R", local = knitr::knit_global())
PPI_gen_sarima <- arima(log(PPI_gen), order = c(2, 1, 0),
                        seasonal = list(order = c(1, 0, 0), period = 12))
PPI_gen_sarima

PPI_gen_resid <- residuals(PPI_gen_sarima)

start.year <- year(as.Date(time(PPI_gen)))[1]
infile <- paste("processed data/flow SF regressor/", PPI_gen_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(PPI_gen)), ] 
h = data.frame(x)


PPI_gen_sarimaSF <- arimax(log(PPI_gen), order = c(2, 1, 0), xreg = h,
                           seasonal = list(order = c(1, 0, 0), period = 12))
PPI_gen_residSF <- residuals(PPI_gen_sarimaSF)
PPI_gen_parsSF <- coefs2poly(PPI_gen_sarimaSF)

plot(PPI_gen_resid, type = "l")
lines(PPI_gen_residSF, col = 2)
abline(0, 0, lty = 6)

PPI_gen_o1SF <- locate.outliers(resid = PPI_gen_residSF, pars = PPI_gen_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(PPI_gen, start_year = 2000)[PPI_gen_o1SF$ind, ], PPI_gen_o1SF)
#The seasonal adjustment shows PPI_gen does not need adjustment and no SF effects.


# 14 Seasonal Adjustment for PPI_dur
PPI_dur <- ts(all_monthly$PPI_dur, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(PPI_dur), name = "log(PPI_dur)")
par(mfrow = c(1, 1))


# PPI_dur_seas <- seas_flow(vara = PPI_dur)
# saveRDS(PPI_dur_seas, "output/seas output/PPI_dur_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
PPI_dur_seas <- readRDS("output/seas output/PPI_dur_seas.RDS")

summary(PPI_dur_seas$x13res)
PPI_dur_seas$SF_ind
head(PPI_dur_seas$x13res$data)

plot(PPI_dur_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_dur: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(PPI_dur_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "PPI_dur: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

PPI_dur_final <- PPI_dur_seas$x13res$data[, 1]
PPI_dur_season <- PPI_dur_seas$x13res$data[, 2]
PPI_dur_trend <- PPI_dur_seas$x13res$data[, 4]


#The seasonal adjustment shows PPI does not need adjustment and no SF effects. Similarly, the sanme as PPI_pro, PPI_ex, PPI_raw, 
#      PPI_fac, PPI_liv, PPI_food, PPI_gen and PPI_dur, they all do not need adjustment and no SF effects.


# 15 Seasonal Adjustment for IPI

IPI <- ts(all_monthly$IPI, start = c(2000, 01), frequency = 12)
plot_acf_diff(vara = log(IPI), name = "log(IPI)")
par(mfrow = c(1, 1))

# IPI_seas <- seas_flow(vara = IPI)
# saveRDS(IPI_seas, "output/seas output/IPI_seas.RDS")
# seas_flow is time consuming, save output so that we do not need to run it everytime
IPI_seas <- readRDS("output/seas output/IPI_seas.RDS")
IPI_seas$SF_ind
summary(IPI_seas$x13res)
head(IPI_seas$x13res$data)
plot(IPI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "IPI: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)
plot(IPI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "IPI: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

IPI_final <- IPI_seas$x13res$data[, 1]
IPI_season <- IPI_seas$x13res$data[, 2]
IPI_trend <- IPI_seas$x13res$data[, 4]


source("code/SARIMAX funcs.R", local = knitr::knit_global())
IPI_sarima <- arima(log(IPI), order = c(2, 1, 0), 
                    seasonal = list(order = c(1, 0, 1), period = 12))

IPI_sarima
IPI_resid <- residuals(IPI_sarima)

infile <- paste("processed data/flow SF regressor/", IPI_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x <- x[1:length(IPI), ]
h = data.frame(x)

IPI_sarimaSF <- arimax(log(IPI), order = c(2, 1, 0), xreg = h,
                       seasonal = list(order = c(1, 0, 1), period = 12))
IPI_residSF <- residuals(IPI_sarimaSF)
IPI_parsSF <- coefs2poly(IPI_sarimaSF)

plot(IPI_resid, type = "l")
lines(IPI_residSF, col = 2)
abline(0, 0, lty = 6)

IPI_o1SF <- locate.outliers(resid = IPI_residSF, pars = IPI_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(IPI, start_year = 2000)[IPI_o1SF$ind, ], IPI_o1SF)
#no seasonal adjustment for IPI

# 16 Seasonal Adjustment for Loan
Loan <- ts(all_monthly$Loan, start = c(2000, 01), frequency = 12)
plot_acf_diff(Loan, name = "Loan")
par(mfrow = c(1, 1))
# Loan_seas <- seas_stock(vara = Loan)
# saveRDS(Loan_seas, "output/seas output/Loan_seas.RDS")
Loan_seas <- readRDS("output/seas output/Loan_seas.RDS")

summary(Loan_seas$x13res)
Loan_seas$SF_ind
head(Loan_seas$x13res$data)
plot(Loan_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Loan: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(Loan_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Loan: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)


Loan_sarima <- arima(Loan, order = c(0,2,1),
                     seasonal = list(order = c(0, 1, 1), period = 12))
Loan_sarima
Loan_resid <- residuals(Loan_sarima)

start.year <- year(as.Date(time(Loan)))[1]
infile <- paste("processed data/stock SF regressor/", Loan_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(Loan)), ] 
h = data.frame(x)


Loan_sarimaSF <- arimax(Loan, order = c(0, 2, 1), xreg = h,
                        seasonal = list(order = c(0, 1, 1), period = 12))
Loan_residSF <- residuals(Loan_sarimaSF)

plot(Loan_resid, type = "l")
lines(Loan_residSF, col = 2)
abline(0, 0, lty = 6)
# From the result, no seasonal adjustment is needed for Loan. Similarly, 
# no seasonal adjustment for Deposit.

# 17 Seasonal Adjustment for import
import <- ts(all_monthly$Imp, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(import), name = "log(Import)")
par(mfrow = c(1, 1))

# import_seas <- seas_flow(vara = import)
# saveRDS(import_seas, "output/seas output/import_seas.RDS")
import_seas <- readRDS("output/seas output/import_seas.RDS")

import_seas
summary(import_seas$x13res)
import_seas$SF_ind
head(import_seas$x13res$data)
# final is the seasonal adjusted series, which equals to seasonaladj
plot(import_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "import: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(import_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Import: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

import_final <- import_seas$x13res$data[, 1]
import_season <- import_seas$x13res$data[, 2]
import_trend <- import_seas$x13res$data[, 4]

import_sarima <- arima(log(import), order = c(3,1,1),
                       seasonal = list(order = c(0, 1, 1), period = 12))
import_resid <- residuals(import_sarima)

start.year <- year(as.Date(time(import)))[1]
infile <- paste("processed data/flow SF regressor/", import_seas$SF_ind , ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(import)), ] 
h = data.frame(x)


import_sarimaSF <- arimax(log(import), order = c(3,1,1), xreg = h,
                          seasonal = list(order = c(0, 1, 1), period = 12))
import_residSF <- residuals(import_sarimaSF)
import_parsSF <- coefs2poly(import_sarimaSF)

plot(import_resid, type = "l")
lines(import_residSF, col = 2)
abline(0, 0, lty = 6)

import_o1SF <- locate.outliers(resid = import_residSF, pars = import_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(import, start_year = 2000)[import_o1SF$ind, ], import_o1SF)


# 18 Seasonal Adjustment for export
export <- ts(all_monthly$Exp, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(export), name = "log(export)")
par(mfrow = c(1, 1))

# export_seas <- seas_flow(vara = export)
# saveRDS(export_seas, "output/seas output/export_seas.RDS")
export_seas <- readRDS("output/seas output/export_seas.RDS")

export_seas
summary(export_seas$x13res)
export_seas$SF_ind
head(export_seas$x13res$data)
# final is the seasonal adjusted series, which equals to seasonaladj
plot(export_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Export: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(export_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Export: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

export_final <- export_seas$x13res$data[, 1]
export_season <- export_seas$x13res$data[, 2]
export_trend <- export_seas$x13res$data[, 4]

export_sarima <- arima(log(export), order = c(0,1,1),
                       seasonal = list(order = c(0, 1, 1), period = 12))
export_sarima
export_resid <- residuals(export_sarima)

start.year <- year(as.Date(time(export)))[1]
infile <- paste("processed data/flow SF regressor/", export_seas$SF_ind , ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(export)), ] 
h = data.frame(x)


export_sarimaSF <- arimax(log(export), order = c(0, 1, 1), xreg = h,
                          seasonal = list(order = c(0, 1, 1), period = 12))
export_residSF <- residuals(export_sarimaSF)
export_parsSF <- coefs2poly(export_sarimaSF)

plot(export_resid, type = "l")
lines(export_residSF, col = 2)
abline(0, 0, lty = 6)

export_o1SF <- locate.outliers(resid = export_residSF, pars = export_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(export, start_year = 2000)[export_o1SF$ind, ], export_o1SF)

# 19 Seasonal Adjustment for Commodity_imputed

Commodity_imputed <- ts(all_monthly$Commodity_imputed[121 : nrow(all_monthly)], start = c(2010, 01), frequency = 12)
plot_acf_diff(Commodity_imputed, name = "Commodity_imputed")
par(mfrow = c(1, 1))

# Commodity_imputed_seas <- seas_flow(vara = Commodity_imputed)
# saveRDS(Commodity_imputed_seas, "output/seas output/Commodity_imputed_seas.RDS")
Commodity_imputed_seas <- readRDS("output/seas output/Commodity_imputed_seas.RDS")

Commodity_imputed_seas
summary(Commodity_imputed_seas$x13res)
Commodity_imputed_seas$SF_ind
head(Commodity_imputed_seas$x13res$data)
# final is the seasonal adjusted series, which equals to seasonaladj
plot(Commodity_imputed_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Commodity_imputed: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(Commodity_imputed_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Commodity_imputed: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

Commodity_imputed_final <- Commodity_imputed_seas$x13res$data[, 1]
Commodity_imputed_season <- Commodity_imputed_seas$x13res$data[, 2]
Commodity_imputed_trend <- Commodity_imputed_seas$x13res$data[, 4]

Commodity_imputed_sarima <- arima(log(Commodity_imputed), order = c(0, 1, 2),
                       seasonal = list(order = c(1, 1, 0), period = 12))
Commodity_imputed_sarima
Commodity_imputed_resid <- residuals(Commodity_imputed_sarima)

start.year <- year(as.Date(time(Commodity_imputed)))[1]
infile <- paste("processed data/flow SF regressor/", Commodity_imputed_seas$SF_ind , ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(Commodity_imputed)), ] 
h = data.frame(x)


Commodity_imputed_sarimaSF <- arimax(log(Commodity_imputed), order = c(0, 1, 2), xreg = h,
                          seasonal = list(order = c(1, 1, 0), period = 12))
Commodity_imputed_residSF <- residuals(Commodity_imputed_sarimaSF)
Commodity_imputed_parsSF <- coefs2poly(Commodity_imputed_sarimaSF)

plot(Commodity_imputed_resid, type = "l")
lines(Commodity_imputed_residSF, col = 2)
abline(0, 0, lty = 6)

Commodity_imputed_o1SF <- locate.outliers(resid = Commodity_imputed_residSF, pars = Commodity_imputed_parsSF, types = c("AO", "LS", "TC", "IO"))
cbind(obtain_time(Commodity_imputed, start_year = 2010)[Commodity_imputed_o1SF$ind, ], Commodity_imputed_o1SF)


# 20 Seasonal Adjustment for Consumption_imputed

Consumption_imputed <- ts(all_monthly$Consumption_imputed, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(Consumption_imputed), name = "log(Consumption_imputed)")
par(mfrow = c(1, 1))

# Consumption_imputed_seas <- seas_flow(vara = Consumption_imputed)
# saveRDS(Consumption_imputed_seas, "output/seas output/Consumption_imputed_seas.RDS")
Consumption_imputed_seas <- readRDS("output/seas output/Consumption_imputed_seas.RDS")
Consumption_imputed_seas$SF_ind
summary(Consumption_imputed_seas$x13res)
# final is the seasonal adjusted series, which equals to seasonaladj
plot(Consumption_imputed_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Consumption_imputed: Index Data", transform = c("none"))
legend("topleft", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(Consumption_imputed_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "Consumption_imputed: Percentage Change of Year", transform = c("PCY"))
legend("topright", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

Consumption_imputed_final <- Consumption_imputed_seas$x13res$data[, 1]
Consumption_imputed_season <- Consumption_imputed_seas$x13res$data[, 2]
Consumption_imputed_trend <- Consumption_imputed_seas$x13res$data[, 4]
#SF_ind = 1，no SF effects, only seasonal adjustment.

# 21 Seasonal Adjustment for M0
M0 <- ts(all_monthly$M0, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(M0), name = "log(M0)")
par(mfrow = c(1, 1))


# M0_seas <- seas_flow(vara = M0)
# saveRDS(M0_seas, "output/seas output/M0_seas.RDS")
M0_seas <- readRDS("output/seas output/M0_seas.RDS")

summary(M0_seas$x13res)
M0_seas$SF_ind
head(M0_seas$x13res$data)
plot(M0_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M0: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(M0_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M0: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

M0_final <- M0_seas$x13res$data[, 1]
M0_season <- M0_seas$x13res$data[, 2]
M0_trend <- M0_seas$x13res$data[, 4]

M0_sarima <- arima(log(M0), order = c(1, 1, 1),
                   seasonal = list(order = c(1, 1, 1), period = 12))
M0_sarima
M0_resid <- residuals(M0_sarima)

start.year <- year(as.Date(time(M0)))[1]
infile <- paste("processed data/flow SF regressor/", M0_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(M0)), ] 
h = data.frame(x)


M0_sarimaSF <- arimax(log(M0), order = c(1, 1, 1), xreg = h,
                      seasonal = list(order = c(1, 1, 1), period = 12))
M0_residSF <- residuals(M0_sarimaSF)

plot(M0_resid, type = "l")
lines(M0_residSF, col = 2)
abline(0, 0, lty = 6)

#I tried both SF effect for flow and stock data, and find the flow-type SF effect 
#gives better model fitting. The seasonal adjustment use flow-type SF effects.


# 22 Seasonal Adjustment for M1
M1 <- ts(all_monthly$M1, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(M1), name = "log(M1)")
par(mfrow = c(1, 1))


# M1_seas <- seas_stock(vara = M1)
# saveRDS(M1_seas, "output/seas output/M1_seas.RDS")
M1_seas <- readRDS("output/seas output/M1_seas.RDS")

summary(M1_seas$x13res)
M1_seas$SF_ind
head(M1_seas$x13res$data)
plot(M1_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M1: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(M1_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M1: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

M1_final <- M1_seas$x13res$data[, 1]
M1_season <- M1_seas$x13res$data[, 2]
M1_trend <- M1_seas$x13res$data[, 4]

M1_sarima <- arima(log(M1), order = c(0, 1, 0),
                   seasonal = list(order = c(1, 1, 0), period = 12))
M1_sarima
M1_resid <- residuals(M1_sarima)

start.year <- year(as.Date(time(M1)))[1]
infile <- paste("processed data/flow SF regressor/", M1_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(M1)), ] 
h = data.frame(x)


M1_sarimaSF <- arimax(log(M1), order = c(0, 1, 0), xreg = h,
                      seasonal = list(order = c(1, 1, 0), period = 12))
M1_residSF <- residuals(M1_sarimaSF)

plot(M1_resid, type = "l")
lines(M1_residSF, col = 2)
abline(0, 0, lty = 6)


# 23 Seasonal Adjustment for M2
M2 <- ts(all_monthly$M2, start = c(2000, 01), frequency = 12)
plot_acf_diff(log(M2), name = "log(M2)")
par(mfrow = c(1, 1))


# M2_seas <- seas_stock(vara = M2)
# saveRDS(M1_seas, "output/seas output/M2_seas.RDS")
M2_seas <- readRDS("output/seas output/M2_seas.RDS")

summary(M2_seas$x13res)
M2_seas$SF_ind
head(M2_seas$x13res$data)
plot(M2_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M2: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(M2_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "M2: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

M2_final <- M2_seas$x13res$data[, 1]
M2_season <- M2_seas$x13res$data[, 2]
M2_trend <- M2_seas$x13res$data[, 4]

M2_sarima <- arima(log(M2), order = c(0, 1, 0),
                   seasonal = list(order = c(1, 1, 0), period = 12))
M2_sarima
M2_resid <- residuals(M2_sarima)

start.year <- year(as.Date(time(M2)))[1]
infile <- paste("processed data/flow SF regressor/", M2_seas$SF_ind, ".txt", sep = "")
x <- read.table(infile, header = FALSE)
x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + length(M2)), ] 
h = data.frame(x)


M2_sarimaSF <- arimax(log(M2), order = c(0, 1, 0), xreg = h,
                      seasonal = list(order = c(1, 1, 0), period = 12))
M2_residSF <- residuals(M2_sarimaSF)

plot(M2_resid, type = "l")
lines(M2_residSF, col = 2)
abline(0, 0, lty = 6)

#The seasonal adjustment shows M2 does not need adjustment and no SF effects.


# 24 Seasonal Adjustment for NEERI
NEERI <- ts(all_monthly$NEERI[!is.na(all_monthly$NEERI)], start = c(2000, 01), frequency = 12)
plot_acf_diff(NEERI, name = "NEERI")
par(mfrow = c(1, 1))


# NEERI_seas <- seas_stock(vara = NEERI)
# saveRDS(NEERI_seas, "output/seas output/NEERI_seas.RDS")
NEERI_seas <- readRDS("output/seas output/NEERI_seas.RDS")

summary(NEERI_seas$x13res)
NEERI_seas$SF_ind
head(NEERI_seas$x13res$data)
plot(NEERI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "NEERI: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(NEERI_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "NEERI: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)

#The seasonal adjustment shows NEERI does not need adjustment and no SF effects.


# 24 Seasonal Adjustment for SHPE
SHPE <- ts(all_monthly$SHPE[25 : nrow(all_monthly)], start = c(2002, 01), frequency = 12)
plot_acf_diff(SHPE, name = "SHPE")
par(mfrow = c(1, 1))


# SHPE_seas <- seas_stock(vara = SHPE)
# saveRDS(SHPE_seas, "output/seas output/SHPE_seas.RDS")
SHPE_seas <- readRDS("output/seas output/SHPE_seas.RDS")

summary(SHPE_seas$x13res)
SHPE_seas$SF_ind
head(SHPE_seas$x13res$data)
plot(SHPE_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "SHPE: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(SHPE_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "SHPE: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
#The seasonal adjustment shows SHPE does not need adjustment and no SF effects,the same as SZPE.

# 25 Seasonal Adjustment for SZPE
SZPE <- ts(all_monthly$SZPE[25 : nrow(all_monthly)], start = c(2002, 01), frequency = 12)
plot_acf_diff(log(SZPE), name = "log(SZPE)")
par(mfrow = c(1, 1))


# SZPE_seas <- seas_stock(vara = SZPE)
# saveRDS(SZPE_seas, "output/seas output/SZPE_seas.RDS")
SZPE_seas <- readRDS("output/seas output/SZPE_seas.RDS")

summary(SZPE_seas$x13res)
SZPE_seas$SF_ind
head(SZPE_seas$x13res$data)
plot(SZPE_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "SZPE: Index Data", transform = c("none"))
legend("topright", c("original","seasonal adjusted"),
       lty=c(1, 1), col = c(1, 2), lw = 2)

plot(SZPE_seas$x13res, outliers = TRUE, trend = FALSE,
     main = "SZPE: Percentage Change of Year", transform = c("PCY"))
legend("topleft", c("original", "seasonal adjusted"),
       lty = c(1, 1), col = c(1, 2),lw = 2)
#The seasonal adjustment shows SZPE does not need adjustment and no SF effects.

# Output Seasonal Adjusted Data

library(dplyr)

## Note, need to check all if last several obs are missing or not!! for example,
## RPI data in 2023 is not published yet, so c(rpi_final, rep(NA, 3))
all_monthly_seasadj <- all_monthly %>% select(-c(CPI, CPI_food, CPI_health, CPI_trans, Commodity, Commodity_imputed, Consumption, Consumption_imputed,
                                                Imp, Exp, Exp_YoY, Imp_YoY, Trade, M0, M1, M0_YoY, M1_YoY, M2_YoY)) %>%
  mutate(CPI_seasadj = cpi_final,
         CPI_food_seasadj = cpi_food_final,
         CPI_health_seasadj = cpi_health_final,
         CPI_trans_seasadj = cpi_trans_final,
         Commodity_imputed_seasadj = c(rep(NA, 120), Commodity_imputed_final),
         Consumption_imputed_seasadj = Consumption_imputed_final,
         Import_seasadj = import_final,
         Export_seasadj = export_final,
         M0_seasadj = M0_final,
         M1_seasadj = M1_final) %>%
  mutate(Trade_balance_seaadj = Export_seasadj - Import_seasadj)
colnames(all_monthly_seasadj)

all_monthly_season <- data.frame("CPI_season" = cpi_season,
                                 "CPI_food_season" = cpi_food_season,
                                 "CPI_health_season" = cpi_health_season,
                                 "CPI_trans_season" = cpi_trans_season,
                                 "Commodity_imputed_season" = c(rep(NA, 120), Commodity_imputed_season),
                                 "Consumption_imputed_season" = Consumption_imputed_season,
                                 "Import_season" = import_season,
                                 "Export_season" = export_season,
                                 "M0_season" = M0_season,
                                 "M1_season" = M1_season)

write.csv(all_monthly_seasadj, 
          file = "processed data/all_monthly_seasadj.csv",
          row.names = FALSE)
write.csv(all_monthly_season, 
          file = "processed data/all_monthly_season.csv",
          row.names = FALSE)

