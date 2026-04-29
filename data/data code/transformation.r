all_monthly_seasadj <- read.csv("processed data/all_monthly_seasadj.csv", header = TRUE)

varas <- colnames(all_monthly_seasadj)
varas

par(mfrow = c(1, 2))
for(i in 2:ncol(all_monthly_seasadj)){
  plot(ts(all_monthly_seasadj[i], start = c(2000, 1), frequency = 12), 
       main = varas[i], type = "o")
  plot(ts(log(all_monthly_seasadj[i]), start = c(2000, 1), frequency = 12), 
       main = paste("log(", varas[i], ")", sep=""), type = "o")
}

library(tseries)
adf.test(diff(all_monthly_seasadj$MEC[!is.na(all_monthly_seasadj$MEC)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CFCI[!is.na(all_monthly_seasadj$CFCI)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CMI[!is.na(all_monthly_seasadj$CMI)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$RECI[!is.na(all_monthly_seasadj$RECI)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PMI_man[!is.na(all_monthly_seasadj$PMI_man)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PMI[!is.na(all_monthly_seasadj$PMI)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$NEERI[!is.na(all_monthly_seasadj$NEERI)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$SHPE[!is.na(all_monthly_seasadj$SHPE)]), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$SZPE[!is.na(all_monthly_seasadj$SZPE)]), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Deposit)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Loan)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$M2)), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_pro), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_ex), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_raw), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_fac), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_liv), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_food), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_gen), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$PPI_dur), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$IPI), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$RECI_imputed), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CPI_seasadj), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CPI_food_seasadj), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CPI_health_seasadj), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$CPI_trans_seasadj), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Commodity_imputed_seasadj[!is.na(all_monthly_seasadj$Commodity_imputed_seasadj)])), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Consumption_imputed_seasadj)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Import_seasadj)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$Export_seasadj)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$M0_seasadj)), alternative = "stationary")

adf.test(diff(log(all_monthly_seasadj$M1_seasadj)), alternative = "stationary")

adf.test(diff(all_monthly_seasadj$Trade_balance_seaadj[!is.na(all_monthly_seasadj$Trade_balance_seaadj)]), alternative = "stationary")
