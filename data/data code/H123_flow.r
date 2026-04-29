
######### create SF regressors for 2000-01 to 2100-12  ###########
rm(list=ls())
library(dplyr)
library(lubridate)
setwd("/Users/liumeishao/paper")
SFdate <- read.csv("raw data/1956 to 2101 Chinese New Year dates.csv", header = TRUE)


# The monthly data start from 2000-01, so we only need SFdate start from 2000
# date of SF is shown in numeric value: 1.22 means Jan 22, 2.05 means Feb 05

SFdate <- SFdate[SFdate$year >= 1999,]

year.start <- 2000
year.end <-2101 # SF in 2101 may cause SF effect on 2100-12
HolidayTime <- data.frame("year" = rep(c(year.start:year.end), each = 12), 
                          "month" = rep(c(1:12), (year.end - year.start + 1)),
                          "monthdate"=rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), 
                                          (year.end - year.start + 1)))
HolidayTime$monthdate[HolidayTime$month == 2 & 
                        leap_year(HolidayTime$year)] <- 29
# set 1999/12/31 as time = 0 in days
# mt[i] is the time of the last day of month [i - 1], 
# i.e. how many days the last day of month [i - 1] is away from 1999/12/31
# ht[i] is the time of SF in the year of month i, 
# i.e. how many days the SF date in that year is away from 1999/12/31
# For example of month i = 2000/01: mt[i] = 0 since last day of month [i - 1] is 1999/12/31, 
# ht[i] = 36 since SF of 2000 is on 2000/02/05, 2000/02/05 - 1999/12/31 = 36 
HolidayTime$mt <- cumsum(c(0, HolidayTime$monthdate))[1 : nrow(HolidayTime)]
HolidayTime <- left_join(HolidayTime, SFdate, by = c("year")) %>% 
  mutate(SFmonth = floor(date), 
         SFday = (date - SFmonth) * 100, 
         ht_tmp = mt + SFday, 
         ht = NA)
for(yr in unique(HolidayTime$year)){
  HolidayTime$ht[HolidayTime$year == yr] <- HolidayTime$ht_tmp[HolidayTime$year==yr & HolidayTime$month == HolidayTime$SFmonth]
}

############## generate regressor files #################
w <-  data.frame("w1" = rep(c(0, 4, 8, 12, 16, 20, 24), 21),
                 "w2" = rep(rep(c(0, 4, 8), each = 7), 7),
                 "w3" = rep(c(0, 4, 8, 12, 16, 20, 24), each = 21))
# (w1, w2, w3) are (tau1, tau2, tau3) in the "Seasonal Adjustment for Monthly Data"
ht <- HolidayTime$ht
mt <- HolidayTime$mt   
J <- nrow(w)
Nmonth <- (year.end - year.start) * 12 # number of months from 1990-01 to 2100-12

for(j in 2 : J) { # 1st row of w is (w1, w2, w3) = (0, 0, 0), i.e. no SF effects
  w1 <- w[j, 1]
  w2 <- w[j, 2]
  w3 <- w[j, 3]
  before_start <- ht - w1 # time of the first day of before period of SF in the year of month i 
  before_end <- ht - 1 # time of the last day of before period of SF in the year of month i 
  during_start <- ht # time of the first day of during period of SF in the year of month i 
  during_end <- ht + w2 - 1 # time of the last day of during period of SF in the year of month i 
  after_start <- ht + w2 # time of the first day of after period of SF in the year of month i 
  after_end <- ht + w2 + w3 - 1 # time of the last day of after period of SF in the year of month i 
  
  # h1, h2, h3 are vectors of SF regressors H_{it} (i = 1, 2, 3)
  h1 <- rep(NA, Nmonth)
  h2 <- h1
  h3 <- h1
 
  ## calculate regressor for before period, i.e. h1
  if(w1 != 0) {
    for(i in 1 : Nmonth) {
      if((before_start[i] > mt[i]) & (before_start[i] <= mt[(i + 1)])) { # start day of before period falls in month i
        if(before_end[i] <= mt[(i + 1)]) h1[i] <- 1 # end of before period falls in month i
        else h1[i] <- (mt[(i + 1)] - before_start[i] + 1) / w1
      } else if ((before_end[i] > mt[i]) & (before_end[i] <= mt[(i + 1)])) { # end of before period falls in month i, but start day of before period is not in month i 
        h1[i] <- (before_end[i] - mt[i]) / w1
      } else if ((before_start[(i + 1)] > mt[i]) & (before_start[(i + 1)] <= mt[(i + 1)])) { 
        # month i is Dec and month i+1 is Jan of next yeasr. Then before_start[i] != before_start[(i + 1)]
        # before_start of next year's SF falls into month i, e.g. before_start of 2010 SF falls in to 2009-12 (i = 2009-12)
        # then before_end of next year's SF cannot fall into month i, so only one scenario as follow
        h1[i] <- (mt[(i + 1)] - before_start[(i + 1)] + 1) / w1
      } else h1[i] <- 0
    } # end of for loop
  } # end of if(w1 != 0)
  
  ## calculate regressor for during period, i.e. h2
  if (w2 != 0) {
    for(i in 1 : Nmonth) {
      if ((during_start[i] > mt[i]) & (during_start[i] <= mt[(i + 1)])) {
        if ((during_end[i]) <= mt[(i+1)]) h2[i] <- 1 
        else h2[i] <- (mt[(i + 1)] - during_start[i] + 1) / w2
      } else if ((during_end[i] > mt[i]) & ((during_end[i]) <= mt[(i + 1)])) {
        h2[i] = (during_end[i] - mt[i]) / w2
      } else h2[i] <- 0
    }
  }
  
  ## calculate regressor for after period, i.e. h3
  if (w3 != 0) { 
    for(i in 1 : Nmonth) {
      if((after_start[i] > mt[i]) & (after_start[i] <= mt[(i + 1)])) {
        if(after_end[i] <= mt[(i + 1)]) h3[i] <- 1 
        else h3[i] <- (mt[(i + 1)] - after_start[i] + 1) / w3
      }
      else if((after_end[i] > mt[i]) & (after_end[i] <= mt[(i + 1)])) {
        h3[i] <- (after_end[i] - mt[i]) / w3
      } else h3[i] <- 0
    }
  }
  
  # since i != 1, w1, w2, w3 cannot be all 0.
  if(w1 == 0) {
    if (w2 == 0) SF_regressors <- h3
    else if (w3 == 0) SF_regressors <- h2
    else SF_regressors <- cbind(h2, h3) 
  } else if (w2 == 0) {
    if (w3 == 0) SF_regressors <- h1
    else SF_regressors <- cbind(h1, h3)
  } else if (w3 == 0) {
    SF_regressors <- cbind(h1, h2)
  } else SF_regressors <- cbind(h1, h2, h3)
  
  SF_regressors <- round(SF_regressors, digits = 4)
  outfile  <- paste("processed data/flow SF regressor/", j, ".txt", sep = "")
  write.table(SF_regressors, outfile, col.names = FALSE, row.names = FALSE, sep = " ")
}



