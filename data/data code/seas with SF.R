library(dplyr)
library(seasonal)
library(seasonalview)
library(zoo)
library(lubridate)

#### n is the sample size of vara
plot_acf_diff <- function(vara, name) {
  par(mfrow = c(3,2), mgp=c(1.5, 0.5, 0), mar=c(3, 3, 3, 1), cex = 0.8)
  plot(vara, main = paste("(a)", name, ' Series'))
  acf(as.vector(vara), ci.type = 'ma',main = paste('(b) SACFs of ', name), lag = 48)
  plot(diff(vara), main = paste('(c)',  name, 'after 1st Order Difference'))
  acf(as.vector(diff(vara)), lag = 48, ci.type = 'ma', 
      main = paste('(d) SACFs of', name, 'after 1st Order Difference'))
  plot(diff(diff(vara),lag = 12),
       main = paste('(e)', name, 'after 1st and 12th Order Difference'))
  acf(as.vector(diff(diff(vara),lag = 12)),ci.type = 'ma',
      main = paste('(f) SACFs of', name, 'after 1st and 12th Order Difference'), lag = 48)
}

seas_flow = function(vara, outlier = NULL, transform.function = "auto") {
  if (!is.ts(vara) | (month(as.Date(time(vara)))[1] != 1)) stop("vara must be an ts object and start with Jan of a year")
  n <- length(vara)
  start.year <- year(as.Date(time(vara)))[1]
  aicc <- rep(NA, 147)
  
  x13out <- seas(vara, outlier = outlier, 
                 transform.function = transform.function, 
                 x11 = '', regression.aictest = NULL)
  aicc[1] <- udg(x13out, stat = "aicc")
  
  for (i in 2 : 147) {
    infile <- paste("processed data/flow SF regressor/", i, ".txt", sep = "")
    x <- read.table(infile, header = FALSE)
    x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + n + 12), ] 
    # total n + 12 rows of SF regressors, 12 month ahead of vara since seas will do 12-month forecast
    h = ts(x, frequency = 12, start = c(start.year, 01))
    
    x13out <- NULL
    try(x13out <- seas(vara, xreg = h, regression.usertype = "holiday",
                       outlier = outlier, 
                       transform.function = transform.function, 
                   x11 = '', regression.aictest = NULL))
    if(!is.null(x13out)) aicc[i] <- udg(x13out, stat = "aicc")
  }
  
  SF_ind <- which.min(aicc)
  
  if(SF_ind >= 2){
    infile <- paste("processed data/flow SF regressor/", SF_ind, ".txt", sep = "")
    x <- read.table(infile, header = FALSE)
    x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + n + 12), ] # total n + 12 rows of SF regressors
    h = ts(x, frequency = 12, start = c(start.year, 01))
    
    x13out <- seas(vara, xreg = h, regression.usertype = "holiday",
                   outlier = outlier, 
                   transform.function = transform.function, 
                   x11 = '', regression.aictest = NULL)
  }
  
  if(SF_ind == 1){
    x13out <- seas(vara, outlier = outlier, 
                   transform.function = transform.function, 
                   x11 = '', regression.aictest = NULL)
  }
  
  # outfilepng1<-paste("./plots/",colnames(vara)," Index Data",".png",sep="")
  # plotname1<-paste(colnames(vara),": Plot of Index data",sep="")
  # png(filename=outfilepng1,width=900,height=600)
  # plot(x13out, outliers = TRUE, trend = TRUE,
  #      main = plotname1, transform = c("none"))
  # legend("topleft",c("original","seasonal adjusted","trend"),
  #        lty=c(1,1,2),col=c(1,2,4),lw=2)
  # dev.off()
  # 
  # outfilepng2<-paste("./plots/",colnames(vara)," Percentage of Year",".png",sep="")
  # plotname2<-paste(colnames(vara),": Plot of Percentage Change of Year",sep="")
  # png(filename=outfilepng2,width=900,height=600)
  # plot(x13out, outliers = TRUE, trend = FALSE,
  #      main = plotname2, transform = c("PCY"))
  # legend("topleft",c("original","seasonal adjusted","trend"),
  #        lty=c(1,1,2),col=c(1,2,4),lw=2)
  # dev.off()
  
  return(list(aicc = aicc, SF_ind = SF_ind, x13res = x13out))
}

seas_stock = function(vara, outlier = NULL, transform.function = "auto") {
  if (!is.ts(vara) | (month(as.Date(time(vara)))[1] != 1)) stop("vara must be an ts object and start with Jan of a year")
  n <- length(vara)
  start.year <- year(as.Date(time(vara)))[1]
  aicc <- rep(NA, 147)
  
  x13out <- seas(vara, outlier = outlier, 
                 transform.function = transform.function, 
                 x11 = '', regression.aictest = NULL)
  aicc[1] <- udg(x13out, stat = "aicc")
  
  for (i in 2 : 147) {
    infile <- paste("processed data/stock SF regressor/", i, ".txt", sep = "")
    x <- read.table(infile, header = FALSE)
    x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + n + 12), ] 
    # total n + 12 rows of SF regressors, 12 month ahead of vara since seas will do 12-month forecast
    h = ts(x, frequency = 12, start = c(start.year, 01))
    
    x13out <- NULL
    try(x13out <- seas(vara, xreg = h, regression.usertype = "holiday",
                       outlier = outlier, 
                       transform.function = transform.function, 
                       x11 = '', regression.aictest = NULL))
    if(!is.null(x13out)) aicc[i] <- udg(x13out, stat = "aicc")
  }
  
  SF_ind <- which.min(aicc)
  
  if(SF_ind >= 2){
    infile <- paste("processed data/stock SF regressor/", SF_ind, ".txt", sep = "")
    x <- read.table(infile, header = FALSE)
    x = x[((start.year - 2000) * 12 + 1) : ((start.year - 2000) * 12 + n + 12), ] # total n + 12 rows of SF regressors
    h = ts(x, frequency = 12, start = c(start.year, 01))
    
    x13out <- seas(vara, xreg = h, regression.usertype = "holiday",
                   outlier = outlier, 
                   transform.function = transform.function, 
                   x11 = '', regression.aictest = NULL)
  }
  
  if(SF_ind == 1){
    x13out <- seas(vara, outlier = outlier, 
                   transform.function = transform.function, 
                   x11 = '', regression.aictest = NULL)
  }
  
  # outfilepng1<-paste("./plots/",colnames(vara)," Index Data",".png",sep="")
  # plotname1<-paste(colnames(vara),": Plot of Index data",sep="")
  # png(filename=outfilepng1,width=900,height=600)
  # plot(x13out, outliers = TRUE, trend = TRUE,
  #      main = plotname1, transform = c("none"))
  # legend("topleft",c("original","seasonal adjusted","trend"),
  #        lty=c(1,1,2),col=c(1,2,4),lw=2)
  # dev.off()
  # 
  # outfilepng2<-paste("./plots/",colnames(vara)," Percentage of Year",".png",sep="")
  # plotname2<-paste(colnames(vara),": Plot of Percentage Change of Year",sep="")
  # png(filename=outfilepng2,width=900,height=600)
  # plot(x13out, outliers = TRUE, trend = FALSE,
  #      main = plotname2, transform = c("PCY"))
  # legend("topleft",c("original","seasonal adjusted","trend"),
  #        lty=c(1,1,2),col=c(1,2,4),lw=2)
  # dev.off()
  
  return(list(aicc = aicc, SF_ind = SF_ind, x13res = x13out))
}
  