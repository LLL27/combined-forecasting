get_index <- function(Time, YoY, MoM, index_month = "2001-01"){
  if(length(Time) != length(YoY) |
     length(Time) != length(MoM) ) stop("Time, YoY, MoM must have same length!")
  if(!(index_month %in% Time)) stop("index_month must be within Time!")
  n <- length(YoY)
  index <- rep(NA, n)
  ind <- c(1:n)[Time == index_month]
  
  # get Jan value using YoY
  index[ind] <- 100
  i <- ind - 12
  j <- ind + 12
  while(i >= 1){
    index[i] <- index[(i + 12)] / YoY[(i + 12)] * 100
    i <- i - 12
  }
  while(j <= n){
    index[j] <- index[(j - 12)] * YoY[j] / 100
    j <- j + 12
  }
  
  # get other month's value from MoM and YoY
  for(l in 1:11){
    index[(ind + l)] <- index[(ind + l - 1)] * MoM[(ind + l)] / 100
    i <- ind + l - 12
    j <- ind + l + 12
    while(i >= 1){
      index[i] <- index[(i + 12)] / YoY[(i + 12)] * 100
      i <- i - 12
    }
    while(j <= n){
      index[j] <- index[(j - 12)] * YoY[j] / 100
      j <- j + 12
    }
  }
  
  return(index)
}


get_quarter_cum <- function(Time, Cum, start_quarter){
  print("get_quarter_cum(): start_quarter must be the first quarter of a year, and the first entry of Time must the first quarter of a year. Otherwise the output will be wrong. ")
  if(length(Time) != length(Cum)) stop("Time, Cum must have same length!")
  if(!(start_quarter %in% Time)) stop("start_quarter must be within Time!")
  n <- length(Cum)
  Cum <- c(Cum, rep(NA, (ceiling(n / 4) * 4 - n)))
  Quarter <- sapply(1 : ceiling(n / 4), FUN = function(k){
    c(Cum[(4 * k - 3)], diff(Cum[(4 * k - 3):(4 * k)]))
  })
  Quarter <- as.vector(Quarter)[1:n]
  return(Quarter)
}

get_GDPquarter_YoY <- function(Time, gdp, YoY){
  if(length(Time) != length(YoY) |
     length(Time) != length(gdp) ) stop("Time, gdp, MoM must have same length!")
  n <- length(gdp)
  i <- c(1:n)[Time == "2006-12"]
  while(i >= 1){
    gdp[i] <- gdp[(i + 4)] / (100 + YoY[(i + 4)]) * 100
    i <- i - 1
  }
  return(gdp)
}

