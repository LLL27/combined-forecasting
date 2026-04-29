obtain_time <- function(x, start_year){
  print("function obtain_time(): x must be a monthly time series start with Jan of a year.")
  nyear <- ceiling(length(x) / 12)
  yr_mon <- data.frame("year" = c(start_year : (start_year + nyear - 1)),
                       "month" = rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), nyear))
  yr_mon <- yr_mon[1 : length(x), ]
  return(yr_mon)
}

single_arima <- function(x, p, d, q, P, D, Q){
  mdl <- arima(x, order = c(p, d, q), 
               seasonal = list(order = c(P, D, Q), period = 12))
  npara <- p + q + P + Q
  bic <- npara * log(length(x)) - 2 * mdl$loglik
  aicc <- 2 * npara / (1 - (npara + 1) / (length(x) - 12 * D - d)) - 2 * mdl$loglik
  crtr <- c(bic, mdl$aic, aicc, mdl$loglik, mdl$sigma2)
  names(crtr) <- c("bic", "aic", "aicc", "loglik", "sigma2")
  return(crtr)
}


select_sarima <- function(x, d, D, order_mat){
  # Each row of order_mat is a vector of (p, q, P, Q)
  # to be noticed, colnames of order_mat must be p, q, P, Q
  if(! all.equal(colnames(order_mat), c("p", "q", "P", "Q"))) stop("Error: Colnames of order_mat must be p, q, P, Q!")
  slct_criteria <- sapply(c(1 : nrow(order_mat)), FUN = function(i){
    criteria <- rep(NA, 5)
    p = order_mat[i, 1]
    q = order_mat[i, 2]
    P = order_mat[i, 3]
    Q = order_mat[i, 4]
    tmp <- tryCatch({criteria <- single_arima(x = x, p = p, d = d, q = q,
                                              P = P, D = D, Q = Q)},
             error = function(cnd){"error"},
             warning = function(cnd){"warning"}
    )
    if(! length(tmp) == 5){
      if(tmp == "error") print(paste("Error: Row", i, "with order (", p, d, q, ") x (", P, D, Q, ")_12 does not converge"))
      if(tmp == "warning"){
        criteria <- rep(NA, 5)
        print(paste("Warning: Row", i, "with order (", p, d, q, ") x (", P, D, Q, ")_12 has warning in covergence"))
        }
    }
    return(criteria)
  }) # some orders results in convergence issue, tryCatch() function can make sure the execution does not stop and
     # show which row and order has error or warning
  slct_criteria <- t(slct_criteria)
  slct_criteria <- cbind(order_mat, slct_criteria)
  
  return(slct_criteria)
}


# select_SF<- function(x, p, d, q, P, D, Q){
#   # Each row of order_mat is a vector of (p, q, P, Q)
#   # to be noticed, colnames of order_mat must be p, q, P, Q
#   
#   if (!is.ts(vara) | (month(as.Date(time(vara)))[1] != 1)) stop("vara must be an ts object and start with Jan of a year")
#   n <- length(vara)
#   start.year <- year(as.Date(time(vara)))[1]
#   aicc <- rep(NA, 147)
#   
#   x13out <- seas(vara, transform.function = "auto", 
#                  x11 = '', regression.aictest = NULL)
#   aicc[1] <- udg(x13out, stat = "aicc")
#   
#   for (i in 2 : 147) {
#     infile <- paste("../processed data/flow SF regressor/", i, ".txt", sep = "")
#     x <- read.table(infile, header = FALSE)
#     x = x[((start.year - 1990) * 12 + 1) : ((start.year - 1990) * 12 + n + 12), ] 
#     # total n + 12 rows of SF regressors, 12 month ahead of vara since seas will do 12-month forecast
#     h = ts(x, frequency = 12, start = c(start.year, 01))
#     
#     x13out <- NULL
#     try(x13out <- seas(vara, xreg = h, regression.usertype = "holiday",
#                        transform.function = "auto", 
#                        x11 = '', regression.aictest = NULL))
#     if(!is.null(x13out)) aicc[i] <- udg(x13out, stat = "aicc")
#   }
#   
#   
# }

