# ==============================================================================
# functions_new.R — 不确定性组合预测方法：函数定义文件（Rolling Window版）
# ==============================================================================
# 包含: 数据加载、MBB、五模型预测、权重计算、评估指标
#
# ========================= 方法设计说明 =========================
#
# 【数据划分】
# - 验证集: Rolling window，训练窗口固定长度，逐期滚动
# - 测试集: 固定训练集为最大可用（验证集之前全部数据）
# - 参数: H_val=24 (验证集), N_test=12 (测试集)
#
# 【MBB Bootstrap】
# - 块长度: l = round(2 * n^(1/3))
# - Bootstrap次数: B = 50
# - 最后一块固定为真实数据（保证尾部不受截断影响）
#
# 【五模型】
# 1. AR: 递推预测，OLS估计，BIC选阶(最大12阶)
# 2. FAVAR: 因子数K=4，滞后p=2，Bayesian估计
# 3. DFM: 动态因子模型，r=1因子，p=2滞后
# 4. XGBoost: Direct Forecasting，embed_lag=4
# 5. NNAR: 神经网络自回归，size=24，repeats=20
#
# 【组合预测】
# - 权重: 逆MSPE加权
# - 文献: Adhikari & Agrawal (2014), Baumeister & Kilian (2015)
#
# ==============================================================================

# ========================= 加载依赖包 =========================

library(readr)
library(FAVAR)
library(dfms)
library(xgboost)
library(forecast)   # for nnetar 和 AR

# ========================= 1.1 数据加载与变换 =========================
# 不变

load_and_transform <- function(filepath) {
  DATA <- read.csv(filepath)
  DATA <- na.omit(DATA)
  DATE <- DATA[, 1]
  Y <- DATA[, -1]
  Y <- na.omit(Y)

  Y_adf <- data.frame(
    "MEC" = diff(Y$MEC),
    "NEERI" = diff(Y$NEERI),
    "SHPE" = diff(Y$SHPE),
    "SZPE" = diff(Y$SZPE),
    "Deposit" = diff(log(Y$Deposit)),
    "Loan" = diff(log(Y$Loan)),
    "M2" = diff(log(Y$M2)),
    "PPI" = diff(Y$PPI),
    "PPI_pro" = diff(Y$PPI_pro),
    "PPI_ex" = diff(Y$PPI_ex),
    "PPI_raw" = diff(Y$PPI_raw),
    "PPI_fac" = diff(Y$PPI_fac),
    "PPI_liv" = diff(Y$PPI_liv),
    "PPI_food" = diff(Y$PPI_food),
    "PPI_clo" = diff(Y$PPI_clo),
    "PPI_gen" = diff(Y$PPI_gen),
    "PPI_dur" = diff(Y$PPI_dur),
    "IPI" = diff(Y$IPI),
    "RECI_imputed" = diff(Y$RECI_imputed),
    "CPI_food_seasadj" = diff(Y$CPI_food_seasadj),
    "CPI_health_seasadj" = diff(Y$CPI_health_seasadj),
    "CPI_trans_seasadj" = diff(Y$CPI_trans_seasadj),
    "Consumption_imputed_seasadj" = diff(log(Y$Consumption_imputed_seasadj)),
    "Import_seasadj" = diff(log(Y$Import_seasadj)),
    "Export_seasadj" = diff(log(Y$Export_seasadj)),
    "M0_seasadj" = diff(log(Y$M0_seasadj)),
    "M1_seasadj" = diff(log(Y$M1_seasadj)),
    "Trade_balance_seasadj" = diff(Y$Trade_balance_seasadj),
    "CPI_seasadj" = diff(Y$CPI_seasadj)
  )

  DATE_adf <- DATE[-1]
  return(list(Y_adf = Y_adf, DATE = DATE_adf))
}

# ========================= 1.2 数据划分（Rolling Window版） =========================
# 验证集: rolling window，训练窗口固定长度，逐期滚动
# 测试集: 固定训练集为最大可用（验证集之前全部数据）

split_data <- function(Y_adf, H_val = 24, N_test = 12) {
  n <- nrow(Y_adf)

  # 固定训练窗口大小
  train_size <- n - H_val - N_test

  # 验证集范围
  val_start <- train_size + 1
  val_end <- n - N_test

  # 测试集范围
  test_start <- val_end + 1
  test_end <- n

  # 测试集使用的训练集（最大可用 = 验证集最后一行）
  test_train_end <- val_end

  cat(sprintf("数据总行数: %d\n", n))
  cat(sprintf("训练窗口大小: %d\n", train_size))
  cat(sprintf("验证集: %d ~ %d (%d期)\n", val_start, val_end, H_val))
  cat(sprintf("测试集: %d ~ %d (%d期)\n", test_start, test_end, N_test))
  cat(sprintf("测试集训练集: 1 ~ %d (固定)\n", test_train_end))

  return(list(
    train_size     = train_size,
    val_start      = val_start,
    val_end        = val_end,
    test_start     = test_start,
    test_end       = test_end,
    test_train_end = test_train_end
  ))
}

# ========================= 1.3 MBB函数 =========================
# 注意：滚动窗口下，每次循环的 train_data 长度不同 → l 不同 → MBB 样本自然不同
# 这是正确行为，不需额外处理

get_mbb_sample <- function(data, l) {
  n <- nrow(data)
  num_blocks <- n - l + 1
  k <- ceiling(n / l)

  # 前 (k-1) 块随机有放回抽样，最后一块固定为真实的最后块 (n-l+1)
  # 保证尾部始终是真实数据，不受截断影响
  sampled_blocks <- c(
    sample(1:(num_blocks - 1), size = k - 1, replace = TRUE),
    num_blocks
  )

  # 前 k-1 块拼接后，截断到 (n - l) 行（给最后一块留出完整空间）
  mbb_concat <- do.call(rbind, lapply(sampled_blocks[-length(sampled_blocks)], function(block_start) {
    data[block_start:(block_start + l - 1), ]
  }))
  max_concat <- n - l
  mbb_concat <- if (nrow(mbb_concat) > max_concat) mbb_concat[1:max_concat, , drop = FALSE] else mbb_concat

  # 最后一块固定为真实尾部，不截断
  last_block <- data[(n - l + 1):n, , drop = FALSE]

  mbb_data <- rbind(mbb_concat, last_block)

  return(mbb_data)
}

calc_block_length <- function(n) {
  round(2 * n^(1/3))
}

# ========================= 1.4 模型预测函数 =========================

# ---------- AR模型 ----------

predict_AR_h_step <- function(train_data, target_var, h, params = list()) {
  y <- train_data[, target_var]
  order_max <- ifelse(is.null(params$order_max), 12, params$order_max)

  fit <- ar(y, order.max = order_max, method = "ols", demean = TRUE, intercept = TRUE)
  pred <- predict(fit, n.ahead = h)

  return(as.numeric(pred$pred))
}

# ---------- FAVAR模型 ----------

predict_FAVAR_h_step <- function(train_data, target_var, h, params = list()) {
  K <- ifelse(is.null(params$K), 4, params$K)
  p_lag <- ifelse(is.null(params$p_lag), 2, params$p_lag)
  nrep <- ifelse(is.null(params$nrep), 1500, params$nrep)
  nburn <- ifelse(is.null(params$nburn), 1000, params$nburn)

  Y_train <- train_data[, target_var]
  X_train <- train_data[, !names(train_data) %in% target_var]

  fit <- FAVAR(
    Y = Y_train,
    X = X_train,
    fctmethod = 'BBE',
    factorprior = list(b0 = 0, vb0 = NULL, c0 = 0.01, d0 = 0.01),
    varprior = list(b0 = 0, vb0 = 10, nu0 = 0, s0 = 0),
    nrep = nrep, nburn = nburn,
    K = K, plag = p_lag
  )

  n_vars <- K + 1
  varcoef <- coef(fit)$varcoef

  A_list <- list()
  for (i in 1:p_lag) {
    start_col <- (i - 1) * n_vars + 1
    end_col <- i * n_vars
    A_list[[i]] <- varcoef[, start_col:end_col]
  }

  if (ncol(varcoef) > p_lag * n_vars) {
    const <- varcoef[, ncol(varcoef)]
  } else {
    const <- rep(0, n_vars)
  }

  F_hist <- fit$factorx
  Z_hist <- cbind(Y_train, F_hist)
  n_obs <- nrow(Z_hist)

  predictions <- numeric(h)
  Z_buffer <- Z_hist[(n_obs - p_lag + 1):n_obs, , drop = FALSE]

  for (step in 1:h) {
    Z_pred <- const
    n_buf <- nrow(Z_buffer)
    for (i in 1:p_lag) {
      Z_pred <- Z_pred + A_list[[i]] %*% Z_buffer[n_buf - i + 1, ]
    }
    predictions[step] <- Z_pred[1]
    Z_buffer <- rbind(Z_buffer, as.numeric(Z_pred))
  }

  return(predictions)
}

# ---------- DFM模型 ----------

predict_DFM_h_step <- function(train_data, target_var, h, params = list()) {
  r <- ifelse(is.null(params$r), 1, params$r)
  p <- ifelse(is.null(params$p), 2, params$p)

  if (any(is.na(train_data))) {
    train_data <- tsnarmimp(train_data)
  }

  dfm_fit <- DFM(train_data, r = r, p = p)
  forecast_result <- predict(dfm_fit, h = h)
  forecast_df <- as.data.frame(forecast_result$X_fcst)

  predictions <- forecast_df[1:h, target_var]
  return(as.numeric(predictions))
}

# ---------- XGBoost模型 (Direct Forecasting) ----------

predict_XGBoost_h_step <- function(train_data, target_var, h, params = list()) {
  embed_lag <- ifelse(is.null(params$embed_lag), 4, params$embed_lag)
  nrounds   <- ifelse(is.null(params$nrounds),   1000, params$nrounds)
  eta       <- ifelse(is.null(params$eta),        0.01, params$eta)
  max_depth <- ifelse(is.null(params$max_depth),  4,    params$max_depth)
  colsample <- ifelse(is.null(params$colsample),  2/3,  params$colsample)
  subsample <- ifelse(is.null(params$subsample),  1,    params$subsample)

  target_col <- which(colnames(train_data) == target_var)
  data_matrix <- as.matrix(train_data)
  X0 <- cbind(data_matrix[, target_col], data_matrix[, -target_col])

  aux <- embed(X0, embed_lag + h)
  y <- aux[, 1]
  X <- aux[, -c(1:(ncol(X0) * h))]

  if (h == 1) {
    X.out <- tail(aux, 1)[1:ncol(X)]
  } else {
    X.out <- aux[, -c(1:(ncol(X0) * (h - 1)))]
    X.out <- tail(X.out, 1)[1:ncol(X)]
  }

  y <- y[1:(length(y) - h + 1)]
  X <- X[1:(nrow(X) - h + 1), ]

  dtrain <- xgb.DMatrix(data = X, label = y)
  xgb_params <- xgb.params(
    learning_rate     = eta,
    max_depth         = max_depth,
    colsample_bylevel = colsample,
    subsample         = subsample,
    nthread           = 1,
    objective         = "reg:squarederror"
  )
  model <- xgb.train(
    params  = xgb_params,
    data    = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  pred <- predict(model, xgb.DMatrix(matrix(X.out, nrow = 1)))

  predictions <- rep(NA, h)
  predictions[h] <- pred
  return(predictions)
}

# ---------- NNAR模型 ----------

predict_NNAR_h_step <- function(train_data, target_var, h, params = list()) {
  p_order <- if (is.null(params$p)) NULL else params$p
  size    <- ifelse(is.null(params$size), 24, params$size)
  repeats <- ifelse(is.null(params$repeats), 20, params$repeats)

  y <- ts(train_data[, target_var], frequency = 1)

  if (is.null(p_order)) {
    fit <- nnetar(y, P = 0, size = size, repeats = repeats)
  } else {
    fit <- nnetar(y, p = p_order, P = 0, size = size, repeats = repeats)
  }

  fcast <- forecast::forecast(fit, h = h)
  predictions <- as.numeric(fcast$mean)

  return(predictions)
}

# ========================= 1.5 权重计算 =========================

compute_inverse_mspe_weights <- function(mspe_vec) {
  model_names_vec <- names(mspe_vec)

  inv_mspe <- 1 / mspe_vec
  inv_mspe[is.infinite(inv_mspe) | is.na(inv_mspe)] <- 0

  total <- sum(inv_mspe, na.rm = TRUE)

  if (total == 0 || is.na(total)) {
    stop("所有模型的MSPE均无效（全为NA/NaN），请检查模型预测是否正常运行。")
  }

  weights <- inv_mspe / total
  names(weights) <- model_names_vec
  return(weights)
}

# ========================= 1.6 评估指标 =========================

calc_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

calc_mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}

calc_smape <- function(actual, predicted) {
  denominator <- abs(actual) + abs(predicted)
  valid <- denominator > 0
  if (sum(valid) == 0) return(NA)
  mean(2 * abs(actual[valid] - predicted[valid]) / denominator[valid], na.rm = TRUE) * 100
}

calc_theils_u1 <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2)) / (sqrt(mean(actual^2)) + sqrt(mean(predicted^2)))
}

calc_mdape <- function(actual, predicted) {
  valid <- actual != 0
  if (sum(valid) == 0) return(NA)
  median(abs((actual[valid] - predicted[valid]) / actual[valid]), na.rm = TRUE) * 100
}

evaluate_all <- function(actual, predicted) {
  results <- c(
    RMSE = calc_rmse(actual, predicted),
    MAE = calc_mae(actual, predicted),
    SMAPE = calc_smape(actual, predicted),
    Theils_U1 = calc_theils_u1(actual, predicted),
    MDAPE = calc_mdape(actual, predicted)
  )
  return(results)
}

# ========================= 辅助函数 =========================

safe_predict <- function(predict_fn, train_data, target_var, h, params, model_name) {
  tryCatch({
    pred <- predict_fn(train_data, target_var, h, params)
    return(pred)
  }, error = function(e) {
    cat(sprintf("  %s 预测失败: %s\n", model_name, e$message))
    return(rep(NA, h))
  })
}

# ========================= MSPE 保存（用于后续权重对比） =========================
# 保存验证集的MSPE矩阵，供后续新权重方法对比使用

save_mspe_for_comparison <- function(mspe_by_horizon, model_names, horizons, filepath_prefix = "mspe_data") {
  # mspe_by_horizon: list，每个元素是 named vector (模型名 -> MSPE)
  # model_names: 模型名向量
  # horizons: 预测步长向量

  for (h in horizons) {
    mspe_h <- mspe_by_horizon[[as.character(h)]]
    if (!is.null(mspe_h)) {
      df <- data.frame(
        Horizon = h,
        Model = names(mspe_h),
        MSPE = as.numeric(mspe_h),
        stringsAsFactors = FALSE
      )
      filename <- sprintf("%s_h%d.csv", filepath_prefix, h)
      write.csv(df, filename, row.names = FALSE)
      cat(sprintf("已保存: %s\n", filename))
    }
  }

  # 保存汇总表（所有h的平均MSPE）
  summary_df <- data.frame(Horizon = horizons, stringsAsFactors = FALSE)
  for (m in model_names) {
    summary_df[[m]] <- sapply(horizons, function(h) {
      mspe_h <- mspe_by_horizon[[as.character(h)]]
      if (is.null(mspe_h)) NA else mspe_h[m]
    })
  }
  summary_file <- sprintf("%s_summary.csv", filepath_prefix)
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat(sprintf("已保存: %s\n", summary_file))
}
