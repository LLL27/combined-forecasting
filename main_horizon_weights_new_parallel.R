# ==============================================================================
# main_horizon_weights_new_parallel.R — MBB组合预测（Windows兼容版）
# ==============================================================================

rm(list = ls())
set.seed(123)

# ========================= 加载函数 =========================
source("/Users/liumeishao/Documents/claude/paper writing/forecast_code/functions_new.R")

# ========================= 参数设置 =========================

data_path <- "/Users/liumeishao/Documents/claude/paper writing/data/processed data/all_monthly_seasadj_1208_var_adjust.csv"
target_var <- "CPI_seasadj"

H_val <- 24
N_test <- 12
B <- 50

horizons <- c(1, 2, 3, 6, 9, 12)

use_models <- c(
  AR      = TRUE,
  FAVAR   = TRUE,
  DFM     = TRUE,
  XGBoost = TRUE,
  NNAR    = TRUE
)

model_names <- names(use_models[use_models])
n_models <- length(model_names)

params_list <- list(
  AR      = list(order_max = 12),
  FAVAR   = list(K = 4, p_lag = 2, nrep = 300, nburn = 200),
  DFM     = list(r = 1, p = 2),
  XGBoost = list(embed_lag = 3, nrounds = 1000, eta = 0.01, max_depth = 4),
  NNAR    = list(p = NULL, size = 24, repeats = 20)
)

predict_fns <- list(
  AR      = predict_AR_h_step,
  FAVAR   = predict_FAVAR_h_step,
  DFM     = predict_DFM_h_step,
  XGBoost = predict_XGBoost_h_step,
  NNAR    = predict_NNAR_h_step
)

# ========================= 加载数据 =========================

cat("=== 加载数据 ===\n")
data_result <- load_and_transform(data_path)
Y_adf <- data_result$Y_adf
DATE <- data_result$DATE

cat(sprintf("数据维度: %d 行 × %d 列\n", nrow(Y_adf), ncol(Y_adf)))
cat(sprintf("启用的模型: %s\n", paste(model_names, collapse = ", ")))

# ========================= 数据划分 =========================

cat("\n=== 数据划分 ===\n")
splits <- split_data(Y_adf, H_val = H_val, N_test = N_test)

# ==============================================================================
# 第一阶段: 验证集 — 计算MBB权重
# ==============================================================================

cat("\n============================================\n")
cat("第一阶段: 验证集 — 按步长分别计算MBB权重\n")
cat("============================================\n")

weights_by_horizon <- list()
mspe_by_horizon <- list()

# Helper: 单次bootstrap预测所有模型
run_one_boot <- function(b, train_data, l, h, target_var, params_list, model_names, predict_fns) {
  mbb_sample <- get_mbb_sample(train_data, l)

  preds <- numeric(length(model_names))
  names(preds) <- model_names

  for (m in model_names) {
    preds[m] <- tryCatch({
      pred <- predict_fns[[m]](mbb_sample, target_var, h = h, params_list[[m]])
      pred[h]
    }, error = function(e) NA)
  }

  preds
}

for (h in horizons) {
  cat(sprintf("\n====== 计算 h=%d 步权重 ======\n", h))

  n_val_points <- H_val - h + 1
  cat(sprintf("  可用验证点数: %d, B=%d\n", n_val_points, B))

  MSPE_val_h <- matrix(NA, nrow = n_val_points, ncol = n_models)
  colnames(MSPE_val_h) <- model_names

  for (v in 1:n_val_points) {
    train_start <- v
    train_end <- train_start + splits$train_size - 1
    actual_idx <- train_end + h
    actual_val <- Y_adf[[target_var]][actual_idx]

    train_data <- Y_adf[train_start:train_end, ]
    l <- calc_block_length(nrow(train_data))

    cat(sprintf("  验证点 %d/%d (训练 %d:%d, 预测第%d行)\n",
                v, n_val_points, train_start, train_end, actual_idx))

    # Bootstrap循环（Windows用lapply）
    seed_vec <- 123L + seq_len(B) * 1000

    if (.Platform$OS.type == "windows") {
      boot_list <- lapply(1:B, function(b) {
        set.seed(seed_vec[b])
        run_one_boot(b, train_data, l, h, target_var, params_list, model_names, predict_fns)
      })
    } else {
      boot_list <- parallel::mclapply(1:B, function(b) {
        set.seed(seed_vec[b])
        run_one_boot(b, train_data, l, h, target_var, params_list, model_names, predict_fns)
      }, mc.cores = min(8, B))
    }

    boot_preds <- do.call(rbind, boot_list)

    for (m in model_names) {
      preds_m <- boot_preds[, m]
      MSPE_val_h[v, m] <- mean((preds_m - actual_val)^2, na.rm = TRUE)
    }

    cat(sprintf("    MSPE: "))
    for (m in model_names) {
      cat(sprintf("%s=%.4f ", m, MSPE_val_h[v, m]))
    }
    cat("\n")
  }

  avg_MSPE_h <- colMeans(MSPE_val_h, na.rm = TRUE)
  weights_h <- compute_inverse_mspe_weights(avg_MSPE_h)

  weights_by_horizon[[as.character(h)]] <- weights_h
  mspe_by_horizon[[as.character(h)]] <- avg_MSPE_h

  cat(sprintf("\n  h=%d 平均MSPE: ", h))
  print(round(avg_MSPE_h, 6))
  cat(sprintf("  h=%d 权重: ", h))
  print(round(weights_h, 4))
}

# ==============================================================================
# 第二阶段: 测试集 — 版本A（每个h有多个预测，训练数据动态扩展）
# ==============================================================================

cat("\n============================================\n")
cat("第二阶段: 测试集 — 版本A（动态扩展训练集）\n")
cat("============================================\n")

all_results_A <- list()

for (h in horizons) {
  cat(sprintf("\n====== h = %d ======\n", h))

  n_forecasts <- N_test - h + 1
  cat(sprintf("  预测次数: %d\n", n_forecasts))

  weights <- weights_by_horizon[[as.character(h)]]

  preds_matrix <- matrix(NA, nrow = n_forecasts, ncol = n_models + 1)
  colnames(preds_matrix) <- c(model_names, "Combined")
  actuals <- numeric(n_forecasts)

  for (f in 1:n_forecasts) {
    actual_idx <- splits$test_start + f + h - 2
    train_end <- actual_idx - h
    train_data <- Y_adf[1:train_end, ]
    actuals[f] <- Y_adf[[target_var]][actual_idx]

    cat(sprintf("  f=%d/%d: 训练 1:%d, 预测第%d行\n",
                f, n_forecasts, train_end, actual_idx))

    for (m in model_names) {
      pred <- safe_predict(predict_fns[[m]], train_data, target_var, h, params_list[[m]], m)
      preds_matrix[f, m] <- pred[h]
    }

    model_preds <- preds_matrix[f, model_names]
    valid_models <- !is.na(model_preds)
    if (any(valid_models)) {
      w_valid <- weights[valid_models]
      w_valid <- w_valid / sum(w_valid)
      preds_matrix[f, "Combined"] <- sum(w_valid * model_preds[valid_models])
    }
  }

  eval_results <- data.frame(Model = c(model_names, "Combined"), stringsAsFactors = FALSE)
  for (col_name in c(model_names, "Combined")) {
    preds <- preds_matrix[, col_name]
    metrics <- evaluate_all(actuals, preds)
    if (col_name == model_names[1]) {
      for (metric_name in names(metrics)) eval_results[[metric_name]] <- NA
    }
    row_idx <- which(eval_results$Model == col_name)
    for (metric_name in names(metrics)) eval_results[row_idx, metric_name] <- metrics[metric_name]
  }

  eval_results$Horizon <- h
  eval_results$Actual <- mean(actuals)
  eval_results$Actuals <- paste(round(actuals, 4), collapse = "; ")
  eval_results <- cbind(eval_results[, c("Model", "Horizon", "Actual", "Actuals")],
                        eval_results[, c("RMSE", "MAE", "SMAPE", "Theils_U1", "MDAPE")])
  all_results_A[[as.character(h)]] <- eval_results

  # 保存原始预测矩阵（用于后续权重对比）
  preds_raw <- data.frame(
    Forecast = 1:n_forecasts,
    Actual = actuals,
    preds_matrix
  )
  write.csv(preds_raw, sprintf("/Users/liumeishao/Documents/claude/paper writing/forecast_code/result/predictions_raw_A_h%d.csv", h), row.names = FALSE)

  cat(sprintf("\n--- h=%d 评估结果 ---\n", h))
  print(eval_results[, c("Model", "RMSE", "MAE", "SMAPE", "Theils_U1", "MDAPE")])
}

# ==============================================================================
# 第二阶段: 测试集 — 版本B（每个h只有一个预测，训练数据固定）
# ==============================================================================

cat("\n============================================\n")
cat("第二阶段: 测试集 — 版本B（固定训练集）\n")
cat("============================================\n")

all_results_B <- list()

for (h in horizons) {
  cat(sprintf("\n====== h = %d ======\n", h))

  actual_idx <- splits$test_train_end + h
  actual_val <- Y_adf[[target_var]][actual_idx]
  train_data <- Y_adf[1:splits$test_train_end, ]
  weights <- weights_by_horizon[[as.character(h)]]

  cat(sprintf("  训练 1:%d, 预测第%d行\n", splits$test_train_end, actual_idx))

  preds <- numeric(n_models + 1)
  names(preds) <- c(model_names, "Combined")

  for (m in model_names) {
    pred <- safe_predict(predict_fns[[m]], train_data, target_var, h, params_list[[m]], m)
    preds[m] <- pred[h]
  }

  model_preds <- preds[model_names]
  valid_models <- !is.na(model_preds)
  if (any(valid_models)) {
    w_valid <- weights[valid_models]
    w_valid <- w_valid / sum(w_valid)
    preds["Combined"] <- sum(w_valid * model_preds[valid_models])
  }

  eval_results <- data.frame(
    Model = c(model_names, "Combined"),
    Horizon = h,
    Actual = actual_val,
    RMSE = NA, MAE = NA, SMAPE = NA, Theils_U1 = NA, MDAPE = NA,
    stringsAsFactors = FALSE
  )

  for (col_name in c(model_names, "Combined")) {
    metrics <- evaluate_all(actual_val, preds[col_name])
    row_idx <- which(eval_results$Model == col_name)
    for (metric_name in names(metrics)) eval_results[row_idx, metric_name] <- metrics[metric_name]
  }

  all_results_B[[as.character(h)]] <- eval_results

  cat(sprintf("\n--- h=%d 评估结果 ---\n", h))
  print(eval_results[, c("Model", "RMSE", "MAE", "SMAPE", "Theils_U1", "MDAPE")])
}

# ==============================================================================
# 保存结果
# ==============================================================================

cat("\n============================================\n")
cat("保存结果\n")
cat("============================================\n")

# 1. 保存预测值
for (h in horizons) {
  preds_A_h <- all_results_A[[as.character(h)]]
  if (!is.null(preds_A_h)) {
    write.csv(preds_A_h, sprintf("predictions_A_h%d.csv", h), row.names = FALSE)
  }
}

preds_B_all <- do.call(rbind, all_results_B)
rownames(preds_B_all) <- NULL
write.csv(preds_B_all, "predictions_B_all.csv", row.names = FALSE)

# 2. MSPE汇总
cat("\n=== 验证集 MSPE 汇总 ===\n")
mspe_summary <- data.frame(Horizon = horizons, stringsAsFactors = FALSE)
for (m in model_names) {
  mspe_summary[[m]] <- sapply(horizons, function(h) {
    mspe_h <- mspe_by_horizon[[as.character(h)]]
    if (is.null(mspe_h)) NA else mspe_h[m]
  })
}
mspe_summary$Best_Model <- apply(mspe_summary[, model_names, drop=FALSE], 1, function(row) {
  names(which.min(row))
})
print(mspe_summary)
write.csv(mspe_summary, "mspe_summary.csv", row.names = FALSE)

# 3. 最优模型
cat("\n=== 各步长最优模型（版本A）===\n")
for (h in horizons) {
  eval_h <- all_results_A[[as.character(h)]]
  if (!is.null(eval_h)) {
    best_idx <- which.min(eval_h$RMSE)
    combined_idx <- which(eval_h$Model == "Combined")
    cat(sprintf("h=%2d: 最优=%s (RMSE=%.4f), 组合=%s (RMSE=%.4f)\n",
                 h, eval_h$Model[best_idx], eval_h$RMSE[best_idx],
                 eval_h$Model[combined_idx], eval_h$RMSE[combined_idx]))
  }
}

cat("\n=== 各步长最优模型（版本B）===\n")
for (h in horizons) {
  eval_h <- all_results_B[[as.character(h)]]
  if (!is.null(eval_h)) {
    best_idx <- which.min(eval_h$RMSE)
    combined_idx <- which(eval_h$Model == "Combined")
    cat(sprintf("h=%2d: 最优=%s (RMSE=%.4f), 组合=%s (RMSE=%.4f)\n",
                 h, eval_h$Model[best_idx], eval_h$RMSE[best_idx],
                 eval_h$Model[combined_idx], eval_h$RMSE[combined_idx]))
  }
}

# 4. 权重
cat("\n=== 各步长权重 ===\n")
weights_summary <- data.frame(Horizon = horizons, stringsAsFactors = FALSE)
for (m in model_names) {
  weights_summary[[m]] <- sapply(horizons, function(h) {
    w <- weights_by_horizon[[as.character(h)]]
    if (is.null(w)) NA else round(w[m], 4)
  })
}
print(weights_summary)
write.csv(weights_summary, "weights_by_horizon.csv", row.names = FALSE)

# 5. 评估结果
final_results_A <- do.call(rbind, all_results_A)
final_results_B <- do.call(rbind, all_results_B)
rownames(final_results_A) <- NULL
rownames(final_results_B) <- NULL
write.csv(final_results_A, "results_test_A.csv", row.names = FALSE)
write.csv(final_results_B, "results_test_B.csv", row.names = FALSE)

# 6. MSPE详情
save_mspe_for_comparison(mspe_by_horizon, model_names, horizons)

cat("\n=== 完成! ===\n")
cat("结果已保存:\n")
cat("  - predictions_A_h*.csv\n")
cat("  - predictions_B_all.csv\n")
cat("  - mspe_summary.csv\n")
cat("  - weights_by_horizon.csv\n")
cat("  - results_test_A.csv\n")
cat("  - results_test_B.csv\n")
cat("  - mspe_data_h*.csv\n")
