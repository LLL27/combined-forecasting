# ==============================================================================
# weight_comparison.R — 不同权重方法对比
# 基于保存的 MSPE 数据计算权重，对比评估结果
# ==============================================================================

rm(list = ls())

# ========================= 加载数据 =========================

result_path <- "/Users/liumeishao/Documents/claude/paper writing/forecast_code/result"

# 读取MSPE汇总
mspe_summary <- read.csv(file.path(result_path, "mspe_summary.csv"))
horizons <- mspe_summary$Horizon
model_names <- c("AR", "FAVAR", "DFM", "XGBoost", "NNAR")

# 读取各h的MSPE详情
mspe_by_horizon <- list()
for (h in horizons) {
  mspe_by_horizon[[ as.character(h)]] <- read.csv(
    file.path(result_path, sprintf("mspe_data_h%d.csv", h))
  )
}

# 读取测试集预测结果（原始预测值）
predictions_A <- list()
for (h in horizons) {
  predictions_A[[ as.character(h)]] <- read.csv(
    file.path(result_path, sprintf("predictions_raw_A_h%d.csv", h))
  )
}

cat("数据加载完成\n")
cat("模型:", paste(model_names, collapse = ", "), "\n")
cat("Horizons:", paste(horizons, collapse = ", "), "\n\n")

# ========================= 权重计算函数 =========================

compute_inverse_mspe_weights <- function(mspe_vec) {
  inv_mspe <- 1 / mspe_vec
  inv_mspe[is.infinite(inv_mspe) | is.na(inv_mspe)] <- 0
  total <- sum(inv_mspe, na.rm = TRUE)
  if (total == 0 || is.na(total)) stop("所有模型的MSPE均无效")
  weights <- inv_mspe / total
  names(weights) <- names(mspe_vec)
  return(weights)
}

compute_sqrt_inverse_mspe_weights <- function(mspe_vec) {
  inv_mspe <- 1 / sqrt(mspe_vec)
  inv_mspe[is.infinite(inv_mspe) | is.na(inv_mspe)] <- 0
  total <- sum(inv_mspe, na.rm = TRUE)
  if (total == 0 || is.na(total)) stop("所有模型的MSPE均无效")
  weights <- inv_mspe / total
  names(weights) <- names(mspe_vec)
  return(weights)
}

compute_geometric_inverse_mspe_weights <- function(mspe_vec, alpha = 2) {
  inv_mspe <- (1 / mspe_vec)^alpha
  inv_mspe[is.infinite(inv_mspe) | is.na(inv_mspe)] <- 0
  total <- sum(inv_mspe, na.rm = TRUE)
  if (total == 0 || is.na(total)) stop("所有模型的MSPE均无效")
  weights <- inv_mspe / total
  names(weights) <- names(mspe_vec)
  return(weights)
}

compute_rank_weights <- function(mspe_vec) {
  ranks <- rank(mspe_vec, ties.method = "average")
  weights <- (length(mspe_vec) - ranks + 1)
  weights <- weights / sum(weights)
  names(weights) <- names(mspe_vec)
  return(weights)
}

compute_offset_penalty_weights <- function(mspe_vec, k = 0.1) {
  adjusted_mspe <- mspe_vec * (1 + k)
  inv_mspe <- 1 / adjusted_mspe
  inv_mspe[is.infinite(inv_mspe) | is.na(inv_mspe)] <- 0
  total <- sum(inv_mspe, na.rm = TRUE)
  if (total == 0 || is.na(total)) stop("所有模型的MSPE均无效")
  weights <- inv_mspe / total
  names(weights) <- names(mspe_vec)
  return(weights)
}

compute_shrinkage_weights <- function(mspe_vec, lambda = 0.5) {
  n <- length(mspe_vec)
  w_uniform <- rep(1/n, n)
  w_inv <- compute_inverse_mspe_weights(mspe_vec)
  weights <- lambda * w_uniform + (1 - lambda) * w_inv
  names(weights) <- names(mspe_vec)
  return(weights)
}

# ========================= 权重方法列表 =========================

weight_methods <- list(
  "Inverse_MSPE" = list(fn = compute_inverse_mspe_weights, params = NULL),
  "Sqrt_Inverse_MSPE" = list(fn = compute_sqrt_inverse_mspe_weights, params = NULL),
  "Geometric_α2" = list(fn = compute_geometric_inverse_mspe_weights, params = list(alpha = 2)),
  "Rank_Weighted" = list(fn = compute_rank_weights, params = NULL),
  "Offset_Penalty_k01" = list(fn = compute_offset_penalty_weights, params = list(k = 0.1)),
  "Shrinkage_λ05" = list(fn = compute_shrinkage_weights, params = list(lambda = 0.5))
)

# ========================= 评估函数 =========================

calc_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

evaluate_combination <- function(preds_df, weights, model_names) {
  # preds_df: 包含 Actual 和各模型预测值的 data.frame
  actuals <- preds_df$Actual
  combined <- rep(0, nrow(preds_df))
  for (m in model_names) {
    combined <- combined + weights[m] * preds_df[[m]]
  }
  rmse <- calc_rmse(actuals, combined)
  return(list(RMSE = rmse, combined = combined))
}

# ========================= 计算各方法权重 =========================

cat("============================================\n")
cat("各方法权重计算结果\n")
cat("============================================\n\n")

all_weights <- list()

for (method_name in names(weight_methods)) {
  method <- weight_methods[[method_name]]
  cat(sprintf("--- %s ---\n", method_name))

  weights_by_h <- list()
  for (h in horizons) {
    mspe_vec <- mspe_by_horizon[[ as.character(h)]]$MSPE
    names(mspe_vec) <- mspe_by_horizon[[ as.character(h)]]$Model

    if (!is.null(method$params)) {
      weights_h <- do.call(method$fn, c(list(mspe_vec), method$params))
    } else {
      weights_h <- method$fn(mspe_vec)
    }
    weights_by_h[[ as.character(h)]] <- weights_h
  }
  all_weights[[method_name]] <- weights_by_h

  # 打印各h权重
  for (h in horizons) {
    w <- weights_by_h[[ as.character(h)]]
    cat(sprintf("h=%2d: ", h))
    for (m in model_names) {
      cat(sprintf("%s=%.3f ", m, w[m]))
    }
    cat("\n")
  }
  cat("\n")
}

# ========================= 对比评估结果 =========================

cat("============================================\n")
cat("各方法评估结果对比\n")
cat("============================================\n\n")

all_evaluations <- list()

for (method_name in names(weight_methods)) {
  cat(sprintf("--- %s ---\n", method_name))

  weights_by_h <- all_weights[[method_name]]
  eval_by_h <- list()

  for (h in horizons) {
    preds_df <- predictions_A[[ as.character(h)]]
    eval_result <- evaluate_combination(preds_df, weights_by_h[[ as.character(h)]], model_names)

    eval_by_h[[ as.character(h)]] <- data.frame(
      Horizon = h,
      Method = method_name,
      RMSE = eval_result$RMSE,
      stringsAsFactors = FALSE
    )
    cat(sprintf("h=%2d: RMSE=%.4f\n", h, eval_result$RMSE))
  }
  cat("\n")
  all_evaluations[[method_name]] <- eval_by_h
}

# ========================= 汇总对比表 =========================

cat("============================================\n")
cat("RMSE 对比汇总\n")
cat("============================================\n\n")

rmse_summary <- data.frame(Horizon = horizons, stringsAsFactors = FALSE)
for (method_name in names(weight_methods)) {
  rmse_summary[[method_name]] <- sapply(horizons, function(h) {
    all_evaluations[[method_name]][[ as.character(h)]]$RMSE
  })
}

# 添加简单平均结果
simple_avg_rmse <- sapply(horizons, function(h) {
  preds_df <- predictions_A[[ as.character(h)]]
  combined <- rowMeans(preds_df[, model_names])
  calc_rmse(preds_df$Actual, combined)
})
rmse_summary$Simple_Average <- simple_avg_rmse

# 每行最优方法
rmse_summary$Best_Method <- apply(rmse_summary[, names(weight_methods), drop=FALSE], 1, function(row) {
  names(which.min(row))
})

print(rmse_summary)
write.csv(rmse_summary, file.path(result_path, "weight_comparison_rmse.csv"), row.names = FALSE)

# ========================= 最优方法统计 =========================

cat("\n============================================\n")
cat("各方法最优次数统计\n")
cat("============================================\n\n")

best_count <- table(rmse_summary$Best_Method)
print(best_count)

# ========================= 保存详细结果 =========================

for (method_name in names(weight_methods)) {
  eval_df <- do.call(rbind, all_evaluations[[method_name]])
  write.csv(eval_df, file.path(result_path, sprintf("eval_%s.csv", method_name)), row.names = FALSE)
}

# 保存权重汇总
weights_df <- data.frame(Horizon = horizons, stringsAsFactors = FALSE)
for (method_name in names(weight_methods)) {
  for (m in model_names) {
    col_name <- sprintf("%s_%s", method_name, m)
    weights_df[[col_name]] <- sapply(horizons, function(h) {
      all_weights[[method_name]][[ as.character(h)]][m]
    })
  }
}
write.csv(weights_df, file.path(result_path, "weights_all_methods.csv"), row.names = FALSE)

cat("\n=== 完成! ===\n")
cat("结果已保存:\n")
cat("  - weight_comparison_rmse.csv (RMSE对比汇总)\n")
cat("  - weights_all_methods.csv (各方法权重)\n")
cat("  - eval_*.csv (各方法详细评估)\n")
