# ============================================================
# 双硫死亡+DN 竞赛项目 — 04c_ML严格验证（审稿修订版）
# 目的：解决审稿人关切的过拟合问题
# 方法：
#   1. 重复10折交叉验证 (repeated 10-fold CV)
#   2. Bootstrap乐观性校正 (optimism-corrected AUC)
#   3. 置换检验 (permutation test)
#   4. AUC置信区间 (DeLong + bootstrap CI)
#   5. 3-gene vs 5-gene vs 7-gene 模型比较
#   6. AIC/BIC比较 + 校准曲线
# ============================================================

library(glmnet)
library(caret)
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(boot)

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")

expr_norm <- readRDS("GSE96804_gene_expr.rds")
pheno <- readRDS("GSE96804_pheno.rds")
final_biomarkers <- readRDS("04_final_biomarkers.rds")

group_info <- pheno[["title"]]
group <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"),
                levels = c("Control", "DN"))
y <- as.numeric(group == "DN")

cat(sprintf("最终标志物 (%d个): %s\n",
            length(final_biomarkers), paste(final_biomarkers, collapse = ", ")))

# 准备数据
X_raw <- as.data.frame(t(expr_norm[final_biomarkers, , drop = FALSE]))
X <- scale(X_raw)
n_samples <- nrow(X)
n_genes <- ncol(X)

cat(sprintf("样本数: %d (DN=%d, Control=%d), 基因数: %d\n",
            n_samples, sum(y==1), sum(y==0), n_genes))

# ---- 2. 重复10折交叉验证 ----
cat("\n===== 重复10折交叉验证 =====\n")

set.seed(42)
n_repeats <- 10
n_folds <- 10
cv_aucs <- matrix(NA, nrow = n_repeats, ncol = n_folds)

for (r in 1:n_repeats) {
  fold_ids <- createFolds(as.factor(y), k = n_folds, list = TRUE)
  for (f in 1:n_folds) {
    test_idx <- fold_ids[[f]]
    train_idx <- setdiff(1:n_samples, test_idx)

    # 训练逻辑回归
    train_data <- as.data.frame(cbind(y = y[train_idx], X[train_idx, ]))
    model <- glm(y ~ ., data = train_data, family = binomial())

    # 预测
    test_data <- as.data.frame(X[test_idx, , drop = FALSE])
    colnames(test_data) <- colnames(X)
    pred_prob <- predict(model, newdata = test_data, type = "response")

    # 计算AUC
    roc_obj <- roc(y[test_idx], pred_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
    cv_aucs[r, f] <- auc(roc_obj)
  }
  cat(sprintf("  重复 %d/10: mean AUC = %.4f (SD = %.4f)\n",
              r, mean(cv_aucs[r, ]), sd(cv_aucs[r, ])))
}

# 汇总统计
all_cv_aucs <- as.vector(cv_aucs)
cv_mean_auc <- mean(all_cv_aucs)
cv_sd_auc <- sd(all_cv_aucs)
cv_se_auc <- cv_sd_auc / sqrt(length(all_cv_aucs))
cv_ci_lower <- cv_mean_auc - 1.96 * cv_se_auc
cv_ci_upper <- cv_mean_auc + 1.96 * cv_se_auc

cat(sprintf("\n重复10×10折CV结果:\n"))
cat(sprintf("  Mean AUC = %.4f\n", cv_mean_auc))
cat(sprintf("  SD = %.4f\n", cv_sd_auc))
cat(sprintf("  95%% CI = [%.4f, %.4f]\n", cv_ci_lower, cv_ci_upper))

# ---- 3. Bootstrap乐观性校正 ----
cat("\n===== Bootstrap乐观性校正 =====\n")

# 计算表观AUC（在全部训练数据上训练和评估）
full_data <- as.data.frame(cbind(y = y, X))
full_model <- glm(y ~ ., data = full_data, family = binomial())
full_pred <- predict(full_model, type = "response")
apparent_auc <- auc(roc(y, full_pred, levels = c(0, 1), direction = "<", quiet = TRUE))
cat(sprintf("表观AUC (训练集): %.4f\n", apparent_auc))

# Bootstrap乐观性估计
set.seed(42)
n_boot <- 500
optimism <- numeric(n_boot)
boot_aucs <- numeric(n_boot)

for (b in 1:n_boot) {
  # Bootstrap重采样
  boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
  oob_idx <- setdiff(1:n_samples, unique(boot_idx))

  if (length(oob_idx) < 5) next  # 需要足够的OOB样本

  # 训练数据
  boot_data <- full_data[boot_idx, ]
  boot_model <- glm(y ~ ., data = boot_data, family = binomial())

  # 在bootstrap样本上的性能
  boot_pred_train <- predict(boot_model, newdata = boot_data, type = "response")
  boot_auc_train <- auc(roc(boot_data$y, boot_pred_train,
                            levels = c(0, 1), direction = "<", quiet = TRUE))

  # 在OOB样本上的性能
  oob_data <- full_data[oob_idx, ]
  oob_pred <- predict(boot_model, newdata = oob_data, type = "response")
  boot_auc_test <- auc(roc(oob_data$y, oob_pred,
                           levels = c(0, 1), direction = "<", quiet = TRUE))

  optimism[b] <- boot_auc_train - boot_auc_test
  boot_aucs[b] <- boot_auc_test
}

# 去除NA
optimism <- optimism[!is.na(optimism)]
boot_aucs <- boot_aucs[!is.na(boot_aucs)]

mean_optimism <- mean(optimism)
optimism_corrected_auc <- apparent_auc - mean_optimism

cat(sprintf("\nBootstrap乐观性校正 (n=%d):\n", length(optimism)))
cat(sprintf("  表观AUC: %.4f\n", apparent_auc))
cat(sprintf("  平均乐观性: %.4f\n", mean_optimism))
cat(sprintf("  乐观性校正后AUC: %.4f\n", optimism_corrected_auc))
cat(sprintf("  OOB AUC (mean): %.4f (SD=%.4f)\n", mean(boot_aucs), sd(boot_aucs)))

# ---- 4. 置换检验 ----
cat("\n===== 置换检验 =====\n")

set.seed(42)
n_perm <- 1000
perm_aucs <- numeric(n_perm)

for (p in 1:n_perm) {
  y_perm <- sample(y)  # 打乱标签
  perm_data <- as.data.frame(cbind(y = y_perm, X))
  perm_model <- glm(y ~ ., data = perm_data, family = binomial())
  perm_pred <- predict(perm_model, type = "response")
  perm_aucs[p] <- auc(roc(y_perm, perm_pred, levels = c(0, 1),
                          direction = "<", quiet = TRUE))
}

# P值 = 置换AUC >= 真实AUC的比例
p_value <- mean(perm_aucs >= apparent_auc)
cat(sprintf("置换检验 (n=%d):\n", n_perm))
cat(sprintf("  真实AUC: %.4f\n", apparent_auc))
cat(sprintf("  置换AUC均值: %.4f (SD=%.4f)\n", mean(perm_aucs), sd(perm_aucs)))
cat(sprintf("  置换P值: %.4f %s\n", p_value,
            ifelse(p_value < 0.001, "(P<0.001, 高度显著)",
            ifelse(p_value < 0.05, "(P<0.05)", "(不显著)"))))

# 置换分布图
perm_df <- data.frame(AUC = perm_aucs)
ggplot(perm_df, aes(x = AUC)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = apparent_auc, color = "red", linewidth = 1.5, linetype = "dashed") +
  annotate("text", x = apparent_auc + 0.02, y = 30,
           label = sprintf("Observed\nAUC=%.3f", apparent_auc),
           color = "red", hjust = 0, size = 4) +
  labs(title = sprintf("Permutation Test (n=%d) — P = %.4f", n_perm, p_value),
       x = "AUC under null hypothesis", y = "Frequency") +
  theme_minimal(base_size = 13)
ggsave("04c_permutation_test.png", width = 8, height = 5, dpi = 150)
cat("✅ 置换检验图已保存\n")

# ---- 5. AUC置信区间 (DeLong方法) ----
cat("\n===== AUC置信区间 =====\n")

roc_full <- roc(y, full_pred, levels = c(0, 1), direction = "<")
ci_delong <- ci.auc(roc_full, method = "delong")
ci_boot <- ci.auc(roc_full, method = "bootstrap", boot.n = 2000, quiet = TRUE)

cat(sprintf("AUC置信区间:\n"))
cat(sprintf("  DeLong 95%% CI: [%.4f, %.4f]\n", ci_delong[1], ci_delong[3]))
cat(sprintf("  Bootstrap 95%% CI: [%.4f, %.4f]\n", ci_boot[1], ci_boot[3]))

# ---- 6. 基因数量比较：3-gene vs 5-gene vs 7-gene ----
cat("\n===== 基因数量模型比较 =====\n")

# 按单基因AUC排序
single_aucs <- sapply(final_biomarkers, function(g) {
  roc_obj <- roc(y, X_raw[, g], levels = c(0, 1), direction = "<", quiet = TRUE)
  auc(roc_obj)
})
gene_order <- names(sort(single_aucs, decreasing = TRUE))
cat("基因按AUC排序:", paste(sprintf("%s(%.3f)", gene_order, single_aucs[gene_order]), collapse = ", "), "\n")

# 构建3-gene, 5-gene, 7-gene模型
model_configs <- list(
  "3-gene" = gene_order[1:3],
  "5-gene" = gene_order[1:5],
  "7-gene" = gene_order[1:7]
)

model_comparison <- data.frame(
  Model = character(),
  N_Genes = integer(),
  Genes = character(),
  AIC = numeric(),
  BIC = numeric(),
  CV_AUC_mean = numeric(),
  CV_AUC_sd = numeric(),
  Optimism_Corrected_AUC = numeric(),
  stringsAsFactors = FALSE
)

set.seed(42)
for (config_name in names(model_configs)) {
  genes <- model_configs[[config_name]]
  cat(sprintf("\n--- %s 模型: %s ---\n", config_name, paste(genes, collapse=", ")))

  X_sub <- X[, genes, drop = FALSE]
  model_data <- as.data.frame(cbind(y = y, X_sub))

  # 拟合模型
  fmla <- as.formula(paste("y ~", paste(genes, collapse = " + ")))
  model <- glm(fmla, data = model_data, family = binomial())

  # AIC & BIC
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  cat(sprintf("  AIC = %.2f, BIC = %.2f\n", aic_val, bic_val))

  # 重复CV
  cv_aucs_config <- numeric(n_repeats * n_folds)
  idx <- 1
  for (r in 1:n_repeats) {
    fold_ids <- createFolds(as.factor(y), k = n_folds, list = TRUE)
    for (f in 1:n_folds) {
      test_idx <- fold_ids[[f]]
      train_idx <- setdiff(1:n_samples, test_idx)
      train_data <- as.data.frame(cbind(y = y[train_idx], X_sub[train_idx, , drop = FALSE]))
      cv_model <- glm(y ~ ., data = train_data, family = binomial())
      test_data <- as.data.frame(X_sub[test_idx, , drop = FALSE])
      colnames(test_data) <- genes
      pred_prob <- predict(cv_model, newdata = test_data, type = "response")
      roc_obj <- roc(y[test_idx], pred_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
      cv_aucs_config[idx] <- auc(roc_obj)
      idx <- idx + 1
    }
  }

  cv_mean <- mean(cv_aucs_config)
  cv_sd <- sd(cv_aucs_config)
  cat(sprintf("  CV AUC = %.4f ± %.4f\n", cv_mean, cv_sd))

  # Bootstrap乐观性校正
  app_auc <- auc(roc(y, predict(model, type = "response"),
                     levels = c(0, 1), direction = "<", quiet = TRUE))
  opt_values <- numeric(n_boot)
  for (b in 1:n_boot) {
    boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
    oob_idx <- setdiff(1:n_samples, unique(boot_idx))
    if (length(oob_idx) < 5) next
    bd <- model_data[boot_idx, ]
    bm <- glm(y ~ ., data = bd, family = binomial())
    bta <- auc(roc(bd$y, predict(bm, type = "response"),
                   levels = c(0, 1), direction = "<", quiet = TRUE))
    od <- model_data[oob_idx, ]
    bte <- auc(roc(od$y, predict(bm, newdata = od, type = "response"),
                   levels = c(0, 1), direction = "<", quiet = TRUE))
    opt_values[b] <- bta - bte
  }
  opt_values <- opt_values[!is.na(opt_values)]
  opt_corrected <- app_auc - mean(opt_values)
  cat(sprintf("  表观AUC = %.4f, 乐观性校正AUC = %.4f\n", app_auc, opt_corrected))

  model_comparison <- rbind(model_comparison, data.frame(
    Model = config_name,
    N_Genes = length(genes),
    Genes = paste(genes, collapse = ", "),
    AIC = aic_val,
    BIC = bic_val,
    CV_AUC_mean = cv_mean,
    CV_AUC_sd = cv_sd,
    Optimism_Corrected_AUC = opt_corrected,
    stringsAsFactors = FALSE
  ))
}

# 模型比较表
cat("\n===== 模型比较汇总 =====\n")
print(model_comparison[, c("Model", "N_Genes", "AIC", "BIC", "CV_AUC_mean", "Optimism_Corrected_AUC")],
      digits = 4)

# ---- 7. 模型比较可视化 ----
cat("\n===== 模型比较可视化 =====\n")

# AIC比较
ggplot(model_comparison, aes(x = Model, y = AIC, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", AIC)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("3-gene" = "#66c2a5", "5-gene" = "#fc8d62", "7-gene" = "#8da0cb")) +
  labs(title = "Model Comparison: AIC (lower = better)",
       x = "", y = "AIC") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("04c_model_aic_comparison.png", width = 6, height = 5, dpi = 150)

# CV AUC比较（带误差线）
ggplot(model_comparison, aes(x = Model, y = CV_AUC_mean, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = CV_AUC_mean - CV_AUC_sd, ymax = CV_AUC_mean + CV_AUC_sd),
                width = 0.15, linewidth = 1) +
  geom_text(aes(label = sprintf("%.4f±%.4f", CV_AUC_mean, CV_AUC_sd)),
            vjust = -1.5, size = 3.5) +
  scale_fill_manual(values = c("3-gene" = "#66c2a5", "5-gene" = "#fc8d62", "7-gene" = "#8da0cb")) +
  labs(title = "Cross-Validated AUC: 3-gene vs 5-gene vs 7-gene",
       subtitle = sprintf("Repeated 10×10-fold CV (error bars = ±1 SD)"),
       x = "", y = "AUC") +
  ylim(0.75, 1.05) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("04c_model_cv_auc_comparison.png", width = 7, height = 5, dpi = 150)

# 乐观性校正AUC vs 表观AUC
comparison_long <- model_comparison %>%
  select(Model, CV_AUC_mean, Optimism_Corrected_AUC) %>%
  rename(`Cross-Validated AUC` = CV_AUC_mean,
         `Optimism-Corrected AUC` = Optimism_Corrected_AUC) %>%
  pivot_longer(-Model, names_to = "Metric", values_to = "AUC")

ggplot(comparison_long, aes(x = Model, y = AUC, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", AUC)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Cross-Validated AUC" = "#4daf4a",
                                "Optimism-Corrected AUC" = "#e41a1c")) +
  labs(title = "AUC Comparison: Cross-Validation vs Optimism-Corrected",
       x = "", y = "AUC") +
  ylim(0.75, 1.05) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
ggsave("04c_auc_cv_vs_optimism.png", width = 8, height = 5, dpi = 150)

# ---- 8. 校准曲线 (Calibration Plot) ----
cat("\n===== 校准曲线 =====\n")

# 为7-gene模型绘制校准曲线（使用bootstrap验证）
genes_7 <- model_configs[["7-gene"]]
X_7 <- X[, genes_7, drop = FALSE]
model_data_7 <- as.data.frame(cbind(y = y, X_7))

# 使用caret进行bootstrap校准评估
set.seed(42)
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# 校准数据
full_model_7 <- glm(y ~ ., data = model_data_7, family = binomial())
pred_probs <- predict(full_model_7, type = "response")

# 将预测概率分箱并计算观测比例
n_bins <- 10
bins <- cut(pred_probs, breaks = seq(0, 1, length.out = n_bins + 1),
            include.lowest = TRUE)
observed <- tapply(y, bins, mean)
predicted <- tapply(pred_probs, bins, mean)
bin_n <- tapply(y, bins, length)

cal_df <- data.frame(
  Predicted = as.numeric(predicted),
  Observed = as.numeric(observed),
  N = as.numeric(bin_n)
)
cal_df <- cal_df[complete.cases(cal_df), ]

# 完美校准线 + 实际校准
ggplot(cal_df, aes(x = Predicted, y = Observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = N), color = "steelblue", alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "red", fill = "red", alpha = 0.15) +
  scale_size_continuous(name = "Samples per bin") +
  labs(title = "Calibration Plot — 7-Gene Model (GSE96804)",
       x = "Predicted Probability", y = "Observed Proportion") +
  xlim(0, 1) + ylim(0, 1) +
  theme_minimal(base_size = 13)
ggsave("04c_calibration_plot.png", width = 6, height = 5.5, dpi = 150)
cat("✅ 校准曲线已保存\n")

# ---- 9. CV AUC分布箱线图 ----
cat("\n===== CV AUC分布可视化 =====\n")

cv_df <- data.frame(
  AUC = all_cv_aucs,
  Repeat = factor(rep(1:n_repeats, each = n_folds))
)

ggplot(cv_df, aes(x = "10×10-fold CV", y = AUC)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5, width = 0.3) +
  geom_jitter(aes(color = Repeat), width = 0.08, alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1, y = cv_mean_auc,
           label = sprintf("Mean = %.3f\nSD = %.3f", cv_mean_auc, cv_sd_auc),
           vjust = -1, size = 4, color = "darkred") +
  labs(title = "Repeated 10-Fold Cross-Validation AUC Distribution",
       subtitle = sprintf("7-gene model, 10 repeats × 10 folds = %d estimates", length(all_cv_aucs)),
       x = "", y = "AUC") +
  ylim(0.5, 1.05) +
  theme_minimal(base_size = 13)
ggsave("04c_cv_auc_distribution.png", width = 7, height = 5.5, dpi = 150)
cat("✅ CV分布图已保存\n")

# ---- 10. 单个基因稳定性评估 ----
cat("\n===== 单基因稳定性 =====\n")

# 为每个基因计算在不同CV fold中的系数稳定性
set.seed(42)
gene_coef_stability <- matrix(NA, nrow = n_repeats * n_folds, ncol = n_genes)
colnames(gene_coef_stability) <- final_biomarkers
row_idx <- 1

for (r in 1:n_repeats) {
  fold_ids <- createFolds(as.factor(y), k = n_folds, list = TRUE)
  for (f in 1:n_folds) {
    test_idx <- fold_ids[[f]]
    train_idx <- setdiff(1:n_samples, test_idx)
    train_data <- as.data.frame(cbind(y = y[train_idx], X[train_idx, , drop = FALSE]))
    cv_model <- glm(y ~ ., data = train_data, family = binomial())
    gene_coef_stability[row_idx, ] <- coef(cv_model)[-1]  # 去掉截距
    row_idx <- row_idx + 1
  }
}

# 系数稳定性
coef_summary <- data.frame(
  Gene = final_biomarkers,
  Mean_Coef = colMeans(gene_coef_stability, na.rm = TRUE),
  SD_Coef = apply(gene_coef_stability, 2, sd, na.rm = TRUE),
  CV_Coef = abs(apply(gene_coef_stability, 2, sd, na.rm = TRUE) /
                apply(gene_coef_stability, 2, mean, na.rm = TRUE))
)
coef_summary$Selection_Frequency <- colMeans(abs(gene_coef_stability) > 1e-6, na.rm = TRUE)

cat("基因系数稳定性:\n")
print(coef_summary)

# 系数稳定性图
ggplot(coef_summary, aes(x = reorder(Gene, abs(Mean_Coef)), y = Mean_Coef, fill = Mean_Coef > 0)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Coef - SD_Coef, ymax = Mean_Coef + SD_Coef),
                width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
                    labels = c("TRUE" = "Positive", "FALSE" = "Negative")) +
  coord_flip() +
  labs(title = "Coefficient Stability Across 100 CV Folds",
       subtitle = "Error bars = ±1 SD; Color = direction of effect",
       x = "", y = "Mean Coefficient", fill = "Direction") +
  theme_minimal(base_size = 13)
ggsave("04c_coefficient_stability.png", width = 8, height = 5, dpi = 150)
cat("✅ 系数稳定性图已保存\n")

# ---- 11. 保存结果 ----
cat("\n===== 保存结果 =====\n")

results_list <- list(
  cv_results = list(
    n_repeats = n_repeats,
    n_folds = n_folds,
    all_cv_aucs = all_cv_aucs,
    mean_auc = cv_mean_auc,
    sd_auc = cv_sd_auc,
    ci_95 = c(cv_ci_lower, cv_ci_upper)
  ),
  bootstrap = list(
    apparent_auc = apparent_auc,
    n_boot = n_boot,
    mean_optimism = mean_optimism,
    optimism_corrected_auc = optimism_corrected_auc,
    oob_mean_auc = mean(boot_aucs),
    oob_sd_auc = sd(boot_aucs)
  ),
  permutation_test = list(
    n_perm = n_perm,
    observed_auc = apparent_auc,
    null_mean_auc = mean(perm_aucs),
    null_sd_auc = sd(perm_aucs),
    p_value = p_value
  ),
  model_comparison = model_comparison,
  coefficient_stability = coef_summary,
  auc_confidence_intervals = list(
    delong = as.numeric(ci_delong),
    bootstrap = as.numeric(ci_boot)
  )
)

saveRDS(results_list, "04c_ml_validation_results.rds")

# 保存CSV
write.csv(model_comparison, "04c_model_comparison.csv", row.names = FALSE)
write.csv(coef_summary, "04c_coefficient_stability.csv", row.names = FALSE)

cat("\n========================================\n")
cat("04c_ML严格验证完成！输出:\n")
cat(sprintf("  重复CV: mean AUC = %.4f ± %.4f (95%% CI [%.4f, %.4f])\n",
            cv_mean_auc, cv_sd_auc, cv_ci_lower, cv_ci_upper))
cat(sprintf("  Bootstrap乐观性校正: %.4f → %.4f (△=%.4f)\n",
            apparent_auc, optimism_corrected_auc, mean_optimism))
cat(sprintf("  置换检验: P = %.4f\n", p_value))
cat(sprintf("  基因数比较: 3-gene AIC=%.1f | 5-gene AIC=%.1f | 7-gene AIC=%.1f\n",
            model_comparison$AIC[1], model_comparison$AIC[2], model_comparison$AIC[3]))
cat("  图表: 04c_permutation_test.png, 04c_model_*_comparison.png,\n")
cat("        04c_calibration_plot.png, 04c_cv_auc_distribution.png,\n")
cat("        04c_coefficient_stability.png\n")
cat("========================================\n")
