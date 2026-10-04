# ============================================================
# 双硫死亡+DN 竞赛项目 — 06_Nomogram+校准+DCA
# 构建列线图诊断模型 + 校准曲线 + DCA决策曲线
# ============================================================

library(rms)          # Nomogram + 校准曲线
library(pROC)         # ROC
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
# 注：DCA决策曲线使用自主实现，不依赖rmda

# ---- 1. 加载数据 ----
cat("\n===== 加载数据和模型 =====\n")

final_genes <- readRDS("04_final_biomarkers.rds")
expr_norm <- readRDS("GSE96804_gene_expr.rds")
pheno <- readRDS("GSE96804_pheno.rds")

# 分组 - 使用title列
group_info <- pheno[["title"]]
group <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"), levels = c("Control", "DN"))

# 准备数据
X <- as.data.frame(t(expr_norm[final_genes, , drop = FALSE]))
X$DN_status <- as.numeric(group == "DN")  # 1=DN, 0=Control
X$DN_factor <- group

cat(sprintf("标志物: %s\n", paste(final_genes, collapse = ", ")))
cat(sprintf("样本: %d DN + %d Control\n",
            sum(group == "DN"), sum(group == "Control")))

# ---- 2. Nomogram (列线图) ----
cat("\n===== 构建Nomogram =====\n")

# 使用rms包构建
dd <- datadist(X)
options(datadist = "dd")

# 构建逻辑回归模型 (rms版本)
fmla_rms <- as.formula(paste("DN_status ~", paste(final_genes, collapse = " + ")))
lrm_model <- lrm(fmla_rms, data = X, x = TRUE, y = TRUE)
cat(sprintf("LRM模型 C-index: %.3f\n", lrm_model$stats["C"]))

# 生成Nomogram
png("06_nomogram.png", width = 1000, height = 700, res = 120)
nom <- nomogram(lrm_model, fun = plogis,
                fun.at = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99),
                funlabel = "Probability of DN",
                lp = FALSE)
plot(nom, xfrac = 0.4, cex.axis = 0.8, cex.var = 0.9,
     main = "Nomogram for Diabetic Nephropathy Diagnosis")
dev.off()
cat("✅ Nomogram已保存\n")

# 如果rms nomogram失败，用ggplot2手绘简化版
if (FALSE) {  # 备用方案
  # 提取模型系数
  coefs <- coef(lrm_model)
  intercept <- coefs[1]
  gene_coefs <- coefs[-1]

  cat("\nNomogram系数:\n")
  cat(sprintf("  Intercept: %.4f\n", intercept))
  for (i in seq_along(gene_coefs)) {
    cat(sprintf("  %s: %.4f\n", names(gene_coefs)[i], gene_coefs[i]))
  }
}

# ---- 3. 校准曲线 (Calibration Curve) ----
cat("\n===== 校准曲线 =====\n")

# 重抽样校准 (Bootstrap=1000)
set.seed(42)
cal <- calibrate(lrm_model, method = "boot", B = 1000)

png("06_calibration.png", width = 800, height = 700, res = 120)
plot(cal,
     xlab = "Predicted Probability",
     ylab = "Observed DN Proportion",
     main = "Calibration Curve (Bootstrap B=1000)",
     subtitles = FALSE)
abline(a = 0, b = 1, lty = 2, col = "grey50")
legend("bottomright",
       c("Apparent", "Bias-corrected", "Ideal"),
       col = c("black", "blue", "grey50"),
       lty = c(1, 1, 2), lwd = 2, bty = "n")
dev.off()
cat("✅ 校准曲线已保存\n")

# ---- 4. DCA 决策曲线 (自主实现，不依赖rmda) ----
cat("\n===== DCA决策曲线 =====\n")

# DCA核心函数：计算净获益
calc_net_benefit <- function(y_true, y_pred_prob, thresholds = seq(0, 1, by = 0.01)) {
  n <- length(y_true)
  nb_model <- numeric(length(thresholds))
  nb_all <- numeric(length(thresholds))

  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    y_pred_class <- as.numeric(y_pred_prob >= pt)

    TP <- sum(y_pred_class == 1 & y_true == 1)
    FP <- sum(y_pred_class == 1 & y_true == 0)

    # Net Benefit = TP/N - FP/N * (pt/(1-pt))
    nb_model[i] <- TP / n - (FP / n) * (pt / (1 - pt))

    # Treat all: NB = prevalence - (1-prevalence) * (pt/(1-pt))
    prevalence <- mean(y_true)
    nb_all[i] <- prevalence - (1 - prevalence) * (pt / (1 - pt))
  }

  return(data.frame(
    threshold = thresholds,
    net_benefit = nb_model,
    treat_all = nb_all,
    treat_none = 0
  ))
}

# 计算组合模型的预测概率
pred_prob <- plogis(predict(lrm_model))

# DCA计算
dca_result <- calc_net_benefit(X$DN_status, pred_prob, thresholds = seq(0, 1, by = 0.01))

# DCA图 (组合模型)
png("06_DCA.png", width = 900, height = 700, res = 120)
plot(dca_result$threshold, dca_result$net_benefit,
     type = "l", lwd = 2.5, col = "darkred",
     xlim = c(0, 0.6), ylim = c(-0.05, 0.5),
     xlab = "Threshold Probability", ylab = "Net Benefit",
     main = "Decision Curve Analysis")
lines(dca_result$threshold, dca_result$treat_all,
      lwd = 2, col = "darkgreen", lty = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topright",
       legend = c("Combined Model", "Treat All", "Treat None"),
       col = c("darkred", "darkgreen", "grey50"),
       lty = c(1, 2, 3), lwd = c(2.5, 2, 1), bty = "n")
dev.off()
cat("✅ DCA决策曲线已保存 (自主实现)\n")

# DCA (含单基因对比)
if (length(final_genes) <= 8) {
  # 单基因逻辑回归DCA
  dca_gene_colors <- c("darkred", RColorBrewer::brewer.pal(min(8, length(final_genes)), "Set1"))

  png("06_DCA_with_genes.png", width = 1000, height = 750, res = 120)
  plot(dca_result$threshold, dca_result$net_benefit,
       type = "l", lwd = 3, col = dca_gene_colors[1],
       xlim = c(0, 0.6), ylim = c(-0.05, 0.5),
       xlab = "Threshold Probability", ylab = "Net Benefit",
       main = "DCA - Model vs Individual Genes")

  legend_labels <- c("Combined Model")
  legend_cols <- c(dca_gene_colors[1])
  legend_lty <- c(1)
  legend_lwd <- c(3)

  for (j in seq_along(final_genes)) {
    gene <- final_genes[j]
    fmla_single <- as.formula(paste("DN_status ~", gene))
    gene_model <- glm(fmla_single, data = X, family = binomial())
    gene_pred <- predict(gene_model, type = "response")
    gene_dca <- calc_net_benefit(X$DN_status, gene_pred, thresholds = seq(0, 1, by = 0.01))

    lines(gene_dca$threshold, gene_dca$net_benefit,
          col = dca_gene_colors[j + 1], lwd = 1.5, lty = 2)

    legend_labels <- c(legend_labels, gene)
    legend_cols <- c(legend_cols, dca_gene_colors[j + 1])
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 1.5)
  }

  lines(dca_result$threshold, dca_result$treat_all,
        lwd = 2, col = "darkgreen", lty = 3)
  abline(h = 0, lty = 3, col = "grey50")

  legend_labels <- c(legend_labels, "Treat All", "Treat None")
  legend_cols <- c(legend_cols, "darkgreen", "grey50")
  legend_lty <- c(legend_lty, 3, 3)
  legend_lwd <- c(legend_lwd, 2, 1)

  legend("topright", legend = legend_labels,
         col = legend_cols, lty = legend_lty, lwd = legend_lwd,
         cex = 0.7, bty = "n")
  dev.off()
  cat("✅ DCA含单基因对比图已保存\n")
}
dca_colors <- c("darkred", "steelblue", "forestgreen", "orange", "purple")

# ---- 5. 预测概率分布 ----
cat("\n===== 预测概率分布 =====\n")

# 计算预测概率
pred_prob <- plogis(predict(lrm_model))

prob_df <- data.frame(
  Sample = seq_along(pred_prob),
  Probability = pred_prob,
  Group = group
)

ggplot(prob_df, aes(x = Group, y = Probability, color = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  stat_compare_means(aes(group = Group), label = "p.signif", label.y = 1.05) +
  scale_color_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Predicted DN Probability Distribution",
       y = "Predicted Probability") +
  theme_minimal(base_size = 13)
ggsave("06_predicted_probability.png", width = 7, height = 6, dpi = 150)
cat("✅ 预测概率分布图已保存\n")

# ---- 6. 混淆矩阵 (最佳截断值) ----
cat("\n===== 混淆矩阵 =====\n")

# Youden指数找最佳截断值
train_roc <- roc(X$DN_status, pred_prob, levels = c(0, 1), direction = "<")
best_coords <- coords(train_roc, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))

cat("最佳截断值 (Youden):\n")
cat(sprintf("  Threshold: %.3f\n", best_coords$threshold))
cat(sprintf("  Sensitivity: %.3f\n", best_coords$sensitivity))
cat(sprintf("  Specificity: %.3f\n", best_coords$specificity))
cat(sprintf("  PPV: %.3f\n", best_coords$ppv))
cat(sprintf("  NPV: %.3f\n", best_coords$npv))

# 预测类别
pred_class <- ifelse(pred_prob > best_coords$threshold, "DN", "Control")
confusion <- table(Predicted = pred_class, Actual = X$DN_factor)
cat("\n混淆矩阵:\n")
print(confusion)

accuracy <- sum(diag(confusion)) / sum(confusion)
cat(sprintf("Accuracy: %.3f\n", accuracy))

# ---- 7. 保存 ----
cat("\n===== 保存结果 =====\n")
saveRDS(lrm_model, "06_lrm_model.rds")
saveRDS(best_coords, "06_best_threshold.rds")

cat("\n========================================\n")
cat("06_Nomogram分析完成！输出:\n")
cat("  06_nomogram.png — 列线图\n")
cat("  06_calibration.png — 校准曲线\n")
cat("  06_DCA.png — DCA决策曲线\n")
cat("  06_predicted_probability.png — 预测概率分布\n")
cat("  06_lrm_model.rds — 逻辑回归模型\n")
cat("  06_best_threshold.rds — 最佳截断值\n")
cat("========================================\n")
