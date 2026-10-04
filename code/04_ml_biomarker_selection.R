# ============================================================
# 双硫死亡+DN 竞赛项目 — 04_机器学习筛选标志物
# 方法：LASSO + SVM-RFE + Random Forest 三法取交集
# 输入：03_candidate_genes (DR-DEGs ∩ WGCNA) 或 02_dr_degs_sig
# 输出：最优诊断标志物基因组合
# ============================================================

library(glmnet)        # LASSO
library(e1071)         # SVM
library(randomForest)  # 随机森林
library(caret)         # 统一ML框架
library(pROC)          # ROC曲线
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")

# 加载基因级表达矩阵和分组
expr_norm <- readRDS("GSE96804_gene_expr.rds")
pheno <- readRDS("GSE96804_pheno.rds")

# 构建分组 (DN=1, Control=0) - 使用title列
group_info <- pheno[["title"]]
group <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"),
                levels = c("Control", "DN"))
cat(sprintf("分组: DN=%d, Control=%d\n", sum(group == "DN"), sum(group == "Control")))

# 加载候选基因
if (file.exists("03_candidate_genes.rds")) {
  candidate_genes <- readRDS("03_candidate_genes.rds")
  cat(sprintf("从WGCNA+DR-DEGs交集加载了 %d 个候选基因\n", length(candidate_genes)))
} else if (file.exists("02_dr_degs_sig.rds")) {
  dr_degs_sig <- readRDS("02_dr_degs_sig.rds")
  candidate_genes <- dr_degs_sig$Gene
  cat(sprintf("从显著DR-DEGs加载了 %d 个候选基因\n", length(candidate_genes)))
} else {
  # 回退：使用disulfidptosis_genes.txt中在表达矩阵中的基因
  drg_list <- readLines("disulfidptosis_genes.txt")
  deg_results <- readRDS("02_deg_results.rds")
  dr_degs <- deg_results[deg_results$Gene %in% drg_list, ]
  dr_degs_sig <- subset(dr_degs, adj.P.Val < 0.05 & abs(logFC) > 0.5)
  candidate_genes <- dr_degs_sig$Gene
  cat(sprintf("从差异分析回退加载了 %d 个候选基因\n", length(candidate_genes)))
}

if (length(candidate_genes) < 5) {
  cat("警告：候选基因不足5个，可能需要放宽筛选标准\n")
  # 扩大候选范围
  deg_results <- readRDS("02_deg_results.rds")
  drg_list <- readLines("disulfidptosis_genes.txt")
  all_dr_degs <- deg_results[deg_results$Gene %in% drg_list, ]
  all_dr_degs <- all_dr_degs[order(all_dr_degs$adj.P.Val), ]
  candidate_genes <- head(all_dr_degs$Gene, 30)
  cat(sprintf("扩大候选基因池至 %d 个\n", length(candidate_genes)))
}

# 准备训练矩阵
train_genes <- intersect(candidate_genes, rownames(expr_norm))
cat(sprintf("最终用于ML的基因数: %d\n", length(train_genes)))

X_raw <- as.data.frame(t(expr_norm[train_genes, , drop = FALSE]))
y <- as.numeric(group == "DN")  # 1=DN, 0=Control

# 标准化
X <- scale(X_raw)

cat(sprintf("训练数据: %d 样本 x %d 基因, DN比例=%.1f%%\n",
            nrow(X), ncol(X), mean(y) * 100))

# ---- 2. LASSO回归 ----
cat("\n===== LASSO回归 =====\n")

set.seed(42)
lasso_cv <- cv.glmnet(X, y, alpha = 1, family = "binomial",
                      nfolds = 10, type.measure = "deviance")

# 最佳lambda
cat(sprintf("lambda.min = %.4f, lambda.1se = %.4f\n",
            lasso_cv$lambda.min, lasso_cv$lambda.1se))

# 提取LASSO系数（用lambda.min）
lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_coef_mat <- as.matrix(lasso_coef)
lasso_genes <- rownames(lasso_coef_mat)[lasso_coef_mat[, 1] != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
cat(sprintf("LASSO筛选出 %d 个基因: %s\n",
            length(lasso_genes), paste(lasso_genes, collapse = ", ")))

# LASSO路径图
png("04_lasso_cv.png", width = 1000, height = 500, res = 120)
par(mfrow = c(1, 2))
plot(lasso_cv, main = "LASSO Cross-Validation")
plot(lasso_cv$glmnet.fit, xvar = "lambda", label = TRUE, main = "LASSO Coefficient Path")
abline(v = log(lasso_cv$lambda.min), lty = 2, col = "red")
abline(v = log(lasso_cv$lambda.1se), lty = 2, col = "blue")
dev.off()
cat("✅ LASSO图已保存\n")

# ---- 3. SVM-RFE（递归特征消除）----
cat("\n===== SVM-RFE =====\n")

# 自定义SVM-RFE函数
svm_rfe <- function(X, y, k = 10) {
  n_features <- ncol(X)
  features <- colnames(X)
  ranks <- numeric(n_features)
  names(ranks) <- features
  remaining <- features

  for (i in seq_len(n_features - 1)) {
    # 训练SVM
    X_sub <- X[, remaining, drop = FALSE]

    # 5折CV评估
    folds <- createFolds(y, k = min(k, floor(nrow(X) / 2)), list = FALSE)
    cv_errors <- numeric(max(folds))
    for (f in seq_len(max(folds))) {
      test_idx <- folds == f
      train_idx <- !test_idx
      svm_model <- svm(X_sub[train_idx, , drop = FALSE], as.factor(y[train_idx]),
                        kernel = "linear", scale = FALSE, probability = TRUE)
      pred <- predict(svm_model, X_sub[test_idx, , drop = FALSE])
      cv_errors[f] <- mean(pred != y[test_idx])
    }

    # 计算权重并移除最小权重特征
    svm_full <- svm(X_sub, as.factor(y), kernel = "linear", scale = FALSE)
    w <- t(svm_full$coefs) %*% svm_full$SV
    w_abs <- abs(as.vector(w))
    names(w_abs) <- remaining

    # 标记被移除的基因
    worst <- names(which.min(w_abs))
    ranks[worst] <- n_features - i + 1
    remaining <- setdiff(remaining, worst)
  }
  # 最后一个基因
  ranks[remaining] <- 1
  return(sort(ranks))
}

set.seed(42)
svm_ranks <- svm_rfe(X, y)
cat(sprintf("SVM-RFE完成，基因排名:\n"))
print(head(svm_ranks, 15))

# 取排名前15的基因
svm_top_n <- min(15, length(svm_ranks))
svm_genes <- names(svm_ranks)[1:svm_top_n]
cat(sprintf("SVM-RFE前%d个基因: %s\n",
            svm_top_n, paste(svm_genes, collapse = ", ")))

# SVM-RFE排名图
svm_rank_df <- data.frame(
  Gene = names(svm_ranks),
  Rank = svm_ranks
)
svm_rank_df$Gene <- factor(svm_rank_df$Gene, levels = rev(names(svm_ranks)))

ggplot(svm_rank_df[1:20, ], aes(x = Gene, y = Rank)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "SVM-RFE Gene Ranking (Top 20)",
       x = "", y = "Rank (lower = more important)") +
  theme_minimal(base_size = 12)
ggsave("04_svm_rfe_rank.png", width = 8, height = 6, dpi = 150)
cat("✅ SVM-RFE排名图已保存\n")

# ---- 4. Random Forest ----
cat("\n===== Random Forest =====\n")

set.seed(42)
rf_model <- randomForest(x = X, y = as.factor(y),
                         ntree = 500, importance = TRUE,
                         mtry = max(1, floor(sqrt(ncol(X)))))

# 按MeanDecreaseGini排序
rf_importance <- importance(rf_model)
rf_importance_df <- data.frame(
  Gene = rownames(rf_importance),
  MeanDecreaseGini = rf_importance[, "MeanDecreaseGini"]
)
rf_importance_df <- rf_importance_df[order(rf_importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# 取重要的基因（MeanDecreaseGini > 中位数的1.5倍）
rf_threshold <- median(rf_importance_df$MeanDecreaseGini) * 1.5
rf_genes <- rf_importance_df$Gene[rf_importance_df$MeanDecreaseGini > rf_threshold]
# 至少取前10个
if (length(rf_genes) < 10) {
  rf_genes <- head(rf_importance_df$Gene, 10)
}
cat(sprintf("RF筛选出 %d 个重要基因 (>%.3f Gini)\n",
            length(rf_genes), rf_threshold))
cat(sprintf("RF Top基因: %s\n",
            paste(head(rf_genes, 10), collapse = ", ")))

# RF重要性图
ggplot(rf_importance_df[1:20, ], aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#4daf4a", alpha = 0.8) +
  coord_flip() +
  labs(title = "Random Forest - Top 20 Important Genes",
       x = "", y = "Mean Decrease Gini") +
  theme_minimal(base_size = 12)
ggsave("04_rf_importance.png", width = 8, height = 6, dpi = 150)
cat("✅ RF重要性图已保存\n")

# ---- 5. 三法取交集 ----
cat("\n===== 三法取交集 =====\n")

final_genes <- Reduce(intersect, list(lasso_genes, svm_genes, rf_genes))
cat(sprintf("LASSO: %d 基因 | SVM-RFE: %d 基因 | RF: %d 基因\n",
            length(lasso_genes), length(svm_genes), length(rf_genes)))
cat(sprintf("三法交集: %d 个基因\n", length(final_genes)))

if (length(final_genes) == 0) {
  cat("三法无重叠，采用两两交集或宽松策略\n")
  # 尝试两两交集
  intersect_lasso_svm <- intersect(lasso_genes, svm_genes)
  intersect_lasso_rf <- intersect(lasso_genes, rf_genes)
  intersect_svm_rf <- intersect(svm_genes, rf_genes)

  # 合并所有两两交集
  final_genes <- unique(c(intersect_lasso_svm, intersect_lasso_rf, intersect_svm_rf))
  cat(sprintf("两两交集合并: %d 个基因\n", length(final_genes)))
}

if (length(final_genes) == 0) {
  cat("仍无交集，取LASSO基因作为最终标志物\n")
  final_genes <- lasso_genes
}

cat(sprintf("最终标志物基因 (%d个): %s\n",
            length(final_genes), paste(final_genes, collapse = ", ")))

# Venn图
library(pheatmap)

# 手动画Venn（不依赖VennDiagram包）
venn_data <- list(
  LASSO = lasso_genes,
  `SVM-RFE` = svm_genes,
  `RF` = rf_genes
)

# 计算交集大小
all_genes <- unique(c(lasso_genes, svm_genes, rf_genes))
intersection_matrix <- sapply(venn_data, function(x) all_genes %in% x)
colnames(intersection_matrix) <- names(venn_data)

# 保存Venn信息
cat("\n基因集重叠信息:\n")
cat(sprintf("  LASSO ∩ SVM-RFE: %d 基因\n",
            length(intersect(lasso_genes, svm_genes))))
cat(sprintf("  LASSO ∩ RF: %d 基因\n",
            length(intersect(lasso_genes, rf_genes))))
cat(sprintf("  SVM-RFE ∩ RF: %d 基因\n",
            length(intersect(svm_genes, rf_genes))))
cat(sprintf("  三者交集: %d 基因\n",
            length(Reduce(intersect, list(lasso_genes, svm_genes, rf_genes)))))

# ---- 6. 标志物评估——ROC曲线 ----
cat("\n===== 标志物ROC评估 =====\n")

if (length(final_genes) > 0) {
  # 单个基因ROC
  roc_list <- list()
  auc_values <- numeric(length(final_genes))
  names(auc_values) <- final_genes

  png("04_ROC_curves.png", width = 1000, height = 800, res = 120)
  plot(NULL, xlim = c(1, 0), ylim = c(0, 1),
       xlab = "Specificity", ylab = "Sensitivity",
       main = "ROC Curves - Final Biomarkers (GSE96804)")
  abline(a = 1, b = -1, lty = 2, col = "grey")

  colors <- brewer.pal(min(8, length(final_genes)), "Set1")
  for (i in seq_along(final_genes)) {
    gene <- final_genes[i]
    roc_obj <- roc(y, X_raw[, gene], levels = c(0, 1), direction = "<")
    roc_list[[gene]] <- roc_obj
    auc_values[i] <- auc(roc_obj)
    lines(1 - roc_obj$specificities, roc_obj$sensitivities,
          col = colors[i], lwd = 2)
  }

  legend("bottomright",
         legend = sprintf("%s (AUC=%.3f)", final_genes, auc_values),
         col = colors, lwd = 2, cex = 0.8)
  dev.off()
  cat("✅ 单基因ROC曲线已保存\n")

  # 多基因逻辑回归模型
  if (length(final_genes) >= 2) {
    model_data <- as.data.frame(X_raw[, final_genes, drop = FALSE])
    model_data$y <- y

    # 构建逻辑回归模型
    fmla <- as.formula(paste("y ~", paste(final_genes, collapse = " + ")))
    logit_model <- glm(fmla, data = model_data, family = binomial())

    # 预测概率
    pred_prob <- predict(logit_model, type = "response")

    # 组合ROC
    combined_roc <- roc(y, pred_prob, levels = c(0, 1), direction = "<")
    combined_auc <- auc(combined_roc)
    cat(sprintf("多基因组合模型 AUC: %.3f\n", combined_auc))

    # 输出每个基因的AUC
    cat("\n各标志物AUC:\n")
    for (i in seq_along(final_genes)) {
      cat(sprintf("  %s: AUC = %.3f\n", final_genes[i], auc_values[i]))
    }
    cat(sprintf("  组合模型: AUC = %.3f\n", combined_auc))
  }
}

# ---- 7. 标志物表达箱线图 ----
cat("\n===== 绘制标志物表达箱线图 =====\n")

if (length(final_genes) > 0) {
  # 准备数据
  plot_data <- data.frame()
  for (gene in final_genes) {
    gene_data <- data.frame(
      Gene = gene,
      Expression = as.numeric(X_raw[, gene]),
      Group = group
    )
    plot_data <- rbind(plot_data, gene_data)
  }

  ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ Gene, scales = "free_y") +
    stat_compare_means(aes(group = Group), label = "p.signif",
                       method = "t.test", label.y.npc = 0.95) +
    scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
    labs(title = "Final Biomarkers Expression (GSE96804)",
         y = "Expression") +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(face = "bold"))
  ggsave("04_biomarkers_boxplot.png", width = max(8, length(final_genes) * 3),
         height = 6, dpi = 150)
  cat("✅ 标志物表达箱线图已保存\n")
}

# ---- 8. 标志物相关性热图 ----
cat("\n===== 标志物相关性 =====\n")

if (length(final_genes) >= 2) {
  cor_matrix <- cor(X_raw[, final_genes], method = "spearman")

  png("04_biomarkers_correlation.png", width = 800, height = 700, res = 120)
  pheatmap(cor_matrix,
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = "Biomarker Correlation (Spearman)",
           fontsize_number = 10)
  dev.off()
  cat("✅ 标志物相关性热图已保存\n")
}

# ---- 9. 保存结果 ----
cat("\n===== 保存结果 =====\n")
saveRDS(lasso_genes, "04_lasso_genes.rds")
saveRDS(svm_genes, "04_svm_genes.rds")
saveRDS(rf_genes, "04_rf_genes.rds")
saveRDS(final_genes, "04_final_biomarkers.rds")
write.csv(data.frame(
  Gene = final_genes,
  AUC = auc_values[final_genes]
), "04_final_biomarkers.csv", row.names = FALSE)

if (exists("combined_auc")) {
  cat(sprintf("\n最终组合模型AUC: %.3f\n", combined_auc))
}

cat("\n========================================\n")
cat("04_机器学习筛选完成！输出:\n")
cat(sprintf("  最终标志物 (%d个): %s\n",
            length(final_genes), paste(final_genes, collapse = ", ")))
cat("  04_lasso_genes.rds / svm_genes.rds / rf_genes.rds\n")
cat("  04_final_biomarkers.rds / .csv\n")
cat("  04_lasso_cv.png — LASSO交叉验证图\n")
cat("  04_svm_rfe_rank.png — SVM-RFE排名图\n")
cat("  04_rf_importance.png — RF重要性图\n")
cat("  04_ROC_curves.png — ROC曲线\n")
cat("  04_biomarkers_boxplot.png — 标志物表达箱线图\n")
cat("  04_biomarkers_correlation.png — 标志物相关性热图\n")
cat("========================================\n")
