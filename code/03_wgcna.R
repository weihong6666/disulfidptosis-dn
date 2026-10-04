# ============================================================
# 双硫死亡+DN 竞赛项目 — 03_WGCNA共表达网络
# 输入：GSE96804标准化表达矩阵
# 输出：关键模块基因
# ============================================================

library(WGCNA)
library(ggplot2)

# 允许WGCNA多线程
allowWGCNAThreads()

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")
# 使用基因级表达矩阵
if (file.exists("GSE96804_gene_expr.rds")) {
  expr_norm <- readRDS("GSE96804_gene_expr.rds")
  cat("加载基因级表达矩阵\n")
} else {
  expr_norm <- readRDS("GSE96804_expr_normalized.rds")
  cat("加载探针级表达矩阵\n")
}
pheno <- readRDS("GSE96804_pheno.rds")

# ---- 2. 准备WGCNA输入 ----
cat("\n===== 准备数据 =====\n")

# 先筛选高变异基因（减少内存压力，标准WGCNA做法）
n_top <- 5000  # 取top 5000高变异基因
gene_vars <- apply(expr_norm, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_top, length(gene_vars))]
expr_wgcna <- expr_norm[top_genes, ]
cat(sprintf("筛选Top %d 高变异基因用于WGCNA\n", length(top_genes)))

# 转置：基因在列，样本在行（WGCNA要求）
datExpr0 <- as.data.frame(t(expr_wgcna))

# 检查缺失值
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  cat("已移除不良基因和样本\n")
}
cat(sprintf("WGCNA输入: %d 样本 x %d 基因\n", nrow(datExpr0), ncol(datExpr0)))

# ---- 3. 样本聚类（检测离群样本） ----
cat("\n===== 样本聚类 =====\n")

sampleTree <- hclust(dist(datExpr0), method = "average")

png("03_sample_clustering.png", width = 1200, height = 600, res = 120)
par(cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
abline(h = 150, col = "red", lty = 2)  # 参考线
dev.off()
cat("✅ 样本聚类图已保存\n")

# ---- 4. 软阈值选择（关键步骤！） ----
cat("\n===== 软阈值(power)优化 =====\n")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# 绘制软阈值曲线
png("03_soft_threshold.png", width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
cex1 <- 0.9

# 1. Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence",
     type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "blue", lty = 2)  # R^2 = 0.85 参考线

# 2. Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity",
     type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()
cat("✅ 软阈值图已保存\n")

# 确定最佳power（第一个达到R^2>0.85的power，且连通性不过低）
best_power <- sft$powerEstimate
if (is.na(best_power)) {
  # 选择R^2>0.8的第一个
  idx <- which(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] > 0.8)
  if (length(idx) > 0) {
    best_power <- min(idx)
  } else {
    # 最后回退：选R^2最高的power
    best_power <- sft$fitIndices$Power[which.max(
      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    )]
    cat(sprintf("⚠ 警告：无power达到R^2>0.8，使用最佳power = %d (R^2=%.3f)\n",
                best_power,
                max(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], na.rm = TRUE)))
  }
}
cat(sprintf("最佳软阈值 power = %d\n", best_power))

# ---- 5. 构建共表达网络 ----
cat("\n===== 构建网络（一步法）=====\n")

net <- blockwiseModules(
  datExpr0,
  power = best_power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "03_GSE96804_TOM",
  verbose = 3
)

# 模块颜色
moduleColors <- labels2colors(net$colors)
cat(sprintf("识别到 %d 个模块\n", length(unique(moduleColors))))

# ---- 6. 模块-性状关联 ----
cat("\n===== 模块-性状关联 =====\n")

# 性状：DN = 1, Control = 0 (使用title列判断)
group_info <- pheno[["title"]]
trait_dn <- as.numeric(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE))
datTraits <- data.frame(DN = trait_dn)
cat(sprintf("性状分布: DN=%d, Control=%d\n", sum(trait_dn), sum(!trait_dn)))

# 计算模块特征基因
MEs0 <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# 模块-性状相关性
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))

# 可视化
png("03_module_trait_heatmap.png", width = 800, height = 1200, res = 120)
par(mar = c(4, 10, 3, 2))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "DN",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(p=",
                     signif(moduleTraitPvalue, 2), ")", sep = ""),
  setStdMargins = FALSE,
  cex.text = 0.7,
  main = "Module-Trait Relationships"
)
dev.off()
cat("✅ 模块-性状热图已保存\n")

# ---- 7. 提取关键模块基因 ----
cat("\n===== 提取关键模块基因 =====\n")

# 找与DN最相关（p<0.05）的模块
sig_modules <- names(which(moduleTraitPvalue[, "DN"] < 0.05))
cat(sprintf("与DN显著相关的模块 (%d个): %s\n",
            length(sig_modules), paste(sig_modules, collapse = ", ")))

# 去掉ME前缀得到颜色名
sig_module_colors <- gsub("^ME", "", sig_modules)
cat(sprintf("对应颜色: %s\n", paste(sig_module_colors, collapse = ", ")))

# 提取关键模块基因
key_module_genes <- colnames(datExpr0)[moduleColors %in% sig_module_colors]
cat(sprintf("关键模块基因总数: %d\n", length(key_module_genes)))

# 按相关性排序
module_gene_cor <- cor(datExpr0[, key_module_genes, drop = FALSE], datTraits$DN)
key_module_genes_sorted <- key_module_genes[order(abs(module_gene_cor), decreasing = TRUE)]

# ---- 8. 与DR-DEGs取交集 ----
cat("\n===== 与DR-DEGs取交集 =====\n")

if (file.exists("02_dr_degs_sig.rds")) {
  dr_degs_sig <- readRDS("02_dr_degs_sig.rds")
} else {
  # 如果还没有差异分析结果，用全部DR基因
  dr_genes <- readLines("disulfidptosis_genes.txt")
  dr_degs_sig <- data.frame(Gene = intersect(dr_genes, rownames(expr_norm)))
}

candidate_genes <- intersect(key_module_genes, dr_degs_sig$Gene)
cat(sprintf("DR-DEGs ∩ WGCNA关键模块 = %d 个候选基因\n", length(candidate_genes)))
if (length(candidate_genes) > 0) {
  cat("候选基因:", paste(candidate_genes, collapse = ", "), "\n")
}

# ---- 9. 保存 ----
cat("\n===== 保存结果 =====\n")
saveRDS(net, "03_wgcna_net.rds")
saveRDS(moduleColors, "03_module_colors.rds")
saveRDS(key_module_genes_sorted, "03_key_module_genes.rds")
saveRDS(candidate_genes, "03_candidate_genes.rds")
write.csv(data.frame(Gene = key_module_genes_sorted), "03_key_module_genes.csv", row.names = FALSE)
if (length(candidate_genes) > 0) {
  write.csv(data.frame(Gene = candidate_genes), "03_candidate_genes.csv", row.names = FALSE)
}

cat("\n========================================\n")
cat("03_WGCNA分析完成！输出:\n")
cat("  03_wgcna_net.rds — 完整网络对象\n")
cat("  03_key_module_genes.rds — 关键模块基因\n")
cat("  03_candidate_genes.csv — DR-DEGs∩WGCNA候选基因\n")
cat("  03_soft_threshold.png — 软阈值曲线\n")
cat("  03_module_trait_heatmap.png — 模块-性状关联\n")
cat("========================================\n")
