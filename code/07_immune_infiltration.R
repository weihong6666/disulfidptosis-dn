# ============================================================
# 双硫死亡+DN 竞赛项目 — 07_免疫浸润分析
# 方法：CIBERSORT + ssGSEA + 免疫检查点差异分析
# ============================================================

library(limma)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(GSVA)         # ssGSEA
library(GSEABase)
# immunedeconv 条件加载（见CIBERSORT分析部分）
# 如果未安装: BiocManager::install("immunedeconv")
library(tidyr)
library(dplyr)

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")

expr_norm <- readRDS("GSE96804_gene_expr.rds")
pheno <- readRDS("GSE96804_pheno.rds")
final_genes <- readRDS("04_final_biomarkers.rds")

group_info <- pheno[["title"]]
group <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"), levels = c("Control", "DN"))

cat(sprintf("样本: %d DN + %d Control\n", sum(group == "DN"), sum(group == "Control")))

# ---- 2. CIBERSORT 免疫细胞浸润 ----
cat("\n===== CIBERSORT分析 =====\n")

# CIBERSORT分析：尝试immunedeconv包
# 注意：CIBERSORT需要线性（非log）TPM值；log2标准化数据结果仅供参考
immune_fractions <- NULL
if (requireNamespace("immunedeconv", quietly = TRUE)) {
  library(immunedeconv)
  tryCatch({
    # immunedeconv::deconvolute 支持的算法: cibersort_abs, mcp_counter, estimate, etc.
    # 对log2数据进行反卷积（算法内部有适应）
    deconv_results <- deconvolute(as.matrix(expr_norm), method = "cibersort_abs")

    if (!is.null(deconv_results) && nrow(deconv_results) > 1) {
      # deconv_results: 第1列是细胞类型，后续列是样本
      immune_fractions <- as.data.frame(t(deconv_results[, -1, drop = FALSE]))
      colnames(immune_fractions) <- deconv_results[, 1]
      rownames(immune_fractions) <- colnames(expr_norm)

      # 移除全0列
      immune_fractions <- immune_fractions[, colSums(immune_fractions) > 0, drop = FALSE]

      cat(sprintf("CIBERSORT识别到 %d 种免疫细胞\n", ncol(immune_fractions)))
      cat("细胞类型:", paste(colnames(immune_fractions), collapse = ", "), "\n")
    } else {
      cat("⚠ CIBERSORT返回结果无效\n")
    }
  }, error = function(e) {
    cat("⚠ CIBERSORT失败:", conditionMessage(e), "\n")
  })
}
if (is.null(immune_fractions)) {
  cat("将使用ssGSEA进行免疫浸润分析\n")
}

# ---- 3. ssGSEA 免疫特征评分 ----
cat("\n===== ssGSEA分析 =====\n")

# 28种免疫细胞基因集 (来自文献：Charoentong et al. 2017)
immune_cell_sets <- list(
  "Activated CD8 T cell" = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "TNF", "GNLY"),
  "Effector memory CD8 T cell" = c("CD8A", "IL7R", "KLRG1", "GZMK", "GZMA", "IFNG", "CCL5", "CX3CR1"),
  "Central memory CD8 T cell" = c("CD8A", "CD27", "CD28", "CCR7", "SELL", "IL7R", "TCF7"),
  "Activated CD4 T cell" = c("CD4", "CD40LG", "IFNG", "IL2", "TNF", "CD69", "ICOS"),
  "Effector memory CD4 T cell" = c("CD4", "IL7R", "KLRG1", "IFNG", "CCR5", "CXCR3", "TBX21"),
  "T follicular helper cell" = c("BCL6", "CXCR5", "ICOS", "IL21", "PDCD1", "CD200", "SLAMF6"),
  "Th1 cell" = c("TBX21", "IFNG", "STAT1", "STAT4", "IL12RB2", "CXCR3", "CCR5"),
  "Th2 cell" = c("GATA3", "IL4", "IL5", "IL13", "STAT6", "CCR3", "CCR4"),
  "Th17 cell" = c("RORC", "IL17A", "IL17F", "IL22", "CCR6", "STAT3", "RORA"),
  "Treg" = c("FOXP3", "CTLA4", "IL2RA", "IL2RB", "TNFRSF18", "IKZF2", "ENTPD1"),
  "Naive B cell" = c("CD19", "MS4A1", "CD22", "CD79A", "CD79B", "FCER2", "PAX5"),
  "Memory B cell" = c("CD19", "CD27", "CD38", "TNFRSF17", "SDC1", "SLAMF7"),
  "Plasma cell" = c("SDC1", "CD38", "TNFRSF17", "SLAMF7", "XBP1", "IRF4", "PRDM1"),
  "NK cell" = c("NCAM1", "KLRD1", "KLRF1", "NCR1", "NCR3", "KLRC1", "KLRK1"),
  "NK T cell" = c("NCAM1", "CD3D", "CD3E", "KLRB1", "IL2RB", "ZBTB16"),
  "Monocyte" = c("CD14", "FCGR3A", "CSF1R", "ITGAM", "CD33", "CD86", "HLA-DRA"),
  "Macrophage M1" = c("CD86", "CD80", "TNF", "IL6", "IL12A", "IL23A", "NOS2", "TLR2", "TLR4"),
  "Macrophage M2" = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "ARG1", "CCL22", "CCL18"),
  "Neutrophil" = c("FCGR3B", "CEACAM8", "CSF3R", "ITGAM", "MMP9", "MPO", "ELANE"),
  "Eosinophil" = c("CCR3", "IL5RA", "SIGLEC8", "RNASE2", "EPX", "PRG2"),
  "Mast cell" = c("KIT", "FCER1A", "TPSAB1", "CPA3", "HDC", "MS4A2"),
  "Myeloid dendritic cell" = c("CD1C", "CLEC10A", "FCER1A", "ITGAX", "CLEC4C", "XCR1"),
  "Plasmacytoid dendritic cell" = c("LILRA4", "CLEC4C", "NRP1", "TCF4", "IRF7", "IRF8"),
  "Immature dendritic cell" = c("CD1A", "FCER1A", "CD207", "CCR6", "TLR2"),
  "Gamma delta T cell" = c("TRGV9", "TRDV2", "KLRB1", "KLRK1", "GNLY"),
  "Follicular helper T cell" = c("BCL6", "CXCR5", "PDCD1", "ICOS", "IL21"),
  "Central memory CD4 T cell" = c("CD4", "CD27", "CD28", "CCR7", "SELL", "IL7R", "TCF7"),
  "Effector memory T cell" = c("CD3D", "CD3E", "PRF1", "GZMB", "IFNG", "TBX21", "EOMES")
)

# 过滤：只保留基因组中包含足够基因的细胞集
immune_sets_filtered <- immune_cell_sets[
  sapply(immune_cell_sets, function(gs) sum(gs %in% rownames(expr_norm))) >= 3
]
cat(sprintf("过滤后保留 %d/%d 个免疫细胞基因集\n",
            length(immune_sets_filtered), length(immune_cell_sets)))

# 运行ssGSEA (GSVA >= 2.0 API)
ssgsea_scores <- gsva(ssgseaParam(as.matrix(expr_norm), immune_sets_filtered))

ssgsea_df <- as.data.frame(t(ssgsea_scores))
cat(sprintf("ssGSEA完成: %d 细胞类型 x %d 样本\n",
            ncol(ssgsea_df), nrow(ssgsea_df)))

# ---- 4. 免疫细胞差异分析 ----
cat("\n===== 免疫细胞差异分析 =====\n")

# 比较DN vs Control
immune_diff <- data.frame(
  CellType = colnames(ssgsea_df),
  Mean_DN = sapply(ssgsea_df[group == "DN", ], mean, na.rm = TRUE),
  Mean_Control = sapply(ssgsea_df[group == "Control", ], mean, na.rm = TRUE),
  log2FC = sapply(ssgsea_df, function(x) {
    mean(x[group == "DN"], na.rm = TRUE) - mean(x[group == "Control"], na.rm = TRUE)
  }),
  P_value = sapply(ssgsea_df, function(x) {
    tryCatch(t.test(x[group == "DN"], x[group == "Control"])$p.value,
             error = function(e) 1)
  })
)
immune_diff$adj_P <- p.adjust(immune_diff$P_value, method = "BH")
immune_diff$Significant <- immune_diff$adj_P < 0.05

cat(sprintf("显著差异的免疫细胞: %d/%d\n",
            sum(immune_diff$Significant), nrow(immune_diff)))
if (sum(immune_diff$Significant) > 0) {
  print(subset(immune_diff, Significant)[, c("CellType", "log2FC", "adj_P")])
}

# ssGSEA热图
png("07_ssGSEA_heatmap.png", width = 1200, height = 900, res = 120)
annotation_col <- data.frame(Group = group, row.names = rownames(ssgsea_df))
pheatmap(t(ssgsea_df),
         annotation_col = annotation_col,
         annotation_colors = list(Group = c(DN = "tomato", Control = "steelblue")),
         show_colnames = FALSE,
         scale = "row",
         main = "ssGSEA Immune Cell Infiltration Scores",
         fontsize_row = 9,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
cat("✅ ssGSEA热图已保存\n")

# ---- 5. 免疫细胞浸润箱线图 ----
cat("\n===== 免疫细胞箱线图 =====\n")

# 选取显著差异的免疫细胞（或Top12）
if (sum(immune_diff$Significant) >= 1) {
  plot_cells <- head(immune_diff$CellType[immune_diff$Significant], 12)
} else {
  plot_cells <- head(immune_diff$CellType[order(immune_diff$P_value)], 12)
}

ssgsea_long <- ssgsea_df
ssgsea_long$Group <- group
ssgsea_long$Sample <- rownames(ssgsea_long)

plot_data <- ssgsea_long %>%
  pivot_longer(cols = all_of(plot_cells),
               names_to = "CellType", values_to = "Score")

ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  stat_compare_means(aes(group = Group), label = "p.signif",
                     method = "t.test", label.y.npc = 0.95, size = 3) +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Immune Cell Infiltration - DN vs Control",
       y = "ssGSEA Score") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("07_immune_boxplot.png", width = 14, height = 10, dpi = 150)
cat("✅ 免疫细胞箱线图已保存\n")

# ---- 6. 标志物-免疫细胞相关性 ----
cat("\n===== 标志物与免疫细胞相关性 =====\n")

gene_expr <- as.data.frame(t(expr_norm[final_genes, , drop = FALSE]))

cor_immune <- matrix(NA, nrow = length(final_genes), ncol = ncol(ssgsea_df),
                     dimnames = list(final_genes, colnames(ssgsea_df)))
cor_pval <- cor_immune

for (g in final_genes) {
  for (c in colnames(ssgsea_df)) {
    ct <- cor.test(gene_expr[, g], ssgsea_df[, c], method = "spearman")
    cor_immune[g, c] <- ct$estimate
    cor_pval[g, c] <- ct$p.value
  }
}

# 热图
png("07_gene_immune_correlation.png", width = 1200, height = 600, res = 120)
pheatmap(cor_immune,
         display_numbers = matrix(ifelse(cor_pval < 0.05, "*", ""),
                                  nrow = nrow(cor_pval)),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Biomarker - Immune Cell Spearman Correlation",
         fontsize_row = 10, fontsize_number = 10,
         number_color = "black")
dev.off()
cat("✅ 标志物-免疫相关性热图已保存\n")

# ---- 7. 免疫检查点基因分析 ----
cat("\n===== 免疫检查点分析 =====\n")

# 常见免疫检查点基因
checkpoint_genes <- c(
  "PDCD1", "CD274", "PDCD1LG2",   # PD-1/PD-L1 axis
  "CTLA4",                        # CTLA-4
  "LAG3",                         # LAG-3
  "HAVCR2",                       # TIM-3
  "TIGIT",                        # TIGIT
  "VSIR",                         # VISTA
  "BTLA",                         # BTLA
  "CD226",                        # DNAM-1
  "ICOS",                         # ICOS
  "CD80", "CD86",                 # B7 family
  "IDO1",                         # IDO
  "TNFRSF4", "TNFSF4",           # OX40/OX40L
  "TNFRSF9", "TNFSF9",           # 4-1BB/4-1BBL
  "TNFRSF18", "TNFSF18",         # GITR/GITRL
  "CD70", "CD27"                  # CD70/CD27
)

checkpoints_in_data <- intersect(checkpoint_genes, rownames(expr_norm))
cat(sprintf("免疫检查点基因: %d/%d 在数据中\n",
            length(checkpoints_in_data), length(checkpoint_genes)))

if (length(checkpoints_in_data) > 0) {
  # 差异分析
  cp_expr <- as.data.frame(t(expr_norm[checkpoints_in_data, , drop = FALSE]))

  cp_diff <- data.frame(
    Gene = checkpoints_in_data,
    Mean_DN = sapply(cp_expr[group == "DN", , drop = FALSE], mean),
    Mean_Control = sapply(cp_expr[group == "Control", , drop = FALSE], mean),
    logFC = sapply(cp_expr, function(x) {
      mean(x[group == "DN"]) - mean(x[group == "Control"])
    }),
    P_value = sapply(cp_expr, function(x) {
      tryCatch(t.test(x[group == "DN"], x[group == "Control"])$p.value,
               error = function(e) 1)
    })
  )
  cp_diff$adj_P <- p.adjust(cp_diff$P_value, method = "BH")
  cp_diff <- cp_diff[order(cp_diff$adj_P), ]

  # 热图
  cp_heatmap <- t(scale(t(expr_norm[cp_diff$Gene, ])))
  colnames(cp_heatmap) <- colnames(expr_norm)

  png("07_checkpoint_heatmap.png", width = 1000, height = 600, res = 120)
  pheatmap(cp_heatmap,
           annotation_col = data.frame(Group = group,
                                       row.names = colnames(expr_norm)),
           annotation_colors = list(Group = c(DN = "tomato", Control = "steelblue")),
           show_colnames = FALSE,
           cluster_cols = TRUE,
           main = "Immune Checkpoint Gene Expression",
           fontsize_row = 10,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
  cat("✅ 免疫检查点热图已保存\n")

  # 显著差异的检查点
  cp_sig <- subset(cp_diff, adj_P < 0.05)
  cat(sprintf("显著差异检查点: %d 个\n", nrow(cp_sig)))
  if (nrow(cp_sig) > 0) print(cp_sig[, c("Gene", "logFC", "adj_P")])
}

# ---- 8. ESTIMATE免疫评分 ----
cat("\n===== ESTIMATE评分 =====\n")

# 简化的Stromal/Immune评分 (基于ESTIMATE方法的核心基因)
estimate_genes <- list(
  Stromal = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "COL6A1",
              "COL6A2", "COL6A3", "LAMA2", "LAMA4", "LAMB1", "LAMC1",
              "FN1", "DCN", "LUM", "BGN", "FBN1", "SPARC", "TNC", "VIM"),
  Immune = c("CD2", "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD19", "CD79A",
             "CD79B", "CD14", "CD68", "CD163", "FCGR3A", "HLA-DRA", "HLA-DRB1",
             "HLA-DPA1", "HLA-DPB1", "TLR2", "TLR4", "CXCL10", "CCL2", "CCL5",
             "GZMA", "GZMB", "PRF1", "IFNG", "TNF", "IL6", "IL10")
)

# 计算富集分数 (GSVA >= 2.0 API)
estimate_scores <- gsva(ssgseaParam(as.matrix(expr_norm), estimate_genes))
estimate_df <- as.data.frame(t(estimate_scores))
estimate_df$Group <- group
estimate_df$ESTIMATE <- estimate_df$Stromal + estimate_df$Immune

# 差异检验
for (score_type in c("Stromal", "Immune", "ESTIMATE")) {
  p <- t.test(estimate_df[group == "DN", score_type],
              estimate_df[group == "Control", score_type])$p.value
  cat(sprintf("  %s Score: DN=%.3f, Control=%.3f, p=%.4f\n",
              score_type,
              mean(estimate_df[group == "DN", score_type]),
              mean(estimate_df[group == "Control", score_type]),
              p))
}

# ESTIMATE箱线图
estimate_long <- estimate_df %>%
  pivot_longer(cols = c("Stromal", "Immune", "ESTIMATE"),
               names_to = "Component", values_to = "Score")

ggplot(estimate_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.5, alpha = 0.5, size = 1) +
  facet_wrap(~ Component) +
  stat_compare_means(aes(group = Group), label = "p.format",
                     method = "t.test") +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "ESTIMATE Immune/Stromal Scores",
       y = "Score") +
  theme_minimal(base_size = 13)
ggsave("07_ESTIMATE_scores.png", width = 10, height = 5, dpi = 150)
cat("✅ ESTIMATE评分图已保存\n")

# ---- 9. 保存 ----
cat("\n===== 保存结果 =====\n")
saveRDS(ssgsea_df, "07_ssgsea_scores.rds")
saveRDS(immune_diff, "07_immune_diff.rds")
saveRDS(cor_immune, "07_gene_immune_cor.rds")
if (exists("cp_diff")) saveRDS(cp_diff, "07_checkpoint_diff.rds")
saveRDS(estimate_df, "07_estimate_scores.rds")

write.csv(immune_diff, "07_immune_diff.csv", row.names = FALSE)
write.csv(as.data.frame(cor_immune), "07_gene_immune_cor.csv")

cat("\n========================================\n")
cat("07_免疫浸润分析完成！输出:\n")
cat("  07_ssGSEA_heatmap.png — 免疫细胞浸润热图\n")
cat("  07_immune_boxplot.png — 免疫细胞箱线图\n")
cat("  07_gene_immune_correlation.png — 标志物-免疫相关性\n")
cat("  07_ESTIMATE_scores.png — ESTIMATE评分\n")
cat("========================================\n")
