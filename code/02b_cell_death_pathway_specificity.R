# ============================================================
# 双硫死亡+DN 竞赛项目 — 02b_细胞死亡通路特异性分析（审稿修订版）
# 目的：验证双硫死亡signature的特异性
# 方法：
#   1. 比较5种细胞死亡通路活性（双硫死亡 vs 铁死亡 vs 凋亡 vs 焦亡 vs 坏死性凋亡）
#   2. ssGSEA通路活性评分
#   3. 氧化应激相关基因深度分析（NADPH代谢、谷胱甘肽通路、ROS）
#   4. 通路活性热图 + 相关性分析
# ============================================================

library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(tidyr)
library(corrplot)

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")

expr_norm <- readRDS("GSE96804_gene_expr.rds")
pheno <- readRDS("GSE96804_pheno.rds")

group_info <- pheno[["title"]]
group <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"),
                levels = c("Control", "DN"))

cat(sprintf("表达矩阵: %d 基因 x %d 样本\n", nrow(expr_norm), ncol(expr_norm)))
cat(sprintf("分组: DN=%d, Control=%d\n", sum(group=="DN"), sum(group=="Control")))

# ---- 2. 定义5种细胞死亡通路基因集 ----
cat("\n===== 构建细胞死亡通路基因集 =====\n")

# 双硫死亡基因集（我们的）
disulfidptosis_genes <- readLines("disulfidptosis_genes.txt")

# 铁死亡基因集（FerrDB核心 + 文献共识）
ferroptosis_genes <- c(
  # 核心调控因子
  "GPX4", "SLC7A11", "ACSL4", "LPCAT3", "FTH1", "FTL",
  "TFRC", "SLC39A8", "SLC39A14", "SLC40A1",
  # 代谢相关
  "GCLC", "GCLM", "GSS", "GSR", "TXN", "TXNRD1",
  # 脂质过氧化
  "ALOX5", "ALOX12", "ALOX15", "ALOXE3", "PTGS2",
  # 调控因子
  "NFE2L2", "KEAP1", "HMOX1", "NQO1", "SQSTM1",
  # 铁代谢
  "IREB2", "FBXL5", "NCOA4", "PCBP1", "PCBP2",
  # 补充
  "VDAC2", "VDAC3", "MAP1LC3A", "ATG5", "ATG7",
  "BECN1", "CISD1", "CISD2", "HSPB1", "NFKB1"
)

# 凋亡基因集（KEGG + GO经典）
apoptosis_genes <- c(
  # 外源通路
  "FAS", "FASLG", "TNF", "TNFRSF1A", "TNFRSF10A", "TNFRSF10B",
  "FADD", "TRADD", "CASP8", "CASP10", "CFLAR",
  # 内源通路（线粒体）
  "BCL2", "BCL2L1", "BAX", "BAK1", "BAD", "BID",
  "BIM", "PUMA", "NOXA", "MCL1",
  # Cytochrome c → apoptosome
  "CYCS", "APAF1", "CASP9",
  # 执行caspase
  "CASP3", "CASP6", "CASP7",
  # IAP家族
  "XIAP", "BIRC2", "BIRC3", "BIRC5",
  # p53通路
  "TP53", "MDM2", "CDKN1A",
  # 其他
  "AIFM1", "ENDOG", "DIABLO", "HTRA2"
)

# 焦亡基因集
pyroptosis_genes <- c(
  # 炎症小体
  "NLRP1", "NLRP3", "NLRC4", "NAIP", "AIM2",
  "PYCARD", "CASP1", "CASP4", "CASP5",
  # Gasdermin家族
  "GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME",
  # 炎症因子
  "IL1B", "IL18", "HMGB1",
  # 调控
  "TLR4", "NFKB1", "RELA", "TNF",
  # 其他
  "NLRP6", "NLRP12", "IFI16", "MEFV", "PSTPIP1"
)

# 坏死性凋亡基因集
necroptosis_genes <- c(
  # 核心通路
  "RIPK1", "RIPK3", "MLKL",
  # 上游激活
  "TNF", "TNFR1", "TRADD", "TRAF2", "TRAF5",
  "CYLD", "CIAP1", "CIAP2",
  # 调控因子
  "FADD", "CASP8", "CFLAR",
  # RIPK1/RIPK3相关
  "TAB1", "TAB2", "TAB3", "TAK1",
  "IKKA", "IKKB", "NEMO",
  # 执行
  "PGAM5", "DRP1",
  # 其他
  "HMGB1", "PARP1", "AIFM1",
  "ZBP1", "TLR3", "TLR4", "IFNAR1", "IFNAR2"
)

# 汇总
cell_death_sets <- list(
  Disulfidptosis = intersect(disulfidptosis_genes, rownames(expr_norm)),
  Ferroptosis = intersect(ferroptosis_genes, rownames(expr_norm)),
  Apoptosis = intersect(apoptosis_genes, rownames(expr_norm)),
  Pyroptosis = intersect(pyroptosis_genes, rownames(expr_norm)),
  Necroptosis = intersect(necroptosis_genes, rownames(expr_norm))
)

cat("\n各通路在表达矩阵中的基因数:\n")
for (nm in names(cell_death_sets)) {
  cat(sprintf("  %s: %d genes\n", nm, length(cell_death_sets[[nm]])))
}

# ---- 氧化应激相关基因集 ----
cat("\n===== 氧化应激相关基因集 =====\n")

oxidative_stress_sets <- list(
  # NADPH代谢
  NADPH_metabolism = intersect(c(
    "G6PD", "PGD", "H6PD", "IDH1", "IDH2", "ME1", "ME2", "ME3",
    "NNT", "GLUD1", "GLUD2"
  ), rownames(expr_norm)),
  # 谷胱甘肽通路
  Glutathione = intersect(c(
    "GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX3", "GPX4",
    "GSTP1", "GSTA1", "GSTM1", "GGT1", "ANPEP", "OPLAH",
    "SLC7A11", "SLC3A2"
  ), rownames(expr_norm)),
  # ROS相关
  ROS = intersect(c(
    "SOD1", "SOD2", "SOD3", "CAT", "PRDX1", "PRDX2", "PRDX3",
    "PRDX4", "PRDX5", "PRDX6", "TXN", "TXNRD1", "TXNRD2",
    "NOX1", "NOX4", "NOX5", "CYBB", "DUOX1", "DUOX2",
    "NFE2L2", "KEAP1", "HMOX1", "NQO1"
  ), rownames(expr_norm)),
  # Pentose Phosphate Pathway
  PPP = intersect(c(
    "G6PD", "PGD", "PGLS", "TKT", "TALDO1", "RPIA", "RPE",
    "PRPS1", "PRPS2"
  ), rownames(expr_norm))
)

cat("\n氧化应激子通路基因数:\n")
for (nm in names(oxidative_stress_sets)) {
  cat(sprintf("  %s: %d genes\n", nm, length(oxidative_stress_sets[[nm]])))
}

# ---- 3. ssGSEA通路活性评分 ----
cat("\n===== ssGSEA通路活性评分 =====\n")

# 合并所有基因集
all_gene_sets <- c(cell_death_sets, oxidative_stress_sets)

# 运行ssGSEA (兼容新旧GSVA版本)
set.seed(42)
# 尝试使用gsvaParam (GSVA >= 1.50)
ssgsea_scores <- tryCatch({
  gsva_param <- gsvaParam(as.matrix(expr_norm), all_gene_sets,
                           kcdf = "Gaussian",
                           minSize = 5,
                           maxSize = 500)
  gsva(gsva_param, verbose = FALSE)
}, error = function(e) {
  # 回退到旧版API
  gsva(as.matrix(expr_norm), all_gene_sets,
       method = "ssgsea",
       kcdf = "Gaussian",
       min.sz = 5,
       max.sz = 500,
       ssgsea.norm = TRUE,
       verbose = FALSE)
})

cat(sprintf("ssGSEA评分矩阵: %d pathways x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

# ---- 4. 通路活性差异比较 ----
cat("\n===== 通路活性差异 (DN vs Control) =====\n")

pathway_diff <- data.frame(
  Pathway = rownames(ssgsea_scores),
  DN_mean = rowMeans(ssgsea_scores[, group == "DN"]),
  Control_mean = rowMeans(ssgsea_scores[, group == "Control"]),
  log2FC = rep(NA, nrow(ssgsea_scores)),
  P_value = rep(NA, nrow(ssgsea_scores)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(ssgsea_scores)) {
  dn_vals <- ssgsea_scores[i, group == "DN"]
  ctrl_vals <- ssgsea_scores[i, group == "Control"]
  test_result <- t.test(dn_vals, ctrl_vals)
  pathway_diff$P_value[i] <- test_result$p.value
  pathway_diff$log2FC[i] <- log2(mean(dn_vals) / mean(ctrl_vals))
}

pathway_diff$adj_P <- p.adjust(pathway_diff$P_value, method = "BH")
pathway_diff$Significance <- ifelse(pathway_diff$adj_P < 0.001, "***",
                             ifelse(pathway_diff$adj_P < 0.01, "**",
                             ifelse(pathway_diff$adj_P < 0.05, "*", "ns")))
pathway_diff <- pathway_diff[order(pathway_diff$adj_P), ]

cat("\n通路活性差异汇总:\n")
print(pathway_diff[, c("Pathway", "log2FC", "P_value", "adj_P", "Significance")])

# ---- 5. 细胞死亡通路热图 ----
cat("\n===== 细胞死亡通路比较热图 =====\n")

# 通路活性热图（仅细胞死亡通路）
death_pathways <- names(cell_death_sets)
death_scores <- ssgsea_scores[death_pathways, ]

# 按分组排序
sample_order <- order(group)
death_scores_ordered <- death_scores[, sample_order]

# 注释条
annotation_col <- data.frame(
  Group = group[sample_order],
  row.names = colnames(death_scores_ordered)
)
ann_colors <- list(Group = c(DN = "tomato", Control = "steelblue"))

png("02b_cell_death_pathway_heatmap.png", width = 1000, height = 600, res = 130)
pheatmap(death_scores_ordered,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         main = "Cell Death Pathway Activity (ssGSEA) — GSE96804",
         fontsize = 11,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
cat("✅ 细胞死亡通路热图已保存\n")

# ---- 6. 通路活性箱线图 ----
cat("\n===== 通路活性箱线图 =====\n")

# 长格式数据
plot_data <- data.frame()
for (i in 1:nrow(ssgsea_scores)) {
  pathway_name <- rownames(ssgsea_scores)[i]
  pathway_data <- data.frame(
    Pathway = pathway_name,
    Score = as.numeric(ssgsea_scores[i, ]),
    Group = group,
    Category = ifelse(pathway_name %in% death_pathways, "Cell Death", "Oxidative Stress")
  )
  plot_data <- rbind(plot_data, pathway_data)
}

# 细胞死亡通路比较
p1 <- ggplot(subset(plot_data, Category == "Cell Death"),
             aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 5) +
  stat_compare_means(aes(group = Group), label = "p.signif",
                     method = "t.test", label.y.npc = 0.9, size = 5) +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Cell Death Pathway Activity Comparison (ssGSEA)",
       y = "ssGSEA Enrichment Score") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
ggsave("02b_cell_death_boxplot.png", width = 14, height = 5, dpi = 150)
cat("✅ 细胞死亡通路箱线图已保存\n")

# 氧化应激通路
p2 <- ggplot(subset(plot_data, Category == "Oxidative Stress"),
             aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
  stat_compare_means(aes(group = Group), label = "p.signif",
                     method = "t.test", label.y.npc = 0.9, size = 5) +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Oxidative Stress Pathway Activity Comparison (ssGSEA)",
       y = "ssGSEA Enrichment Score") +
  theme_minimal(base_size = 12) +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
ggsave("02b_oxidative_stress_boxplot.png", width = 12, height = 5, dpi = 150)
cat("✅ 氧化应激通路箱线图已保存\n")

# ---- 7. 双硫死亡特异性vs其他通路热图 ----
cat("\n===== 通路特异性分析 =====\n")

# 计算各通路效应大小
pathway_effect <- pathway_diff[pathway_diff$Pathway %in% death_pathways, ]
pathway_effect <- pathway_effect[order(abs(pathway_effect$log2FC), decreasing = TRUE), ]

# 效应对比图
ggplot(pathway_effect, aes(x = reorder(Pathway, abs(log2FC)),
                           y = abs(log2FC), fill = log2FC > 0)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("P=%s", Significance)), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
                    labels = c("TRUE" = "Up in DN", "FALSE" = "Down in DN")) +
  coord_flip() +
  labs(title = "Cell Death Pathway Effect Size — DN vs Control",
       subtitle = "Absolute effect size (|log2FC|) by ssGSEA",
       x = "", y = "|log2(Fold Change)|", fill = "Direction") +
  theme_minimal(base_size = 13)
ggsave("02b_pathway_effect_size.png", width = 8, height = 4, dpi = 150)
cat("✅ 通路效应图已保存\n")

# ---- 8. 氧化应激基因深度分析 ----
cat("\n===== 氧化应激关键基因分析 =====\n")

# 关键氧化应激基因
key_ox_genes <- unique(unlist(oxidative_stress_sets))
key_ox_in_data <- intersect(key_ox_genes, rownames(expr_norm))
cat(sprintf("氧化应激关键基因数: %d\n", length(key_ox_in_data)))

# 差异表达分析
ox_expr <- expr_norm[key_ox_in_data, ]

ox_diff <- data.frame(
  Gene = key_ox_in_data,
  row.names = key_ox_in_data
)

# 统计检验
for (gene in key_ox_in_data) {
  dn_vals <- as.numeric(ox_expr[gene, group == "DN"])
  ctrl_vals <- as.numeric(ox_expr[gene, group == "Control"])
  test <- t.test(dn_vals, ctrl_vals)
  ox_diff[gene, "logFC"] <- mean(dn_vals) - mean(ctrl_vals)
  ox_diff[gene, "P_value"] <- test$p.value
}
ox_diff$adj_P <- p.adjust(ox_diff$P_value, method = "BH")
ox_diff$Significant <- ox_diff$adj_P < 0.05
ox_diff <- ox_diff[order(ox_diff$adj_P), ]

cat(sprintf("显著差异氧化应激基因: %d\n", sum(ox_diff$Significant)))
if (sum(ox_diff$Significant) > 0) {
  cat("显著基因:\n")
  sig_ox <- subset(ox_diff, Significant)
  print(sig_ox[, c("logFC", "P_value", "adj_P")])
}

# 氧化应激基因热图
if (sum(ox_diff$Significant) >= 3) {
  sig_genes <- rownames(subset(ox_diff, Significant))
  sig_expr <- ox_expr[sig_genes, sample_order]

  png("02b_oxidative_stress_heatmap.png", width = 900, height = 700, res = 130)
  pheatmap(sig_expr,
           scale = "row",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           main = "Differentially Expressed Oxidative Stress Genes",
           fontsize = 10,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
  cat("✅ 氧化应激基因热图已保存\n")
}

# ---- 9. 双硫死亡基因与氧化应激基因相关性 ----
cat("\n===== Disulfidptosis vs Oxidative Stress 基因相关性 =====\n")

final_biomarkers <- readRDS("04_final_biomarkers.rds")
dr_genes_in_data <- intersect(final_biomarkers, rownames(expr_norm))

# 取显著氧化应激基因
top_ox_genes <- head(rownames(ox_diff), 20)
top_ox_genes <- intersect(top_ox_genes, rownames(expr_norm))

# 合并基因
cor_genes <- unique(c(dr_genes_in_data, top_ox_genes))

if (length(cor_genes) >= 5) {
  cor_matrix <- cor(t(expr_norm[cor_genes, ]), method = "spearman")

  # 只显示双硫死亡基因 vs 氧化应激基因的交叉部分
  ox_genes_only <- intersect(top_ox_genes, cor_genes)
  dr_genes_only <- intersect(dr_genes_in_data, cor_genes)

  cross_cor <- cor_matrix[dr_genes_only, ox_genes_only, drop = FALSE]

  png("02b_dr_vs_oxidative_correlation.png", width = 1000, height = 700, res = 130)
  pheatmap(cross_cor,
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = "Disulfidptosis Biomarkers vs Oxidative Stress Genes (Spearman)",
           fontsize = 10,
           fontsize_number = 8)
  dev.off()
  cat("✅ 双硫死亡vs氧化应激相关性热图已保存\n")
}

# ---- 10. 通路间相关性 ----
cat("\n===== 通路间相关性 =====\n")

pathway_cor <- cor(t(ssgsea_scores), method = "spearman")

png("02b_pathway_correlation.png", width = 900, height = 800, res = 130)
corrplot(pathway_cor, method = "color", type = "upper",
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.cex = 0.9,
         title = "Pathway Activity Correlation (Spearman)",
         mar = c(0, 0, 2, 0))
dev.off()
cat("✅ 通路相关性图已保存\n")

# ---- 11. 双硫死亡特异性指数 ----
cat("\n===== 双硫死亡特异性指数 =====\n")

# 计算每个样本的双硫死亡特异性评分（标准化后与均值的偏差）
death_activity <- t(ssgsea_scores[death_pathways, ])
death_z <- scale(death_activity)

disulfidptosis_specificity <- death_z[, "Disulfidptosis"] -
  rowMeans(death_z[, setdiff(death_pathways, "Disulfidptosis"), drop = FALSE])

specificity_df <- data.frame(
  Sample = colnames(ssgsea_scores),
  Group = group,
  Specificity = disulfidptosis_specificity
)

# 比较特异性
spec_test <- t.test(Specificity ~ Group, data = specificity_df)
cat(sprintf("双硫死亡特异性指数 (DN vs Control):\n"))
cat(sprintf("  DN mean = %.3f, Control mean = %.3f\n",
            mean(specificity_df$Specificity[specificity_df$Group == "DN"]),
            mean(specificity_df$Specificity[specificity_df$Group == "Control"])))
cat(sprintf("  P = %.4f\n", spec_test$p.value))

ggplot(specificity_df, aes(x = Group, y = Specificity, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.format") +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Disulfidptosis Pathway Specificity Index",
       subtitle = "Deviation of disulfidptosis activity from other cell death pathways",
       y = "Specificity Index (z-score)") +
  theme_minimal(base_size = 13)
ggsave("02b_disulfidptosis_specificity.png", width = 6, height = 5, dpi = 150)
cat("✅ 双硫死亡特异性图已保存\n")

# ---- 12. 汇总表格 ----
cat("\n===== 汇总结果 =====\n")

# 各通路显著差异基因统计
summary_table <- data.frame(
  Pathway = names(cell_death_sets),
  N_Genes = sapply(cell_death_sets, length)
)
summary_table$N_Sig_Genes <- sapply(names(cell_death_sets), function(pw) {
  genes <- cell_death_sets[[pw]]
  deg_results <- readRDS("02_deg_results.rds")
  sum(deg_results$Gene %in% genes & deg_results$adj.P.Val < 0.05, na.rm = TRUE)
})
summary_table$Pct_Sig <- round(summary_table$N_Sig_Genes / summary_table$N_Genes * 100, 1)

# 加入通路活性
pathway_diff_death <- pathway_diff[pathway_diff$Pathway %in% names(cell_death_sets), ]
rownames(pathway_diff_death) <- pathway_diff_death$Pathway
summary_table$ssGSEA_log2FC <- pathway_diff_death[summary_table$Pathway, "log2FC"]
summary_table$ssGSEA_P <- pathway_diff_death[summary_table$Pathway, "adj_P"]

cat("\n各细胞死亡通路比较:\n")
print(summary_table)

# ---- 13. 保存结果 ----
cat("\n===== 保存结果 =====\n")

results <- list(
  cell_death_sets = cell_death_sets,
  oxidative_stress_sets = oxidative_stress_sets,
  ssgsea_scores = ssgsea_scores,
  pathway_diff = pathway_diff,
  ox_diff = ox_diff,
  specificity_result = spec_test,
  summary_table = summary_table
)

saveRDS(results, "02b_cell_death_analysis.rds")
write.csv(pathway_diff, "02b_pathway_diff.csv", row.names = FALSE)
write.csv(ox_diff, "02b_oxidative_stress_diff.csv")
write.csv(summary_table, "02b_cell_death_summary.csv", row.names = FALSE)

cat("\n========================================\n")
cat("02b_细胞死亡通路分析完成！关键发现:\n")
cat(sprintf("  - 双硫死亡通路效应: log2FC=%.3f, P=%.2e\n",
            pathway_diff$log2FC[pathway_diff$Pathway=="Disulfidptosis"],
            pathway_diff$adj_P[pathway_diff$Pathway=="Disulfidptosis"]))
cat(sprintf("  - 双硫死亡特异性指数P值: %.4f\n", spec_test$p.value))
if (nrow(sig_ox) > 0) {
  cat(sprintf("  - 显著差异氧化应激基因: %d个 (Top: %s)\n",
              nrow(sig_ox), paste(head(rownames(sig_ox), 5), collapse=", ")))
}
cat("========================================\n")
