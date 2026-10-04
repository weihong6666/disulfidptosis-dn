# ============================================================
# 双硫死亡+DN 竞赛项目 — 02_差异分析+双硫死亡基因交集
# 输入：GSE96804标准化表达矩阵 + 双硫死亡基因集
# 输出：差异双硫死亡基因 (DR-DEGs)
# ============================================================

library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")
if (!exists("expr_norm")) {
  expr_norm <- readRDS("GSE96804_expr_normalized.rds")
}
pheno <- readRDS("GSE96804_pheno.rds")

# ---- 1.5 探针ID转基因符号 ----
cat("\n===== 探针→基因符号转换 =====\n")

# 检查是否已有基因级表达矩阵
if (file.exists("GSE96804_gene_expr.rds")) {
  cat("发现基因级表达矩阵缓存，直接加载...\n")
  expr_norm <- readRDS("GSE96804_gene_expr.rds")
} else {
  cat("读取GPL17586注释文件...\n")

  # 读取注释文件（跳过#注释行，tab分隔）
  annot <- read.delim("GPL17586-45144.txt", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "#",
                       na.strings = c("", "---"))

  cat(sprintf("注释文件: %d 行, 列名: %s\n", nrow(annot),
              paste(colnames(annot)[1:5], collapse = ", ")))

  # 从gene_assignment列提取基因符号
  # 格式: "NM_xxx // GENE_SYMBOL // description /// ENSTxxx // GENE_SYMBOL // ..."
  extract_gene_symbols <- function(ga) {
    if (is.na(ga) || ga == "---" || ga == "") return(NA_character_)
    # 分割多个assignment（用 /// 分隔）
    parts <- strsplit(ga, " /// ")[[1]]
    symbols <- character(0)
    for (part in parts) {
      fields <- strsplit(part, " // ")[[1]]
      if (length(fields) >= 2) {
        symbols <- c(symbols, fields[2])
      }
    }
    symbols <- unique(symbols[symbols != "---" & !is.na(symbols)])
    if (length(symbols) == 0) return(NA_character_)
    return(paste(symbols, collapse = ";"))
  }

  cat("提取基因符号...\n")
  annot$Gene_Symbol <- sapply(annot$gene_assignment, extract_gene_symbols)

  # 统计
  n_mapped <- sum(!is.na(annot$Gene_Symbol))
  cat(sprintf("成功映射: %d/%d (%.1f%%) 探针\n",
              n_mapped, nrow(annot), n_mapped / nrow(annot) * 100))

  # 创建探针→基因映射表
  probe_to_gene <- annot[!is.na(annot$Gene_Symbol), c("ID", "Gene_Symbol")]
  cat(sprintf("唯一基因: %d\n", length(unique(probe_to_gene$Gene_Symbol))))

  # 过滤表达矩阵中的探针
  common_probes <- intersect(rownames(expr_norm), probe_to_gene$ID)
  cat(sprintf("表达矩阵中可映射探针: %d/%d\n",
              length(common_probes), nrow(expr_norm)))

  # 多探针对应同一基因：取均值
  expr_sub <- expr_norm[common_probes, , drop = FALSE]
  probe_genes <- probe_to_gene$Gene_Symbol[match(common_probes, probe_to_gene$ID)]

  # 按基因聚合（取均值）
  gene_list <- unique(probe_genes)
  gene_expr <- matrix(NA, nrow = length(gene_list), ncol = ncol(expr_sub))
  rownames(gene_expr) <- gene_list
  colnames(gene_expr) <- colnames(expr_sub)

  for (g in gene_list) {
    probes <- common_probes[probe_genes == g]
    if (length(probes) == 1) {
      gene_expr[g, ] <- expr_sub[probes, ]
    } else {
      gene_expr[g, ] <- colMeans(expr_sub[probes, , drop = FALSE], na.rm = TRUE)
    }
  }

  cat(sprintf("基因级表达矩阵: %d 基因 x %d 样本\n", nrow(gene_expr), ncol(gene_expr)))

  # 替换表达矩阵
  expr_norm <- gene_expr
  saveRDS(expr_norm, "GSE96804_gene_expr.rds")
  cat("✅ 基因级表达矩阵已保存到 GSE96804_gene_expr.rds\n")
}

# ---- 2. 构建分组信息 ----
cat("\n===== 构建分组 =====\n")
# GSE96804: 41 DN肾小球 vs 20 正常肾小球
# 从title列提取分组（"Human glomeruli DN 01" vs "Human glomeruli Control 01"）
sample_names <- colnames(expr_norm)

# 使用title列判断DN/Control
group_col <- "title"
cat("使用分组列:", group_col, "\n")

# 创建分组向量（DN=1, Control=0）
group_info <- pheno[[group_col]]
group <- factor(ifelse(grepl("DN|diabetic|diabetes|patient|DKD", group_info, ignore.case = TRUE),
                       "DN", "Control"),
                levels = c("Control", "DN"))

cat(sprintf("DN样本: %d, Control样本: %d\n", sum(group == "DN"), sum(group == "Control")))

# ---- 3. 差异表达分析 (limma) ----
cat("\n===== limma差异分析 =====\n")

# 构建设计矩阵
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Control", "DN")

# 构建对比矩阵
contrast <- makeContrasts(DN - Control, levels = design)

# limma流程
fit <- lmFit(expr_norm, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)

# 提取所有差异基因结果
deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)
deg_results$Gene <- rownames(deg_results)

cat(sprintf("总计基因数: %d\n", nrow(deg_results)))
cat(sprintf("显著差异基因 (adj.P.Val<0.05 & |logFC|>0.5): %d\n",
            sum(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 0.5, na.rm = TRUE)))

# ---- 4. 双硫死亡基因集 ----
cat("\n===== 加载双硫死亡基因集 =====\n")

# 从文件读取我们收集的基因集
drg_file <- "双硫死亡基因集汇总.md"
if (!file.exists("disulfidptosis_genes.txt")) {
  # 手动定义组合基因集（来自文献）
  drg_list <- c(
    # 10基因核心集 (Liu et al. 2023)
    "SLC7A11", "SLC3A2", "GYS1", "LRPPRC", "NDUFA11",
    "NDUFS1", "NUBPL", "OXSM", "RPN1", "NCKAP1",
    # 15基因actin集
    "FLNA", "FLNB", "MYH9", "MYH10", "TLN1", "ACTB",
    "MYL6", "CAPZB", "DSTN", "IQGAP1", "ACTN4", "PDLIM1",
    "CD2AP", "INF2",
    # WAVE复合物 + 扩展
    "WASF2", "CYFIP1", "ABI2", "BRK1", "RAC1",
    "SLC2A1", "HK1", "G6PD", "PGD", "PRDX1", "BAK1",
    "NCKAP1L",
    # 从GeneCards补充（Xu 2023方法）
    "TLR4", "NFKB1", "TNF", "IL6", "IL1B",
    "MAPK1", "MAPK3", "AKT1", "TP53",
    "CXCL6", "CD48", "C1QB", "COL6A3",  # Xu的4个基因
    "VEGFA", "MAGI2", "THSD7A", "ANKRD28"  # Wang的4个基因
  )
  drg_list <- unique(drg_list)

  writeLines(drg_list, "disulfidptosis_genes.txt")
  cat("双硫死亡基因集已保存到 disulfidptosis_genes.txt\n")
} else {
  drg_list <- readLines("disulfidptosis_genes.txt")
}

cat(sprintf("双硫死亡基因集共 %d 个基因\n", length(drg_list)))

# 检查哪些基因在表达矩阵中
drg_in_expr <- intersect(drg_list, rownames(deg_results))
cat(sprintf("其中 %d 个在表达矩阵中\n", length(drg_in_expr)))
cat("未在矩阵中的基因:", paste(setdiff(drg_list, rownames(deg_results)), collapse=", "), "\n")

# ---- 5. 筛选差异双硫死亡基因 (DR-DEGs) ----
cat("\n===== 筛选DR-DEGs =====\n")

# 从差异分析结果中提取双硫死亡基因
dr_degs <- deg_results[deg_results$Gene %in% drg_in_expr, ]
dr_degs <- dr_degs[order(dr_degs$adj.P.Val), ]  # 按显著性排序

# 显著DR-DEGs
dr_degs_sig <- subset(dr_degs, adj.P.Val < 0.05 & abs(logFC) > 0.5)
cat(sprintf("显著DR-DEGs数量: %d\n", nrow(dr_degs_sig)))

if (nrow(dr_degs_sig) > 0) {
  cat("\n显著DR-DEGs列表:\n")
  print(dr_degs_sig[, c("Gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")])
}

# ---- 6. 火山图 ----
cat("\n===== 绘制火山图 =====\n")

# 标记：显著DR-DEGs + top标记基因
volcano_data <- deg_results
volcano_data$category <- "NS"
volcano_data$category[volcano_data$adj.P.Val < 0.05 & volcano_data$logFC > 0.5] <- "Up"
volcano_data$category[volcano_data$adj.P.Val < 0.05 & volcano_data$logFC < -0.5] <- "Down"
volcano_data$category[volcano_data$Gene %in% dr_degs_sig$Gene] <- "DR-DEG"

volcano_data$label <- ""
if (nrow(dr_degs_sig) >= 1) {
  top_genes <- head(dr_degs_sig$Gene, 15)
  volcano_data$label[volcano_data$Gene %in% top_genes] <- volcano_data$Gene[volcano_data$Gene %in% top_genes]
}

colors <- c("NS" = "grey70", "Up" = "#e41a1c", "Down" = "#377eb8", "DR-DEG" = "#ff7f00")

ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value), color = category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  color = "black", fontface = "italic") +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "GSE96804: DN vs Control",
       subtitle = sprintf("Significant genes: %d (up) + %d (down) | DR-DEGs: %d",
                          sum(volcano_data$category == "Up"),
                          sum(volcano_data$category == "Down"),
                          nrow(dr_degs_sig)),
       x = "log2(Fold Change)", y = "-log10(P-value)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
ggsave("02_volcano.png", width = 10, height = 8, dpi = 150)
cat("✅ 火山图已保存到 02_volcano.png\n")

# ---- 7. DR-DEGs热图 ----
cat("\n===== 绘制DR-DEGs热图 =====\n")

if (nrow(dr_degs_sig) >= 3) {
  dr_expr <- expr_norm[dr_degs_sig$Gene, ]

  # 按分组排序样本
  sample_order <- order(group)
  dr_expr <- dr_expr[, sample_order]

  # 注释条
  annotation_col <- data.frame(
    Group = group[sample_order],
    row.names = colnames(dr_expr)
  )

  ann_colors <- list(Group = c(DN = "tomato", Control = "steelblue"))

  # 用 pheatmap 自带的 filename 参数管理图形设备：
  # 先前用 png() + dev.off() 的写法在本机曾产出空白图（单色 PNG），改用 filename 更稳
  pheatmap(dr_expr,
           scale = "row",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_colnames = FALSE,
           cluster_cols = TRUE,
           main = "DR-DEGs Expression Heatmap (GSE96804)",
           fontsize_row = 11,
           border_color = NA,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           filename = "02_DRDEGs_heatmap.png",
           width = 8.33, height = 6.67, res = 120)
  cat("✅ 热图已保存到 02_DRDEGs_heatmap.png\n")
}

# ---- 8. 保存结果 ----
cat("\n===== 保存结果 =====\n")
saveRDS(deg_results, "02_deg_results.rds")
saveRDS(dr_degs, "02_dr_degs.rds")
saveRDS(dr_degs_sig, "02_dr_degs_sig.rds")
write.csv(dr_degs, "02_dr_degs.csv")
write.csv(dr_degs_sig, "02_dr_degs_sig.csv")

cat("\n========================================\n")
cat("02_差异分析完成！输出文件:\n")
cat("  02_deg_results.rds — 全部差异表达结果\n")
cat("  02_dr_degs.rds / .csv — 双硫死亡差异基因\n")
cat("  02_dr_degs_sig.rds / .csv — 显著DR-DEGs\n")
cat("  02_volcano.png — 火山图\n")
cat("  02_DRDEGs_heatmap.png — DR-DEGs热图\n")
cat("========================================\n")
