# ============================================================
# 双硫死亡+DN 竞赛项目 — 01_GEO数据下载与质控
# 数据集：GSE96804 (训练集)
# 支持：GEOquery在线下载 / 本地.gz文件手动解析
# ============================================================

# ---- 加载包 ----
library(limma)
library(ggplot2)
library(ggpubr)
library(pheatmap)

# ---- 1. 加载/下载 GSE96804 ----
cat("\n===== 加载 GSE96804 =====\n")
gse_id <- "GSE96804"
local_gz <- "GSE96804_series_matrix.txt.gz"

# 标准化表达矩阵缓存
if (file.exists("GSE96804_expr_normalized.rds") && file.exists("GSE96804_pheno.rds")) {
  cat("发现处理后的数据缓存，直接加载...\n")
  expr_norm <- readRDS("GSE96804_expr_normalized.rds")
  expr_data <- expr_norm  # 已经是标准化后的
  pheno_data <- readRDS("GSE96804_pheno.rds")

} else if (file.exists(local_gz)) {
  # ---- 手动解析本地 .gz 文件 ----
  cat("解析本地GEO文件...\n")

  con <- gzfile(local_gz, "r")
  all_lines <- readLines(con, warn = FALSE)
  close(con)

  # 定位数据表
  table_begin <- grep("^!series_matrix_table_begin", all_lines)
  table_end   <- grep("^!series_matrix_table_end", all_lines)
  stopifnot(length(table_begin) > 0 && length(table_end) > 0)

  # 解析列名（第一列是 ID_REF，后面是样本GSM编号）
  header_line <- all_lines[table_begin + 1]
  columns <- strsplit(header_line, "\t")[[1]]
  columns <- gsub('"', '', columns)
  sample_ids <- columns[-1]
  cat(sprintf("  样本数: %d\n", length(sample_ids)))

  # 解析表达数据
  data_lines <- all_lines[(table_begin + 2):(table_end - 1)]
  n_genes <- length(data_lines)

  # 预分配矩阵
  expr_data <- matrix(NA_real_, nrow = n_genes, ncol = length(sample_ids))
  gene_ids  <- character(n_genes)

  for (i in seq_along(data_lines)) {
    fields <- strsplit(data_lines[i], "\t")[[1]]
    gene_ids[i] <- gsub('"', '', fields[1])
    expr_data[i, ] <- as.numeric(fields[-1])
  }
  rownames(expr_data) <- gene_ids
  colnames(expr_data) <- sample_ids
  cat(sprintf("  表达矩阵: %d 基因 x %d 样本\n", n_genes, length(sample_ids)))

  # 解析样本metadata（从 !Sample_* 行提取）
  meta_lines <- all_lines[1:(table_begin - 1)]
  sample_meta <- meta_lines[grepl("^!Sample_", meta_lines)]

  pheno_data <- data.frame(row.names = sample_ids, stringsAsFactors = FALSE)

  for (ln in sample_meta) {
    # 格式: !Sample_characteristics_ch1  "tissue: glomeruli"  "tissue: glomeruli" ...
    key_val <- strsplit(ln, "\t")[[1]]
    key_raw <- sub("^!Sample_", "", key_val[1])
    vals <- gsub('"', '', key_val[-1])

    if (length(vals) >= length(sample_ids)) {
      pheno_data[[key_raw]] <- vals[1:length(sample_ids)]
    }
  }

  # 提取标题作为补充
  title_lines <- meta_lines[grepl("^!Sample_title", meta_lines)]
  if (length(title_lines) > 0) {
    tv <- strsplit(title_lines[1], "\t")[[1]][-1]
    if (length(tv) >= length(sample_ids)) {
      pheno_data$title <- gsub('"', '', tv)[1:length(sample_ids)]
    }
  }

  cat(sprintf("  表型列: %s\n", paste(colnames(pheno_data), collapse = ", ")))

  # 保存原始数据
  saveRDS(expr_data, "GSE96804_expr_raw.rds")
  saveRDS(pheno_data, "GSE96804_pheno.rds")
  cat("✅ 原始数据已保存 (手动解析)\n")

} else {
  # ---- 在线下载 ----
  cat("从GEO在线下载...\n")
  cat("⚠ 如超时，请浏览器下载:\n")
  cat("  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96804/matrix/GSE96804_series_matrix.txt.gz\n\n")

  library(GEOquery)
  old_timeout <- getOption("timeout")
  options(timeout = 600)
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
  options(timeout = old_timeout)

  expr_data <- exprs(gse[[1]])
  pheno_data <- pData(gse[[1]])
  saveRDS(expr_data, "GSE96804_expr_raw.rds")
  saveRDS(pheno_data, "GSE96804_pheno.rds")
}

# ---- 2. 样本信息 ----
cat("\n===== 样本信息 =====\n")
cat(sprintf("表达矩阵: %d 基因 x %d 样本\n", nrow(expr_data), ncol(expr_data)))

# 找到分组列
meta_cols <- colnames(pheno_data)
cat(sprintf("表型列 (%d): %s\n", length(meta_cols), paste(meta_cols, collapse = ", ")))

# 选择最可能有分组信息的列
group_col <- NA
# 优先检查 title（GSE96804的title含DN/Control）
for (pat in c("title", "characteristics", "source", "disease", "diagnosis")) {
  matches <- grep(pat, meta_cols, ignore.case = TRUE)
  if (length(matches) > 0) {
    group_col <- meta_cols[matches[1]]
    break
  }
}
if (is.na(group_col)) group_col <- meta_cols[1]
cat(sprintf("使用分组列: %s\n", group_col))

# 打印分组分布
group_values <- pheno_data[[group_col]]
cat("前10个样本的", group_col, ":\n")
print(head(group_values, 10))

# ---- 3. 箱线图 ----
cat("\n===== 原始表达箱线图 =====\n")

if (max(expr_data, na.rm = TRUE) > 100) {
  expr_log <- log2(expr_data + 1)
} else {
  expr_log <- expr_data
}

png("GSE96804_boxplot_raw.png", width = 1200, height = 600, res = 120)
boxplot(expr_log,
        main = "GSE96804 - Raw Expression Distribution (log2)",
        xlab = "Samples", ylab = "log2(Expression)",
        las = 2, cex.axis = 0.5, outline = FALSE,
        col = rep(c("steelblue", "tomato"), length.out = ncol(expr_log)))
dev.off()
cat("✅ 箱线图已保存\n")

# ---- 4. PCA图 ----
cat("\n===== PCA图 =====\n")

expr_var <- apply(expr_log, 1, var, na.rm = TRUE)
expr_filtered <- expr_log[expr_var > quantile(expr_var, 0.5, na.rm = TRUE), ]

pca <- prcomp(t(expr_filtered), scale. = TRUE)
pca_var <- round(summary(pca)$importance[2, 1:2] * 100, 1)

# 构建分组
group_val <- pheno_data[[group_col]]
group <- ifelse(grepl("DN|diabetic|diabetes|patient|DKD", group_val, ignore.case = TRUE),
                "DN", "Control")
cat(sprintf("分组结果: DN=%d, Control=%d\n", sum(group == "DN"), sum(group == "Control")))

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = group)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  labs(title = "GSE96804 PCA Plot",
       subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%%", pca_var[1], pca_var[2])) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("DN" = "tomato", "Control" = "steelblue"))
ggsave("GSE96804_PCA.png", width = 8, height = 6, dpi = 150)
cat("✅ PCA图已保存\n")

# ---- 5. 数据标准化 ----
cat("\n===== 数据标准化 (分位数标准化) =====\n")
expr_norm <- normalizeBetweenArrays(expr_data, method = "quantile")
saveRDS(expr_norm, file = "GSE96804_expr_normalized.rds")
cat("✅ 标准化表达矩阵已保存\n")

# ---- 6. 总结 ----
cat("\n========================================\n")
cat("01_GEO数据处理完成！\n")
cat(sprintf("  样本数: %d\n", ncol(expr_data)))
cat(sprintf("  基因数: %d\n", nrow(expr_data)))
cat(sprintf("  DN样本: %d | Control样本: %d\n",
            sum(group == "DN"), sum(group == "Control")))
cat("  输出文件:\n")
cat("    - GSE96804_expr_normalized.rds (标准化表达矩阵)\n")
cat("    - GSE96804_expr_raw.rds (原始表达矩阵)\n")
cat("    - GSE96804_pheno.rds (表型数据)\n")
cat("    - GSE96804_boxplot_raw.png (箱线图)\n")
cat("    - GSE96804_PCA.png (PCA图)\n")
cat("========================================\n")
