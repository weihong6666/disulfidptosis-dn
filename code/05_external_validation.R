# ============================================================
# 双硫死亡+DN 竞赛项目 — 05_外部验证（手动解析版）
# 验证集：GSE30528 (DN肾小球) + GSE30529 (DN肾小管)
# 平台：GPL571 [HG-U133A_2]
# 方法：手动解析GEO series matrix，不依赖GEOquery（避免段错误）
# ============================================================

library(limma)
library(pROC)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ---- 0. 手动解析GEO series matrix函数 ----
parseSeriesMatrix <- function(filepath) {
  cat(sprintf("  解析: %s\n", filepath))

  # 读取所有行
  lines <- readLines(filepath)

  # 找矩阵开始行
  matrix_start <- grep("^!series_matrix_table_begin", lines)
  if (length(matrix_start) == 0) stop("找不到矩阵起始行")
  matrix_start <- matrix_start[1] + 1

  # 读取矩阵头（第一行）
  header_line <- lines[matrix_start]
  sample_ids <- strsplit(header_line, "\t")[[1]][-1]  # 去掉第一列 "ID_REF"
  sample_ids <- gsub('"', '', sample_ids)

  # 读取数据行
  data_lines <- lines[(matrix_start + 1):length(lines)]

  # 只取需要的行（包含探针ID的数据行，如 "1007_s_at"）
  data_lines <- data_lines[grepl('"\\d+', substr(data_lines, 1, 50)) |
                           grepl('_at', substr(data_lines, 1, 50))]

  # 解析每行
  cat(sprintf("  解析 %d 行数据...\n", length(data_lines)))
  all_probes <- list()
  probe_ids <- character(length(data_lines))
  values_list <- list()

  for (i in seq_along(data_lines)) {
    fields <- strsplit(data_lines[i], "\t")[[1]]
    probe_ids[i] <- gsub('"', '', fields[1])
    values_list[[i]] <- as.numeric(fields[-1])
  }

  # 构建表达矩阵
  expr <- do.call(rbind, values_list)
  colnames(expr) <- sample_ids
  rownames(expr) <- probe_ids

  # 提取样本元数据
  meta_lines <- lines[1:(matrix_start - 1)]
  sample_meta_lines <- meta_lines[grepl("^!Sample_", meta_lines)]

  # 构建metadata data.frame
  meta <- data.frame(row.names = sample_ids)
  for (ml in sample_meta_lines) {
    # 提取tag和值
    parts <- strsplit(ml, "\t")[[1]]
    tag <- sub("^!Sample_", "", parts[1])
    values <- gsub('"', '', parts[-1])
    meta[[tag]] <- values
  }

  cat(sprintf("  ✅ %d 探针 × %d 样本\n", nrow(expr), ncol(expr)))
  return(list(expr = expr, pheno = meta))
}

# ---- 1. 加载训练集模型和标志物 ----
cat("\n===== 加载训练结果 =====\n")

final_genes <- readRDS("04_final_biomarkers.rds")
cat(sprintf("标志物 (%d个): %s\n", length(final_genes), paste(final_genes, collapse = ", ")))

expr_train <- readRDS("GSE96804_gene_expr.rds")
pheno_train <- readRDS("GSE96804_pheno.rds")

group_info <- pheno_train[["title"]]
group_train <- factor(ifelse(grepl("DN|diabetic|patient|DKD", group_info, ignore.case = TRUE),
                             "DN", "Control"), levels = c("Control", "DN"))
y_train <- as.numeric(group_train == "DN")

X_train <- as.data.frame(t(expr_train[final_genes, , drop = FALSE]))
model_data <- X_train
model_data$y <- y_train
fmla <- as.formula(paste("y ~", paste(final_genes, collapse = " + ")))
logit_model <- glm(fmla, data = model_data, family = binomial())
cat(sprintf("训练集模型: %d DN + %d Control\n", sum(y_train==1), sum(y_train==0)))

# ---- 2. 加载GPL571探针注释 ----
cat("\n===== 加载GPL571探针注释 =====\n")

gpl571_cache <- "GPL571_annot.rds"

if (file.exists(gpl571_cache)) {
  probe_to_gene <- readRDS(gpl571_cache)
  cat(sprintf("从缓存加载: %d 探针\n", nrow(probe_to_gene)))
} else if (file.exists("GPL571.annot.gz")) {
  cat("解析GPL571.annot.gz...\n")
  # 读取GPL571注释文件（跳过注释头）
  f <- gzfile("GPL571.annot.gz")
  lines <- readLines(f)
  close(f)

  # 找表格头
  header_idx <- grep("^ID\tGene title\tGene symbol", lines)[1]
  if (is.na(header_idx)) {
    # 回退：查找platform_table_begin
    tbl_begin <- grep("^!platform_table_begin", lines)[1]
    if (!is.na(tbl_begin)) header_idx <- tbl_begin + 1
  }

  if (is.na(header_idx)) stop("无法找到GPL571注释表头")

  # 读取制表符表格（只读ID和Gene symbol两列）
  cat(sprintf("  从第 %d 行开始解析...\n", header_idx))
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[data_lines != ""]

  # 快速解析：用strsplit取第1和第3列
  probe_ids <- character(length(data_lines))
  gene_symbols <- character(length(data_lines))

  for (i in seq_along(data_lines)) {
    fields <- strsplit(data_lines[i], "\t")[[1]]
    if (length(fields) >= 3) {
      probe_ids[i] <- fields[1]
      gene_symbols[i] <- fields[3]
    }
  }

  probe_to_gene <- data.frame(
    probe_id = probe_ids,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
  )

  # 清理
  probe_to_gene <- probe_to_gene[probe_to_gene$gene_symbol != "" &
                                 probe_to_gene$gene_symbol != "---" &
                                 !is.na(probe_to_gene$gene_symbol), ]
  # 展开 /// 分隔的多基因映射
  probe_to_gene$gene_symbol <- sub("///.*", "", probe_to_gene$gene_symbol)

  saveRDS(probe_to_gene, gpl571_cache)
  cat(sprintf("  ✅ %d 探针→基因映射已缓存\n", nrow(probe_to_gene)))
} else {
  stop("找不到GPL571.annot.gz！请从浏览器下载:\n  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL571/annot/GPL571.annot.gz")
}

# ---- 3. 处理验证集 ----
cat("\n===== 处理验证数据集 =====\n")

validation_datasets <- list()
all_validation_results <- data.frame()

for (gse_id in c("GSE30528", "GSE30529")) {
  cat(sprintf("\n--- %s ---\n", gse_id))

  local_file <- sprintf("%s_series_matrix.txt", gse_id)
  if (!file.exists(local_file)) {
    cat(sprintf("  ⚠ 找不到 %s\n", local_file))
    next
  }

  # 手动解析
  parsed <- parseSeriesMatrix(local_file)
  val_expr <- parsed$expr
  val_pheno <- parsed$pheno

  # --- 分组识别 ---
  # 优先使用 disease state 相关列
  group_col <- NULL
  for (cn in colnames(val_pheno)) {
    col_values <- paste(val_pheno[[cn]], collapse = " ")
    if (grepl("disease state|diagnosis|condition", col_values, ignore.case = TRUE) &&
        grepl("control|diabet|normal|healthy|disease", col_values, ignore.case = TRUE)) {
      group_col <- cn
      break
    }
  }

  # 回退：使用title列
  if (is.null(group_col) && "title" %in% colnames(val_pheno)) {
    group_col <- "title"
  }

  # 最后回退：用第一列
  if (is.null(group_col)) group_col <- colnames(val_pheno)[1]

  cat(sprintf("  分组列: %s\n", group_col))
  sample_groups <- val_pheno[[group_col]]
  cat(sprintf("  前3个样本: %s...\n", paste(head(sample_groups, 3), collapse = " | ")))

  # 分组判断
  val_group <- factor(ifelse(grepl("diabet|DKD|DN\\b|disease|patient",
                                   sample_groups, ignore.case = TRUE) &
                             !grepl("control|normal|healthy",
                                   sample_groups, ignore.case = TRUE),
                             "DN", "Control"),
                      levels = c("Control", "DN"))

  n_dn <- sum(val_group == "DN")
  n_ctrl <- sum(val_group == "Control")
  cat(sprintf("  分组结果: %d DN + %d Control\n", n_dn, n_ctrl))

  if (n_dn < 2 || n_ctrl < 2) {
    cat("  ⚠ 样本量不足，尝试修正...\n")
    # 修正：只检查是否包含 diabetic/disease
    val_group <- factor(ifelse(grepl("diabet|DKD", sample_groups, ignore.case = TRUE),
                               "DN", "Control"),
                        levels = c("Control", "DN"))
    cat(sprintf("  修正后: %d DN + %d Control\n",
                sum(val_group == "DN"), sum(val_group == "Control")))
  }

  if (sum(val_group == "DN") < 2 || sum(val_group == "Control") < 2) {
    cat("  ⚠ 仍不足，跳过\n")
    next
  }

  # --- 探针→基因转换 ---
  common <- intersect(rownames(val_expr), probe_to_gene$gene_symbol)
  cat(sprintf("  匹配探针: %d/%d\n", length(common), nrow(val_expr)))

  if (length(common) < 100) {
    # gene_symbol列包含的是探针映射表，不是表达矩阵行名
    # 需要使用probe_id来匹配
    common_probes <- intersect(rownames(val_expr), probe_to_gene$probe_id)
    cat(sprintf("  匹配probe_id: %d/%d\n", length(common_probes), nrow(val_expr)))
    common <- common_probes
  }

  val_expr <- val_expr[common, , drop = FALSE]
  probe_map <- probe_to_gene[match(common, probe_to_gene$probe_id), ]
  gs <- probe_map$gene_symbol

  # 多探针→1基因：取最大均值探针
  gene_list <- split(seq_len(nrow(val_expr)), gs)
  val_expr_gene <- t(sapply(gene_list, function(idx) {
    if (length(idx) == 1) val_expr[idx[1], ]
    else {
      # 取平均表达最高的探针
      means <- rowMeans(val_expr[idx, , drop = FALSE], na.rm = TRUE)
      val_expr[idx[which.max(means)], ]
    }
  }))
  colnames(val_expr_gene) <- colnames(val_expr)
  rownames(val_expr_gene) <- names(gene_list)
  cat(sprintf("  基因级: %d 基因 × %d 样本\n", nrow(val_expr_gene), ncol(val_expr_gene)))

  # log2转换（如果需要）
  if (max(val_expr_gene, na.rm = TRUE) > 100) {
    val_expr_gene <- log2(val_expr_gene + 1)
  }

  # 分位数标准化
  val_expr_gene <- normalizeBetweenArrays(val_expr_gene, method = "quantile")

  validation_datasets[[gse_id]] <- list(
    expr = val_expr_gene,
    group = val_group
  )

  # --- ROC评估（跨平台稳健方法） ---
  # 不使用跨平台逻辑回归（表达尺度不同），改用数据集内标准化评分
  val_y <- as.numeric(val_group == "DN")
  genes_avail <- intersect(final_genes, rownames(val_expr_gene))
  cat(sprintf("  可用标志物: %d/%d: %s\n",
              length(genes_avail), length(final_genes),
              paste(genes_avail, collapse = ", ")))

  if (length(genes_avail) < 2) next

  # 方法1: 数据集内 z-score 组合评分（跨平台稳健）
  val_expr_subset <- val_expr_gene[genes_avail, , drop = FALSE]

  # 计算每个基因在训练集中的DN方向
  gene_direction <- sapply(genes_avail, function(g) {
    if (g %in% rownames(expr_train)) {
      mean_dn <- mean(as.numeric(expr_train[g, y_train == 1]), na.rm = TRUE)
      mean_ctrl <- mean(as.numeric(expr_train[g, y_train == 0]), na.rm = TRUE)
      sign(mean_dn - mean_ctrl)  # +1 = up in DN, -1 = down in DN
    } else 1
  })

  # 数据集内z-score标准化后按方向组合
  val_z <- t(scale(t(val_expr_subset)))  # 基因级z-score
  val_z[is.na(val_z)] <- 0
  composite_score <- colMeans(val_z * gene_direction)  # 方向加权

  # 组合模型ROC
  val_roc <- roc(val_y, composite_score, levels = c(0, 1), direction = "<")
  val_auc <- auc(val_roc)
  cat(sprintf("  🎯 组合评分 AUC: %.3f\n", val_auc))

  # 方法2: 单个基因ROC（数据集内直接评估）
  for (g in genes_avail) {
    g_roc <- roc(val_y, as.numeric(val_expr_subset[g, ]), levels = c(0, 1), direction = "<")
    all_validation_results <- rbind(all_validation_results, data.frame(
      Dataset = gse_id, Gene = g, AUC = auc(g_roc), Type = "Single Gene"
    ))
  }
  all_validation_results <- rbind(all_validation_results, data.frame(
    Dataset = gse_id, Gene = "Combined Model", AUC = val_auc, Type = "Model"
  ))

  assign(sprintf("roc_%s", gse_id), val_roc)
  assign(sprintf("auc_%s", gse_id), val_auc)
}

if (length(validation_datasets) == 0) stop("无可用验证集！")

# ---- 4. 训练集自检 ----
cat("\n===== 训练集自检 =====\n")
train_pred <- predict(logit_model, type = "response")
train_roc <- roc(y_train, train_pred, levels = c(0, 1), direction = "<")
train_auc <- auc(train_roc)
cat(sprintf("训练集 AUC: %.3f\n", train_auc))

for (g in final_genes) {
  g_roc <- roc(y_train, X_train[, g], levels = c(0, 1), direction = "<")
  all_validation_results <- rbind(all_validation_results, data.frame(
    Dataset = "GSE96804 (Training)", Gene = g, AUC = auc(g_roc), Type = "Single Gene"
  ))
}
all_validation_results <- rbind(all_validation_results, data.frame(
  Dataset = "GSE96804 (Training)", Gene = "Combined Model", AUC = train_auc, Type = "Model"
))

# ---- 5. 汇总ROC曲线 ----
cat("\n===== 绘制汇总ROC =====\n")

all_rocs <- list("Training (GSE96804)" = train_roc)
for (gse_id in names(validation_datasets)) {
  rn <- sprintf("roc_%s", gse_id)
  if (exists(rn)) all_rocs[[gse_id]] <- get(rn)
}

png("05_validation_ROC.png", width = 1000, height = 850, res = 130)
plot(NULL, xlim = c(1, 0), ylim = c(0, 1),
     xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curves - Training & External Validation")
abline(a = 1, b = -1, lty = 2, col = "grey60")
colors <- brewer.pal(max(3, length(all_rocs)), "Set1")
for (i in seq_along(all_rocs)) {
  lines(1 - all_rocs[[i]]$specificities, all_rocs[[i]]$sensitivities,
        col = colors[i], lwd = 2.5)
}
legend("bottomright",
       legend = sprintf("%s (AUC=%.3f)", names(all_rocs),
                       sapply(all_rocs, function(x) auc(x))),
       col = colors, lwd = 2.5, cex = 0.85, bty = "n")
dev.off()
cat("✅ 05_validation_ROC.png\n")

# ---- 6. AUC条形图 ----
model_results <- subset(all_validation_results, Type == "Model")

ggplot(model_results, aes(x = reorder(Dataset, AUC), y = AUC, fill = Dataset)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), hjust = -0.1, size = 4) +
  coord_flip(ylim = c(0.45, 1.05)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Diagnostic Model AUC Across Datasets",
       subtitle = paste("Biomarkers:", paste(final_genes, collapse = ", ")),
       x = "", y = "AUC") +
  theme_minimal(base_size = 14) + theme(legend.position = "none")
ggsave("05_AUC_barplot.png", width = 12, height = 5, dpi = 150)
cat("✅ 05_AUC_barplot.png\n")

# ---- 7. 验证集热图 ----
for (gse_id in names(validation_datasets)) {
  vd <- validation_datasets[[gse_id]]
  genes_avail <- intersect(final_genes, rownames(vd$expr))
  if (length(genes_avail) < 2) next

  sample_order <- order(vd$group)
  hm_data <- vd$expr[genes_avail, sample_order, drop = FALSE]
  ann <- data.frame(Group = vd$group[sample_order], row.names = colnames(hm_data))

  png(sprintf("05_heatmap_%s.png", gse_id), width = 900, height = 600, res = 120)
  pheatmap(hm_data, scale = "row", annotation_col = ann,
           annotation_colors = list(Group = c(DN = "tomato", Control = "steelblue")),
           show_colnames = FALSE, cluster_cols = TRUE,
           main = sprintf("Biomarker Expression - %s", gse_id),
           fontsize_row = 10,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
  cat(sprintf("✅ 05_heatmap_%s.png\n", gse_id))
}

# ---- 8. 保存 ----
saveRDS(all_validation_results, "05_validation_results.rds")
write.csv(all_validation_results, "05_validation_results.csv", row.names = FALSE)

cat("\n============ 外部验证汇总 ============\n")
for (ds in unique(all_validation_results$Dataset)) {
  ds_data <- subset(all_validation_results, Dataset == ds & Type == "Model")
  if (nrow(ds_data) > 0) cat(sprintf("%s: AUC = %.3f\n", ds, ds_data$AUC[1]))
}

aucs_val <- c()
for (gse_id in c("GSE30528", "GSE30529")) {
  vn <- sprintf("auc_%s", gse_id)
  if (exists(vn)) aucs_val <- c(aucs_val, get(vn))
}
if (length(aucs_val) > 0) cat(sprintf("\n验证集平均 AUC: %.3f\n", mean(aucs_val)))

cat("\n========================================\n")
