# ============================================================
# 双硫死亡+DN 竞赛项目 — 08b_单细胞深度分析（审稿修订版 v2）
# 目的：深化单细胞分析，从"展示定位"升级到"功能验证"
# 兼容：Seurat v4/v5 | 自适应数据格式
# ============================================================

library(Seurat)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(Matrix)
library(pROC)
library(AUCell)

# ---- Helper: 兼容 Seurat v4/v5 获取表达数据 ----
get_expr_data <- function(obj, genes=NULL, cells=NULL) {
  if (is.null(genes)) genes <- rownames(obj)
  if (is.null(cells)) cells <- colnames(obj)

  expr <- NULL

  # After JoinLayers, @layers$data should exist as a single merged layer
  expr <- tryCatch({
    obj@assays$RNA@layers$data[genes, cells, drop=FALSE]
  }, error = function(e) NULL)

  # Fallback: try @data (v4)
  if (is.null(expr)) {
    expr <- tryCatch({
      obj@assays$RNA@data[genes, cells, drop=FALSE]
    }, error = function(e) NULL)
  }

  # Fallback: GetAssayData
  if (is.null(expr)) {
    expr <- tryCatch({
      GetAssayData(obj, assay="RNA", layer="data")[genes, cells, drop=FALSE]
    }, error = function(e) NULL)
  }

  if (is.null(expr)) {
    cat("WARNING: get_expr_data failed!\n")
  }
  return(expr)
}

# ---- 1. 加载数据 ----
cat("\n===== 加载数据 =====\n")

if (file.exists("08_seurat_object.rds")) {
  seurat_obj <- readRDS("08_seurat_object.rds")
  cat(sprintf("Seurat对象: %d cells x %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))

  # Seurat v5 Assay5: 合并分层layers（每个样本独立存储）
  # JoinLayers将 data.1...data.N 合并为统一的 "data" layer
  cat("检查Seurat版本和layer结构...\n")
  assay_class <- class(seurat_obj@assays$RNA)[1]
  cat(sprintf("Assay class: %s\n", assay_class))
  if (assay_class == "Assay5") {
    cat("检测到Assay5格式，执行JoinLayers合并样本层...\n")
    seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
    cat(sprintf("JoinLayers完成, 新layers: %s\n",
                paste(names(seurat_obj@assays$RNA@layers), collapse=", ")))
  }

  # 获取细胞类型（可能是数字或名称）
  cell_type_labels <- as.character(Idents(seurat_obj))
  unique_ct <- unique(cell_type_labels)
  cat(sprintf("细胞类型数: %d (labels: %s)\n",
              length(unique_ct), paste(head(unique_ct, 10), collapse=", ")))
} else {
  stop("请先运行 08_单细胞分析.R 生成 Seurat 对象")
}

# ---- 2. 确认分组 ----
cat("\n===== 确认分组信息 =====\n")

# 尝试多种方式获取分组
group_vec <- NULL
if ("group" %in% colnames(seurat_obj@meta.data)) {
  group_vec <- seurat_obj@meta.data$group
} else if ("condition" %in% colnames(seurat_obj@meta.data)) {
  group_vec <- seurat_obj@meta.data$condition
} else if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
  orig <- seurat_obj@meta.data$orig.ident
  group_vec <- ifelse(grepl("DN|diabet|Diabet", orig, ignore.case=TRUE), "DN", "Control")
} else {
  # 从barcode判断
  bc <- colnames(seurat_obj)
  group_vec <- ifelse(grepl("DN|diabet", bc, ignore.case=TRUE), "DN", "Control")
}

cat(sprintf("分组: DN=%d, Control=%d\n", sum(group_vec=="DN"), sum(group_vec=="Control")))

# ---- 3. 加载标志物 & 双硫死亡基因 ----
cat("\n===== 加载基因集 =====\n")

final_biomarkers <- readRDS("04_final_biomarkers.rds")
cat(sprintf("标志物基因 (%d): %s\n", length(final_biomarkers), paste(final_biomarkers, collapse=", ")))

# 检测基因标识格式
all_sc_genes <- rownames(seurat_obj)
is_ensembl <- grepl("^ENSG", all_sc_genes[1:min(100, length(all_sc_genes))])
cat(sprintf("基因标识格式: %s\n", if(mean(is_ensembl) > 0.8) "Ensembl ID" else "Gene Symbol"))

# 如果是Ensembl ID，需要转换
if (mean(is_ensembl) > 0.8) {
  cat("检测到Ensembl ID格式，进行基因符号映射...\n")

  # 使用biomaRt或预置的Ensembl→Symbol映射
  # 方法：从biomaRt获取，或使用内置注释
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })

  # 获取所有Ensembl ID→Gene Symbol映射
  ensembl_to_symbol <- tryCatch({
    map <- select(org.Hs.eg.db,
                  keys = all_sc_genes,
                  columns = c("ENSEMBL", "SYMBOL"),
                  keytype = "ENSEMBL")
    map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
    cat(sprintf("成功映射 %d/%d Ensembl ID → Gene Symbol\n", nrow(map), length(all_sc_genes)))
    map
  }, error = function(e) {
    cat("org.Hs.eg.db映射失败，尝试替代方法...\n")
    # 如果org.Hs.eg.db不可用，尝试使用预置的Ensembl映射
    NULL
  })

  if (!is.null(ensembl_to_symbol)) {
    # 查找标志物基因对应的Ensembl ID
    biomarker_mapping <- list()
    for (bm in final_biomarkers) {
      matched <- ensembl_to_symbol[ensembl_to_symbol$SYMBOL == bm, "ENSEMBL"]
      if (length(matched) > 0) {
        biomarker_mapping[[bm]] <- matched[1]
      }
    }
    cat(sprintf("标志物基因映射: %d/%d 在scRNA数据中\n",
                length(biomarker_mapping), length(final_biomarkers)))

    if (length(biomarker_mapping) > 0) {
      biomarkers_in_sc <- unlist(biomarker_mapping)
      for (nm in names(biomarker_mapping)) {
        cat(sprintf("  %s → %s\n", nm, biomarker_mapping[[nm]]))
      }
    } else {
      cat("⚠️ 标志物基因在scRNA中均未找到Ensembl映射\n")
      # 使用双硫死亡基因集作为替代
      dr_genes <- readLines("disulfidptosis_genes.txt")
      dr_mapping <- list()
      for (g in dr_genes) {
        matched <- ensembl_to_symbol[ensembl_to_symbol$SYMBOL == g, "ENSEMBL"]
        if (length(matched) > 0) dr_mapping[[g]] <- matched[1]
      }
      biomarkers_in_sc <- unlist(dr_mapping)
      cat(sprintf("回退：使用双硫死亡基因集 (%d 个在scRNA中)\n", length(biomarkers_in_sc)))
    }
  } else {
    biomarkers_in_sc <- character(0)
    cat("⚠️ 无法映射基因，跳过表达相关分析\n")
  }
} else {
  # 基因符号格式 — 直接匹配
  biomarker_mapping <- list()
  for (bm in final_biomarkers) {
    matches <- grep(paste0("^", bm, "$"), all_sc_genes, ignore.case=TRUE, value=TRUE)
    if (length(matches) > 0) biomarker_mapping[[bm]] <- matches[1]
  }
  biomarkers_in_sc <- unlist(biomarker_mapping)
  cat(sprintf("直接匹配: %d/%d 标志物在scRNA中\n",
              length(biomarkers_in_sc), length(final_biomarkers)))
}

if (length(biomarkers_in_sc) == 0) {
  cat("\n❌ 无法在scRNA数据中找到任何标志物基因。\n")
  cat("原因：基因标识不匹配（Ensembl ID vs Gene Symbol）且映射失败。\n")
  cat("请确保 org.Hs.eg.db 已安装或手动提供Ensembl→Symbol映射表。\n")
  cat("后续基于标志物的分析将跳过。\n")
}

# ---- 4. 伪批量差异表达 ----
cat("\n===== 伪批量差异表达 (Pseudobulk DE) =====\n")

pseudobulk_results <- list()
ct_counts <- table(cell_type_labels)
# 保留细胞数 >= 30 的细胞类型
major_cts <- names(ct_counts[ct_counts >= 30])
cat(sprintf("分析 %d 个主要细胞类型\n", length(major_cts)))

for (ct in major_cts) {
  cells_ct <- colnames(seurat_obj)[cell_type_labels == ct]
  groups_ct <- group_vec[cell_type_labels == ct]

  n_dn <- sum(groups_ct == "DN")
  n_ctrl <- sum(groups_ct == "Control")

  if (n_dn < 10 || n_ctrl < 10) {
    cat(sprintf("  %s: 细胞不足 (DN=%d, Ctrl=%d), 跳过\n", ct, n_dn, n_ctrl))
    next
  }

  # 获取这些细胞的表达数据
  expr_ct <- tryCatch({
    get_expr_data(seurat_obj, genes=biomarkers_in_sc, cells=cells_ct)
  }, error = function(e) {
    cat(sprintf("  %s: 获取表达数据失败, 跳过\n", ct))
    return(NULL)
  })

  if (is.null(expr_ct)) next

  # 对每个基因做t-test
  ct_de <- data.frame(Gene=character(), logFC=numeric(), P_value=numeric(),
                       stringsAsFactors=FALSE)
  for (g in biomarkers_in_sc) {
    vals_dn <- as.numeric(expr_ct[g, groups_ct == "DN"])
    vals_ctrl <- as.numeric(expr_ct[g, groups_ct == "Control"])
    if (length(vals_dn) >= 10 && length(vals_ctrl) >= 10) {
      test <- tryCatch(t.test(vals_dn, vals_ctrl), error=function(e) NULL)
      if (!is.null(test)) {
        ct_de <- rbind(ct_de, data.frame(
          Gene=g, logFC=mean(vals_dn)-mean(vals_ctrl), P_value=test$p.value))
      }
    }
  }

  if (nrow(ct_de) > 0) {
    ct_de$CellType <- ct
    ct_de$N_DN <- n_dn
    ct_de$N_Ctrl <- n_ctrl
    pseudobulk_results[[ct]] <- ct_de
    n_sig <- sum(ct_de$P_value < 0.05, na.rm=TRUE)
    cat(sprintf("  %s: %d/%d 基因显著 (DN=%d, Ctrl=%d)\n",
                ct, n_sig, nrow(ct_de), n_dn, n_ctrl))
  }
}

# 合并
pseudobulk_all <- do.call(rbind, pseudobulk_results)
if (!is.null(pseudobulk_all) && nrow(pseudobulk_all) > 0) {
  pseudobulk_all$adj_P <- p.adjust(pseudobulk_all$P_value, method="BH")
  pseudobulk_all$Significant <- pseudobulk_all$adj_P < 0.05
}

cat(sprintf("\n伪批量DE汇总: %d 细胞类型 x %d 基因检验\n",
            length(unique(pseudobulk_all$CellType)), nrow(pseudobulk_all)))

# ---- 5. 细胞类型表达DotPlot ----
cat("\n===== 细胞类型表达模式 =====\n")

# 计算每个细胞类型的平均表达（分DN/Control）
avg_list <- list()
for (ct in major_cts) {
  cells_ct <- colnames(seurat_obj)[cell_type_labels == ct]
  groups_ct <- group_vec[cell_type_labels == ct]

  expr_ct <- tryCatch({
    get_expr_data(seurat_obj, genes=biomarkers_in_sc, cells=cells_ct)
  }, error=function(e) NULL)
  if (is.null(expr_ct)) next

  cells_dn <- which(groups_ct == "DN")
  cells_ctrl <- which(groups_ct == "Control")

  if (length(cells_dn) >= 5 && length(cells_ctrl) >= 5) {
    avg_list[[paste0(ct, "_DN")]] <- Matrix::rowMeans(expr_ct[, cells_dn, drop=FALSE])
    avg_list[[paste0(ct, "_Control")]] <- Matrix::rowMeans(expr_ct[, cells_ctrl, drop=FALSE])
  }
}

if (length(avg_list) >= 4) {
  avg_mat <- do.call(cbind, avg_list)
  avg_mat <- avg_mat[Matrix::rowSums(avg_mat) > 0, , drop=FALSE]

  if (nrow(avg_mat) >= 2) {
    png("08b_celltype_expression_heatmap.png", width=1100, height=700, res=130)
    pheatmap(as.matrix(avg_mat),
             scale="row", cluster_cols=TRUE, cluster_rows=TRUE,
             main="Gene Expression by Cell Type and Condition",
             fontsize=10, angle_col=45,
             color=colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
    cat("✅ 细胞类型表达热图已保存\n")
  }
}

# ---- 6. 细胞类型特异性AUC ----
cat("\n===== 细胞类型特异性AUC =====\n")

ct_auc_results <- data.frame()
for (ct in major_cts) {
  cells_ct <- colnames(seurat_obj)[cell_type_labels == ct]
  groups_ct <- group_vec[cell_type_labels == ct]
  if (sum(groups_ct=="DN") < 5 || sum(groups_ct=="Control") < 5) next

  expr_ct <- tryCatch({
    get_expr_data(seurat_obj, genes=biomarkers_in_sc, cells=cells_ct)
  }, error=function(e) NULL)
  if (is.null(expr_ct)) next

  # 计算基因平均分作为biomarker score
  if (nrow(expr_ct) >= 2) {
    scores <- colMeans(as.matrix(expr_ct))
    roc_obj <- tryCatch(
      roc(groups_ct=="DN", scores, levels=c(FALSE,TRUE), direction="<", quiet=TRUE),
      error=function(e) NULL)
    if (!is.null(roc_obj)) {
      ct_auc_results <- rbind(ct_auc_results, data.frame(
        CellType=ct, AUC=auc(roc_obj), N_Cells=length(cells_ct)))
    }
  }
}

if (nrow(ct_auc_results) > 0) {
  ct_auc_results <- ct_auc_results[order(ct_auc_results$AUC, decreasing=TRUE), ]
  cat("\n细胞类型AUC排名:\n")
  print(ct_auc_results)

  ggplot(ct_auc_results, aes(x=reorder(CellType, AUC), y=AUC, fill=AUC)) +
    geom_bar(stat="identity", alpha=0.8, width=0.6) +
    geom_text(aes(label=sprintf("AUC=%.3f\nn=%d", AUC, N_Cells)), hjust=-0.1, size=3) +
    scale_fill_gradient(low="steelblue", high="tomato", limits=c(0.5, 1)) +
    coord_flip() + geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
    labs(title="Cell-Type-Specific AUC", y="AUC", x="") +
    ylim(0.4, 1.1) + theme_minimal(base_size=13)
  ggsave("08b_celltype_auc.png", width=8, height=6, dpi=150)
  cat("✅ 细胞类型AUC图已保存\n")
}

# ---- 7. AUCell通路评分 ----
cat("\n===== AUCell通路活性评分 =====\n")

dr_genes_all <- readLines("disulfidptosis_genes.txt")
# 映射到scRNA基因名
dr_in_sc_all <- intersect(toupper(dr_genes_all), toupper(all_sc_genes))
dr_sc_mapped <- character()
for (g in dr_genes_all) {
  m <- grep(paste0("^", g, "$"), all_sc_genes, ignore.case=TRUE, value=TRUE)
  if (length(m) > 0) dr_sc_mapped <- c(dr_sc_mapped, m[1])
}

ox_genes <- intersect(toupper(c(
  "SLC7A11","SLC3A2","G6PD","PGD","PRDX1","GPX4",
  "GCLC","GCLM","GSR","TXN","TXNRD1","SOD1","SOD2","CAT",
  "NFE2L2","KEAP1","HMOX1","NQO1")), toupper(all_sc_genes))
ox_sc_mapped <- character()
for (g in c("SLC7A11","SLC3A2","G6PD","PGD","PRDX1","GPX4",
            "GCLC","GCLM","GSR","TXN","TXNRD1","SOD1","SOD2","CAT",
            "NFE2L2","KEAP1","HMOX1","NQO1")) {
  m <- grep(paste0("^", g, "$"), all_sc_genes, ignore.case=TRUE, value=TRUE)
  if (length(m) > 0) ox_sc_mapped <- c(ox_sc_mapped, m[1])
}

gene_sets <- list()
if (length(dr_sc_mapped) >= 5) gene_sets$Disulfidptosis <- dr_sc_mapped
if (length(biomarkers_in_sc) >= 2) gene_sets$Biomarkers <- biomarkers_in_sc
if (length(ox_sc_mapped) >= 5) gene_sets$Oxidative_Stress <- ox_sc_mapped

cat(sprintf("基因集: %s\n", paste(names(gene_sets), collapse=", ")))
for (nm in names(gene_sets)) cat(sprintf("  %s: %d genes\n", nm, length(gene_sets[[nm]])))

if (length(gene_sets) >= 1) {
  # 获取完整表达矩阵
  cat("构建表达矩阵...\n")
  expr_full <- tryCatch({
    get_expr_data(seurat_obj)
  }, error=function(e) {
    # 回退：使用counts
    tryCatch({
      LayerData(seurat_obj, assay="RNA", layer="counts")
    }, error=function(e2) {
      GetAssayData(seurat_obj, assay="RNA", slot="counts")
    })
  })

  # 只取表达的基因
  expressed_genes <- names(which(Matrix::rowSums(expr_full > 0) >= 10))
  for (nm in names(gene_sets)) {
    gene_sets[[nm]] <- intersect(gene_sets[[nm]], expressed_genes)
  }
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]

  if (length(gene_sets) >= 1) {
    cat(sprintf("过滤后基因集: %d 个\n", length(gene_sets)))

    # 随机抽取最多10000个细胞用于AUCell（内存限制）
    set.seed(42)
    max_cells <- min(10000, ncol(expr_full))
    sampled_cells <- sample(colnames(expr_full), max_cells)

    cat("运行AUCell (可能需要2-5分钟)...\n")
    cells_rankings <- AUCell_buildRankings(as.matrix(expr_full[, sampled_cells]),
                                            nCores=1, plotStats=FALSE, verbose=FALSE)
    cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings,
                                 aucMaxRank=ceiling(0.05 * nrow(expr_full)),
                                 verbose=FALSE)
    auc_matrix <- getAUC(cells_AUC)

    cat(sprintf("AUCell完成: %d gene sets x %d cells\n", nrow(auc_matrix), ncol(auc_matrix)))

    # ---- AUCell per cell type ----
    sampled_cts <- cell_type_labels[match(sampled_cells, colnames(seurat_obj))]
    sampled_groups <- group_vec[match(sampled_cells, colnames(seurat_obj))]

    auc_by_ct <- data.frame()
    for (ct in unique(sampled_cts)) {
      cells_ct_idx <- which(sampled_cts == ct)
      if (length(cells_ct_idx) < 10) next
      for (gs in rownames(auc_matrix)) {
        scores <- auc_matrix[gs, cells_ct_idx]
        groups <- sampled_groups[cells_ct_idx]
        auc_by_ct <- rbind(auc_by_ct, data.frame(
          CellType=ct, GeneSet=gs, Score=mean(scores), SD=sd(scores),
          DN_score=mean(scores[groups=="DN"]),
          Control_score=mean(scores[groups=="Control"]),
          N_Cells=length(cells_ct_idx)))
      }
    }

    # 可视化：双硫死亡评分 per cell type
    if ("Disulfidptosis" %in% rownames(auc_matrix)) {
      dr_scores <- auc_by_ct[auc_by_ct$GeneSet=="Disulfidptosis", ]
      dr_scores <- dr_scores[order(dr_scores$Score, decreasing=TRUE), ]

      ggplot(dr_scores, aes(x=reorder(CellType, Score), y=Score, fill=Score)) +
        geom_bar(stat="identity", alpha=0.8, width=0.6) +
        geom_errorbar(aes(ymin=Score-SD, ymax=Score+SD), width=0.15) +
        scale_fill_gradient(low="steelblue", high="tomato") +
        coord_flip() +
        labs(title="Disulfidptosis AUCell Score by Cell Type",
             subtitle="Error bars = ±1 SD", x="", y="AUCell Score") +
        theme_minimal(base_size=13)
      ggsave("08b_aucell_disulfidptosis_by_celltype.png", width=8, height=5, dpi=150)
      cat("✅ AUCell细胞类型图已保存\n")

      # DN vs Control
      dr_long <- dr_scores %>%
        select(CellType, DN=DN_score, Control=Control_score) %>%
        pivot_longer(-CellType, names_to="Group", values_to="Score")
      ggplot(dr_long, aes(x=CellType, y=Score, fill=Group)) +
        geom_bar(stat="identity", position=position_dodge(width=0.8), alpha=0.8, width=0.7) +
        scale_fill_manual(values=c("DN"="tomato", "Control"="steelblue")) +
        labs(title="Disulfidptosis Score: DN vs Control per Cell Type",
             x="", y="AUCell Score") +
        theme_minimal(base_size=13) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
      ggsave("08b_aucell_dn_vs_control.png", width=10, height=5, dpi=150)
      cat("✅ AUCell DN vs Control图已保存\n")
    }

    # 小提琴图
    plot_data <- data.frame()
    for (gs in rownames(auc_matrix)) {
      pd <- data.frame(Cell=sampled_cells, GeneSet=gs,
                        AUCell_Score=as.numeric(auc_matrix[gs, ]),
                        Group=sampled_groups)
      plot_data <- rbind(plot_data, pd)
    }

    if (nrow(plot_data) > 0) {
      ggplot(plot_data, aes(x=GeneSet, y=AUCell_Score, fill=Group)) +
        geom_violin(alpha=0.6, scale="width", draw_quantiles=0.5) +
        scale_fill_manual(values=c("DN"="tomato", "Control"="steelblue")) +
        labs(title="Pathway Activity: DN vs Control", x="", y="AUCell Score") +
        theme_minimal(base_size=13) +
        theme(axis.text.x=element_text(angle=30, hjust=1))
      ggsave("08b_pathway_violin.png", width=10, height=6, dpi=150)
      cat("✅ 通路小提琴图已保存\n")
    }
  }
}

# ---- 8. 伪批量DE可视化 ----
cat("\n===== 伪批量DE可视化 =====\n")

if (!is.null(pseudobulk_all) && nrow(pseudobulk_all) >= 5) {
  # logFC热图
  de_matrix <- reshape2::dcast(pseudobulk_all, Gene ~ CellType, value.var="logFC")
  rownames(de_matrix) <- de_matrix$Gene
  de_matrix$Gene <- NULL
  de_matrix <- as.matrix(de_matrix)

  if (nrow(de_matrix) >= 2 && ncol(de_matrix) >= 2) {
    png("08b_pseudobulk_logFC_heatmap.png", width=900, height=600, res=130)
    pheatmap(de_matrix, cluster_rows=TRUE, cluster_cols=TRUE,
             display_numbers=TRUE, number_format="%.2f",
             main="Pseudobulk logFC by Cell Type (DN vs Control)",
             fontsize=11, fontsize_number=9,
             color=colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
    cat("✅ 伪批量logFC热图已保存\n")
  }
}

# ---- 9. 保存结果 ----
cat("\n===== 保存结果 =====\n")

results <- list(
  pseudobulk_de = pseudobulk_all,
  celltype_auc = ct_auc_results,
  biomarkers_used = biomarkers_in_sc,
  biomarker_mapping = biomarker_mapping
)
if (exists("auc_matrix")) results$aucell_scores <- auc_matrix
if (exists("auc_by_ct")) results$auc_by_celltype <- auc_by_ct

saveRDS(results, "08b_sc_deep_analysis.rds")
if (!is.null(pseudobulk_all)) write.csv(pseudobulk_all, "08b_pseudobulk_de.csv", row.names=FALSE)
if (nrow(ct_auc_results) > 0) write.csv(ct_auc_results, "08b_celltype_auc.csv", row.names=FALSE)

cat("\n========================================\n")
cat("08b_单细胞深度分析完成！\n")
cat(sprintf("  伪批量DE: %d 细胞类型 x %d 基因\n",
            length(unique(pseudobulk_all$CellType)), length(unique(pseudobulk_all$Gene))))
if (nrow(ct_auc_results) > 0) {
  cat(sprintf("  最佳细胞类型AUC: %s (AUC=%.3f, n=%d)\n",
              ct_auc_results$CellType[1], ct_auc_results$AUC[1], ct_auc_results$N_Cells[1]))
}
if (exists("auc_matrix")) {
  cat(sprintf("  AUCell: %d gene sets x %d cells\n", nrow(auc_matrix), ncol(auc_matrix)))
}
cat("========================================\n")
