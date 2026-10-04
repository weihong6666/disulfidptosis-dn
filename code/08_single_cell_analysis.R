# ============================================================
# 双硫死亡+DN 竞赛项目 — 08_单细胞分析（RDS版）
# 数据集：GSE131882 (3 DN + 3 Control 肾组织)
# 方法：Seurat + 双硫死亡评分 + CellChat + monocle3
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)

set.seed(42)

# ---- 1. 读取单细胞RDS数据 ----
cat("\n===== 读取GSE131882单细胞数据 =====\n")

data_dir <- "GSE131882_data"
# 读取未压缩的.rds文件（忽略.gz，已在外部用gunzip解压）
all_rds <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files <- all_rds[!grepl("\\.rds\\.gz$", all_rds)]
cat(sprintf("找到 %d 个RDS文件\n", length(rds_files)))
print(basename(rds_files))

# 识别样本名和分组
sample_names <- gsub("\\.rds$", "", basename(rds_files))
sample_group <- ifelse(grepl("diabet", sample_names, ignore.case = TRUE),
                       "DN", "Control")
cat(sprintf("样本分组: %s\n", paste(sprintf("%s=%s", sample_names, sample_group), collapse=", ")))

# ---- 2. 读取每个样本并创建Seurat对象列表 ----
cat("\n===== 读取各样本数据 =====\n")

seurat_list <- list()

for (i in seq_along(rds_files)) {
  cat(sprintf("\n[%d/%d] %s (%s)\n", i, length(rds_files),
              sample_names[i], sample_group[i]))

          cat("  读取中...
")
          raw_data <- readRDS(rds_files[i])
  cat(sprintf("  数据大小: %.1f MB\n", as.numeric(object.size(raw_data)) / 1e6))

  cat(sprintf("  原始对象类型: %s\n", class(raw_data)[1]))

  # 自动检测数据结构并提取计数矩阵
  counts <- NULL

          # 递归查找矩阵（自动钻入多层嵌套，如 umicount$exon$all）
          find_matrix <- function(obj, path = "", depth = 0) {
            if (depth > 5) return(list(mat = NULL, path = path))
            if (inherits(obj, "dgCMatrix") || inherits(obj, "dgeMatrix") ||
                inherits(obj, "Matrix") || inherits(obj, "matrix")) {
              return(list(mat = obj, path = path))
            }
            if (is.list(obj) && !is.data.frame(obj)) {
              for (nm in names(obj)) {
                result <- find_matrix(obj[[nm]], paste0(path, "$", nm), depth + 1)
                if (!is.null(result$mat)) return(result)
              }
            }
            return(list(mat = NULL, path = path))
          }

          if (inherits(raw_data, "dgCMatrix") || inherits(raw_data, "matrix") ||
              inherits(raw_data, "dgeMatrix")) {
            counts <- raw_data
            cat("  直接是矩阵格式\n")
          } else if (is.list(raw_data)) {
            cat(sprintf("  列表元素: %s\n", paste(names(raw_data), collapse=", ")))
            result <- find_matrix(raw_data)
            if (!is.null(result$mat)) {
              counts <- result$mat
              cat(sprintf("  找到计数矩阵: %s (%d x %d)\n",
                          result$path, nrow(counts), ncol(counts)))
            }
          }

          if (is.null(counts)) {
            cat("  ⚠ 无法提取计数矩阵，跳过此样本\n")
            next
          }
  colnames(counts) <- paste0(sample_names[i], "_", colnames(counts))
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_group[i],
    min.cells = 3,
    min.features = 200
  )
  seurat_obj$sample <- sample_names[i]
  seurat_obj$group <- sample_group[i]
  seurat_obj$orig.ident <- sample_names[i]

  cat(sprintf("  ✅ Seurat对象: %d 基因 x %d 细胞\n", nrow(seurat_obj), ncol(seurat_obj)))
  seurat_list[[sample_names[i]]] <- seurat_obj
}

if (length(seurat_list) == 0) stop("没有成功读取任何样本！")

cat(sprintf("\n共读取 %d 个样本\n", length(seurat_list)))

# ---- 3. 合并所有样本 ----
cat("\n===== 合并样本 =====\n")

if (length(seurat_list) == 1) {
  sc_obj <- seurat_list[[1]]
} else {
  sc_obj <- merge(seurat_list[[1]],
                  y = seurat_list[-1],
                  add.cell.ids = names(seurat_list))
}

cat(sprintf("合并后: %d 基因 x %d 细胞\n", nrow(sc_obj), ncol(sc_obj)))

# ---- 4. 质控 ----
cat("\n===== 质控过滤 =====\n")

sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj[["percent.ribo"]] <- PercentageFeatureSet(sc_obj, pattern = "^RP[SL]")

# 质控小提琴图
VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        group.by = "group", ncol = 4, pt.size = 0.1)
ggsave("08_QC_violin.png", width = 16, height = 5, dpi = 150)

# 过滤
sc_obj <- subset(sc_obj,
                 nFeature_RNA > 200 & nFeature_RNA < 6000 &
                 nCount_RNA > 500 & percent.mt < 25)
cat(sprintf("QC后: %d 基因 x %d 细胞\n", nrow(sc_obj), ncol(sc_obj)))

# ---- 5. 标准化和降维 ----
cat("\n===== 标准化+降维 =====\n")

sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize",
                        scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj, vars.to.regress = c("percent.mt", "nCount_RNA"))
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(sc_obj))

ElbowPlot(sc_obj, ndims = 50)
ggsave("08_elbow.png", width = 8, height = 5, dpi = 150)

# 选择PCA维度
ndim <- min(30, max(which(sc_obj@reductions$pca@stdev > 1)))
cat(sprintf("使用 %d 个PCA维度\n", ndim))

sc_obj <- FindNeighbors(sc_obj, dims = 1:ndim)
sc_obj <- FindClusters(sc_obj, resolution = 0.4)
sc_obj <- RunUMAP(sc_obj, dims = 1:ndim)

# UMAP - 按聚类
DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  labs(title = sprintf("DN Kidney scRNA-seq (%d clusters)", length(unique(Idents(sc_obj)))))
ggsave("08_UMAP_clusters.png", width = 9, height = 7, dpi = 150)

# UMAP - 按分组
DimPlot(sc_obj, reduction = "umap", group.by = "group", pt.size = 0.5,
        cols = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "UMAP by Condition")
ggsave("08_UMAP_group.png", width = 9, height = 7, dpi = 150)

# UMAP - 按样本
DimPlot(sc_obj, reduction = "umap", group.by = "sample", pt.size = 0.3) +
  labs(title = "UMAP by Sample")
ggsave("08_UMAP_sample.png", width = 10, height = 7, dpi = 150)
cat("✅ UMAP已保存\n")


# ---- 6. 细胞注释（Ensembl ID版本） ----
cat("\n===== 细胞注释 =====\n")

kidney_markers <- list(
  "Podocyte" = c("ENSG00000161270", "ENSG00000116218", "ENSG00000184937",
                 "ENSG00000171992", "ENSG00000128567"),
  "Proximal_Tubule" = c("ENSG00000131183", "ENSG00000081479", "ENSG00000107611",
                         "ENSG00000197901", "ENSG00000140675"),
  "Loop_of_Henle" = c("ENSG00000169344", "ENSG00000074803", "ENSG00000113946"),
  "Distal_Tubule" = c("ENSG00000070915", "ENSG00000134873", "ENSG00000104327"),
  "Collecting_Duct" = c("ENSG00000167580", "ENSG00000126895", "ENSG00000111319",
                         "ENSG00000004939"),
  "Mesangial" = c("ENSG00000113721", "ENSG00000077943", "ENSG00000107485"),
  "Endothelial" = c("ENSG00000261371", "ENSG00000174059", "ENSG00000110799",
                     "ENSG00000102755", "ENSG00000106991"),
  "Fibroblast" = c("ENSG00000108821", "ENSG00000011465", "ENSG00000139329",
                    "ENSG00000134853"),
  "Immune" = c("ENSG00000081237", "ENSG00000167286", "ENSG00000170458",
                "ENSG00000129226", "ENSG00000156738"),
  "Macula_Densa" = c("ENSG00000089250", "ENSG00000073756")
)

for (ct in names(kidney_markers)) {
  genes_avail <- intersect(kidney_markers[[ct]], rownames(sc_obj))
  if (length(genes_avail) >= 2) {
    sc_obj <- AddModuleScore(sc_obj, features = list(genes_avail),
                             name = paste0("score_", ct))
  }
}

score_cols <- grep("^score_", colnames(sc_obj@meta.data), value = TRUE)
if (length(score_cols) >= 4) {
  FeaturePlot(sc_obj, features = head(score_cols, 4), ncol = 2,
              pt.size = 0.3, order = TRUE)
  ggsave("08_celltype_scores.png", width = 12, height = 10, dpi = 150)
  cat("Cell type score plots saved\n")
}

if (length(score_cols) > 0) {
  score_mat <- sc_obj@meta.data[, score_cols, drop = FALSE]
  assigned <- apply(score_mat, 1, function(r) {
    if (max(r) > 0.1) return(gsub("^score_", "", names(which.max(r))))
    return("Unassigned")
  })
  sc_obj$cell_type <- assigned
  cat(sprintf("Cell types: %s\n", paste(names(table(sc_obj$cell_type)), collapse = ", ")))
  DimPlot(sc_obj, reduction = "umap", group.by = "cell_type",
          label = TRUE, pt.size = 0.3, repel = TRUE) +
    labs(title = "Cell Type Annotation")
  ggsave("08_celltype_umap.png", width = 11, height = 8, dpi = 150)
  cat("Cell type UMAP saved\n")
}

# ---- 7. 双硫死亡基因评分（Ensembl ID版本） ----
cat("\n===== Disulfidptosis Score =====\n")

drg_ensembl <- c(
  "ENSG00000151012", "ENSG00000168003", "ENSG00000104812", "ENSG00000138095",
  "ENSG00000174886", "ENSG00000023228", "ENSG00000151413", "ENSG00000151093",
  "ENSG00000163902", "ENSG00000061676", "ENSG00000196924", "ENSG00000136068",
  "ENSG00000100345", "ENSG00000133026", "ENSG00000137076", "ENSG00000075624",
  "ENSG00000092841", "ENSG00000077549", "ENSG00000125868", "ENSG00000140575",
  "ENSG00000130402"
)

biomarker_ensembl <- c(
  NCKAP1L = "ENSG00000123338", THSD7A = "ENSG00000005108",
  MYH10 = "ENSG00000133026", PDLIM1 = "ENSG00000107438",
  IL1B = "ENSG00000125538", FLNB = "ENSG00000136068",
  SLC3A2 = "ENSG00000168003"
)

all_score_ens <- unique(c(drg_ensembl, biomarker_ensembl))
drg_avail <- intersect(all_score_ens, rownames(sc_obj))
cat(sprintf("Available disulfidptosis/biomarker genes: %d/%d\n",
            length(drg_avail), length(all_score_ens)))

if (length(drg_avail) >= 5) {
  sc_obj <- AddModuleScore(sc_obj, features = list(drg_avail),
                           name = "Disulfidptosis_Score")

  FeaturePlot(sc_obj, features = "Disulfidptosis_Score1",
              pt.size = 0.5, order = TRUE) +
    scale_color_viridis_c(option = "inferno") +
    labs(title = "Disulfidptosis + Biomarker Score in DN Kidney")
  ggsave("08_disulfidptosis_UMAP.png", width = 10, height = 8, dpi = 150)
  cat("Disulfidptosis UMAP saved\n")

  VlnPlot(sc_obj, features = "Disulfidptosis_Score1",
          group.by = "group", pt.size = 0.1,
          cols = c("Control" = "steelblue", "DN" = "tomato")) +
    labs(title = "Disulfidptosis Score: DN vs Control", y = "Score")
  ggsave("08_disulfidptosis_violin.png", width = 8, height = 5, dpi = 150)
  cat("Disulfidptosis violin saved\n")
}

biomarker_avail <- intersect(biomarker_ensembl, rownames(sc_obj))
if (length(biomarker_avail) >= 3) {
  FeaturePlot(sc_obj, features = head(biomarker_avail, 4),
              ncol = 2, pt.size = 0.3, order = TRUE)
  ggsave("08_biomarker_features.png", width = 12, height = 10, dpi = 150)
  cat("Biomarker feature plots saved\n")

  tryCatch({
    VlnPlot(sc_obj, features = head(biomarker_avail, 6),
            group.by = "group", pt.size = 0, ncol = 3,
            cols = c("Control" = "steelblue", "DN" = "tomato"))
    ggsave("08_biomarker_violin.png", width = 14, height = 8, dpi = 150)
    cat("Biomarker violin plots saved\n")
  }, error = function(e) {
    cat("Biomarker violin skipped:", conditionMessage(e), "\n")
  })
}

# ---- 8. DN vs Control差异 ----
cat("\n===== DN vs Control =====\n")

cluster_comp <- table(sc_obj$seurat_clusters, sc_obj$group)
cluster_prop <- prop.table(cluster_comp, margin = 2)
print(round(cluster_prop, 3))

cluster_counts <- as.data.frame(table(sc_obj$seurat_clusters, sc_obj$group))
colnames(cluster_counts) <- c("Cluster", "Group", "Count")

ggplot(cluster_counts, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Control" = "steelblue", "DN" = "tomato")) +
  labs(title = "Cell Count per Cluster: DN vs Control", y = "Number of Cells") +
  theme_minimal()
ggsave("08_cluster_composition.png", width = 12, height = 5, dpi = 150)
cat("Cluster composition saved\n")

Idents(sc_obj) <- "seurat_clusters"
dn_markers_all <- tryCatch({
  FindMarkers(sc_obj, ident.1 = "DN", group.by = "group",
              min.pct = 0.1, logfc.threshold = 0.25)
}, error = function(e) NULL)

if (!is.null(dn_markers_all) && nrow(dn_markers_all) > 0) {
  write.csv(dn_markers_all, "08_DN_vs_Ctrl_DEGs.csv")
  cat(sprintf("DN vs Control DEGs: %d\n", nrow(dn_markers_all)))
}

# ---- 9. 拟时序分析 (monocle3) ----
cat("\n===== Pseudotime (monocle3) =====\n")

if (requireNamespace("monocle3", quietly = TRUE)) {
  library(monocle3)
  tryCatch({
    cds <- as.cell_data_set(sc_obj)
    cds <- cluster_cells(cds, reduction_method = "UMAP")
    cds <- learn_graph(cds)
    cds <- order_cells(cds)

    plot_cells(cds, color_cells_by = "pseudotime",
               label_cell_groups = FALSE, label_leaves = FALSE,
               label_branch_points = FALSE, cell_size = 1)
    ggsave("08_pseudotime.png", width = 10, height = 8, dpi = 150)

    plot_cells(cds, color_cells_by = "group",
               label_cell_groups = FALSE, cell_size = 1)
    ggsave("08_pseudotime_group.png", width = 10, height = 8, dpi = 150)

    pseudotime_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 2)
    pseudotime_sig <- subset(pseudotime_genes, q_value < 0.05)
    cat(sprintf("Pseudotime significant genes: %d\n", nrow(pseudotime_sig)))

    drg_pseudo <- intersect(drg_avail, rownames(pseudotime_sig))
    if (length(drg_pseudo) > 0) {
      cat(sprintf("DRG genes significant in pseudotime: %s\n",
                  paste(drg_pseudo, collapse = ", ")))
    }

    saveRDS(cds, "08_monocle_cds.rds")
    cat("monocle3 pseudotime completed\n")
  }, error = function(e) {
    cat("monocle3 error:", conditionMessage(e), "\n")
  })
} else {
  cat("monocle3 not installed, skipping\n")
}

# ---- 10. 保存 ----
cat("\n===== Saving =====\n")
saveRDS(sc_obj, "08_seurat_object.rds")

cat("\n========================================\n")
cat("08 single-cell analysis complete!\n")
cat("========================================\n")
