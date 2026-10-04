# ============================================================
# 双硫死亡+DN — 11_GO_KEGG富集分析（补充分析）
# 目的：对DR-DEGs和最终标志物进行功能和通路富集
# 补充当前管线的缺失环节
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(DOSE)
library(httr)

# 增加KEGG API超时时间（国内访问KEGG较慢）
httr::set_config(httr::timeout(300))

set.seed(42)

# ---- 1. 加载基因列表 ----
cat("\n===== 加载基因 =====\n")

# 9个DR-DEGs
dr_degs <- c("SLC3A2", "IL1B", "FLNB", "COL6A3", "MYH10",
             "PDLIM1", "MAGI2", "NCKAP1L", "THSD7A")

# 7个最终标志物
final_biomarkers <- c("NCKAP1L", "THSD7A", "MYH10", "PDLIM1",
                      "IL1B", "FLNB", "SLC3A2")

# 差异分析中上调 vs 下调的基因
up_genes <- c("FLNB", "COL6A3", "MYH10", "PDLIM1", "NCKAP1L")
down_genes <- c("SLC3A2", "IL1B", "MAGI2", "THSD7A")

cat(sprintf("DR-DEGs: %d genes\n", length(dr_degs)))
cat(sprintf("Up-regulated in DN: %s\n", paste(up_genes, collapse=", ")))
cat(sprintf("Down-regulated in DN: %s\n", paste(down_genes, collapse=", ")))

# ---- 2. 基因Symbol → Entrez ID ----
dr_degs_entrez <- tryCatch({
  bitr(dr_degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("Entrez ID转换失败，尝试其他方式...\n")
  NULL
})

final_entrez <- tryCatch({
  bitr(final_biomarkers, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) NULL)

if (!is.null(dr_degs_entrez)) {
  cat(sprintf("成功转换 %d/%d DR-DEGs 为 Entrez ID\n",
              nrow(dr_degs_entrez), length(dr_degs)))
}

# ---- 3. GO富集分析 ----
cat("\n===== GO富集分析 =====\n")

if (!is.null(dr_degs_entrez)) {
  # BP (Biological Process)
  go_bp <- enrichGO(
    gene = dr_degs_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  cat(sprintf("GO BP: %d enriched terms\n", ifelse(is.null(go_bp), 0, nrow(go_bp))))

  # CC (Cellular Component)
  go_cc <- enrichGO(
    gene = dr_degs_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  cat(sprintf("GO CC: %d enriched terms\n", ifelse(is.null(go_cc), 0, nrow(go_cc))))

  # MF (Molecular Function)
  go_mf <- enrichGO(
    gene = dr_degs_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  cat(sprintf("GO MF: %d enriched terms\n", ifelse(is.null(go_mf), 0, nrow(go_mf))))
}

# ---- 4. KEGG通路富集 ----
cat("\n===== KEGG富集分析 =====\n")

if (!is.null(dr_degs_entrez)) {
  kegg <- tryCatch({
    enrichKEGG(
      gene = dr_degs_entrez$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat("KEGG enrichment failed (API timeout/network):", conditionMessage(e), "\n")
    cat("Will try alternative: gseKEGG or skip KEGG for now\n")
    return(NULL)
  })
  if (!is.null(kegg) && nrow(kegg) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  cat(sprintf("KEGG: %d enriched pathways\n", ifelse(is.null(kegg), 0, nrow(kegg))))
}

# ---- 5. 可视化 ----
cat("\n===== 生成富集分析图 =====\n")

# 5a. GO BP dotplot
if (exists("go_bp") && !is.null(go_bp) && nrow(go_bp) > 0) {
  p1 <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("11_GO_BP_dotplot.png", p1, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: 11_GO_BP_dotplot.png\n")
}

# 5b. GO CC dotplot
if (exists("go_cc") && !is.null(go_cc) && nrow(go_cc) > 0) {
  p2 <- dotplot(go_cc, showCategory = 10, title = "GO Cellular Component") +
    theme_minimal(base_size = 11)
  ggsave("11_GO_CC_dotplot.png", p2, width = 9, height = 5, dpi = 300, bg = "white")
  cat("Saved: 11_GO_CC_dotplot.png\n")
}

# 5c. GO MF dotplot
if (exists("go_mf") && !is.null(go_mf) && nrow(go_mf) > 0) {
  p3 <- dotplot(go_mf, showCategory = 10, title = "GO Molecular Function") +
    theme_minimal(base_size = 11)
  ggsave("11_GO_MF_dotplot.png", p3, width = 9, height = 5, dpi = 300, bg = "white")
  cat("Saved: 11_GO_MF_dotplot.png\n")
}

# 5d. KEGG dotplot
if (exists("kegg") && !is.null(kegg) && nrow(kegg) > 0) {
  p4 <- dotplot(kegg, showCategory = 15, title = "KEGG Pathways") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("11_KEGG_dotplot.png", p4, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: 11_KEGG_dotplot.png\n")
}

# 5e. GO多维度组合图（Top5 BP + Top5 CC + Top5 MF）
if (exists("go_bp") && !is.null(go_bp) && nrow(go_bp) >= 3 &&
    exists("go_cc") && !is.null(go_cc) && nrow(go_cc) >= 3 &&
    exists("go_mf") && !is.null(go_mf) && nrow(go_mf) >= 3) {

  # 简化：选取top terms并合并
  go_bp_top <- go_bp[1:min(5, nrow(go_bp)), ]
  go_bp_top$Category <- "BP"

  go_cc_top <- go_cc[1:min(5, nrow(go_cc)), ]
  go_cc_top$Category <- "CC"

  go_mf_top <- go_mf[1:min(5, nrow(go_mf)), ]
  go_mf_top$Category <- "MF"

  go_combined <- rbind(
    go_bp_top[, c("Description", "p.adjust", "Count", "Category")],
    go_cc_top[, c("Description", "p.adjust", "Count", "Category")],
    go_mf_top[, c("Description", "p.adjust", "Count", "Category")]
  )

  p5 <- ggplot(go_combined,
               aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)),
                   size = Count, color = Category)) +
    geom_point() +
    facet_wrap(~Category, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("BP" = "#3b82f6", "CC" = "#f59e0b", "MF" = "#10b981")) +
    labs(
      title = "GO Enrichment Analysis of Disulfidptosis-Related DEGs in DN",
      x = expression(-log[10](P.adj)),
      y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 12))

  ggsave("11_GO_combined.png", p5, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved: 11_GO_combined.png\n")
}

# ---- 6. 保存结果表 ----
cat("\n===== 保存富集结果 =====\n")

enrichment_results <- list()

if (exists("go_bp") && !is.null(go_bp)) enrichment_results$GO_BP <- as.data.frame(go_bp)
if (exists("go_cc") && !is.null(go_cc)) enrichment_results$GO_CC <- as.data.frame(go_cc)
if (exists("go_mf") && !is.null(go_mf)) enrichment_results$GO_MF <- as.data.frame(go_mf)
if (exists("kegg") && !is.null(kegg)) enrichment_results$KEGG <- as.data.frame(kegg)

saveRDS(enrichment_results, "11_enrichment_results.rds")

# 输出核心结果到CSV
if (!is.null(enrichment_results$GO_BP) && nrow(enrichment_results$GO_BP) > 0) {
  write.csv(enrichment_results$GO_BP[1:min(20, nrow(enrichment_results$GO_BP)),
            c("Description", "p.adjust", "Count", "geneID")],
            "11_GO_BP_top20.csv", row.names = FALSE)
}
if (!is.null(enrichment_results$KEGG) && nrow(enrichment_results$KEGG) > 0) {
  write.csv(enrichment_results$KEGG[, c("Description", "p.adjust", "Count", "geneID")],
            "11_KEGG_results.csv", row.names = FALSE)
}

cat("\n===== 富集分析完成 =====\n")

# ---- 7. 打印核心结果（供论文使用） ----
cat("\n========== 论文可引用的关键富集结果 ==========\n")

if (exists("go_bp") && !is.null(go_bp) && nrow(go_bp) >= 3) {
  cat("\nTop 5 GO Biological Process:\n")
  for (i in 1:min(5, nrow(go_bp))) {
    cat(sprintf("  %d. %s (P.adj=%.2e, %d genes)\n",
                i, go_bp$Description[i], go_bp$p.adjust[i], go_bp$Count[i]))
  }
}

if (exists("go_cc") && !is.null(go_cc) && nrow(go_cc) >= 3) {
  cat("\nTop 5 GO Cellular Component:\n")
  for (i in 1:min(5, nrow(go_cc))) {
    cat(sprintf("  %d. %s (P.adj=%.2e, %d genes)\n",
                i, go_cc$Description[i], go_cc$p.adjust[i], go_cc$Count[i]))
  }
}

if (exists("kegg") && !is.null(kegg) && nrow(kegg) >= 3) {
  cat("\nTop 5 KEGG Pathways:\n")
  for (i in 1:min(5, nrow(kegg))) {
    cat(sprintf("  %d. %s (P.adj=%.2e, %d genes)\n",
                i, kegg$Description[i], kegg$p.adjust[i], kegg$Count[i]))
  }
}

# ============================================================
# ---- 8. 全部显著DEGs的GO+KEGG富集（补充分析） ----
# 目的：不仅限于DR-DEGs，对所有显著DEGs做富集，
#       为论文提供更全面的生物学上下文
# ============================================================
cat("\n\n===== 全部显著DEGs富集分析 =====\n")

# 从CSV加载全部DEG结果
deg_all <- read.csv("02_dr_degs.csv", row.names = 1)
sig_deg_all <- deg_all[deg_all$adj.P.Val < 0.05, ]
all_sig_genes <- sig_deg_all$Gene
cat(sprintf("全部显著DEGs: %d genes\n", length(all_sig_genes)))

# Entrez ID转换
all_entrez <- tryCatch({
  bitr(all_sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("Entrez ID转换失败\n")
  NULL
})

if (!is.null(all_entrez) && nrow(all_entrez) >= 5) {
  cat(sprintf("成功转换 %d/%d 全部DEGs\n", nrow(all_entrez), length(all_sig_genes)))

  # GO BP
  all_go_bp <- enrichGO(gene = all_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
    qvalueCutoff = 0.2, readable = TRUE)
  cat(sprintf("All-DEGs GO BP: %d terms\n", ifelse(is.null(all_go_bp), 0, nrow(all_go_bp))))

  # GO CC
  all_go_cc <- enrichGO(gene = all_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
    ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
    qvalueCutoff = 0.2, readable = TRUE)
  cat(sprintf("All-DEGs GO CC: %d terms\n", ifelse(is.null(all_go_cc), 0, nrow(all_go_cc))))

  # GO MF
  all_go_mf <- enrichGO(gene = all_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
    ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
    qvalueCutoff = 0.2, readable = TRUE)
  cat(sprintf("All-DEGs GO MF: %d terms\n", ifelse(is.null(all_go_mf), 0, nrow(all_go_mf))))

  # KEGG
  all_kegg <- tryCatch({
    enrichKEGG(gene = all_entrez$ENTREZID, organism = "hsa",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) {
    cat("All-DEGs KEGG failed:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(all_kegg) && nrow(all_kegg) > 0) {
    all_kegg <- setReadable(all_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  cat(sprintf("All-DEGs KEGG: %d pathways\n", ifelse(is.null(all_kegg), 0, nrow(all_kegg))))

  # ---- 可视化 ----
  # 8a. GO BP dotplot
  if (!is.null(all_go_bp) && nrow(all_go_bp) >= 5) {
    p_a1 <- dotplot(all_go_bp, showCategory = 15, title = "GO BP (All Significant DEGs)") +
      theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
    ggsave("11_GO_BP_all_DEGs.png", p_a1, width = 11, height = 7, dpi = 300, bg = "white")
    cat("Saved: 11_GO_BP_all_DEGs.png\n")
  }

  # 8b. KEGG dotplot
  if (!is.null(all_kegg) && nrow(all_kegg) >= 3) {
    p_a2 <- dotplot(all_kegg, showCategory = 15, title = "KEGG (All Significant DEGs)") +
      theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
    ggsave("11_KEGG_all_DEGs.png", p_a2, width = 11, height = 7, dpi = 300, bg = "white")
    cat("Saved: 11_KEGG_all_DEGs.png\n")
  }

  # 8c. GO 三合一组合图 (BP+CC+MF top5 each)
  has_bp <- !is.null(all_go_bp) && nrow(all_go_bp) >= 3
  has_cc <- !is.null(all_go_cc) && nrow(all_go_cc) >= 3
  has_mf <- !is.null(all_go_mf) && nrow(all_go_mf) >= 3

  if (has_bp || has_cc || has_mf) {
    go_parts <- list()
    if (has_bp) {
      tmp <- all_go_bp[1:min(5, nrow(all_go_bp)), ]
      tmp$Category <- "BP"; go_parts[[length(go_parts) + 1]] <- tmp
    }
    if (has_cc) {
      tmp <- all_go_cc[1:min(5, nrow(all_go_cc)), ]
      tmp$Category <- "CC"; go_parts[[length(go_parts) + 1]] <- tmp
    }
    if (has_mf) {
      tmp <- all_go_mf[1:min(5, nrow(all_go_mf)), ]
      tmp$Category <- "MF"; go_parts[[length(go_parts) + 1]] <- tmp
    }

    go_all_combined <- do.call(rbind, lapply(go_parts, function(x) {
      x[, c("Description", "p.adjust", "Count", "Category")]
    }))

    p_a3 <- ggplot(go_all_combined,
      aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)),
          size = Count, color = Category)) +
      geom_point() +
      facet_wrap(~Category, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c("BP" = "#3b82f6", "CC" = "#f59e0b", "MF" = "#10b981")) +
      labs(title = "GO Enrichment — All Significant DEGs in DN",
           x = expression(-log[10](P.adj)), y = "") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 12))
    ggsave("11_GO_combined_all_DEGs.png", p_a3, width = 11, height = 9, dpi = 300, bg = "white")
    cat("Saved: 11_GO_combined_all_DEGs.png\n")
  }

  # 保存结果
  all_enrich <- list()
  if (exists("all_go_bp") && !is.null(all_go_bp)) all_enrich$GO_BP <- as.data.frame(all_go_bp)
  if (exists("all_go_cc") && !is.null(all_go_cc)) all_enrich$GO_CC <- as.data.frame(all_go_cc)
  if (exists("all_go_mf") && !is.null(all_go_mf)) all_enrich$GO_MF <- as.data.frame(all_go_mf)
  if (exists("all_kegg") && !is.null(all_kegg)) all_enrich$KEGG <- as.data.frame(all_kegg)
  saveRDS(all_enrich, "11_all_DEGs_enrichment.rds")

  # CSV输出
  if (!is.null(all_go_bp) && nrow(all_go_bp) > 0) {
    write.csv(all_go_bp[1:min(20, nrow(all_go_bp)), c("Description","p.adjust","Count","geneID")],
              "11_GO_BP_all_DEGs_top20.csv", row.names = FALSE)
  }
  if (!is.null(all_kegg) && nrow(all_kegg) > 0) {
    write.csv(all_kegg[, c("Description","p.adjust","Count","geneID")],
              "11_KEGG_all_DEGs.csv", row.names = FALSE)
  }

  # 打印关键结果
  cat("\n----- 全部DEGs Top GO BP -----\n")
  if (has_bp) {
    for (i in 1:min(10, nrow(all_go_bp))) {
      cat(sprintf("  %d. %s (P.adj=%.2e, %d genes)\n",
        i, all_go_bp$Description[i], all_go_bp$p.adjust[i], all_go_bp$Count[i]))
    }
  }
  cat("\n----- 全部DEGs Top KEGG -----\n")
  if (!is.null(all_kegg) && nrow(all_kegg) > 0) {
    for (i in 1:min(10, nrow(all_kegg))) {
      cat(sprintf("  %d. %s (P.adj=%.2e, %d genes)\n",
        i, all_kegg$Description[i], all_kegg$p.adjust[i], all_kegg$Count[i]))
    }
  }
} else {
  cat("全部DEGs Entrez ID转换失败或基因数不足 (<5)，跳过大基因集富集\n")
}

cat("\n===== 全部富集分析完成 =====\n")
