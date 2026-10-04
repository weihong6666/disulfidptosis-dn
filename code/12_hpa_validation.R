# ============================================================
# HPA蛋白表达验证 Figure
# 目的：可视化7个标志物基因在肾脏组织中的蛋白水平表达
# 数据来源：Human Protein Atlas (proteinatlas.org)
# ============================================================

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# ---- 1. HPA肾脏IHC数据 ----
# 半定量评分: 0=Negative, 1=Low, 2=Medium, 3=High

hpa_data <- data.frame(
  Gene = c("NCKAP1L", "MYH10", "FLNB", "PDLIM1", "THSD7A", "SLC3A2", "IL1B"),
  Glomeruli = c(0, 0.5, 0.5, 1, 1, 0, 0),
  Tubules = c(3, 3, 3, 0, 0, 3, 0.5),
  Endothelium = c(0, 0.5, 2, 0, 0, 0.5, 0),
  check.names = FALSE
)

# 每个基因的蛋白功能注释
hpa_data$Function <- c(
  "WAVE complex\n(actin polymerization)",
  "Non-muscle myosin IIB\n(actin binding)",
  "Filamin B\n(actin crosslinking)",
  "PDZ-LIM adaptor\n(actin cytoskeleton)",
  "Thrombospondin domain\n(podocyte slit diaphragm)",
  "CD98 heavy chain\n(cystine transporter subunit)",
  "Pro-inflammatory\ncytokine (secreted)"
)

# 每个基因的HPA可靠性
hpa_data$Reliability <- c(
  "2 antibodies",
  "3 antibodies",
  "3 antibodies",
  "2 antibodies",
  "1 antibody",
  "2 antibodies",
  "Secreted (IHC N/A)"
)

# ---- 2. 转长格式 ----
hpa_long <- melt(hpa_data[, c("Gene", "Glomeruli", "Tubules", "Endothelium")],
                 id.vars = "Gene", variable.name = "Compartment", value.name = "Score")

hpa_long$Gene <- factor(hpa_long$Gene, levels = rev(hpa_data$Gene))
hpa_long$Compartment <- factor(hpa_long$Compartment,
                                levels = c("Tubules", "Glomeruli", "Endothelium"))

# ---- 3. 主图：热图 ----
p1 <- ggplot(hpa_long, aes(x = Compartment, y = Gene, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +
  scale_fill_gradientn(
    colors = c("#f1f5f9", "#93c5fd", "#3b82f6", "#1d4ed8"),
    values = c(0, 0.33, 0.66, 1),
    limits = c(0, 3),
    breaks = c(0, 1, 2, 3),
    labels = c("Negative", "Low", "Medium", "High"),
    name = "IHC Staining\nIntensity"
  ) +
  geom_text(aes(label = ifelse(Score == 0, "−",
                                ifelse(Score == 0.5, "±",
                                       ifelse(Score == 1, "+",
                                              ifelse(Score == 2, "++", "+++")))),
                color = ifelse(Score >= 2, "white", "#334155")),
            size = 4.5, fontface = "bold") +
  scale_color_identity() +
  labs(
    title = "Protein-level Validation: IHC Staining in Kidney Tissue",
    subtitle = "Human Protein Atlas (proteinatlas.org) — Normal Human Kidney",
    x = "Renal Compartment",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#1e293b"),
    plot.subtitle = element_text(size = 10, color = "#64748b"),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 11, color = "#475569"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

# ---- 4. 补充图：基因功能条形图 ----
hpa_data$Gene_f <- factor(hpa_data$Gene, levels = rev(hpa_data$Gene))

p2 <- ggplot(hpa_data, aes(x = Gene_f, y = 1, fill = Reliability)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.85, height = 0.85) +
  scale_fill_manual(
    values = c("2 antibodies" = "#22c55e", "3 antibodies" = "#16a34a",
               "1 antibody" = "#fbbf24", "Secreted (IHC N/A)" = "#94a3b8"),
    name = "HPA Validation\nReliability"
  ) +
  geom_text(aes(label = Reliability), size = 3.2, color = "#1e293b") +
  labs(y = "", x = "") +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0))

# ---- 5. 保存 ----
ggsave("HPA_validation_heatmap.png", p1, width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: HPA_validation_heatmap.png\n")

# ---- 6. 打印汇总表 ----
cat("\n========== HPA肾脏IHC汇总 ==========\n")
cat(sprintf("%-10s %-12s %-12s %-15s %s\n",
            "Gene", "Glomeruli", "Tubules", "Endothelium", "Reliability"))
cat(strrep("-", 70), "\n")
for (i in 1:nrow(hpa_data)) {
  cat(sprintf("%-10s %-12s %-12s %-15s %s\n",
              hpa_data$Gene[i],
              ifelse(hpa_data$Glomeruli[i]==0, "Negative",
                     ifelse(hpa_data$Glomeruli[i]<=1, "Low", "Medium/High")),
              ifelse(hpa_data$Tubules[i]==0, "Negative",
                     ifelse(hpa_data$Tubules[i]<=1, "Low",
                            ifelse(hpa_data$Tubules[i]==2, "Medium", "High"))),
              ifelse(hpa_data$Endothelium[i]==0, "Negative",
                     ifelse(hpa_data$Endothelium[i]<=1, "Low", "Medium/High")),
              hpa_data$Reliability[i]))
}

cat("\n✅ 关键发现:\n")
cat("  - MYH10, FLNB, NCKAP1L, SLC3A2 在肾小管中强表达 (actin骨架+胱氨酸转运)\n")
cat("  - THSD7A, PDLIM1 主要在肾小球/足细胞水平表达\n")
cat("  - IL1B为分泌蛋白，正常肾脏IHC不适用；DN中由炎症细胞分泌\n")
cat("  - 多个基因与actin细胞骨架相关 → 支持双硫死亡机制假说\n")
