# 加载必要的包
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

# 设置随机种子以确保可重复性
set.seed(123)

# 创建模拟数据 ------------------------------------------------------------

# 定义细胞类型和基因组
cell_types <- c("T_NK cell", "Myeloid cell", "B cell", 
                "Endothelial cell", "Mesenchymal cell", "HSC", "Malignant cell")

gene_groups <- list(
  group1 = c("CCL5", "GZMB", "GNLY", "NKG7", "GZMA", "GZMB"),
  group2 = c("CD74", "CST3", "TXN", "S100A8", "S100A9", "CCL4"),
  group3 = c("JCHAIN", "IGHA2", "MZB1", "IGHG2", "CD79A", "IGHA1"),
  group4 = c("CLDN5", "A2M", "GBP4", "FCN3", "DEPP1", "GNG11"),
  group5 = c("MGP", "ACTA2", "RGS5", "MYL9", "IGFBP7", "TAGLN"),
  group6 = c("KRT18", "KRT8", "SCGB3A1", "CLDN4", "KRT19", "ANXA4"),
  group7 = c("APOH", "APOC3", "APOA1", "HP", "ALB", "TTR")
)

# 设置每组基因的颜色
gene_colors <- c(
  "#F9A23D", "#F47C40", "#E7A56F", "#8CC37A",
  "#79BC5C", "#F7968C", "#ED6564"
)

# 创建GO富集条目
go_terms <- list(
  group1 = c("Positive regulation of lymphocyte chemotaxis", 
             "Natural killer cell chemotaxis", 
             "Positive regulation of T cell chemotaxis"),
  group2 = c("Antigen processing and presentation", 
             "Interferon-gamma-mediated signaling pathway", 
             "Neutrophil degranulation"),
  group3 = c("Regulation of B cell activation", 
             "B cell receptor signaling pathway", 
             "B cell activation"),
  group4 = c("Positive regulation of protein kinase A signaling", 
             "Endothelial cell activation", 
             "Endothelial cell differentiation"),
  group5 = c("Muscle contraction", 
             "Actin-myosin filament sliding", 
             "Myofibril assembly"),
  group6 = c("Cornification", 
             "Hepatocyte apoptotic process", 
             "Skin development"),
  group7 = c("Platelet degranulation", 
             "High-density lipoprotein particle remodeling", 
             "Protein activation cascade")
)

# 创建表达数据矩阵
n_genes <- sum(sapply(gene_groups, length))
n_cells <- length(cell_types)

expr_matrix <- matrix(rnorm(n_genes * n_cells, mean = 0, sd = 0.5), nrow = n_genes, ncol = n_cells)
rownames(expr_matrix) <- unlist(gene_groups)
colnames(expr_matrix) <- cell_types

# 设置特定的表达模式
expr_matrix[gene_groups$group1, "T_NK cell"] <- rnorm(length(gene_groups$group1), mean = 2, sd = 0.3)
expr_matrix[gene_groups$group2, "Myeloid cell"] <- rnorm(length(gene_groups$group2), mean = 2, sd = 0.3)
expr_matrix[gene_groups$group3, "B cell"] <- rnorm(length(gene_groups$group3), mean = 2, sd = 0.3)
expr_matrix[gene_groups$group4, "Endothelial cell"] <- rnorm(length(gene_groups$group4), mean = 1.5, sd = 0.3)
expr_matrix[gene_groups$group5, "Mesenchymal cell"] <- rnorm(length(gene_groups$group5), mean = 1.7, sd = 0.3)
expr_matrix[gene_groups$group6, "HSC"] <- rnorm(length(gene_groups$group6), mean = 1.2, sd = 0.3)
expr_matrix[gene_groups$group7, "Malignant cell"] <- rnorm(length(gene_groups$group7), mean = 2, sd = 0.3)

# 创建GO富集p值数据
go_pvalues <- data.frame(
  GO_term = unlist(go_terms),
  log10_adj_pvalue = runif(sum(sapply(go_terms, length)), min = 3, max = 15),
  gene_group = rep(names(go_terms), sapply(go_terms, length))
)

# 调整显著p值
go_pvalues$log10_adj_pvalue[go_pvalues$gene_group == "group3"] <- runif(3, min = 12, max = 15)
go_pvalues$log10_adj_pvalue[go_pvalues$gene_group == "group7"] <- runif(3, min = 11, max = 14)
go_pvalues$log10_adj_pvalue[go_pvalues$gene_group == "group2"] <- runif(3, min = 9, max = 12)
go_pvalues$log10_adj_pvalue[go_pvalues$gene_group == "group5"] <- runif(3, min = 5, max = 8)

# 为可视化准备数据 ------------------------------------------------------

# 行分组和颜色
row_groups <- rep(names(gene_groups), sapply(gene_groups, length))
names(row_groups) <- unlist(gene_groups)
row_colors <- rep(gene_colors, sapply(gene_groups, length))

# 创建热图注释
gene_group_anno <- rowAnnotation(
  foo = anno_block(
    gp = gpar(fill = gene_colors),
    labels = names(gene_groups),
    labels_gp = gpar(col = "black", fontsize = 8)
  )
)

# 计算热图高度
n_total_genes <- sum(sapply(gene_groups, length))
heatmap_height <- unit(n_total_genes * 4.5, "mm")  # 关键调整参数

# 创建热图对象
heatmap_obj <- Heatmap(
  expr_matrix,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(row_groups, levels = names(gene_groups)),
  row_gap = unit(0, "mm"),
  show_row_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(col = row_colors, fontsize = 8),
  column_names_gp = gpar(
    col = c("#F8CB58", "#B65C22", "#8560A8", "#00A896", "#4A8F3D", "#F7A3A6", "#DE281B"),
    fontsize = 8
  ),
  height = heatmap_height,
  left_annotation = gene_group_anno
)

# 创建GO富集柱状图 --------------------------------------------------------
go_pvalues$gene_group_factor <- factor(go_pvalues$gene_group, levels = names(gene_groups))
go_pvalues$GO_term_factor <- factor(go_pvalues$GO_term, levels = rev(unlist(go_terms)))

go_barplot <- ggplot(go_pvalues, aes(x = log10_adj_pvalue, y = GO_term_factor, fill = gene_group_factor)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = gene_colors) +
  scale_x_continuous(limits = c(0, 16), breaks = seq(0, 15, 5)) +
  coord_fixed(ratio = 2.8) +  # 关键调整参数
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8, margin = margin(r = -10)),
    plot.margin = margin(l = 0, r = 0),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# 创建图例 --------------------------------------------------------------
count_legend <- ggplot(data.frame(x=1,y=1)) + 
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FF7F50", "#FF4500"),
    breaks = c(2.5, 7.5, 12.5),
    limits = c(2, 13)
  ) +
  theme_void() +
  theme(legend.key.size = unit(4, "mm"))

expr_legend <- ggplot(data.frame(x=1,y=1)) +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    breaks = c(-2, 0, 2)
  ) +
  theme_void() +
  theme(legend.key.size = unit(4, "mm"))

# 组合图形 --------------------------------------------------------------
heatmap_plot <- grid.grabExpr(draw(heatmap_obj))
count_leg <- get_legend(count_legend)
expr_leg <- get_legend(expr_legend)

combined_plot <- plot_grid(
  plot_grid(
    heatmap_plot,
    go_barplot,
    count_leg,
    expr_leg,
    ncol = 4,
    rel_widths = c(4.2, 3, 0.5, 0.5),  # 关键宽度比例
    align = "h",
    axis = "tb"
  ),
  nrow = 1
)

ggsave("final_plot.pdf", combined_plot, width = 12, height = n_total_genes * 0.25 + 1)  # 动态调整高度
print(combined_plot)



