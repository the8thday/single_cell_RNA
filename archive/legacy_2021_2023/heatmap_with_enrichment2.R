library(ComplexHeatmap)
library(circlize)

# 设置随机种子
set.seed(123)

# 定义细胞类型和基因组 ---------------------------------------------------------
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

gene_colors <- c(
  "#F9A23D", "#F47C40", "#E7A56F", "#8CC37A",
  "#79BC5C", "#F7968C", "#ED6564"
)

# 生成表达矩阵 -----------------------------------------------------------------
n_genes <- sum(sapply(gene_groups, length))
expr_matrix <- matrix(
  rnorm(n_genes * length(cell_types)),
  nrow = n_genes,
  ncol = length(cell_types),
  dimnames = list(unlist(gene_groups), cell_types)
)

# 模拟高表达模式
for (i in seq_along(gene_groups)) {
  expr_matrix[gene_groups[[i]], cell_types[i]] <- 
    rnorm(length(gene_groups[[i]]), mean = 2, sd = 0.3)
}

# 生成GO注释数据 ---------------------------------------------------------------
n_genes_per_group <- sapply(gene_groups, length)
row_groups <- rep(names(gene_groups), n_genes_per_group)

# 为每个基因生成对应的GO值（总长度 = n_genes）
go_values <- unlist(lapply(1:7, function(i) {
  rep(runif(3, min = 5, max = 15),  # 每组3个GO值
      each = ceiling(n_genes_per_group[i] / 3))[1:n_genes_per_group[i]]
}))

# 构建右侧注释（确保长度 = n_genes）
go_anno <- rowAnnotation(
  GO_enrich = anno_barplot(
    x = go_values,
    gp = gpar(fill = gene_colors[row_groups]),  # 按分组着色
    bar_width = 0.6,
    axis_param = list(side = "top", labels_rot = 0),
    split = factor(row_groups, levels = names(gene_groups)),
    width = unit(3, "cm"),
    show_annotation_name = FALSE
  )
)

# 创建主热图 ------------------------------------------------------------------
main_heatmap <- Heatmap(
  expr_matrix,
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_split = row_groups,
  row_gap = unit(0, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  left_annotation = rowAnnotation(
    Group = anno_block(
      gp = gpar(fill = gene_colors),
      labels = names(gene_groups),
      labels_gp = gpar(col = "black", fontsize = 8)
    )
  ),
  right_annotation = go_anno,
  height = unit(n_genes * 5, "mm")
)

# 绘制并保存结果 ---------------------------------------------------------------
pdf("Integrated_Heatmap_GO.pdf", width = 10, height = 8)
draw(main_heatmap, heatmap_legend_side = "bottom")
dev.off()


