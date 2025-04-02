# 10x Genomics Visium空间转录组数据分析完整流程
# 作者：Claude 3.7 Sonnet
# 日期：2025-03-24

# 0. 加载所需的R包
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(cowplot)
library(future)
plan("multiprocess", workers = 4) # 并行计算设置，根据您的计算机性能调整

# 设置工作目录（请替换为您的数据目录）
setwd("/path/to/your/data/directory")

# 1. 数据导入和预处理
# =====================

# 1.1 导入Visium数据
# 您需要替换以下路径为您自己的样本路径
visium_dir <- "/path/to/your/10x_visium_data/"
sample_names <- c("sample1", "sample2", "sample3") # 替换为您的样本名称
tissue_img_list <- list()
spatial_obj_list <- list()

for (sample in sample_names) {
  sample_dir <- file.path(visium_dir, sample)
  spatial_obj <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample,
    filter.matrix = TRUE,
    image = NULL
  )
  # 添加样本信息
  spatial_obj$sample <- sample
  spatial_obj$orig.ident <- sample
  spatial_obj_list[[sample]] <- spatial_obj
  
  # 保存组织图像以便后续使用
  tissue_img_list[[sample]] <- spatial_obj@images[[sample]]
  
  cat(sprintf("已加载样本: %s\n", sample))
  cat(sprintf("  Spots数量: %d\n", ncol(spatial_obj)))
  cat(sprintf("  基因数量: %d\n", nrow(spatial_obj)))
}

# 1.2 质量控制
for (sample in sample_names) {
  # 计算质控指标
  spatial_obj_list[[sample]][["percent.mt"]] <- PercentageFeatureSet(spatial_obj_list[[sample]], pattern = "^MT-")
  
  # 可视化质控参数
  p1 <- VlnPlot(spatial_obj_list[[sample]], features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
                pt.size = 0, ncol = 3) + NoLegend()
  p2 <- SpatialFeaturePlot(spatial_obj_list[[sample]], 
                           features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
                           ncol = 3, alpha = c(0.1, 1))
  print(p1 / p2)
  
  # 保存质控图
  ggsave(paste0("QC_", sample, ".pdf"), plot = p1 / p2, width = 12, height = 8)
}

# 1.3 过滤低质量的spots
# 根据您的数据特性调整这些过滤阈值
for (sample in sample_names) {
  # 提取当前的spatial对象
  spatial <- spatial_obj_list[[sample]]
  
  # 根据基因计数数量和线粒体比例过滤spots
  # 请根据您的数据特性调整这些参数
  spatial_filtered <- subset(spatial, 
                             subset = nFeature_Spatial > 200 & 
                               nCount_Spatial > 500 & 
                               percent.mt < 20)
  
  cat(sprintf("样本 %s: 过滤前 %d spots, 过滤后 %d spots\n", 
              sample, ncol(spatial), ncol(spatial_filtered)))
  
  # 更新已过滤的对象
  spatial_obj_list[[sample]] <- spatial_filtered
}

# 2. 数据整合和标准化
# =====================

# 2.1 整合多个样本 (如果您有多个样本)
if (length(spatial_obj_list) > 1) {
  # 标准化每个样本
  for (sample in sample_names) {
    spatial_obj_list[[sample]] <- SCTransform(spatial_obj_list[[sample]], 
                                              assay = "Spatial", 
                                              verbose = FALSE)
  }
  
  # 选择用于整合的特征
  features <- SelectIntegrationFeatures(object.list = spatial_obj_list, 
                                        nfeatures = 3000)
  
  # 准备整合
  spatial_obj_list <- PrepSCTIntegration(object.list = spatial_obj_list, 
                                         anchor.features = features)
  
  # 寻找锚点
  anchors <- FindIntegrationAnchors(object.list = spatial_obj_list, 
                                    normalization.method = "SCT", 
                                    anchor.features = features)
  
  # 整合数据
  combined_spatial <- IntegrateData(anchorset = anchors, 
                                    normalization.method = "SCT")
  
  # 设置默认assay为integrated
  DefaultAssay(combined_spatial) <- "integrated"
  
} else {
  # 如果只有一个样本，则只需标准化
  combined_spatial <- SCTransform(spatial_obj_list[[1]], 
                                  assay = "Spatial", 
                                  verbose = FALSE)
}

# 3. 降维、聚类和UMAP可视化
# =====================

# 3.1 PCA降维
combined_spatial <- RunPCA(combined_spatial, assay = "SCT", verbose = FALSE)

# 画出PCA Elbow Plot帮助确定合适的PC数量
ElbowPlot(combined_spatial, ndims = 30)
ggsave("PCA_ElbowPlot.pdf", width = 8, height = 6)

# 3.2 基于PC进行聚类
# 根据elbow plot选择合适的PC数量，这里使用20个
n_pcs <- 20
combined_spatial <- FindNeighbors(combined_spatial, dims = 1:n_pcs)

# 尝试不同的分辨率参数
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  combined_spatial <- FindClusters(combined_spatial, resolution = res, verbose = FALSE)
}

# 选择一个适当的分辨率作为默认聚类
combined_spatial$seurat_clusters <- combined_spatial$SCT_snn_res.0.6
Idents(combined_spatial) <- "seurat_clusters"

# 3.3 UMAP降维可视化
combined_spatial <- RunUMAP(combined_spatial, dims = 1:n_pcs)

# 在UMAP上可视化聚类结果
umap_cluster_plot <- DimPlot(combined_spatial, reduction = "umap", group.by = "seurat_clusters", 
                             label = TRUE, label.size = 4) + NoLegend()

# 如果有多个样本，可视化样本信息
if (length(spatial_obj_list) > 1) {
  umap_sample_plot <- DimPlot(combined_spatial, reduction = "umap", group.by = "sample") + NoLegend()
  combined_umap_plot <- umap_cluster_plot + umap_sample_plot
  print(combined_umap_plot)
  ggsave("UMAP_Clustering.pdf", plot = combined_umap_plot, width = 12, height = 6)
} else {
  print(umap_cluster_plot)
  ggsave("UMAP_Clustering.pdf", plot = umap_cluster_plot, width = 8, height = 6)
}

# 4. 寻找空间差异表达基因
# =====================

# 4.1 为每个聚类寻找标记基因
markers <- FindAllMarkers(combined_spatial, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

# 保存标记基因列表
write.csv(markers, "All_Cluster_Markers.csv", row.names = FALSE)

# 查看每个聚类的前5个标记基因
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 打印前5个标记基因
print(top5_markers)

# 可视化热图
top10_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# 热图可视化
marker_heatmap <- DoHeatmap(combined_spatial, features = top10_markers$gene, size = 4, 
                            angle = 45, hjust = 0) + theme(axis.text.y = element_text(size = 7))
print(marker_heatmap)
ggsave("Top10_Markers_Heatmap.pdf", plot = marker_heatmap, width = 14, height = 10)

# 4.2 通过气泡图可视化一些标记基因
selected_markers <- top5_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

bubble_plot <- DotPlot(combined_spatial, features = unique(selected_markers$gene), 
                       cols = c("lightgrey", "blue"), dot.scale = 5) + 
  RotatedAxis() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
print(bubble_plot)
ggsave("Selected_Markers_BubblePlot.pdf", plot = bubble_plot, width = 12, height = 8)

# 5. 空间特征可视化
# =====================

# 5.1 在空间坐标上可视化聚类结果
for (sample in sample_names) {
  if (length(spatial_obj_list) > 1) {
    # 如果是整合后的数据，需要提取对应样本的subset
    sample_subset <- subset(combined_spatial, subset = sample == sample)
  } else {
    sample_subset <- combined_spatial
  }
  
  # 可视化聚类
  cluster_plot <- SpatialDimPlot(sample_subset, label = TRUE, label.size = 3)
  print(cluster_plot)
  ggsave(paste0("Spatial_Clusters_", sample, ".pdf"), plot = cluster_plot, width = 8, height = 7)
  
  # 可视化组织图像
  image_plot <- SpatialDimPlot(sample_subset, images = sample, crop = TRUE, pt.size.factor = 1.2)
  print(image_plot)
  ggsave(paste0("Tissue_Image_", sample, ".pdf"), plot = image_plot, width = 8, height = 7)
}

# 5.2 在空间中可视化特定基因的表达
# 为每个样本可视化前10个显著标记基因
for (sample in sample_names) {
  if (length(spatial_obj_list) > 1) {
    # 如果是整合后的数据，需要提取对应样本的subset
    sample_subset <- subset(combined_spatial, subset = sample == sample)
  } else {
    sample_subset <- combined_spatial
  }
  
  DefaultAssay(sample_subset) <- "SCT"
  
  # 选择一些标记基因进行可视化
  for (i in 1:min(10, nrow(top5_markers))) {
    gene <- top5_markers$gene[i]
    if (gene %in% rownames(sample_subset)) {
      gene_plot <- SpatialFeaturePlot(sample_subset, features = gene, 
                                      alpha = c(0.1, 1), pt.size.factor = 1.5)
      print(gene_plot)
      ggsave(paste0("Spatial_Gene_", gene, "_", sample, ".pdf"), 
             plot = gene_plot, width = 8, height = 7)
    }
  }
}

# 6. 空间领域分析
# =====================

# 6.1 空间自相关分析
# 识别具有相似表达模式的基因
# 这里以标记基因作为例子
DefaultAssay(combined_spatial) <- "SCT"
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

# 计算基因间的相关性矩阵
# 提取标记基因的表达矩阵
marker_exp <- GetAssayData(combined_spatial, slot = "data")[top_markers, ]

# 计算基因间的相关性
spatial_corr <- cor(t(as.matrix(marker_exp)), method = "pearson")
write.csv(spatial_corr, "Spatial_Gene_Correlation.csv", row.names = TRUE)

# 创建热图展示空间相关性
library(circlize)  # 用于创建颜色映射
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

heatmap_plot <- Heatmap(spatial_corr, 
                        name = "Gene\nCorrelation", 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        clustering_distance_rows = "euclidean",
                        clustering_distance_columns = "euclidean",
                        clustering_method_rows = "complete",
                        clustering_method_columns = "complete",
                        col = col_fun)

pdf("Spatial_Correlation_Heatmap.pdf", width = 10, height = 10)
draw(heatmap_plot)
dev.off()

# 额外：计算基因表达的空间自相关性（Moran's I）
# 注意：需要安装spdep包
# if (!requireNamespace("spdep", quietly = TRUE))
#   install.packages("spdep")
# library(spdep)
# 
# # 获取空间坐标
# spatial_coords <- GetTissueCoordinates(combined_spatial)
# 
# # 创建空间权重矩阵
# coords_dist <- as.matrix(dist(spatial_coords))
# # 定义邻近距离阈值（调整为适合您数据的值）
# threshold <- 100
# adj_matrix <- coords_dist <= threshold
# diag(adj_matrix) <- FALSE
# 
# # 将邻接矩阵转换为spdep权重列表对象
# listw <- mat2listw(adj_matrix, style = "W")
# 
# # 计算每个标记基因的Moran's I统计量
# morans_i_results <- data.frame(gene = character(),
#                               morans_i = numeric(), 
#                               p_value = numeric(),
#                               stringsAsFactors = FALSE)
# 
# for (gene in top_markers) {
#   # 获取基因表达值
#   gene_exp <- GetAssayData(combined_spatial, slot = "data")[gene, ]
#   # 计算Moran's I
#   moran_test <- moran.test(gene_exp, listw)
#   # 保存结果
#   morans_i_results <- rbind(morans_i_results, 
#                           data.frame(gene = gene,
#                                    morans_i = moran_test$estimate[1],
#                                    p_value = moran_test$p.value))
# }
# 
# # 保存空间自相关分析结果
# write.csv(morans_i_results, "Spatial_Morans_I_Results.csv", row.names = FALSE)
# 
# # 可视化空间自相关结果
# morans_i_results <- morans_i_results[order(-morans_i_results$morans_i), ]
# morans_plot <- ggplot(morans_i_results, aes(x = reorder(gene, morans_i), y = morans_i)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   theme_minimal() +
#   coord_flip() +
#   labs(x = "Gene", y = "Moran's I", title = "Spatial Autocorrelation of Gene Expression")
# 
# print(morans_plot)
# ggsave("Spatial_Morans_I_Plot.pdf", plot = morans_plot, width = 10, height = 8)

# 6.2 空间聚类的共同定位分析
# 检测聚类之间的空间关系
nClusters <- length(unique(Idents(combined_spatial)))
cluster_co_occur <- matrix(0, nrow = nClusters, ncol = nClusters)
rownames(cluster_co_occur) <- colnames(cluster_co_occur) <- sort(unique(as.character(Idents(combined_spatial))))

# 定义邻近范围（以像素为单位）
proximity_range <- 50

for (sample in sample_names) {
  if (length(spatial_obj_list) > 1) {
    # 如果是整合后的数据，需要提取对应样本的subset
    sample_subset <- subset(combined_spatial, subset = sample == sample)
  } else {
    sample_subset <- combined_spatial
  }
  
  # 获取空间坐标
  coords <- GetTissueCoordinates(sample_subset)
  clusters <- as.character(Idents(sample_subset))
  
  # 计算每个点对之间的距离
  for (i in 1:(nrow(coords)-1)) {
    for (j in (i+1):nrow(coords)) {
      dist <- sqrt((coords[i, "x"] - coords[j, "x"])^2 + (coords[i, "y"] - coords[j, "y"])^2)
      if (dist <= proximity_range) {
        # 如果两个点在邻近范围内，增加相应聚类的共现计数
        cluster_i <- clusters[i]
        cluster_j <- clusters[j]
        cluster_co_occur[cluster_i, cluster_j] <- cluster_co_occur[cluster_i, cluster_j] + 1
        cluster_co_occur[cluster_j, cluster_i] <- cluster_co_occur[cluster_j, cluster_i] + 1
      }
    }
  }
}

# 正则化共现矩阵
cluster_counts <- table(Idents(combined_spatial))
for (i in 1:nClusters) {
  for (j in 1:nClusters) {
    if (i != j) {
      expected <- cluster_counts[i] * cluster_counts[j] / sum(cluster_counts)
      if (expected > 0) {
        cluster_co_occur[i, j] <- log2(cluster_co_occur[i, j] / expected + 0.01)
      }
    }
  }
}
diag(cluster_co_occur) <- 0

# 可视化空间共定位关系
co_occur_heatmap <- Heatmap(cluster_co_occur,
                            name = "Log2(Observed/Expected)",
                            column_title = "Spatial Cluster Co-localization",
                            col = colorRampPalette(c("blue", "white", "red"))(100),
                            cluster_rows = TRUE,
                            cluster_columns = TRUE)

pdf("Cluster_Colocalization_Heatmap.pdf", width = 8, height = 8)
draw(co_occur_heatmap)
dev.off()

# 7. 空间领域识别和空间轨迹分析
# ===============================

# 7.1 SPOTlight空间领域解析（如果您需要更详细的空间领域分析）
# 注意：需要先安装SPOTlight包
# if (!requireNamespace("SPOTlight", quietly = TRUE))
#   remotes::install_github("MarcElosua/SPOTlight")
# library(SPOTlight)

# 7.2 空间轨迹分析 (使用slingshot包，需要先安装)
# if (!requireNamespace("slingshot", quietly = TRUE))
#   BiocManager::install("slingshot")
# library(slingshot)

# 8. 功能富集和通路分析
# =====================

# 8.1 对标记基因进行GO和KEGG富集分析
# 注意：需要先安装clusterProfiler包
# if (!requireNamespace("clusterProfiler", quietly = TRUE))
#   BiocManager::install("clusterProfiler")
# library(clusterProfiler)
# library(org.Hs.eg.db)  # 人类基因注释数据库，如果是小鼠则使用org.Mm.eg.db

# # 为每个聚类进行富集分析
# for (cluster_id in unique(markers$cluster)) {
#   # 获取当前聚类的标记基因
#   cluster_markers <- markers %>% 
#     filter(cluster == cluster_id & p_val_adj < 0.05) %>%
#     arrange(desc(avg_log2FC))
#   
#   if (nrow(cluster_markers) > 10) {  # 只有当有足够多的基因时才进行富集分析
#     # 准备基因列表
#     gene_list <- cluster_markers$gene
#     
#     # 将基因名转换为Entrez ID
#     gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#     
#     # GO富集分析
#     go_bp <- enrichGO(gene = gene_entrez$ENTREZID,
#                       OrgDb = org.Hs.eg.db,
#                       ont = "BP",  # 生物过程
#                       pAdjustMethod = "BH",
#                       pvalueCutoff = 0.05,
#                       qvalueCutoff = 0.05)
#     
#     # KEGG通路富集分析
#     kegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
#                        organism = "hsa",  # 人类，小鼠为"mmu"
#                        pAdjustMethod = "BH",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.05)
#     
#     # 保存结果
#     if (!is.null(go_bp) && nrow(go_bp) > 0) {
#       go_results <- as.data.frame(go_bp)
#       write.csv(go_results, paste0("Cluster", cluster_id, "_GO_enrichment.csv"), row.names = FALSE)
#       
#       # 可视化前10个富集通路
#       if (nrow(go_results) >= 5) {
#         p_go <- barplot(go_bp, showCategory = 10, title = paste("Cluster", cluster_id, "GO-BP Enrichment"))
#         print(p_go)
#         ggsave(paste0("Cluster", cluster_id, "_GO_barplot.pdf"), plot = p_go, width = 10, height = 8)
#         
#         p_go_dot <- dotplot(go_bp, showCategory = 10, title = paste("Cluster", cluster_id, "GO-BP Enrichment"))
#         print(p_go_dot)
#         ggsave(paste0("Cluster", cluster_id, "_GO_dotplot.pdf"), plot = p_go_dot, width = 10, height = 8)
#       }
#     }
#     
#     if (!is.null(kegg) && nrow(kegg) > 0) {
#       kegg_results <- as.data.frame(kegg)
#       write.csv(kegg_results, paste0("Cluster", cluster_id, "_KEGG_enrichment.csv"), row.names = FALSE)
#       
#       # 可视化前10个富集通路
#       if (nrow(kegg_results) >= 5) {
#         p_kegg <- barplot(kegg, showCategory = 10, title = paste("Cluster", cluster_id, "KEGG Pathway Enrichment"))
#         print(p_kegg)
#         ggsave(paste0("Cluster", cluster_id, "_KEGG_barplot.pdf"), plot = p_kegg, width = 10, height = 8)
#         
#         p_kegg_dot <- dotplot(kegg, showCategory = 10, title = paste("Cluster", cluster_id, "KEGG Pathway Enrichment"))
#         print(p_kegg_dot)
#         ggsave(paste0("Cluster", cluster_id, "_KEGG_dotplot.pdf"), plot = p_kegg_dot, width = 10, height = 8)
#       }
#     }
#   }
# }

# 9. 保存结果
# ============

# 9.1 保存Seurat对象
saveRDS(combined_spatial, "Visium_integrated_analysis.rds")

# 如果需要，可以保存每个样本的单独对象
for (sample in sample_names) {
  saveRDS(spatial_obj_list[[sample]], paste0("Visium_", sample, ".rds"))
}

# 9.2 打印分析摘要
cat("\n======== 分析摘要 ========\n")
cat(sprintf("样本数量: %d\n", length(sample_names)))
cat(sprintf("总spot数量: %d\n", ncol(combined_spatial)))
cat(sprintf("聚类数量 (分辨率0.6): %d\n", length(unique(combined_spatial$seurat_clusters))))
cat(sprintf("标记基因总数: %d\n", nrow(markers)))
cat(sprintf("差异表达基因数量: %d\n", sum(markers$p_val_adj < 0.05)))
cat("========================\n")

# 10. 可视化综合结果
# ==================

# 创建一个总结性的多面板图
for (sample in sample_names) {
  if (length(spatial_obj_list) > 1) {
    # 如果是整合后的数据，需要提取对应样本的subset
    sample_subset <- subset(combined_spatial, subset = sample == sample)
  } else {
    sample_subset <- combined_spatial
  }
  
  # Plot 1: 组织图像与聚类
  p1 <- SpatialDimPlot(sample_subset, label = TRUE, label.size = 3) + 
    ggtitle(paste("Sample", sample, "- Spatial Clusters"))
  
  # Plot 2: 基因计数空间分布
  p2 <- SpatialFeaturePlot(sample_subset, features = "nCount_Spatial") + 
    ggtitle("Spatial Gene Counts")
  
  # Plot 3: 选择一个重要的标记基因
  # 选择第一个聚类的top标记基因
  if(nrow(top5_markers) > 0) {
    selected_gene <- top5_markers$gene[1]
    p3 <- SpatialFeaturePlot(sample_subset, features = selected_gene) + 
      ggtitle(paste("Expression of", selected_gene))
    
    # 结合所有图
    combined_plot <- plot_grid(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"))
    print(combined_plot)
    ggsave(paste0("Summary_", sample, ".pdf"), plot = combined_plot, width = 18, height = 6)
  }
}

# 打印完成消息
cat("空间转录组分析完成！结果已保存到工作目录。\n")


