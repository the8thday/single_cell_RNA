## ---------------------------
##
## Script name: Monocle3.R
##
## Purpose of script:
##
## Author: LiuCong
##
## Date Created: 2022-03-29
##
## Copyright (c) cliu, 2022
## Email: dibatian@live.com
##
## ---------------------------
##
## Notes: Monocle3的使用记录
##   
##
## ---------------------------


library(monocle3)


# basic at a glance-----------------------------------------------------------

# read 10X data
# Provide the path to the Cell Ranger output.0x_data/outs/filtered_feature_bc_matrix/
# cds <- load_cellranger_data("~/Downloads/10x_data", umi_cutoff = 100)
# Monocle 3 is designed for use with absolute transcript counts
cds <- load_mm_data(mat_path = "~/Downloads/Rawdata/SRR7722937/matrix.mtx.gz", 
                    feature_anno_path = "~/Downloads/Rawdata/SRR7722937/features.tsv.gz", 
                    cell_anno_path = "~/Downloads/Rawdata/SRR7722937/barcodes.tsv.gz")

big_cds <- combine_cds(list(cds, cds2))
colData(cds)


# 
# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, 
                      method = 'PCA',
                      num_dim = 100)

## Step 2: Remove batch effects with cell alignment(Optional)
# 可以在colData中加一列，10X的数据怎么加呢
# plot_cells(cds, color_cells_by="batch", label_cell_groups=FALSE)
# cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds = cds)
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))

## Step 4: Cluster the cells
# 此处也注意resolution的选择
cds <- cluster_cells(cds)
plot_cells(cds) # it colors the cells by cluster according to default
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") # partitions是较粗略的分类

# Order cells in pseudotime along a trajectory
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)
plot_cells(cds)

# 寻找轨迹相关基因（图形自相关算法）
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))

# Find marker genes expressed by each cluster
colData(cds)
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)


# Perform differential expression analysis
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time", # 怎么加入分类信息
                        expression_family = "quasipoisson"
                        )


# Annotate your cells according to type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
cds_subset <- choose_cells(cds) # 进一步筛选亚组，交互式操作
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
# Automated annotation with Garnett



# single cells along pseudotime -------------------------------------------
# 轨迹分析常用于阐释细胞发育等细胞变化的过程
# 在完成上述step4后，再展开trajectory analysis

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)










