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
cds <- load_mm_data(mat_path = "./datasets/ctrl_raw_feature_bc_matrix/matrix.mtx.gz", 
                    feature_anno_path = "./datasets/ctrl_raw_feature_bc_matrix/features.tsv.gz", 
                    cell_anno_path = "./datasets/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz")

cds2 <- load_mm_data(mat_path = "./datasets/stim_raw_feature_bc_matrix/matrix.mtx.gz", 
                     feature_anno_path = "./datasets/stim_raw_feature_bc_matrix/features.tsv.gz", 
                     cell_anno_path = "./datasets/stim_raw_feature_bc_matrix/barcodes.tsv.gz")

big_cds <- combine_cds(list(cds, cds2))

cds <- big_cds
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
# 可以在colData中加一列
# plot_cells(cds, color_cells_by="batch", label_cell_groups=FALSE)
cds <- align_cds(cds, alignment_group = "sample")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method="UMAP", 
                        preprocess_method = 'Aligned')
plot_cells(cds = cds, label_groups_by_cluster=FALSE)
plot_cells(cds, genes=c(''),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE
           )

## Step 4: Cluster the cells
# 此处也注意resolution的选择
cds <- cluster_cells(cds)
plot_cells(cds) # it colors the cells by cluster according to default
plot_cells(cds, color_cells_by="partition", group_cells_by="partition") # partitions是较粗略的分类
#  When you are learning trajectories, each partition will eventually become a separate trajectory



# Order cells in pseudotime along a trajectory
## Step 5: Learn a graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

## Step 6: Order cells
# ordering each cell according to its progress along a learned trajectory
# In general, any cell on a parition that lacks a root node will be assigned an infinite pseudotime.
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# 手动选择root节点不比自动选择，下述函数在用于自己的数据集时注意修改对应的列



# Subset cells by branch
cds_sub <- choose_graph_segments(cds)



# 寻找轨迹相关基因（图形自相关算法）
# graph_test可以帮助计算和拟时序有关的基因
pr_test_res <- graph_test(cds,
                          neighbor_graph="principal_graph", #tells it to test whether cells at similar positions on the trajectory have correlated expression
                          cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))




# Find marker genes expressed by each cluster
colData(cds)
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)


# Perform differential expression analysis
# 通过拟合quasipoisson分布，依据其coefficient进行判断, 可以对细胞所处的多个变量batch、group进行分析
# Quasipoisson is a a bit less accurate than the negative binomial but much faster to fit, 
# making it well suited to datasets with thousands of cells
gene_fits <- fit_models(cds_subset, 
                        model_formula_str = "~group + cluster", # 怎么加入分类信息
                        expression_family = "quasipoisson"
                        )
fit_coefs <- coefficient_table(gene_fits)


# Annotate your cells according to type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
cds_subset <- choose_cells(cds) # 进一步筛选亚组，交互式操作
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
# Automated annotation with Garnett













