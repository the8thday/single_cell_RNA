# Seurat的一些用法

library(tidyverse)
library(patchwork)
library(Seurat)
library(clustree)


# read 10X data -----------------------------------------------------------

# 读取单样本counts 数据
# pbmc.data <- Read10X(data.dir = "/Users/congliu/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
# pbmc <- CreateSeuratObject(counts = pbmc.data, 
#                            project = "pbmc3k", 
#                            assay = 'RNA',
#                            min.cells = 1, min.features = 200)

## 读取多个文件的counts数据并merge在一起
for (file in c("SRR7722939", "SRR7722940","SRR7722941","SRR7722942")){
  seurat_data <- Read10X(data.dir = paste0("/Users/congliu/Downloads/10X_Rawdata/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells = 3,
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj)
}

SRR7722939@meta.data$group <- "PBMC_Pre"
SRR7722940@meta.data$group <- "PBMC_EarlyD27"
SRR7722941@meta.data$group <- "PBMC_RespD376"
SRR7722942@meta.data$group <- "PBMC_ARD614"


# 对于不同10X的count数据是采用merge还是integrate来合并数据呢
# merge将不同数据集中的cell相加，合并基因
merged_seurat <- merge(x = SRR7722939,
                       y = c(SRR7722940, SRR7722941,SRR7722942),
                       add.cell.id = c("Pre", "EarlyD27",
                                       "RespD376","ARD614"),
                       project = 'merged_seurat'
                       )
merged_seurat
head(merged_seurat@meta.data)

pbmc <- merged_seurat
slotNames(pbmc)
dim(pbmc)
# pbmc[["RNA"]]@counts # counts矩阵
# head(as.data.frame(pbmc@assays$RNA@counts))[1:5]

# 直接读入https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469下的数据
# foo <- read.csv('/Users/congliu/Downloads/GSE115469_Data.csv')
# raw_sce <- CreateSeuratObject(counts = foo)



# metadata是一个灵活的数据，可以增加多种meta信息


# Standard pre-processing workflow ----------------------------------------
# 计算MT基因的占比
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, 
                                             pattern = "^MT-") # add new column to pbmc object
head(pbmc@meta.data, 5)
# nFeature_RNA为：基因
# nCount_RNA为：应该就是count数
# The number above each plot is a Pearson correlation coefficient
FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.mt')

# 还可以对核糖体基因、红血细胞基因、管家基因等
pbmc <- PercentageFeatureSet(pbmc, "^RP[SL]",col.name = "percent.ribo")
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent.hb")
head(pbmc@meta.data, 5)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.ribo', 'percent.hb'), 
        ncol = 5)


# unique feature counts over 2,500 or less than 200 & cells that have >5% mitochondrial counts
# 此处只是为了演示筛选, 筛选满足条件的细胞
dim(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
dim(pbmc)
# 筛选基因
pbmc.qc <- pbmc
if(T){
  # Extract counts
  counts <- GetAssayData(object = pbmc.qc, slot = "counts")
  # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
  nonzero <- counts > 0
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  # Only keeping those genes expressed in more than 10 cells
  pbmc.qc_counts <- counts[keep_genes,]
  # Reassign to filtered Seurat object
  pbmc.qc <- CreateSeuratObject(pbmc.qc_counts, meta.data = pbmc.qc@meta.data)
}
dim(pbmc.qc)
table(Idents(pbmc.qc))
pbmc <- pbmc.qc

# normalizing the data
# normalization 对测序深度和基因长度进行均一化后求log, 此处是否需要对每一个样本apply?
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
head(pbmc[['RNA']]@data) # Normalized values
GetAssay(pbmc, assay = 'RNA')

head(Idents(pbmc))


# feature selection -------------------------------------------------------

# 在不同细胞群体中表达差异较大的基因往往最具有研究价值
# highly variable features 指的是在一群细胞中高表达，在另一群中低表达的基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000 # nfeature 依据具体的情况进行选择，此处选择2000个
                             )
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
CombinePlots(plots = list(plot1, plot2))


# scaling the data --------------------------------------------------------

# a linear transformation
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# 似乎是做了个Z-SCORE

all.genes <- rownames(pbmc)
head(all.genes)
pbmc <- ScaleData(pbmc, features = all.genes) # feature 默认为上一步的2000个,此处为全部基因

pbmc[["RNA"]]@scale.data[1:5, 1:5]

# scaling这步是完全必要的，默认的feature参数对下游的PCA和聚类没有影响，但可能会影响DoHeatmap()的展示

# we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
# pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), 
#                           features = rownames(seurat_phase))


# 细胞周期归类
# This has to be done after normalization and scaling
pbmc <- CellCycleScoring(object = pbmc,
                         # g2m.features = cc.genes.updated.2019$g2m.genes,
                         # s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes$g2m.genes,
                         s.features = cc.genes$s.genes
                         )
table(pbmc[[]]$Phase)
DimPlot(pbmc,reduction = "umap",label = TRUE,group.by="Phase",pt.size = .5)
head(pbmc@meta.data)
VlnPlot(seurat_phase,features =c("S.Score","G2M.Score"))


# run sctransform, 这步似乎需要在计算MT基因后进行，似乎不需要
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

write_rds(pbmc.qc,file = "./Output/Step2.after.qc_merged_seurat.rds")

# Perform linear dimensional reduction ------------------------------------

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc)) # By default, only the previously determined variable features are used as input
# VariableFeatures 返回的是FindVariableFeatures所选择的2000个基因作为特征
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:5, reduction = "pca") &
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
# it looks for UMAP, then (if not available) tSNE, then PCA
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, 
           dims = 1:4, # PC
           nfeatures = 30,
           cells = 500, 
           balanced = TRUE)
# 如何选择 PC 的数目？文档中用到两种方法，此处选择Elbowplot,不过不够精准定量, err on the higher side
# 官网教程中所推荐的三点建议值得关注
ElbowPlot(pbmc)
# 另外一种是JackStraw, 但是相对比较占资源
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)


# cluster the cell --------------------------------------------------------

# cluster the cell by previously identified PC
# seurat依据PCA的结果进行聚类
pbmc <- FindNeighbors(pbmc, dims = 1:15) # KNN+Jaccard similarity, dims为上一步确定的PC数目,此处刻意选择20
pbmc <- FindClusters(pbmc, 
                     resolution = 0.5, #0.4-1.2 typically returns good results for single-cell datasets of around 3K 
                     algorithm = 1
                     ) # Louvain algorithm
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
table(pbmc@meta.data$seurat_clusters) # identity在此被设定为clusters

# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis
pbmc <- RunUMAP(pbmc, dims = 1:15) # dims参数同上
DimPlot(pbmc, 
        reduction = "umap",
        label = TRUE
        )
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")

pbmc <- RunTSNE(object = pbmc, dims = 1:20, do.fast = TRUE)
DimPlot(pbmc,reduction = "tsne",label=TRUE)

FeaturePlot(pbmc, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

# saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


### 检查批次效应
colnames(pbmc@meta.data)
p1.compare <- wrap_plots(ncol = 3,
                      DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("Before_tSNE"),
                      DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare
# ggsave(plot=p1.compare,filename="./Output/Step3.Before_inter_sum.pdf", width = 20,height = 9)

DimPlot(pbmc,reduction = "umap", # pca, umap, tsne       
        group.by = "orig.ident",       
        label = F,       
        split.by = "orig.ident") + ggtitle("Before_intergrate_UMAP")
# 判断批次效应是通过不同样本在不同聚集处的分布情况


# 整合orig.ident的样本信息
phe <- data.frame(cell=rownames(pbmc@meta.data),
               cluster =pbmc@meta.data$seurat_clusters)
head(phe)
table(phe$cluster)
tsne_pos <- Embeddings(pbmc,'tsne')
DimPlot(pbmc,reduction = "umap",
        group.by  ='orig.ident') # 查看样本的聚类信息
DimPlot(pbmc,reduction = "umap",label=TRUE,
        split.by ='orig.ident')

# 在此处进行各种探索是合适的, 不同细胞类别中几个参数的分布
# 这里的分布也可以提供批次效应的不同
VlnPlot(pbmc,features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(pbmc,features = "percent.rb",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
FeaturePlot(pbmc,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

# 做完这步之后，active assay现在是SCT
# SCT后再次对细胞进行聚类
pbmc <- SCTransform(pbmc, method = "glmGamPoi", 
                    # ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
pbmc
srat <- pbmc
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = F)
table(srat[[]]$seurat_clusters)
DimPlot(srat, label = T)

FeaturePlot(srat,"PPBP") & 
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))


# 下一步可以通过数据整合消除批次效应


# CCA RPCA  ----------------------------------------------------------

# 将不同的数据集整合在一起, 可以消除批次效应
# RPCA 整合的更快，更适合数据集间细胞种类差别大，或者数据集分lane测的时候，或者数据集较多的情况
# CCA 适合多个样本具有大致相似的细胞类型分布

# integrate datasets by RPCA

## split
pbmc_list <- SplitObject(pbmc, split.by = "orig.ident")
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pbmc_list)
# 找到整合的锚点anchors, 默认为CCA
pbmc_anchors <- FindIntegrationAnchors(object.list = pbmc_list, anchor.features = features,
                                       reduction = 'rpca'
                                       )
# 利用anchors进行整合
pbmc_seurat <- IntegrateData(anchorset = pbmc_anchors)
names(pbmc_seurat@assays)
pbmc_seurat@active.assay


# 对于以上过程的integrate过程，有很多整合方法的选择
# harmony 替换上述的整合过程
# library(harmony)
# pbmc_seurat <- pbmc_seurat %>% RunHarmony("orig.ident", plot_convergence = T)
# harmony_embeddings <- Embeddings(pbmc_harmony, 'harmony')
# harmony_embeddings[1:5, 1:5]


# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = FALSE) # 整合后仍需scale
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = FALSE)
ElbowPlot(pbmc_seurat)
# 另外一种是JackStraw, 但是相对比较占资源
pbmc_seurat <- JackStraw(pbmc_seurat, num.replicate = 100)
pbmc_seurat <- ScoreJackStraw(pbmc_seurat, dims = 1:20)
JackStrawPlot(pbmc_seurat, dims = 1:20) # 对此图的解释

print(x = pbmc_seurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# 数据整合后同上选择PCA的维度用于后面分析
pbmc_seurat <- RunTSNE(pbmc_seurat, reduction = "pca", dims = 1:15)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:15)
pbmc_seurat <- FindNeighbors(pbmc_seurat, reduction = "pca", dims = 1:15)
# pbmc_seurat <- FindClusters(pbmc_seurat, resolution = 0.5)


p2.compare <- wrap_plots(ncol = 3,
                      DimPlot(pbmc_seurat, reduction = "pca", group.by = "group")+NoAxes()+ggtitle("After_PCA"),
                      DimPlot(pbmc_seurat, reduction = "tsne", group.by = "group")+NoAxes()+ggtitle("After_tSNE"),
                      DimPlot(pbmc_seurat, reduction = "umap", group.by = "group")+NoAxes()+ggtitle("After_UMAP"),
                      guides = "collect"
)
p1.compare / p2.compare

# 寻找最佳的分辨率，resolution参数。
for (res in c(0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2)) {
  # res=0.01
  print(res)
  pbmc_seurat <- FindClusters(pbmc_seurat, graph.name = "integrated_snn", resolution = res, algorithm = 1)
}
# pbmc_seurat <- FindClusters(pbmc_seurat, resolution = 0.5)
pbmc_seurat.res <- FindClusters(
  object = pbmc_seurat,
  resolution = c(seq(0,1.6,.2))
)
clustree::clustree(pbmc_seurat.res, prefix = "integrated_snn_res.",node_size = 10,
                   node_alpha = 0.8)

pbmc_seurat@meta.data %>% select(starts_with("integrated_snn_res"))  %>% 
  mutate(integrated_snn_res.0.0=0) %>% clustree( prefix = "integrated_snn_res.",layout = "sugiyama")   
####画全部的resolution, 依据XX筛选合适的分辨率


#umap可视化
cluster_umap <- cowplot::plot_grid(ncol = 4,
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.05", label = T)& NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.1", label = T) & NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.2", label = T)& NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.3", label = T)& NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.5", label = T) & NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T) & NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.1", label = T) & NoAxes(),
                          DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.1.2", label = T) & NoAxes()
)
cluster_umap

## 选择合适的分辨率（resolution，此处结合后续的marker gene在反反复复的进行确定
DimPlot(pbmc_seurat, reduction = "umap", group.by = "integrated_snn_res.0.5", label = T)
pbmc_seurat <- SetIdent(pbmc_seurat,value = "integrated_snn_res.0.5")
DimPlot(pbmc_seurat,label = T)

write_rds(pbmc_seurat,file = "./Output/Step3.Intergration_seurat.rds")


## marker genes，用于进行初步的细胞注释
pbmc <- pbmc_seurat
DimPlot(pbmc, reduction = "tsne",label = T) + DimPlot(pbmc, reduction = "umap",label = T)


DefaultAssay(pbmc)
DefaultAssay(pbmc) <- "RNA"
markerGenes <- c("CD3D", # 定位T细胞
                 "CD3E", # 定位T细胞
                 "TRAC", # 定位T细胞
                 "IL7R", # CD4 T cells
                 "GZMA", # NK T /效应T
                 "NKG7", # NK T cells
                 "CD8B", # CD8 T cells
                 "FCGR3A", # CD16(FCGR3A)+ Mono / NK / 效应 CD8+ T
                 "CD14", # CD14+ monocyte
                 "LYZ", #Mono
                 "MS4A1", #B细胞
                 "FCER1A","LILRA4","TPM2", #DC
                 "PPBP","GP1BB"# platelets
)
VlnPlot(pbmc, features = markerGenes)
DotPlot(pbmc, features = markerGenes, dot.scale = 8) + RotatedAxis()
p3 = FeaturePlot(pbmc,
                 features = markerGenes,
                 label.size = 4,
                 repel = T,label = T)&
  theme(plot.title = element_text(size=10),
        legend.position = 'none'
        )&
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
p3


## 5.3 根据markers重新命名clusters
celltype <- data.frame(ClusterID=0:13,celltype='NA')
celltype[celltype$ClusterID %in% c(0,13),2]='B cells'
celltype[celltype$ClusterID %in% c(1),2]='NK cells'
celltype[celltype$ClusterID %in% c(2),2]='CD8+ Effector T cells'
celltype[celltype$ClusterID %in% c(12),2]='CD8+ cytotoxic T cells'
celltype[celltype$ClusterID %in% c(4,6),2]='Naive/memory T cells'
celltype[celltype$ClusterID %in% c(3),2]='CD4+ T cells'
celltype[celltype$ClusterID %in% c(10),2]='Dendritic cells'
celltype[celltype$ClusterID %in% c(5,7),2]='CD14+ monocyte'
celltype[celltype$ClusterID %in% c(8),2]='CD16+ monocyte'  
celltype[celltype$ClusterID %in% c(9),2]='Platelets'
celltype[celltype$ClusterID %in% c(11),2]='RBC'
celltype

pbmc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@active.ident == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(pbmc@meta.data$celltype)

## 可视化
pbmc = SetIdent(pbmc, value = "celltype")
wrap_plots(ncol = 2,
           DimPlot(pbmc, reduction = "umap", label = T,repel = T),
           DimPlot(pbmc, reduction = "tsne", label = T,repel = T),
           guides = "collect"
)

ggsave(filename = "Output/Step5.annotation_plot.pdf",
       width = 9,height = 17)

## 保存
write_rds(pbmc,file = "Output/Step5.annotation_pbmc.rds")


# Finding differentially expressed features (cluster biomarkers) ----------

# find cluster's marker gene
# find all markers of cluster 2, 选择cluster2的marker 基因, 如果ident.2为NULL，则默认比全部其他
DefaultAssay(pbmc) <- "RNA"
unique(pbmc@meta.data$seurat_clusters) # list all clusters

cluster2.markers <- FindMarkers(pbmc, 
                                slot = 'data',
                                ident.1 = 'CD8+ Effector T cells', # which cluster, 由SetIdent决定
                                min.pct = 0.25, # a minimum percentage in either of the two groups of cells
                                test.use = "wilcox"
)
head(cluster2.markers, n = 5)
# p_val : p_val (unadjusted)
# avg_log2FC : log fold-change of the average expression between the two groups. 
# Positive values indicate that the feature is more highly expressed in the first group.
# pct.1 : The percentage of cells where the feature is detected in the first group
# pct.2 : The percentage of cells where the feature is detected in the second group
# p_val_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.

# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster之间的差异marker, 此步骤为DEGs.
# 此处只是找到一个样本内不同cell type间的DEGs，不同样本间同样的cell type间的差异？

# pbmc <- SetIdent(pbmc, value = "seurat_clusters")
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, 
                                ident.2 = c(0, 3), 
                                min.pct = 0.25,
                                test.use = "wilcox"
)
head(cluster5.markers, n = 5)

cluster5.markers2 <- FindAllMarkers(pbmc,
                                    ident.1 = 5,
                                    ident.2 = NULL,
                                    test.use = 'DESeq2'
)
head(cluster5.markers2, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.5)
print(head(pbmc.markers))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) # 每个cluster的top2 marker gene

# DT::datatable(pbmc.markers)

# some plot methods
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", 'GZMK'),
            slot = 'data'
            ) # 此处可以依据已知的marker gene或者上一步得到的marker gene做进一步的筛选
RidgePlot(pbmc, features = c("MS4A1", "CD79A"))
CellScatter(pbmc, cell1 = '', cell2 = '')
DotPlot(pbmc)

# 依据metadata中的分组信息进行差异分析
pbmc.diff <- pbmc
pbmc.diff <- SetIdent(pbmc.diff, value = "orig.ident")
diff.markers <- FindMarkers(pbmc.diff, 
                                ident.1 = 'SRR7722939', 
                                ident.2 = 'SRR7722940', 
                                min.pct = 0.25,
                                test.use = "wilcox"
)
head(diff.markers, n = 5)


## Identification of conserved markers in all conditions
## 分组查找marker gene
DefaultAssay(seurat_integrated) <- "RNA"
FindConservedMarkers(seurat_integrated,
                     ident.1 = cluster, # 欲查找DE的cluster
                     grouping.var = "sample",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.5)



# Assigning cell type identity to clusters --------------------------------

# 对于已有marker的细胞分群，可以进行命名
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.3) + 
  NoLegend()


# 细胞聚类分群后，如何定义每一类的细胞分群
# cellmarker，可以查找组织的细胞构成
# 可参考SingleR.Rmd脚本中关于SingleR的用法。



# ‘SubsetData’ then re-run the Seurat -------------------------------------
# use SubsetData to select particular clusters 
dermal.subset <- SubsetData(object = pbmc, ident.use = '')



# 细胞重命名 -------------------------------------------------------------------

nCoV.integrated1 <- RenameIdents(object = seurat_integrated,
                                 "0" = "NK",
                                 "1" = "T_cell",
                                 "2" = "T_cell",
                                 "3" = "T_cell",
                                 "4" = "NK",
                                 "5" = "T_cell",
                                 "6" = "T_cell",
                                 "7" = "Fibroblast",
                                 "8" = "T_cell",
                                 "9" = "T_cell",
                                 "10" = "Macrophage",
                                 "11" = "NK",
                                 "12" = "Malignant",
                                 "13" = "Fibroblast",
                                 "14" = "B_cell",
                                 "15" = "Malignant",
                                 "16" = "Endothelial",
                                 "17" = "Macrophage",
                                 "18" = "Endothelial",
                                 "19" = "Endothelial",
                                 "20" = "T_cell",
                                 "21" = "Fibroblast",
                                 "22" = "Fibroblast",
                                 "23" = "Fibroblast",
                                 "24" = "Mast",
                                 "25" = "T_cell",
                                 "26" = "Malignant",
                                 "27" = "NK",
                                 "28" = "Cholangiocyte",
                                 "29" = "B_cell",)
table(Idents(object =nCoV.integrated1))



