# Seurat的一些用法

library(tidyverse)
library(Seurat)


# read 10X data -----------------------------------------------------------

# counts 数据
pbmc.data <- Read10X(data.dir = "/Users/congliu/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           assay = 'RNA',
                           min.cells = 1, min.features = 200)

slotNames(pbmc)
dim(pbmc)
# pbmc[["RNA"]]@counts # counts矩阵

# 直接读入https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469下的数据
# a <- read.csv('GSE115469_Data.csv')
# raw_sce <- CreateSeuratObject(counts = a)


# Standard pre-processing workflow ----------------------------------------
# 计算MT基因的占比
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, 
                                             pattern = "^MT-") # add new column to pbmc object
head(pbmc@meta.data, 5)
# nFeature_RNA为：基因
# nCount_RNA为：应该就是count数
# 
FeatureScatter(object = pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# unique feature counts over 2,500 or less than 200 & cells that have >5% mitochondrial counts
# 此处只是为了演示筛选, 筛选满足条件的细胞
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dim(pbmc)

# normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
head(pbmc[['RNA']]@data) # Normalized values
GetAssay(pbmc, assay = 'RNA')


# feature selection -------------------------------------------------------

# 在不同细胞群体中差异较大的基因往往最具有研究价值
# highly variable features 指的是在一群细胞中高表达，在另一群中低表达的基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000 # nfeature 依据具体的情况进行选择，此处选择2000个
                             )
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


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


# run sctransform, 这步需要在计算MT基因后进行
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)



# Perform linear dimensional reduction ------------------------------------

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc)) # By default, only the previously determined variable features are used as input
# VariableFeatures 返回的是FindVariableFeatures所选择的2000个基因作为特征
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, 
           dims = 1:4, # PC
           cells = 500, 
           balanced = TRUE)
# 如何选择 PC 的数目？文档中用到两种方法，此处选择Elbowplot,不过不够精准定量, err on the higher side
# 官网教程中所推荐的三点建议值得关注
ElbowPlot(pbmc)
# 另外一种是JackStraw, 但是相对比较占资源
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)


# cluster the cell --------------------------------------------------------

# cluster the cell by previously identified PC
# seurat依据PCA的结果进行聚类
pbmc <- FindNeighbors(pbmc, dims = 1:20) # KNN+Jaccard similarity, dims为上一步确定的PC数目,此处刻意选择20
pbmc <- FindClusters(pbmc, 
                     resolution = 0.5 #0.4-1.2 typically returns good results for single-cell datasets of around 3K 
                     ) # Louvain algorithm
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, 
        reduction = "umap",
        label = TRUE
        )

pbmc <- RunTSNE(object = pbmc, dims = 1:20, do.fast = TRUE)
DimPlot(pbmc,reduction = "tsne",label=TRUE)

# saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# 整合orig.ident的样本信息
phe <- data.frame(cell=rownames(pbmc@meta.data),
               cluster =pbmc@meta.data$seurat_clusters)
head(phe)
table(phe$cluster)
tsne_pos <- Embeddings(pbmc,'tsne') 
DimPlot(pbmc,reduction = "tsne",
        group.by  ='orig.ident') # 查看样本的聚类信息
DimPlot(pbmc,reduction = "tsne",label=TRUE,
        split.by ='orig.ident')


# Finding differentially expressed features (cluster biomarkers) ----------

# find cluster's marker gene
# find all markers of cluster 2, 选择cluster2的marker 基因
unique(sce@meta.data$seurat_clusters) # list all clusters

cluster2.markers <- FindMarkers(pbmc, 
                                ident.1 = 2, # which cluster
                                min.pct = 0.25, # a minimum percentage in either of the two groups of cells
                                test.use = "wilcox"
                                )
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster之间的差异marker
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(head(pbmc.markers))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) # 没个cluster的top2

# DT::datatable(pbmc.markers)

# some plot methods
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E"),
            slot = 'data'
            )
RidgePlot(pbmc, features = c("MS4A1", "CD79A"))
CellScatter(pbmc)
DotPlot(pbmc)

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







