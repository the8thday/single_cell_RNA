# Seurat的一些用法

library(tidyverse)
library(Seurat)


# read 10X data -----------------------------------------------------------

pbmc.data <- Read10X(data.dir = "/Users/congliu/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           assay = 'RNA',
                           min.cells = 1, min.features = 200)

slotNames(pbmc)



# Standard pre-processing workflow ----------------------------------------

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, 
                                             pattern = "^MT-") # add new column to pbmc object
head(pbmc@meta.data, 5)

# unique feature counts over 2,500 or less than 200 & cells that have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
head(pbmc[['RNA']]@data)



# feature selection -------------------------------------------------------

# 在不同细胞群体中差异较大的基因往往最具有研究价值
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000 # nfeature
                             )
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


# scaling the data --------------------------------------------------------

# a linear transformation
# the mean expression across cells is 0
# the variance across cells is 1

all.genes <- rownames(pbmc)
head(all.genes)
pbmc <- ScaleData(pbmc, features = all.genes) # feature 默认为上一步的2000个,此处为全部基因

SCTransform


# Perform linear dimensional reduction ------------------------------------

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc)) # By default, only the previously determined variable features are used as input
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, 
           dims = 1, # PC
           cells = 500, 
           balanced = TRUE)
# 如何选择 PC 的数目？文档中用到两种方法，此处选择Elbowplot,不过不够精准定量, err on the higher side
ElbowPlot(pbmc)


# cluster the cell --------------------------------------------------------

# cluster the cell by previously identified PC
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

# saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


# Finding differentially expressed features (cluster biomarkers) ----------

# find cluster's marker gene
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, 
                                ident.1 = 2, # which cluster
                                min.pct = 0.25, # a minimum percentage in either of the two groups of cells
                                test.use = "wilcox"
                                )
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

# some plot methods
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
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







