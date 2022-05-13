
#seurat结构
#counts：主要是 counts或者TPKM的raw data，未经normalized  ###https://www.jianshu.com/p/8d9fcef27b05
#data：是经过normalized的表达矩阵seurat_integrated[["RNA"]]@data[1:6,1:5] 
#scale.data：是已经scaled out的表达矩阵seurat_integrated[["RNA"]]@scale.data[1:6,1:5] 

#####加载包
rm(list = ls())
gc()
options(repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")################# 设置R包安装镜像清华镜像
packageVersion("Seurat") ###查看包的版本‘3.2.3’ 4.0.3
options(stringsAsFactors = F)
####setwd("/Node22Data/home/wangting/project/icc/single_cell_RNA_seq/cellranger_count/Rproj")
suppressPackageStartupMessages(library(future)) 
suppressPackageStartupMessages(library(future.apply)) 
plan("multicore", workers = 10) #多线程
options(future.globals.maxSize = 40000 * 1024^2)
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(RCurl))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
#suppressMessages(library(SingleCellExperiment))
#suppressMessages(library(AnnotationHub))
#suppressMessages(library(ensembldb))
library("ggsci")


setwd('/Users/congliu/Downloads/10X_Rawdata/')

### 9种 细胞的颜色
my_colors=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF",
            "#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF")

#sample_colors= pal_uchicago("light", alpha =1)(5)[2:6]
sample_colors= c("#D6D6CEFF" ,"#FFB547FF", "#ADB17DFF" ,"#5B8FA8FF","#800000FF")
group_colors=c("#7AA6DCCC","#003C67CC","#A73030CC" ) #分组颜色

###### #### #####生信技能树的单细胞合并merge教程：使用seurat3的merge功能整合8个10X单细胞转录组样本 
## https://mp.weixin.qq.com/s/nHyijzvonEadXEiO9SCy8A

##读入GEO的矩阵
#new_counts <- read.table(file="GSM4008624_Adult-Liver1-2_dge.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t",row.names=1)
#liver1.2 <- CreateSeuratObject(counts = new_counts, min.cells = 3, project= "liver1.2")
#多个object合并
#sce.big <- merge(liver1.1, y = c(liver1.2,liver2,liver4.1,liver4.2), add.cell.ids = c("liver1.1", "liver1.2", "liver2","liver4.1","liver4.2"), project = "HCLliver")
##### B748 <- Read10X(data.dir = "./B748/")
#sce.big <- CreateSeuratObject(counts = B748, project = "B748", min.cells = 3, min.features = 200) 单独一个的读入

####多个文件读入合并
folders=list.files('/Users/congliu/Downloads/10X_Rawdata/')   #########当前路径下的文件
folders ########[] "N368" "N725" "T368" "T725" ######每个文件包含3个gz：barcodes.tsv.gz features.tsv.gz matrix.mtx.gz
sceList <- lapply(folders, function(folder) {
  CreateSeuratObject(counts = Read10X(folders), min.cells = 3, min.features = 200, project = folder)
}) ######## 创建Seurat对象#lapply函数,可以循环处理列表中的每一个元素,#lapply(参数)：lapply(列表,函数/函数名,其他参数)#总是返回一个列表
####在数据导入的时候，数据集中测到的少于200个基因的细胞（min.features = 200），和少于3个细胞覆盖的基因（min.cells = 3），就已经被过滤掉了。

sce.big <- merge(sceList[[1]], 
                 y = c(sceList[[2]], sceList[[3]], sceList[[4]], sceList[[5]], sceList[[6]]), 
                 add.cell.ids = folders, project = "icc") ######## merge全合并
head(colnames(sce.big))
tail(colnames(sce.big))
unique(sapply(X = strsplit(colnames(sce.big), split = "_"), FUN = "[", 1))
table(sce.big$orig.ident) ## 每个样本有多少个cells
##########merge  结束#################################################################################################################       

#############################哈佛大学单细胞课程：笔记汇总前篇https://mp.weixin.qq.com/s/GTYMsgb_CzCV-K1X6MPfQw 质控，过滤，更全面

##########接着上一步已经合并后的大的merged的sce.big

sce.big
##An object of class Seurat
###22034 features across 24495 samples within 1 assay
####Active assay: RNA (22034 features, 0 variable features)
#####行代表基因(22034)，列代表细胞(24495)
head(sce.big@meta.data)  ######查看信息

##1、#####################将每个细胞的UMI数量进行log10转换并加入到metadata中log10GenesPerUMI这一列
sce.big$log10GenesPerUMI <- log10(sce.big$nFeature_RNA)/log10(sce.big$nCount_RNA)

##2、#######计算线粒体相关基因比例，线粒体基因以MT-开头，注意！("^MT-")只限于人哦。。。
sce.big$mitoRatio <- PercentageFeatureSet(object = sce.big, pattern = "^MT-")
sce.big$mitoRatio <- sce.big@meta.data$mitoRatio/100   #####线粒体相关基因百分比

###3、####### Create metadata dataframe，把merged_seurat中的@meta.data取出来给 metadata
metadata <- sce.big@meta.data 

###4、################ Add cell IDs to metadata在metadata中添加一列叫cells，行名等于metadata的行名
##
metadata$cells <- rownames(metadata)
#metadata$cells <- paste ("B748",metadata$cells, sep = "_") ##单个，给细胞加上sampleID
###5、############################# Rename columns，修改列名，改名，把 nCount_RNA改成nUMI, nFeature_RNA改成nGENE; orig.ident改成seq_folder
metadata <- metadata %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA, seq_folder = orig.ident,)
head(metadata)                      
#####6、########## Create sample column在meta中增加sample这一列
library(stringr)
metadata$sample <- metadata$seq_folder

metadata$batch <- NA ##增加一列批次信息 batch1 batch2
metadata$batch[which(str_detect(metadata$sample, "^N"))] <- "batch1" #匹配 ^行首
metadata$batch[which(str_detect(metadata$sample, "^T"))] <- "batch1"
metadata$batch[which(str_detect(metadata$sample, "P"))] <- "batch2"
metadata$batch[which(str_detect(metadata$sample, "T$"))] <- "batch2"
metadata$batch[which(str_detect(metadata$sample, "T1$"))] <- "batch2"  ## $行尾
metadata$batch[which(str_detect(metadata$sample, "T2$"))] <- "batch2"

metadata$group <-NA ##增加一列批次信息分组P T
metadata$group[which(str_detect(metadata$sample, "^N"))] <- "P" #匹配
metadata$group[which(str_detect(metadata$sample, "P"))] <- "P"
metadata$group[which(str_detect(metadata$sample, "T"))] <- "T"

###7、####################将新的metadata重新添加回seurat，把metadata赋值给merged_seurat中的@meta.data
head(metadata)
sce.big@meta.data <- metadata  
head(sce.big@meta.data)
######################################################################保存原始合并后的
save(sce.big, file="merged_seurat.RData")   #####保存

      
#########所有样本合并质控画图
pdf("1.8.qualitlyQC.pdf",width= 15, height = 10)
VlnPlot(sce.big,features = c("nGene", "nUMI", "mitoRatio"),ncol = 3)
dev.off()

 ###############画相关性
pdf("1.9.corla.pdf",width= 15, height = 10)
plot1 <- FeatureScatter(sce.big,feature1 = "nUMI",feature2 = "mitoRatio")
plot2 <- FeatureScatter(sce.big,feature1 = "nUMI",feature2 = "nGene")
####plot2 <- FeatureScatter(sce.big,feature1 = "nUMI",feature2 = "nGene",group.by= "sample")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#######细胞计数
library(ggplot2)
pdf("1.1.0.numbers_of_cells.pdf",width= 15, height = 10)
ggplot(metadata,aes(x=sample, fill=sample)) + geom_bar() +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) +  ggtitle("Numbers of Cells")######细胞计数
dev.off()

#####UMI计数
pdf("1.2.number_of_UMIs.pdf",width= 15, height = 10)
ggplot(metadata,aes(color=sample, x=nUMI, fill= sample)) + geom_density(alpha = 0.2) +
    scale_x_log10() +
   theme_classic() +
    ylab("Cell density") +
   geom_vline(xintercept = 700)  +
   theme(plot.title = element_text(hjust=0.5, face="bold")) +
   ggtitle("UMI counts (transcripts) per cell")
 dev.off()
   
#####每个细胞观察到的基因 genes 
pdf("1.3.number_of_genes.pdf",width= 15, height = 10)
   ggplot(metadata,aes(color=sample, x=nGene, fill= sample)) + geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
   geom_vline(xintercept = 300)  +
   theme(plot.title = element_text(hjust=0.5, face="bold")) +
   ggtitle("Genes detected per cell")
dev.off()

###### 一起评估UMI的数量和每个细胞检测到的基因数量genes
pdf("1.4.number_of_genes_boxplot.pdf",width= 15, height = 10)
  ggplot(metadata,aes(x=sample, y=log10(nGene), fill=sample)) + geom_boxplot() +
      theme_classic() +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      ggtitle("Num of Cells vs Num of Genes") 
dev.off()

#######UMIs vs. genes评估线的斜率以及该图右下象限中数据点的分布
pdf("1.5.gene_umi_cor.pdf",width= 15, height = 10)
ggplot(metadata,aes(x=nUMI, y=nGene, color=mitoRatio)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + stat_smooth(method=lm) +
 scale_x_log10() +
scale_y_log10() +
  theme_classic() +
    geom_vline(xintercept = 2000) +
    geom_hline(yintercept = 500) +
    facet_wrap(~sample)  
dev.off()

 #####线粒体counts比例Mitochondrial counts ratio
 pdf("1.6.mtio.counts.ratio.pdf",width= 15, height = 10)
ggplot(metadata,aes(color=sample, x=mitoRatio, fill=sample)) + geom_density(alpha = 0.2) + scale_x_log10() +
      theme_classic() +
      geom_vline(xintercept = 0.2)  +
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      ggtitle("Mitochondrial counts ratio") 
dev.off()

 #######复杂度：GenesPerUMI：每个UMI检测到更多的基因，我们的数据更复杂
pdf("1.7.complexity_gene_expression.pdf",width= 15, height = 10)
ggplot(metadata,aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8) +
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      ggtitle("complexity of the gene expression") 
dev.off()

     
#############细胞层面过滤过滤   Filter out low quality reads using selected thresholds - these will change with experiment  过滤条件包括：nUMI  + nGene +  log10GenesPerUMI + mitoRatio

filtered_seurat <- subset(x = sce.big, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) &            ####10的0.8次方=6.1*****
                           (mitoRatio < 0.20))  
                           
# filtered_seurat <- subset(x = sce.big, subset = (nUMI >= 2001) & (nGene >= 500) & (nGene < 6000)  & (mitoRatio < 0.20)  & (log10GenesPerUMI > 0.80))
str(filtered_seurat)    
######################查看对象每一个@开头的都是对应一个插槽(slot)，里面已经存放或者未来会存放分析所需的信息，比如我们的原始数据就放在了这个对象里面。


################基因层面过滤，删除0表达值的基因，和在少于10个细胞中表达的基因。

# Output a logical vector for every gene on whether the more than zero counts per cell# Extract counts

counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# #####Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
#####Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
######## Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

############################################   细胞周期
######Seurat亮点之细胞周期评分和回归  https://www.jianshu.com/p/5473dcf4480d参考

######NormalizeData ,关于Seurat归一化原理，可以看这一篇：https://www.biorxiv.org/content/biorxiv/early/2019/03/18/576827.full.pdf,默认是进行一个全局的LogNormalize操作：log1p(value/colSums[cell-idx] *scale_factor)，其中log1p指的是log(x + 1)
######使用默认参数默认的方法LogNormalize，数据标准化。
###结果在这里seurat_phase[["RNA"]]@data

seurat_phase <- NormalizeData(filtered_seurat) 
seurat_phase[['RNA']][1:3,1:3]  ##查看NormalizeData的结果
# Load cell cycle markers下载细胞周期markers
load("/Node22Data/home/wangting/project/icc/single_cell_RNA_seq/cellranger_count/Rproj_N368_N725/cycle.rda")

######计算细胞周期评分
seurat_phase <- CellCycleScoring(seurat_phase, g2m.features = g2m_genes, s.features = s_genes)
head(seurat_phase[[]])

####利用细胞周期基因进行PCA分析，不出所料，细胞完全按 phase 分离开 phaseS.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期

seurat_phase <- FindVariableFeatures(seurat_phase,selection.method = "vst",nfeatures = 2000,verbose = FALSE)
# Scale the counts                     
seurat_phase <- ScaleData(seurat_phase)
set.seed(123)
seurat_phase <- RunPCA(seurat_phase)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30)
seurat_phase <- FindClusters(seurat_phase, resolution = 0.4)
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30)

pdf("2.0.1.phase.PCA.pdf",width= 15, height = 10)
DimPlot(seurat_phase, reduction = "pca",group.by= "Phase", split.by = "Phase")
dev.off()
pdf("2.0.1.phase.umap.pdf",width= 15, height = 10)
DimPlot(seurat_phase,reduction = "umap",group.by = "Phase")
dev.off()
pdf("2.0.1.sample.umap.pdf",width= 15, height = 10)
DimPlot(seurat_phase,reduction = "umap",group.by = "sample")
dev.off()
####看看细胞周期是否有影响，如果有的话，要移除
#####细胞周期评分分布图
pdf("2.0.2.phase.vlnplot.pdf",width= 15, height = 10)
VlnPlot(seurat_phase,features =c("S.Score","G2M.Score"))###总体
dev.off()
###s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
pdf("2.0.2.phase.MKI67.umap.pdf",width= 15, height = 10)
FeaturePlot(seurat_phase,features = c("PCNA","MKI67"),reduction = "umap")
dev.off()
##head(seurat_phase[[]])####查看评分


####在数据标归一化时去除细胞周期影响，这一步耗时非常长
seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))
####再次跑流程做PCA时，就看不到细胞周期相关基因对主成分的贡献了
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), nfeatures.print = 10)
seurat_phase <- RunPCA(seurat_phase, npcs = 50, verbose = FALSE)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30)
seurat_phase <- FindClusters(seurat_phase, resolution = 0.4)
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30)
seurat_phase <- RunTSNE(seurat_phase, dims = 1:30)
pdf("2.0.1.after.phase.umap.pdf",width= 15, height = 10)
DimPlot(seurat_phase,reduction = "umap",group.by = "Phase")
dev.off()

#### 整合前, 先普通的VST标准化，使用 PCA看有无批次效应 ################################################

options(future.globals.maxSize = 4000 * 1024^2) 
before_seurat_integrated <- CellCycleScoring(filtered_seurat,g2m.features=g2m_genes,s.features=s_genes)
DefaultAssay(before_seurat_integrated) <- "RNA" ####选择过滤后 $RNA assay，
before_seurat_integrated <- NormalizeData(object = filtered_seurat, normalization.method = "LogNormalize") #scale.factor = 10000)  ####标准化[["RNA"]]@data
before_seurat_integrated[['RNA']][1:3,1:3] ###结果在这里before_seurat_integrated[["RNA"]]@data

#################################查找高变基因
before_seurat_integrated <- FindVariableFeatures(before_seurat_integrated, selection.method = "vst") ###查找高变基因
length(VariableFeatures(before_seurat_integrated)) ###查找到的高变基因

###可视化展示HVGs，火山图可以给top差异基因加上标签
top10 <- head(VariableFeatures(before_seurat_integrated), 10);
top10
pdf("2.0.3.top10.HGC.pdf",width= 15, height = 10)
plot1 <- VariableFeaturePlot(before_seurat_integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

################## 如果最终不整合，就在这里对before_seurat_integrated进行 scale 标准化。如果整合成功，就不需要before_seurat_integrated 就不在这里对 all.genes scale.
all.genes <- rownames(before_seurat_integrated)
DefaultAssay(before_seurat_integrated) <- "RNA"  #原始RNA数据修改
before_seurat_integrated <- ScaleData(before_seurat_integrated, 
                                      features = VariableFeatures(before_seurat_integrated), 
                                      vars.to.regress = c("mitoRatio","nUMI")) 
####features = all.genes，features =VariableFeatures(before_seurat_integrated)之前鉴定的HVGs进行标准化。
#scale函数默认是只针对之前鉴定的HVGs进行标准化（版本3中默认得到2000个HVGs）#这样操作的结果中降维和聚类不会被影响，但是只对HVGs基因进行标准化，下面画的热图可能会受到全部基因中某些极值的影响。所以为了热图的结果，还是对所有基因进行归一化比较好
before_seurat_integrated[["RNA"]]@scale.data[30:34,1:3]  ##Scale 后的结果在这里scale.data
length(rownames(before_seurat_integrated))

###############对高变基因做PCA（scale后才能PCA）
set.seed(123)
before_seurat_integrated <- RunPCA(before_seurat_integrated, 
                                   features = VariableFeatures(object = before_seurat_integrated),
                                   ndims.print = 1:5, nfeatures.print = 5)  ####### features = VariableFeatures(object = before_seurat_integrated)意思是对高变基因 run PCA
head(before_seurat_integrated@reductions$pca@feature.loadings) #### 高变基因 run PCA,这50个主成分的全部结果在这里

###############  plot PCA  按照批次检查PCA，期望是不同的批次batch不分开
pdf(file="2.0.3.before_integrated.PCA.sample.pdf",width= 15, height = 10)
PCAPlot(before_seurat_integrated,group.by = "sample") ####分开
dev.off()
#################################    Run UMAP
set.seed(123)
before_seurat_integrated <- RunUMAP(before_seurat_integrated,dims = 1:30,reduction = "pca")
######plot Umap，整体画
pdf(file="2.0.5.before_integrated.umap.pdf",width= 15, height = 10)
DimPlot(before_seurat_integrated, label=T)
dev.off()
######按照sample分开画umap
pdf(file="2.0.6.before_integrated.UMP.sample.pdf", width= 30, height = 10)
DimPlot(before_seurat_integrated, group.by = "sample",label=T)
dev.off()
save(before_seurat_integrated,file="before_seurat_integrated.RData") ####保存未整合，但VST标准化后的



###########   整合分析.     VST标准化，再整合 #######

########https://satijalab.org/seurat/articles/integration_rpca.html

options(future.globals.maxSize = 8000 * 1024^2)  #####调整允许的对象大小的限制，因为太大了，所以改改默认设置,默认值为(500 * 1024 ^ 2 = 500 Mb）
library(patchwork)
table(filtered_seurat$sample)
ifnb.list <- SplitObject(filtered_seurat, split.by = "sample")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- CellCycleScoring(x, g2m.features = g2m_genes , s.features = s_genes )
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
ifnb.list

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures  =  2000) ##择选择用于数据整合的一些features，默认2000

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, vars.to.regress = c("mitoRatio","nUMI","S.Score","G2M.Score"), verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

#### 执行整合

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) ## ####找共同的锚点
table( filtered_seurat$sample)
#immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features,dims = 1:30, k.filter = 40)  ## ####k.filter = 40,细胞少，找共同的锚点 
seurat_integrated <- IntegrateData(anchorset = immune.anchors) #k.weight = 40 #####进行整合分析
dim(seurat_integrated)  ###整合后只有2000个基因

DefaultAssay(seurat_integrated)<- "integrated"
seurat_integrated <- ScaleData(seurat_integrated, features = VariableFeatures(seurat_integrated), verbose = T, vars.to.regress = c("mitoRatio","nUMI","S.Score","G2M.Score"))
dim(seurat_integrated[["integrated"]]@scale.data)
save(seurat_integrated, file="integrated_seurat.RData")    ########### 保存seurat_integrated##saveRDS(seurat_integrated, "integrated_seurat.rds ")

################################  对整合后的数据进行下游的降维可视化  ################################
DefaultAssay(seurat_integrated) <-  "integrated"
################    run PCA
set.seed(123)
seurat_integrated <- RunPCA(object  =  seurat_integrated, verbose = FALSE)
##################   plot PCA
pdf(file="2.0.7.integrated.sample.PCA.pdf",width= 30, height = 10)
PCAPlot(seurat_integrated,split.by = "sample") ####分开
dev.off()

pdf(file="2.0.8.integrated.PCA.pdf",width= 15, height = 10)
DimPlot(seurat_integrated, reduction = "pca")  ####一起
dev.off()

###############      Run UMAP
set.seed(123)
seurat_integrated <- RunUMAP(seurat_integrated,dims = 1:30,reduction = "pca")
###### plot Umap，整体画 按照标本名标注
pdf(file="2.0.9.integrated.umap.pdf",width= 15, height = 10)
DimPlot(seurat_integrated,label=T)
dev.off()
######按照sample分开画umap
pdf(file="2.1.0.integrated.sample.UMP.pdf", width= 15, height = 10)
DimPlot(seurat_integrated, group.by = "sample") #
dev.off()

####### 参考教程简书：https://www.jianshu.com/p/78c650345c5f

############################Seurat可以使用VizDimReduction, DimPlot, 和DimHeatmap函数对PCA的结果进行可视化

###PCA主成分热图
pdf(file="2.1.1.integrated.PCA.Heatmap.pdf",width= 15, height = 20)
DimHeatmap(seurat_integrated,dims = 1:50,cells = 500,balanced = TRUE,ncol = 4)
dev.off()

#######################选择PCA降维后的维数，用于后续的分析

####PCA肘拐点图
pdf(file="2.1.2.ElbowPlot.pdf", width= 15, height = 10)
ElbowPlot(object = seurat_integrated, reduction = "pca",
          ndims = 50)
dev.off()
head(seurat_integrated@reductions$pca@feature.loadings)##PC储存在这里

##################使用JackStraw函数计算每个PC的P值的分布，显著的PC会有较低的p-value
seurat_integrated <- JackStraw(seurat_integrated,num.replicate = 100)# 重复一百次
seurat_integrated <- ScoreJackStraw(seurat_integrated, dims = 1:20)##选择20个PC
## 可视化每个PC的P value分布
pdf(file="2.1.3.JackStrawPlot.pdf",width= 15, height = 10)
JackStrawPlot(seurat_integrated, dims = 1:20)
dev.off()

######打印出每个PCA成分基因
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
##VizDimLoadings(seurat_integrated, dims = 1:2, reduction = "pca") ####前2个PCs基因点图


#####################  Cluster the cells 细胞聚类#################################

seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:30)  ####选择前30个PCs
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = seq(from=0,by=.1,length=10))

####seurat_clusters是最后一个值（这里是1.4）产生的值.#########使用FindClusters函数对数据进行优化，并包含一个分辨率参数，该参数设置下游集群的“granularity”，增加的值将导致更多的集群。将该参数设置在0.4-1.4之间，对于3K左右的单细胞数据集通常会得到良好的结果。对于较大的数据集，最佳分辨率通常会增加。可以使用Idents函数找到clusters。分类结果存在seurat_integrated@active.ident 中，也可以通过Idents(seurat_integrated)调出查看。

##############画不同resolution = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4))  时的分群的关系树状图,参考 clustree ： 聚类可视化利器https://www.jianshu.com/p/f997c2f41c48
library(clustree)
library(dplyr)
pdf("3.0.0.cluster_tree.pdf",width= 25, height = 20)
seurat_integrated@meta.data %>% select(starts_with("integrated_snn_res"))  %>% mutate(integrated_snn_res.0.0=0) %>% clustree( prefix = "integrated_snn_res.",layout = "sugiyama")   ####画全部的resolution
dev.off()
                               

####要选择开始的分辨率，我们通常会在范围的中间进行选择，例如0.6或0.8。我们将使用Idents()功能分配聚类的标识，以0.8的分辨率开始。
table(seurat_integrated@meta.data$integrated_snn_res.1)   #####################查看有多少个分群，每个群有多少个数

levels(seu_obj@active.ident) ##查看活跃的Idents
table(seurat_integrated$integrated_snn_res.0.2)
Idents(object = seurat_integrated) <- "integrated_snn_res.1"
length(levels(Idents(object = seurat_integrated))) #####################查看一共产生了多少个群

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
        dplyr::count(ident, sample) %>%
        tidyr::spread(ident, n)
write.table(n_cells,file="n_cells.xls",se= "\t", row.names =F , col.names= T, quote = F)   ##保存     
head(n_cells) #######按sample查看有多少个分群，每个群有多少个细胞

DefaultAssay(seurat_integrated)
DefaultAssay(seurat_integrated)<- "RNA"
all.genes <- rownames(seurat_integrated)
seurat_integrated <- ScaleData(seurat_integrated, features =all.genes, vars.to.regress = c("mitoRatio","nUMI","S.Score","G2M.Score"))  ####features = all.genes，是对所有基因标准化。#scale函数默认是只针对之前鉴定的HVGs进行标准化（版本3中默认得到2000个HVGs）#这样操作的结果中降维和聚类不会被影响，但是只对HVGs基因进行标准化，下面画的热图seurat_integrated[["RNA"]]@scale.data# 可能会受到全部基因中某些极值的影响。所以为了热图的结果，还是对所有基因进行归一化比较好
dim(seurat_integrated[["RNA"]]@scale.data) ##18450  1826
save(seurat_integrated, file="integrated_seurat.findcluster.RData")   ########### 保存findcluster后的seurat_integrated)

DefaultAssay(seurat_integrated)
DefaultAssay(seurat_integrated)<-  "integrated"
# Plot the UMAP 全部的UMAP

###################   cluster 后的 总的 Umap, 没有标签

pdf(file="3.0.9.umap.pdf", width = 15, height = 15)
pp = DimPlot(object = seurat_integrated, reduction = 'umap',pt.size = 0.8, label = TRUE, label.size = 5)  ###+ NoLegend() 不要侧面标签
pp = pp + theme(axis.title = element_text(size = 15),axis.text =  element_text(size = 15,family = 'sans'),legend.text = element_text(size = 15),axis.line = element_line(size = 0.8))
print(pp)
dev.off()


##########################   按照sample 分开画tsne/umap
dpi = 300
pdf(file="3.0.2.umap.sample.pdf", width = 15, height = 15)
#png(file="3.0.2.umap.sample.png", width = dpi*8, height = dpi*8, units = "px",res = dpi,type='cairo')
pp = DimPlot(object = seurat_integrated, reduction = 'umap',pt.size = 0.5, label = TRUE, label.size = 5, group.by= "sample")  ###+ NoLegend() 不要侧面标签
pp = pp + theme(axis.title = element_text(size = 15),axis.text =  element_text(size = 15,family = 'sans'),
                legend.text = element_text(size = 15),axis.line = element_line(size = 0.8))
print(pp)
dev.off()


############################  ######    探索细胞周期的TSNE
pdf("3.0.3.TSNE.Phase.pdf",width= 15, height = 10)
DimPlot(seurat_integrated, label = TRUE, label.size = 5, pt.size = 0.5)  + NoLegend()
dev.off()

############################探索其他指标，例如每个细胞的UMI和基因数量，S期和G2M期标记，以及通过UMAP进行的线粒体基因表达
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
pdf("3.0.4.TSNE.metrics.pdf",width= 15, height = 10)
FeaturePlot(seurat_integrated, reduction = "umap", features = metrics, sort.cell = TRUE,min.cutoff = 'q10',label = TRUE,label.size = 5, pt.size = 0.5)
dev.off()

###############################################     探索前 16PC在cluster中的比例
library(cowplot)
columns <- c(paste0("PC_", 1:16), "ident",  "UMAP_1", "UMAP_2")
pc_data <- FetchData(seurat_integrated,  vars = columns)   ############################## The FetchData() function just allows us to extract the data more easily. FetchData()提取数据
#########################Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, vars = c("ident", "UMAP_1", "UMAP_2"))  %>% group_by(ident) %>% summarise(x=mean(UMAP_1), y=mean(UMAP_2))
########   画前 16PC在cluster中的比例
pdf("3.0.5.umap.16PCA.pdf",width= 15, height = 10)
map(paste0("PC_", 1:16), function(pc){ ggplot(pc_data,  aes(UMAP_1, UMAP_2)) + geom_point(aes_string(color=pc), alpha = 0.7) + scale_color_gradient(guide = FALSE,  low = "grey90",  high = "blue")  + geom_text(data=umap_label,  aes(label=ident, x, y)) + ggtitle(pc)}) %>% plot_grid(plotlist = .)
dev.off()

#############################打印PCA的基因       
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)


########################### ##############   FeaturePlot探索细胞类型分布：   marker   gene，染色画图
## FeaturePlot()功能染色感兴趣的标记基因，要访问所有基因的表达水平，而不仅仅是3000个高度可变的基因，我们使用存储在RNA slot中的标准化计数数据。

seurat_integrated@assays ######查看有多少个个slot
DefaultAssay(seurat_integrated)   ###查看目前在哪个Assays，共有三个，RNA +SCT + integrated
DefaultAssay(seurat_integrated)<-  "RNA"
seurat_integrated  <- NormalizeData(seurat_integrated,verbose=FALSE)
seurat_integrated[["RNA"]]@data[1:6,1:5] #标准化后的结果在这

################  marker gene  图
Idents(seurat_integrated)  <- seurat_integrated$integrated_snn_res.0.9
new_order = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16")
Idents(seurat_integrated) <- factor(Idents(seurat_integrated),levels = new_order,labels =
 new_order)
markers = c("CD3D","KLRF1","CD79A","CD14","CD1C","TPSB2","ACTA2","VWF","KRT19","FXYD2")
#markers = c("CD14","LILRA4","CD1C","CD79A","VWF","CD3D","KLRF1","ASGR1","KRT19","FXYD2","ACTA2","KIT")
#markers = c("CD3D","CD79A","FGFBP2","GNLY","CSF1R","CLEC9A","CD1C","LILRA4","TPSB2","KRT19","TM4SF4","APOA1","VWF","ACTA2")

pdf(file="3.0.7.umap_marker.pdf", width = 16, height = 12)
pp_temp = FeaturePlot(object = seurat_integrated, features = markers,cols = c("lightgrey","#ff0000"), sort.cell = TRUE, label = TRUE, min.cutoff = 'q10', combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 18),axis.text= element_text(size = 18),
                        plot.title = element_text(family = 'sans',face='italic',size=20),
                        legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 3,legend = 'right')
print(pp)
dev.off()

############################################  marker gene泡泡图

pdf(file="3.0.8.marker_heatmap.pdf", width = 8, height = 7)
pp = DotPlot(seurat_integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) +labs(x='',y='') + 
    guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
    theme(axis.line = element_line(size = 0.6))
print(pp)
dev.off()

############################################  marker gene 小提琴图https://mp.weixin.qq.com/s/DKJnXEAjB0F7jimXQ2_bZw

#Idents(nCoV.integrated1) <- "celltype"
markers <- CaseMatch(markers, rownames(seurat_integrated))
#group_list <- factor(seurat_integrated.aggregate@meta.data$celltype,levels = c('IGG+ B cells','IGA+ B cells','GC B cells','Memory B cells','Naive B cells'))
#seurat_integrated@meta.data$celltype2 <- group_list 

pdf("3.0.8.marker_VlnPlot.pdf", width = 15,height = 15)
VlnPlot(seurat_integrated, features = markers, pt.size = 0, stack = T)+ NoLegend() + ggsci::scale_fill_jco()
dev.off()
 
###########   肝癌 Marker gene画图

#############################################      findallmarkers 差异分析,要用RNA里的原始数据做差异分析。

DefaultAssay(seurat_integrated) <-  "RNA"
annotations <- read.csv("/Node22Data/home/wangting/project/icc/single_cell_RNA_seq/cellranger_count/Rproj_N368_N725/annotation.csv") 
#cluster0_8.markers <- FindMarkers(seurat_integrated, ident.1 = 0, ident.2 = 8, min.pct = 0.25) %>% arrange(avg_logFC) 

get_findallmarker <- function(cluster){
FindMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
               by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
  }

###对所有群执行
#Idents(object = seurat_integrated) <- "integrated_snn_res.1"
table(Idents(object = seurat_integrated))###多少个群
markers <- map_dfr(c(0:23), get_findallmarker) 
 ##按照avg_logFC降序排
#markers <- findallmarkers[,c(3,4,5,6,7,1,2,8)] ###只单独做几个cluster
##markers <- FindAllMarkers(object = seurat_integrated,  only.pos = TRUE, logfc.threshold = 0.25) %>% arrange(avg_logFC) ###所有的cluster
write.table(markers,file="findallmarkers.xls",col.names=T, row.names=F, quote=F,sep='\t') ##保存marker gene##保存
#markers <- read.table(file="findallmarkers.xls", header=T , sep='\t', quote = "\"")
markers_p_0.05 <- filter(markers, p_val < 0.05) ###过滤出p_val < 0.05
write.table(markers_p_0.05,file="findallmarkers.markers_p_0.05.xls",col.names=T, row.names=F, quote=F,sep='\t') ##保存p_val < 0.05 marker gene

DefaultAssay(seurat_integrated)   ###查看目前在哪个Assays，共有三个，RNA +SCT + integrated
DefaultAssay(seurat_integrated)<-  "RNA"
top10 <- markers %>% group_by(cluster_id) %>% top_n(n = 10, wt = avg_log2FC) #desc() ,seurat 4.0 是avg_log2FC
write.table(top10,file="top10_findallmarkers.xls",col.names=T, row.names=T,quote=F,sep='\t') ##保存top10
dim(top10)

##########   每个cluster的前10的markerhttps://mp.weixin.qq.com/s/B2PI19sdNkiE_ZGeguqcuQ,https://www.jianshu.com/p/db844d570e5a  基因画差异基因热图，用的是RNA assay slot scale.data中的数据

dim(seurat_integrated[["RNA"]]@scale.data)
head(seurat_integrated[["RNA"]]@scale.data[,1:6] ) ###画之前查看scale数据，确保RNA 做了scale

### 改热图 https://www.it610.com/article/1211113345411813376.htm；https://mp.weixin.qq.com/s/B2PI19sdNkiE_ZGeguqcuQ

pdf("3.0.6.top10.allmarkergenes.DoHeatmap.pdf",width =30, height = 25)
#seurat_integrated$integrated_snn_res.0.6  <- factor(x = seurat_integrated$integrated_snn_res.0.6, levels = c('0','1','2','3','4','6','5','7','8')) #更改顺序
DoHeatmap(seurat_integrated, features = top10$gene, assay = "RNA")+ NoLegend()
dev.off()

######### 平均的热图
markers = c()
temp <- DotPlot(seurat_integrated, features = markers)
data= temp$data
data= data[,c(3,4,5)]
head(data)
a= spread(data=data, key=id, value=avg.exp.scaled) ##宽数据变长数据
head(a)
rownames(a) = a$features.plot; a$features.plot=NULL;head(a)
pheatmap(a, cluster_rows = F, cluster_cols = F, 
            show_rownames = T, show_colnames =F, 
            color =colorRampPalette(c("DarkBlue", "white","Orange1"))(100), 
            cellwidth = 30, cellheight = 10, fontsize = 10,
            border_color="white", #
            annotation_col=annotation_col, # annotation_colors = 
            filename="heatmap.2.pdf",width=10, height=20)

dpi=400
png(file="3.0.6.top10.allmarkergenes.DoHeatmap.png", width = dpi*15, height = dpi*5, units = "px",res = dpi,type='cairo')
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()
dev.off()

###top2小提琴图
top2 <- markers %>% group_by(cluster_id) %>% top_n(n = 2, wt = avg_log2FC)
pdf(file="top2_Vlnplot.pdf",width= 30, height = 20)
VlnPlot(seurat_integrated, features =top2$gene,pt.size = 0)+NoLegend()
dev.off()

###top2UMP图
pdf(file="top2_FeaturePlot.pdf",width= 30, height = 50)
FeaturePlot(seurat_integrated, reduction = "umap",features =top2$gene, sort.cell = TRUE, label = TRUE)+NoLegend()
dev.off()

############################# 统计cluster占比条形图

# Idents(object = seurat_integrated) <- "integrated_snn_res.1"
length(levels(Idents(object = seurat_integrated))) #####################查看一共产生了多少个群
table(seurat_integrated@meta.data$integrated_snn_res.1) 
nCoV.integrated <- seurat_integrated

DefaultAssay(nCoV.integrated) <- "RNA"
table(Idents(object =nCoV.integrated))         
head(Idents(object = nCoV.integrated)) #####查看每个细胞属于哪个cluster
nCoV.integrated[["cluster"]] <- Idents(object = nCoV.integrated) #####增加一列，每个细胞对应哪个cluster
###write.table(nCoV.integrated[["cluster"]],file="cell_cluster.txt",row.names=T,quote = TRUE)
big.cluster = nCoV.integrated@meta.data  #####提出meta.data
organ.summary = table(big.cluster$sample,big.cluster$cluster)   ###行为sample，列为对应属于的cluster
head(organ.summary)
class(organ.summary)#检查数据格式
organ.summary <- as.data.frame(unclass(organ.summary))  ##变成data.frame
class(organ.summary)
organ.summary$group <- NA ####是data.frame才能加，加分组信息，内容为空
organ.summary$group[which(str_detect(rownames(organ.summary), "N"))] <- "P" #####，匹配以N开头，则为P组
organ.summary$group[which(str_detect(rownames(organ.summary), "T"))] <- "T"  #####，匹配T，则为P组
organ.summary$group[which(str_detect(rownames(organ.summary), "B"))] <- "N"   #####，匹配P，则为P组
organ.summary$batch <- NA ####是data.frame才能加，加分组信息，内容为空
organ.summary$batch[which(str_detect(rownames(organ.summary), "368"))] <- "batch1" #####，匹配以N开头，则为P组
organ.summary$batch[which(str_detect(rownames(organ.summary), "725"))] <- "batch2"
organ.summary$batch[which(str_detect(rownames(organ.summary), "748"))] <- "batch3"
head(organ.summary) ### 列为cluster，行为sample
write.table(organ.summary,file = 'sample_cluster.xls',quote = FALSE, sep = '\t',row.names=T) ###保存文件，行为细胞名，列为对应属于的cluster  与前面的n_cells.xls一样的。

library(reshape2)
organ.summary = read.delim2("sample_cluster.xls",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t",row.names=1)
organ.summary$sample <- rownames(organ.summary)
organ.summary.dataframe = melt(organ.summary)   ######宽数据变成长数据
head(organ.summary.dataframe)
colnames(organ.summary.dataframe) = c('group','batch','sample','cluster','cell') ####修改列名为group','cluster','cell'
head(organ.summary.dataframe)
organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)
samples_name_new = c("B748","N368","T368","N725","T725")  
organ.summary.dataframe$sample = factor(organ.summary.dataframe$sample,labels = samples_name_new,levels = samples_name_new) ##更改固定样本名顺序
head(organ.summary.dataframe)
### new_order = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15")##更改固定cluster顺序## Idents(seurat_integrated) = factor(Idents(seurat_integrated),levels = new_order,labels = new_order)
# organ.summary.dataframe$cluster = factor(organ.summary.dataframe$cluster,levels = new_order,labels = new_order)  ###修改cluster列的顺序= new_order

######################################  挑选出我们的数据, 分cluster，看P T 的对比

organ.summary.dataframe_own <- organ.summary.dataframe %>% filter(.,sample %in% c('B748','N368','N725','T368','T725'))###挑出 'N368','N725','T368','T725'

#cols = c('#32b8ec','#60c3f0','#8ccdf1','#cae5f7','#92519c','#b878b0','#d7b1d2','#e7262a','#e94746','#eb666d','#ee838f','#f4abac','#fad9d9')

cols = c('#fad9d9','#32b8ec','#92519c')#粉、蓝、紫
#104E8B

sample_colors= c("#D6D6CEFF" ,"#FFB547FF", "#ADB17DFF" ,"#5B8FA8FF","#800000FF")
group_colors=c("#7AA6DCCC","#003C67CC","#A73030CC" ) #分组颜色
##### 按照分组展示cluster +coord_flip()翻转90度
dpi = 300
png(file="1_own_sample_cluster_group.png", width = dpi*18, height = dpi*9, units = "px",res = dpi,type='cairo')
ggplot(data=organ.summary.dataframe_own, aes(x=cluster, y=cell, fill=group)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of group per cluster (%)') +
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
    scale_fill_manual(values = group_colors)
dev.off()

#### 按照sample填充颜色
sample_colors = c('#32b8ec','#cae5f7','#d7b1d2','#f4abac','#e7262a')
png(file="1_own_sample_cluster_sample.png", width = dpi*18, height = dpi*9, units = "px",res = dpi,type='cairo')
ggplot(data=organ.summary.dataframe_own, aes(x=cluster, y=cell, fill=sample)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Fraction of sample per cluster (%)') +
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),
          legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) +
    scale_fill_manual(values = sample_colors)
dev.off()
table(Idents(seurat_integrated))

######################################细胞重命名

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
table(Idents(object =nCoV.integrated1)) ####################查看重名后的细胞种类个数
 # nCoV.integrated1@meta.data <- nCoV.integrated1@meta.data %>%  dplyr::rename(old_celltype = celltype)
nCoV.integrated1$celltype = Idents(nCoV.integrated1) 
table(nCoV.integrated1@meta.data$celltype)####################查看重名后的每种细胞种类个数
nCoV_groups = c('NK','T_cell','B_cell','Malignant','Cholangiocyte','Macrophage','Fibroblast','Dendritic','Endothelial','Mast') 

###修改细胞顺序
nCoV.integrated1$celltype = factor(nCoV.integrated1$celltype,levels = nCoV_groups,labels = nCoV_groups) ###factor()修改细胞顺序
Idents(nCoV.integrated1) = nCoV.integrated1$celltype 
celltype <- as.data.frame.matrix(table(nCoV.integrated1$celltype,nCoV.integrated1$group))####################查看重名后的每种细胞来源于NPT的个数，as.data.frame.matrix把table变data.frame
celltype$all <- rowSums(celltype);celltype <- celltype[,c(4,1:3)]
write.table(celltype,file="celltype_group.xls",quote = FALSE,sep = '\t', col.names=T, row.names=T) ##保存每种细胞的数量
##celltype <- read.csv(file="celltype_group.xls",header = TRUE, row.names=1, check.names = FALSE, sep='\t')
#saveRDS(nCoV.integrated1,file="seurat_labelled.rds") 
save(nCoV.integrated1,file="seurat_labelled.RData")####保存修改名、修改细胞顺序后的RDatat

################# 总的，画重命名后的UMAP

library("ggsci")
#sample_colors= pal_uchicago("light", alpha =1)(5)[2:6]
sample_colors= c("#D6D6CEFF" ,"#FFB547FF", "#ADB17DFF" ,"#5B8FA8FF","#800000FF")
group_colors=c("#7AA6DCCC","#003C67CC","#A73030CC" ) #分组颜色
### 10种 细胞的颜色
celltype_colors=c("#374E55FF", "#DF8F44FF" ,"#00A1D5FF" ,"#B24745FF", "#698B22", "#6A6599FF","#80796BFF","#c5a1fc","#91D1C2FF","#CD853F")

######### 平均的热图每种细胞3个marker gene； nCoV_groups = c('NK','T_cell','B_cell','Malignant','Cholangiocyte','Macrophage','Fibroblast','Dendritic','Endothelial','Mast')

markers = c("KLRF1","FCGR3A","FGFBP2","CD3D","CD3E","CD2","MS4A1","CD79A","IGHG4","EPCAM","KRT19","KRT13","PIGR","TM4SF4","ANXA4","CD14","C1QA","C1QC","ACTA2","COL18A1","SOD3","CD1C","CLEC9A","CLEC4A","ENG","VWF","PECAM1","TPSB2","TPSAB1","CPA3")
temp <- DotPlot(nCoV.integrated1, features = markers)
data= temp$data;data= data[,c(3,4,5)]
head(data)##平均表达值
a= spread(data=data, key=id, value=avg.exp.scaled) ##宽数据变长数据
head(a)
rownames(a) = a$features.plot; a$features.plot=NULL; head(a)

# 注释
annotation_col = data.frame(Cell_types = factor(colnames(a)));rownames(annotation_col) = annotation_col$Cell_types
head(annotation_col) ##行名，列都是celltype
# 注释颜色需要是个list
ann_colors = list( Cell_types= c(NK= "#374E55FF", T_cell= "#DF8F44FF" ,B_cell= "#00A1D5FF" ,  Malignant="#B24745FF", Cholangiocyte= "#698B22", Macrophage="#6A6599FF", Fibroblast="#80796BFF",Dendritic="#c5a1fc",Endothelial= "#91D1C2FF", Mast= "#CD853F")
);head(ann_colors)

#### pheatmap
pheatmap::pheatmap(a, cluster_rows = F, cluster_cols = F, 
            show_rownames = T, show_colnames =F, 
            color =colorRampPalette(c("DarkBlue", "white","Orange1"))(100), 
            cellwidth = 40, cellheight = 20, fontsize = 10,
            border_color="NA",#white", #
            annotation_col= annotation_col, annotation_colors = ann_colors ,
            filename= "/Node22Data/home/wangting/project/icc/single_cell_RNA_seq/bile_duct_normal_tissue_20210115/Rproject/merge_ICC_IBD/SCT.without_integrated/test/celltypes.marker.heatmap.2.pdf", width=9, height=10)


#### ComplexHeatmap https://blog.csdn.net/weixin_45822007/article/details/115587526;https://blog.csdn.net/kMD8d5R/article/details/81008892/

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
col <- celltype_colors;names(col) <- levels(factor(rownames(annotation_col)))
#col_fun = colorRamp2(c(-2, 0, 2), c("DarkBlue", "white","Orange1"))
#col_fun = clorRamp2(c(-2, 0, 2), brewer.pl(n=3,name='RdBu'")))
col_fun = colorRamp2(c (-2,0,2),  c("blue", "white", "orange"))
annotation = HeatmapAnnotation(
                 block = anno_block ( gp = gpar( fill = col ), # 框图设置填充色
                 labels = rownames(annotation_col), 
                 labels_gp = gpar(cex = 0.6, col = "white") ))
                 
pdf("/Node22Data/home/wangting/project/icc/single_cell_RNA_seq/bile_duct_normal_tissue_20210115/Rproject/merge_ICC_IBD/SCT.without_integrated/test/celltypes.marker.pdf", width=8, height=10)
p1 = Heatmap(a, name="Scaled Exp", color= col_fun,  cluster_rows = FALSE,cluster_columns = FALSE, show_column_names = FALSE, show_row_names = T, row_names_side = "left", top_annotation = annotation,  column_split = levels(factor(rownames(annotation_col))), column_title = NULL, )
p1
dev.off()

DefaultAssay(nCoV.integrated1)<-  "RNA"
#nCoV.integrated1  <- NormalizeData(nCoV.integrated1,verbose=FALSE)
nCoV.integrated1[["RNA"]]@data[1:6,1:5] #标准化后的结果在这

markers = c("KLRF1","CD3D","CD79A","KRT19","FXYD2","CD14","ACTA2","CD1C","VWF","TPSB2")
#### marker gene 染图
pdf(file="3.0.7.umap_marker.celltype.pdf", width = 16, height = 12)
pp_temp = FeaturePlot(object = nCoV.integrated1, features = markers,cols = c("lightgrey","#ff0000"), sort.cell = TRUE, label = TRUE, min.cutoff = 'q10', combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + theme(axis.title = element_text(size = 18), axis.text= element_text(size = 18),
plot.title = element_text(family = 'sans',face='italic',size=20),
legend.text = element_text(size = 20),legend.key.height = unit(1.8,"line"),
legend.key.width = unit(1.2,"line") ))
pp = CombinePlots(plots = plots,ncol = 3,legend = 'right')
print(pp)
dev.off()

############ 小提琴图 rev() 字符逆序

pdf("3.0.8.marker_VlnPlot.celltype.pdf", width =10, height = 10)
VlnPlot(nCoV.integrated1, features = rev(markers), pt.size = 0, stack = T)+ NoLegend() + scale_fill_manual(values= rev(celltype_colors))
dev.off()

####  marker gene 泡泡图 红灰点
pdf(file="3.0.8.marker_heatmap.celltype.pdf", width = 8, height = 7)
pp = DotPlot(nCoV.integrated1, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + labs(x='',y='') +
    guides(color = guide_colorbar(title = 'Scaled Expression'),size = guide_legend(title = 'Percent Expressed')) +
    theme(axis.line = element_line(size = 0.6))
print(pp)
dev.off()



################# 总的，画重命名后的UMAP

pdf(file="1-umap_all_1.pdf", width = 15, height = 15)
pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = TRUE, label.size = 6, repel = TRUE,combine = TRUE,pt.size = 1.2, cols= celltype_colors)
pp_temp = pp_temp + theme(axis.title = element_text(size = 17),axis.text = element_text(size = 17),strip.text = element_text(family = 'arial',face='plain',size=17),legend.text = element_text(size = 17),axis.line = element_line(size = 1),axis.ticks = element_line(size = 0.8),legend.key.height = unit(1.4,"line"))
print(pp_temp)
dev.off()

#################  group (癌症、癌旁) , 画重命名后的UMAP


pdf(file="1-umap-split-group_1.pdf", width =45, height = 15)
pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',pt.size =1.2, label = TRUE, label.size = 4,repel = T,  split.by = 'group', ncol = 3, combine = TRUE, cols= celltype_colors)
pp_temp = pp_temp + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22), strip.text = element_text(family = 'sans',face='plain',size=22),legend.text = element_text(size =17),axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
print(pp_temp)
dev.off()
  
##################################  按照sample 画 分开的UMAP
pdf(file="1-umap-split-sample.pdf", width =15, height = 15)
pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',pt.size = 1.2, label = F, group.by = 'sample',cols=sample_colors)
pp_temp = pp_temp + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22), strip.text = element_text(family = 'sans',face='plain',size=22),legend.text = element_text(size =17),axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
print(pp_temp)
dev.off()

##################################  按照batch 画 分开的UMAP
dpi=400
png(file="1-umap-split-batch.png", width = dpi*18, height = dpi*9, units = "px",res = dpi,type='cairo')
pp_temp = DimPlot(object = nCoV.integrated1, reduction = 'umap',label = TRUE, split.by = 'batch', label.size = 4, ncol = 2, repel = TRUE, combine = TRUE,pt.size = 1.2, cols=celltype_colors )
pp_temp = pp_temp + theme(axis.title = element_text(size = 22),axis.text = element_text(size = 22), strip.text = element_text(family = 'sans',face='plain',size=22),legend.text = element_text(size =17),axis.line = element_line(size = 1.2),axis.ticks = element_line(size = 1.2))
print(pp_temp)
dev.off()
  

############################## #### 统计所有细胞类型比例图   
library(ggplot2)
library(dplyr)
library(reshape2)
###生成all——cell.txt
all_cell <- (nCoV.integrated1@meta.data %>% select(c(sample, batch, group, celltype)) %>% group_by(sample,batch,group,celltype) %>% summarise(n=n()))
head(all_cell)
write.table(all_cell, file="all_cell.txt",quote = FALSE,sep = '\t', col.names=T, row.names=F)
############################################
organ.summary = read.delim2("all_cell.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
class(organ.summary)
head(organ.summary)
colnames(organ.summary)=c('sample','batch','group','celltype','cell')
head(organ.summary)          ####所有画图需要的长数据
organ.summary$cell = as.numeric(organ.summary$cell)  ######转换成数字


#################################### 按照sample 选出我们的数据 每个sample 的所有的细胞类型比例图
organ.summary.dataframe_own <- organ.summary %>% filter(.,sample %in% c('N368','N725','T368','T725','B748'))###挑出 'N368','N725','T368','T725'
organ.summary.dataframe_own$celltype <- factor(organ.summary.dataframe_own$celltype,levels = c('NK','T_cell','B_cell','Malignant','Cholangiocyte','Macrophage','Fibroblast','Dendritic','Endothelial','Mast'))   ####固定细胞顺序
organ.summary.dataframe_own$sample <- factor(organ.summary.dataframe_own$sample,levels = c("B748","N368","T368","N725","T725"))###固定样本顺序
head(organ.summary.dataframe_own)

#每个sample 的所有的细胞类型比例图 Stacked barplot with multiple groups
pdf("1-own_sample-percentage.pdf", width = 10, height = 10) 
ggplot(data=organ.summary.dataframe_own, aes(x=sample, y=cell, fill= celltype)) + geom_bar(stat="identity",width = 0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Cell type (Percentage)') + 
    theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'),legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) + 
    scale_fill_manual(values = celltype_colors) 
dev.off()


############## ############################################################################## 我们的数据，过滤出需要的细胞类型（免疫细胞），按照分组 sample 或者按照分组 P T 展示 展示######https://github.com/zhangzlab/covid_balf/blob/master/total_fig.R ##第 120行代码

organ.summary.dataframe_own_immue = organ.summary %>% filter(.,celltype %in% c('T_cell','NK','B_cell','Macrophage','Dendritic','Mast'))  ### organ.summary.dataframe_own  的  cluster这一列只剩下过滤出来的免疫细胞类型了
####免疫细胞绘制比例图
  
 library(ggplot2)
 nCoV.integrated1[["cluster"]] <- Idents(object = nCoV.integrated1)
 big.cluster = nCoV.integrated1@meta.data
 organ.summary = table(big.cluster$group,big.cluster$cluster)
 head(organ.summary)
 write.table(organ.summary,file = '1-nCoV-percentage-group.txt',quote = FALSE,sep = '\t')
 
 library(ggplot2)
 library(dplyr)
 library(reshape2)
 organ.summary = read.delim2("1-nCoV-percentage-group.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
 organ.summary$group = rownames(organ.summary)
 organ.summary.dataframe = melt(organ.summary) #长数据
 colnames(organ.summary.dataframe) = c('group','cluster','cell') ###列名等于
 organ.summary.dataframe = organ.summary.dataframe %>% filter(.,cluster %in% c('T_cell','NK','B_cell','Macrophage','Dendritic','Mast'))
 organ.summary.dataframe$cell = as.numeric(organ.summary.dataframe$cell)
 
 #### 不同NPT组中免疫细胞比例比例图

pdf(file="1-nCoV-group-percentage.pdf", width = 30, height =15)
# Stacked barplot with multiple groups
pp = ggplot(data=organ.summary.dataframe, aes(x=cluster,y=cell, fill=group)) + geom_bar(stat="identity",width =0.6,position=position_fill(reverse = TRUE),size = 0.3,colour = '#222222') + labs(x='',y='Ratio of immune cells (%)') + theme(axis.title.x = element_blank(),axis.text = element_text(size = 16),axis.title.y =element_text(size = 16), legend.text = element_text(size = 15),legend.key.height = unit(5,'mm'), legend.title = element_blank(),panel.grid.minor = element_blank()) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(size = 10)) + scale_fill_manual(values = group_colors)
 print(pp)
 dev.off()

conda activate R.4.0

library(ggpubr)  ###ggpubr包中的两个函数: compare_means(): 在R.4.0上跑
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)

sample.percent = read.delim2("1-nCoV-percentage-sample-percent.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE,  sep = "\t") ###各种细胞的百分比表
sample.percent
####
     T_cell    NK B_cell Macrophage Dendritic  Mast
B748  0.764 0.051  0.013      0.109         0 0.064
N368  0.605 0.307  0.034      0.045     0.009     0
N725  0.572 0.371  0.015      0.042         0 0.001
T368  0.833 0.037  0.067      0.056         0 0.006
T725   0.75 0.028  0.017      0.171     0.004 0.029
########
sample.percent$sample <- rownames(sample.percent) ####修改名,新名=旧名
sample.percent$group <-NA ##增加一列分组P、T
library(stringr)
sample.percent$group[which(str_detect(sample.percent$sample, "^N"))] <- "P" #匹配
sample.percent$group[which(str_detect(sample.percent$sample, "T"))] <- "T"
sample.percent$group[which(str_detect(sample.percent$sample, "B"))] <- "N"

sample.percent$batch <-NA ##增加一列批次 batch 
sample.percent$batch[which(str_detect(sample.percent$sample, "^B"))] <- "batch3" #匹配
sample.percent$batch[which(str_detect(sample.percent$sample, "368"))] <- "batch1"
sample.percent$batch[which(str_detect(sample.percent$sample, "725"))] <- "batch2"
sample.percent 
##### 行，为sample, 列为细胞类型， 分组
     T_cell    NK B_cell Macrophage Dendritic  Mast sample group  batch
B748  0.764 0.051  0.013      0.109         0 0.064   B748     N batch3
N368  0.605 0.307  0.034      0.045     0.009     0   N368     P batch1
N725  0.572 0.371  0.015      0.042         0 0.001   N725     P batch2
T368  0.833 0.037  0.067      0.056         0 0.006   T368     T batch1
T725   0.75 0.028  0.017      0.171     0.004 0.029   T725     T batch2
dim(sample.percent )
sample.percent <- sample.percent[,c(7,8,9,1:6)]
write.table(sample.percent,file = 'immune_cell_sample-percent.txt',row.names = TRUE,quote = FALSE,sep='\t')  ####自己加的，详细版的保存免疫细胞百分比
#sample.percent <- read.table(file = '1.immune_cell_sample-percent.txt',header=T,sep="\t")
##################################################################   按照分组P T 查看细胞分布情况
nCoV_groups = c('NK','T_cell','B_cell','Macrophage','Dendritic','Mast') #  
pplist = list()
for(group_ in nCoV_groups){
sample.percent_  = sample.percent %>% select(one_of(c('sample','group',group_)))
colnames(sample.percent_) = c('sample','group','percent') #group
sample.percent_$percent = as.numeric(sample.percent_$percent)
sample.percent_ <- sample.percent_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), #该样本中所有数值由小到大排列后第75%的数字
                                              lower = quantile(percent, 0.25), ##该样本中所有数值由小到大排列后第25%的数字
                                              mean = mean(percent), ##样本的均值
                                              median = median(percent))  ###数据的中位数
print(group_)
print(sample.percent_$median)
#sample.percent_ %>% arrange(group)  ####按group 排序，查看算出来的 mean/ median , upper  lower
pp1 = ggplot(sample.percent_,aes(x=group,y=percent)) + geom_jitter(shape = 21,aes(fill=group),width = 0.25) + stat_summary(fun.y=mean, geom="point", color="grey60") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 6),axis.title = element_text(size = 6),legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),plot.title = element_text(size = 6,face = 'plain'),legend.position = 'none') + labs(title = group_,y='Percentage') +
geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  0.25)
###perform statistics analysis
labely = max(sample.percent_$percent)
compare_means(percent ~ sample,  data = sample.percent_) #percent是数值型变量，sample是因子
my_comparisons <- list( c("P", "T")) ## P T 两两比较
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 1,method = "t.test")
  pplist[[group_]] = pp1
}
#pplist
pdf(file="group-percentage.pdf", width = 6, height = 3)
print(plot_grid(pplist[['T_cell']],pplist[['NK']],
                pplist[['B_cell']],pplist[['Macrophage']],
                pplist[['Dendritic']],pplist[['Mast']], ncol = 3, nrow = 2))
dev.off()
write.table(sample.percent_,file = '1-all-immunecell-mean.median.xls',row.names = TRUE,quote = FALSE,sep='\t') ##保存中位值结果

[1] "NK"
[1] 0.0510 0.3390 0.3390 0.0325 0.0325
Adding missing grouping variables: `group`
[1] "T_cell"
[1] 0.7640 0.5885 0.5885 0.7915 0.7915
Adding missing grouping variables: `group`
[1] "B_cell"
[1] 0.0130 0.0245 0.0245 0.0420 0.0420
Adding missing grouping variables: `group`
[1] "Macrophage"
[1] 0.1090 0.0435 0.0435 0.1135 0.1135
Adding missing grouping variables: `group`
[1] "Dendritic"
[1] 0.0000 0.0045 0.0045 0.0020 0.0020
Adding missing grouping variables: `group`
[1] "Mast"
[1] 0.0640 0.0005 0.0005 0.0175 0.0175
Adding missing grouping variables: `group`

### 不同cluster
print(plot_grid(pplist[['C0']],pplist[['C1']],
                pplist[['C2']],pplist[['C3']],
                pplist[['C4']],pplist[['C5']],
                pplist[['C6']],pplist[['C7']],
                pplist[['C8']],pplist[['C9']],
                pplist[['C10']],pplist[['C11']],
                pplist[['C12']],pplist[['C13']],
                pplist[['C14']],pplist[['C15']],
                ncol =4, nrow = 4))
dev.off()



 
## 保存不同组中的细胞百分比表格
organ.summary = read.delim2("1-nCoV-percentage-group.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
organ.summary1 = organ.summary %>% select(c('T_cell','NK','B_cell','Macrophage','Dendritic','Mast'))
organ.summary2 <- round(organ.summary1 / rowSums(organ.summary1) * 100, 2)
write.table(organ.summary2,file = '1-nCoV-percentage-group-percent.txt',row.names = TRUE,quote =FALSE,sep='\t')  
  ####### 保存不同sample中细胞百分比
library(ggplot2)
nCoV.integrated1[["cluster"]] <- Idents(object = nCoV.integrated1)
big.cluster = nCoV.integrated1@meta.data
organ.summary = as.data.frame.matrix(table(big.cluster$sample,big.cluster$cluster))
write.table(organ.summary,file = '1-nCoV-percentage-sample.txt',quote = FALSE,sep = '\t')
organ.summary1 = organ.summary %>% select(c('T_cell','NK','B_cell','Macrophage','Dendritic','Mast'))
organ.summary2 <- round(organ.summary1 / rowSums(organ.summary1),3) ##round四舍五入函数
write.table(organ.summary2,file = '1-nCoV-percentage-sample-percent.txt',row.names = TRUE,quote =FALSE,sep='\t')  ####保存免疫细胞百分比
  


######把过滤后的，重新分配到Seurat对象中，Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
dim(sce.big) #####过滤前
dim(filtered_seurat)####过滤后
#####画过滤后的细胞数
pdf("1.9.1.filterd.numbers.of.cells.pdf",width= 15, height = 10)
ggplot(filtered_seurat@meta.data,aes(x=sample, fill=sample)) + geom_bar() +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) +  ggtitle("Number of filtered cells")######细胞计数
dev.off()

save(filtered_seurat, file="seurat_filtered.RData")#############################################保存filtered_seurat，Create .RData object to load at any time保存数据
rm("sce.big") ##删除不用的变量，清内存
gc()

