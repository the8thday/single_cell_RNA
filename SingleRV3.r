# ## 1.human
# ref.BE = readRDS(file = './BlueprintEncodeData.rds')
# ref.BE$label.main %>% unique()
# ref.BE$label.fine %>% unique()
# 
# ref.Monaco = readRDS(file = './MonacoImmuneData.rds')
# ref.Monaco$label.main %>% unique()
# ref.Monaco$label.fine %>% unique()
# 
# ref.hpca = readRDS(file = './HumanPrimaryCellAtlasData.rds')
# ref.hpca$label.main %>% unique()
# 
# ref.diced = readRDS(file = './DatabaseImmuneCellExpressionData.rds')
# ref.diced$label.main %>% unique()
# ref.diced$label.fine %>% unique()
# 
# ref.nhd = readRDS(file = './NovershternHematopoieticData.rds')
# ref.nhd$label.main %>% unique()
# 
# ref_human = list(ref.BE,ref.Monaco,ref.diced,ref.hpca,ref.nhd)
# names(ref_human) = c("BE","Monaco","diced","hpca","nhd")
# save(ref_human,file = "./singleR_ref_human.Rdata")
# 
# ## 2.mouse
# ref.igd = readRDS(file = './ImmGenData.rds') #mouse
# ref.mrd = readRDS(file = './MouseRNAseqData.rds') #mouse
# 
# ref_mus = list(ref.igd,ref.mrd)
# names(ref_mus) = c("igd","mrd")
# save(ref_mus,file = "./singleR_ref_mouse.Rdata")
load("~/Downloads/singleR_ref_mouse.Rdata")

# 此文件为整理好的上述五个文件，也可以。。。
load("~/Downloads/singleR_ref_human.Rdata")

################# 自定义函数
singleR_vis = function(input_sce = sce_singleR,
                       seurat_Data = pbmc,
                       clusters = sce_singleR$integrated_snn_res.0.5,
                       ref = ref.nhd,
                       labels = ref.nhd$label.fine,
                       title = "nhd",legend.position ="none",
                       reduction = "tsne",
                       SetIdent = "integrated_snn_res.0.5",
                       label = TRUE, pt.size = 0.5,label.box = T,repel = T){
  hpca.fine <- SingleR(test = input_sce,
                       # assay.type.test = 1,
                       ref = ref,
                       assay.type.test = "logcounts", 
                       assay.type.ref = "logcounts",
                       clusters = clusters,
                       # method= "cluster",
                       labels = labels)
  # 分群注释
  new.cluster.ids <- hpca.fine$pruned.labels
  seurat_Data = SetIdent(seurat_Data, value = SetIdent)
  names(new.cluster.ids) <- levels(seurat_Data)
  
  # 细胞分群的重命名
  seurat_Data <- RenameIdents(seurat_Data, new.cluster.ids)
  p.nhc = DimPlot(seurat_Data, reduction = reduction, label = label, 
                  pt.size = pt.size,label.box = label.box,repel = repel) + 
   ggtitle(title)+ theme_bw() +
    theme(panel.grid=element_blank(), # 去网格线
          plot.title = element_text(size = 15,color="black",hjust = 0.5),
          axis.text.x = element_text(size = 12, color = 'black'),
          axis.text.y = element_text(size = 12, color = 'black'),
          axis.title.x = element_text(size = 12, color = 'black'),
          axis.title.y = element_text(size = 12, color = 'black'),
          axis.ticks = element_line(color = 'black', lineend = 'round'),
          legend.position = legend.position,
          legend.text = element_text(size = 12, color = 'black'),
          legend.title = element_text(size = 12, color = 'black'),
          panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
    scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
    scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))
  
  return(p.nhc)
}
## 运行
# library(SingleR)
# for (i in 1:5) {
#   p.tmp = singleR_vis(input_sce = sce_singleR,
#                       seurat_Data = seurat.data,
#                       clusters = sce_singleR$RNA_snn_res.0.5,
#                       SetIdent = "RNA_snn_res.0.5",
#                       title =names(ref_human)[i],
#                       ref = ref_human[[i]],
#                       labels = ref_human[[i]]$label.fine,
#                       reduction = "umap")
#   assign(paste0("p.",names(ref_human)[i]),p.tmp)
# }
# 
# p.singleR = (p0 | p.BE | p.Monaco) / (p.diced | p.hpca | p.nhd)
# p.singleR
# ggsave(p.singleR, filename = "Outplot/Step5.annotation_singleR.pdf",
#        width = 13,height = 9)
