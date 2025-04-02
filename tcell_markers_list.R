

# T细胞层级注释标记基因列表
t_cell_markers <- list(
  # 第一层级：T细胞总体识别
  "T_cells" = list(
    markers = c("CD3D", "CD3E", "CD3G", "CD2", "CD7"),
    description = "T细胞总体标记"
  ),
  
  # 第二层级：主要T细胞亚群
  "CD4_T" = list(
    markers = c("CD4", "CD3D", "CD3E", "CD3G"),
    description = "CD4+ T细胞",
    
    # 第三层级：CD4+ T细胞亚型
    subtypes = list(
      "Th" = list(
        markers = c("CD4", "CD3D", "CD3E", "CD3G"),
        description = "辅助T细胞",
        
        # 第四层级：辅助T细胞细分
        subtypes = list(
          "Th1" = list(
            markers = c("TBX21", "IFNG", "CXCR3", "IL12RB2", "STAT4"),
            description = "Th1细胞"
          ),
          "Th2" = list(
            markers = c("GATA3", "IL4", "IL13", "IL5", "CCR4", "STAT6"),
            description = "Th2细胞"
          ),
          "Th17" = list(
            markers = c("RORC", "IL17A", "IL17F", "IL23R", "CCR6", "STAT3"),
            description = "Th17细胞"
          ),
          "Tfh" = list(
            markers = c("CXCR5", "PDCD1", "BCL6", "ICOS", "IL21"),
            description = "滤泡辅助T细胞"
          ),
          "Th9" = list(
            markers = c("SPI1", "IL9", "IRF4"),
            description = "Th9细胞"
          ),
          "Th22" = list(
            markers = c("AHR", "IL22", "CCR10", "CCR4", "CCR6"),
            description = "Th22细胞"
          )
        )
      ),
      "Treg" = list(
        markers = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
        description = "调节性T细胞",
        
        # 调节性T细胞亚型
        subtypes = list(
          "tTreg" = list(
            markers = c("FOXP3", "IKZF2", "TNFRSF18"),
            description = "胸腺源性Treg"
          ),
          "pTreg" = list(
            markers = c("FOXP3", "RORC", "IL10"),
            description = "周围诱导Treg"
          ),
          "Fr_I" = list(
            markers = c("CD45RA", "IL2RA", "FOXP3"),
            description = "初始调节性T细胞"
          ),
          "Fr_II" = list(
            markers = c("CD45RO", "IL2RA", "FOXP3", "CTLA4", "ICOS"),
            description = "效应调节性T细胞"
          ),
          "Fr_III" = list(
            markers = c("CD45RO", "IL2RA", "FOXP3lo"),
            description = "非抑制性T细胞"
          )
        )
      )
    )
  ),
  
  "CD8_T" = list(
    markers = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G"),
    description = "CD8+ T细胞",
    
    # 第三层级：CD8+ T细胞亚型
    subtypes = list(
      "CTL" = list(
        markers = c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY"),
        description = "细胞毒性T细胞",
        
        # 细胞毒性T细胞亚型
        subtypes = list(
          "CTL_Eff" = list(
            markers = c("GZMB", "PRF1", "CX3CR1", "FGFBP2", "FCGR3A", "KLRG1"),
            description = "效应CTL"
          ),
          "CTL_Naive" = list(
            markers = c("CCR7", "SELL", "TCF7", "LEF1", "CD45RA", "IL7R"),
            description = "初始CTL"
          )
        )
      ),
      "CD8_Mem" = list(
        markers = c("CD44", "CD45RO", "CXCR3"),
        description = "CD8记忆T细胞",
        
        # CD8记忆T细胞亚型
        subtypes = list(
          "CD8_Tcm" = list(
            markers = c("CCR7", "SELL", "CD45RO", "CD27", "IL7R"),
            description = "中枢记忆T细胞"
          ),
          "CD8_Tem" = list(
            markers = c("CCR7lo", "CD45RO", "GZMB", "PRF1"),
            description = "效应记忆T细胞"
          ),
          "CD8_Trm" = list(
            markers = c("ITGAE", "ITGA1", "CXCR6", "CD69", "CD103"),
            description = "组织驻留记忆T细胞"
          ),
          "CD8_Tscm" = list(
            markers = c("CCR7", "CD45RA", "IL7R", "CD95", "CD122", "CXCR3"),
            description = "干细胞记忆T细胞"
          )
        )
      ),
      "CD8_Ex" = list(
        markers = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4", "ENTPD1"),
        description = "CD8耗竭T细胞"
      )
    )
  ),
  
  "gdT" = list(
    markers = c("TRGC1", "TRGC2", "TRDC", "CD3D", "CD3E", "CD3G"),
    description = "γδT细胞",
    
    # γδT细胞亚型
    subtypes = list(
      "Vd1" = list(
        markers = c("TRDV1", "TRGV9neg", "KLRC1"),
        description = "Vδ1+ γδT细胞"
      ),
      "Vd2" = list(
        markers = c("TRDV2", "TRGV9", "IL17A", "KLRC1"),
        description = "Vδ2+ γδT细胞"
      )
    )
  ),
  
  "MAIT" = list(
    markers = c("SLC4A10", "TRAV1-2", "CD3D", "CD3E", "CD3G", "KLRB1", "DPP4"),
    description = "黏膜相关不变性T细胞"
  ),
  
  "NKT" = list(
    markers = c("CD3D", "CD3E", "CD3G", "CD56", "KLRB1", "ZBTB16", "CD161"),
    description = "自然杀伤T细胞"
  ),
  
  # 细胞状态标记物
  "T_cell_states" = list(
    "Activation" = c("CD69", "HLA-DRA", "CD25", "CD38", "CD40LG", "TNFRSF4", "CD27"),
    "Proliferation" = c("MKI67", "PCNA", "TOP2A", "CDK1", "TYMS"),
    "Exhaustion" = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1"),
    "Cytotoxicity" = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "FASLG"),
    "Migration" = c("SELL", "CCR7", "CXCR3", "CCR5", "CCR4", "CCR6", "CXCR5", "CCR10")
  )
)

# 使用方法示例
# 获取CD8 T细胞标记物
# cd8_markers <- t_cell_markers$CD8_T$markers

# 获取Th17细胞标记物
# th17_markers <- t_cell_markers$CD4_T$subtypes$Th$subtypes$Th17$markers

# 获取耗竭状态标记物
# exhaustion_markers <- t_cell_markers$T_cell_states$Exhaustion


