#！/urs/bin/env Rscript
# Orchestrating Single-Cell Analysis with Bioconductor
# https://osca.bioconductor.org/index.html


# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
library(scRNAseq)
library(SingleCellExperiment)
library(scater)


# a_test ------------------------------------------------------------------

counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix) # must be a matrix object!

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
sce <- scater::logNormCounts(sce)
cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

sce <- SingleCellExperiment(assays = list(counts = counts_matrix),
                            colData = cell_metadata)# 还是一次性建好的好






