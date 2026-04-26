# Spatial Transcriptomics Workflow Templates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize `spatial/` into 10 topic `.qmd` template files + 4 shared helpers covering Visium / Visium HD / Xenium / MERFISH / CosMx and the spatial-specific analyses (deconvolution, SVG, domain, niche, cell talk, ST × scRNA integration).

**Architecture:** Each `.qmd` is a self-contained, runnable Quarto document following the same 7-section structure (intro / setup / data load / analysis / viz / hand-off / refs). Reusable cross-file logic lives in `singlecell_utiles.R` as namespaced helpers. Three legacy R/qmd files are merged into the new `spatial_workflow.qmd` and archived.

**Tech Stack:** R (Seurat v5, qs2, ggplot2), Quarto, key spatial packages (spacexr, CARD, SPOTlight, BayesSpace, hoodscanR, MISTy, CellChat v2, ReactomePA). Public data: `SeuratData::stxBrain`, `SeuratData::allen_reference`, 10x public Xenium / Visium HD samples.

**Repo conventions enforced everywhere:**
- Seurat v5: `JoinLayers()` before counts-consuming tools
- Serialization via `qs2`, never `qs`
- No `ggpubr`; use base R stats + `ggplot2::annotate()`
- Reactome included in any GO/KEGG enrichment
- Comparison direction (numerator/denominator) consistent across variable names, titles, comments
- Datasets, Results/, figures/ stay gitignored — never commit

**Conventions for this plan:**
- "Verification" for helpers = source + smoke test on tiny input.
- "Verification" for `.qmd` = syntax check via `knitr::purl()` + `parse()`. Full render is a manual final step the user runs separately.
- All commits are made from the repo root: `/Users/CongLiu/Library/CloudStorage/OneDrive-Personal/kintor/Script/single_cell_RNA`.

---

## Task 1: Add `LoadVisiumMulti()` helper

**Files:**
- Modify: `singlecell_utiles.R` (append at end)

- [ ] **Step 1: Append helper function**

Append to `singlecell_utiles.R`:

```r

#' Batch-load multi-sample 10x Visium data
#'
#' @param dirs Character vector of spaceranger output directories.
#' @param sample_names Character vector of sample names (same length as dirs).
#' @param filename H5 filename inside each dir. Default "filtered_feature_bc_matrix.h5".
#' @param mt_pattern Regex for mitochondrial genes. Default "^MT-" (human). Use "^mt-" for mouse.
#' @return A named list of Seurat objects, each with `orig.ident`, `sample`, `percent.mt` set.
LoadVisiumMulti <- function(dirs,
                            sample_names,
                            filename = "filtered_feature_bc_matrix.h5",
                            mt_pattern = "^MT-") {
  stopifnot(length(dirs) == length(sample_names))
  out <- vector("list", length(dirs))
  names(out) <- sample_names
  for (i in seq_along(dirs)) {
    obj <- Seurat::Load10X_Spatial(
      data.dir = dirs[i],
      filename = filename,
      assay = "Spatial",
      slice = sample_names[i],
      filter.matrix = TRUE
    )
    obj$orig.ident <- sample_names[i]
    obj$sample <- sample_names[i]
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
    out[[i]] <- obj
  }
  out
}
```

- [ ] **Step 2: Smoke test**

Run from repo root:

```bash
Rscript -e '
source("singlecell_utiles.R")
stopifnot(exists("LoadVisiumMulti"))
stopifnot(is.function(LoadVisiumMulti))
fmls <- names(formals(LoadVisiumMulti))
stopifnot(identical(fmls, c("dirs", "sample_names", "filename", "mt_pattern")))
cat("OK\n")
'
```

Expected output: `OK`

- [ ] **Step 3: Commit**

```bash
git add singlecell_utiles.R
git commit -m "add LoadVisiumMulti helper for multi-sample Visium loading"
```

---

## Task 2: Add `SpatialQCPlot()` helper

**Files:**
- Modify: `singlecell_utiles.R` (append at end)

- [ ] **Step 1: Append helper function**

Append to `singlecell_utiles.R`:

```r

#' Standard QC panel: VlnPlot stacked above SpatialFeaturePlot
#'
#' @param obj Seurat object with spatial assay.
#' @param features Character vector of QC feature names.
#'   Default c("nCount_Spatial", "nFeature_Spatial", "percent.mt").
#' @param pt.size.factor Passed to SpatialFeaturePlot.
#' @return A patchwork object combining the two panels.
SpatialQCPlot <- function(obj,
                          features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                          pt.size.factor = 1.6) {
  miss <- setdiff(features, c(colnames(obj@meta.data), rownames(obj)))
  if (length(miss)) {
    stop("Features not in object: ", paste(miss, collapse = ", "))
  }
  p_vln <- Seurat::VlnPlot(obj, features = features, pt.size = 0, ncol = length(features)) +
    patchwork::plot_layout(guides = "collect")
  p_sp <- Seurat::SpatialFeaturePlot(obj, features = features,
                                     pt.size.factor = pt.size.factor,
                                     ncol = length(features))
  p_vln / p_sp
}
```

- [ ] **Step 2: Smoke test**

```bash
Rscript -e '
source("singlecell_utiles.R")
stopifnot(exists("SpatialQCPlot"))
stopifnot(is.function(SpatialQCPlot))
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add singlecell_utiles.R
git commit -m "add SpatialQCPlot helper for combined VlnPlot + SpatialFeaturePlot QC panel"
```

---

## Task 3: Add `PlotDeconvProportion()` helper

**Files:**
- Modify: `singlecell_utiles.R` (append at end)

- [ ] **Step 1: Append helper function**

Append to `singlecell_utiles.R`:

```r

#' Plot deconvolution proportions: stacked bar (per cluster/region) + spatial overlay
#'
#' @param obj Seurat object.
#' @param decon_assay Assay name holding the cell-type proportion matrix
#'   (rows = cell types, columns = spots, values in [0,1]).
#' @param group_col meta.data column to summarise proportions by (e.g. seurat_clusters or Region).
#' @param top_n Show the top N cell types by overall mean proportion. Default 12.
#' @return A list with elements `bar` and `spatial` (patchwork of per-celltype maps).
PlotDeconvProportion <- function(obj, decon_assay, group_col, top_n = 12) {
  stopifnot(decon_assay %in% Seurat::Assays(obj))
  stopifnot(group_col %in% colnames(obj@meta.data))
  m <- Seurat::GetAssayData(obj, assay = decon_assay, slot = "data")
  m <- as.matrix(m)
  overall <- sort(rowMeans(m), decreasing = TRUE)
  keep <- names(overall)[seq_len(min(top_n, length(overall)))]
  m_keep <- m[keep, , drop = FALSE]

  grp <- obj@meta.data[[group_col]]
  agg <- t(apply(m_keep, 1, function(x) tapply(x, grp, mean)))
  agg_df <- as.data.frame(agg)
  agg_df$celltype <- rownames(agg_df)
  long <- tidyr::pivot_longer(agg_df, cols = -celltype,
                              names_to = "group", values_to = "prop")

  bar <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = prop, fill = celltype)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::labs(x = group_col, y = "Proportion", fill = "Cell type") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  Seurat::DefaultAssay(obj) <- decon_assay
  spatial <- Seurat::SpatialFeaturePlot(obj, features = keep, ncol = 4,
                                        pt.size.factor = 1.6, alpha = c(0.1, 1))

  list(bar = bar, spatial = spatial)
}
```

- [ ] **Step 2: Smoke test**

```bash
Rscript -e '
source("singlecell_utiles.R")
stopifnot(exists("PlotDeconvProportion"))
stopifnot(is.function(PlotDeconvProportion))
fmls <- names(formals(PlotDeconvProportion))
stopifnot(identical(fmls, c("obj", "decon_assay", "group_col", "top_n")))
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add singlecell_utiles.R
git commit -m "add PlotDeconvProportion helper for deconvolution result visualization"
```

---

## Task 4: Add `NeighborComposition()` helper

**Files:**
- Modify: `singlecell_utiles.R` (append at end)

- [ ] **Step 1: Append helper function**

Append to `singlecell_utiles.R`:

```r

#' For each spot, compute the k-NN neighborhood composition over a categorical column.
#'
#' Uses the spatial coordinates from obj@images[[1]] (first image slice).
#'
#' @param obj Seurat object with at least one spatial image slot.
#' @param group_col meta.data column with the categorical label per spot
#'   (e.g. cell-type assignment or seurat_clusters).
#' @param k Number of nearest neighbors. Default 6 (Visium hex grid).
#' @return A long-format data.frame with columns: spot, group, neighbor_group, fraction.
NeighborComposition <- function(obj, group_col, k = 6) {
  stopifnot(group_col %in% colnames(obj@meta.data))
  if (length(obj@images) == 0) {
    stop("No image slot found on Seurat object.")
  }
  coords <- Seurat::GetTissueCoordinates(obj@images[[1]])
  coords <- as.matrix(coords[, c("imagerow", "imagecol")])
  nn <- FNN::get.knn(coords, k = k)$nn.index
  labels <- obj@meta.data[[group_col]]
  spots <- rownames(coords)
  out <- vector("list", length(spots))
  for (i in seq_along(spots)) {
    nb_lab <- labels[nn[i, ]]
    tab <- table(nb_lab) / k
    out[[i]] <- data.frame(
      spot = spots[i],
      group = labels[i],
      neighbor_group = names(tab),
      fraction = as.numeric(tab),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
```

- [ ] **Step 2: Smoke test**

```bash
Rscript -e '
source("singlecell_utiles.R")
stopifnot(exists("NeighborComposition"))
stopifnot(is.function(NeighborComposition))
fmls <- names(formals(NeighborComposition))
stopifnot(identical(fmls, c("obj", "group_col", "k")))
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add singlecell_utiles.R
git commit -m "add NeighborComposition helper for k-NN spatial neighborhood composition"
```

---

## Task 5: Build `spatial_workflow.qmd` (end-to-end Visium template)

This file consolidates the QC / SCT / clustering / marker logic from `Spatial_Transcripts.qmd`, `ST_test.qmd`, and `10x_Genomics_Visium.R` into one canonical Visium workflow.

**Files:**
- Create: `spatial/spatial_workflow.qmd`

- [ ] **Step 1: Create file with header + setup section**

Create `spatial/spatial_workflow.qmd` with this content:

````markdown
---
title: "10x Visium 端到端工作流"
format: html
---

## 简介

本文件是 10x Genomics Visium 空转数据的端到端模板：从数据加载、QC、标准化、聚类、空间区域差异分析到 marker 可视化。后续的反卷积、空间可变基因、空间 domain、niche、空间通讯、ST × scRNA 整合等分析，分别在 `spatial/` 下的同名 `.qmd` 文件中展开。

适用平台：标准 Visium（55µm spot）。Visium HD 见 `visium_hd.qmd`，亚细胞分辨率成像平台见 `xenium_workflow.qmd` / `merfish_cosmx.qmd`。

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs2)
source(here::here("singlecell_utiles.R"))
```

## 数据加载

使用 `SeuratData::stxBrain` 作演示样本（小鼠脑前后冠状切片各一）。第一次运行需要 `InstallData("stxBrain")`。

```{r}
# InstallData("stxBrain")  # 首次运行解开
brain_a <- LoadData("stxBrain", type = "anterior1")
brain_p <- LoadData("stxBrain", type = "posterior1")
brain_a$sample <- "anterior"
brain_p$sample <- "posterior"
```

实际项目里把样本路径放进 `LoadVisiumMulti()`：

```{r}
#| eval: false
sample_dirs <- c("/path/to/sample1/outs",
                 "/path/to/sample2/outs")
sample_names <- c("sample1", "sample2")
obj_list <- LoadVisiumMulti(sample_dirs, sample_names, mt_pattern = "^mt-")
```

## QC

```{r}
brain_a[["percent.mt"]] <- PercentageFeatureSet(brain_a, pattern = "^mt-")
SpatialQCPlot(brain_a)
```

按 nCount / nFeature / 线粒体比例过滤。阈值需按组织类型调整：

```{r}
brain_a <- subset(brain_a,
                  subset = nCount_Spatial > 500 &
                           nFeature_Spatial > 200 &
                           percent.mt < 20)
```

## 标准化与降维

```{r}
brain_a <- SCTransform(brain_a, assay = "Spatial", verbose = FALSE)
brain_a <- RunPCA(brain_a, assay = "SCT", verbose = FALSE)
ElbowPlot(brain_a, ndims = 30)
```

## 聚类与可视化

```{r}
brain_a <- FindNeighbors(brain_a, reduction = "pca", dims = 1:20)
brain_a <- FindClusters(brain_a, resolution = 0.4, verbose = FALSE)
brain_a <- RunUMAP(brain_a, reduction = "pca", dims = 1:20)

p1 <- SpatialDimPlot(brain_a, label = TRUE, label.size = 4)
p2 <- DimPlot(brain_a, reduction = "umap", label = TRUE) + NoLegend()
p1 + p2
```

## 空间区域定义与差异分析

通过 H&E 形态学 + 聚类位置手工定义解剖区域，然后做区域间差异：

```{r}
brain_a$Region <- dplyr::case_when(
  brain_a$seurat_clusters %in% c("0", "3") ~ "Cortex",
  brain_a$seurat_clusters %in% c("1", "2") ~ "Striatum",
  TRUE ~ "Other"
)
SpatialDimPlot(brain_a, group.by = "Region")
```

区域间差异基因（Cortex vs Striatum，方向写明保持一致）：

```{r}
Idents(brain_a) <- "Region"
de_cortex_vs_striatum <- FindMarkers(brain_a,
                                     ident.1 = "Cortex",
                                     ident.2 = "Striatum",
                                     test.use = "wilcox")
head(de_cortex_vs_striatum, 10)
```

## Marker 基因

```{r}
Idents(brain_a) <- "seurat_clusters"
all_markers <- FindAllMarkers(brain_a,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
top5 <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(avg_log2FC, n = 5)
SpatialFeaturePlot(brain_a, features = head(top5$gene, 4))
```

## 多样本整合

```{r}
#| eval: false
obj_list <- list(anterior = brain_a, posterior = brain_p)
obj_list <- lapply(obj_list, SCTransform, assay = "Spatial", verbose = FALSE)
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
anchors  <- FindIntegrationAnchors(obj_list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
brain    <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
```

## 保存

```{r}
#| eval: false
qs2::qs_save(brain_a, "Results/brain_anterior.qs")
```

## 下游分析

- 反卷积： `spatial/spatial_deconvolution.qmd`
- 空间可变基因： `spatial/spatial_SVG.qmd`
- 空间 domain 检测： `spatial/spatial_domain.qmd`
- niche / 邻域： `spatial/spatial_niche.qmd`
- 空间通讯： `spatial/spatial_celltalk.qmd`
- ST × scRNA 整合： `spatial/spatial_scRNA_integration.qmd`

## 注意事项

- Seurat v5 多样本整合后调用反卷积 / CytoTRACE2 等之前需 `JoinLayers()`
- 序列化用 `qs2::qs_save`，不要用 `qs::`
- 区域差异分析比较方向（如 Cortex vs Striatum）变量名 / 标题 / 注释三者保持一致
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_workflow.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_workflow.qmd
git commit -m "add spatial_workflow.qmd: end-to-end Visium template (replaces 3 legacy files)"
```

---

## Task 6: Archive superseded files

**Files:**
- Move: `spatial/Spatial_Transcripts.qmd` → `archive/superseded_rmd/`
- Move: `spatial/ST_test.qmd` → `archive/superseded_rmd/`
- Move: `spatial/10x_Genomics_Visium.R` → `archive/superseded_rmd/`
- Delete: `spatial/spacexr_STalign.qmd`

- [ ] **Step 1: Verify destination exists**

```bash
ls archive/superseded_rmd/ | head -5
```

Expected: directory exists with at least some files.

- [ ] **Step 2: Move three superseded files**

```bash
git mv spatial/Spatial_Transcripts.qmd archive/superseded_rmd/Spatial_Transcripts.qmd
git mv spatial/ST_test.qmd archive/superseded_rmd/ST_test.qmd
git mv spatial/10x_Genomics_Visium.R archive/superseded_rmd/10x_Genomics_Visium.R
```

- [ ] **Step 3: Delete the placeholder file**

```bash
git rm spatial/spacexr_STalign.qmd
```

- [ ] **Step 4: Verify**

```bash
ls spatial/
```

Expected: only `spatial_workflow.qmd` and `spatial_scRNA_workflow.ipynb` are present.

- [ ] **Step 5: Commit**

```bash
git commit -m "archive superseded spatial files; folded into spatial_workflow.qmd"
```

---

## Task 7: Build `spatial_SVG.qmd`

**Files:**
- Create: `spatial/spatial_SVG.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_SVG.qmd`:

````markdown
---
title: "空间可变基因 (SVG)"
format: html
---

## 简介

空间可变基因（spatially variable genes, SVG）是指表达模式在组织空间上**非随机**的基因——既包括聚集分布的，也包括沿梯度变化的。比单纯的高变基因（HVG）更能捕捉空间生物学。

本文件覆盖三种主流方法：
- Seurat `FindSpatiallyVariableFeatures`（Moran's I / `markvariogram`）
- SPARK-X（基于核的快速 SVG 测试，适合大样本）
- `Voyager` 包计算的 Moran's I（更细粒度的空间统计）

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SPARK)        # SPARK-X
library(Voyager)
library(SpatialFeatureExperiment)
library(clusterProfiler)
library(ReactomePA)   # 富集必须含 Reactome
source(here::here("singlecell_utiles.R"))
```

## 数据加载

```{r}
brain <- LoadData("stxBrain", type = "anterior1")
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")
brain <- subset(brain, subset = nFeature_Spatial > 200 & percent.mt < 20)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
```

## 方法 1：Seurat `FindSpatiallyVariableFeatures`

```{r}
brain <- FindSpatiallyVariableFeatures(
  brain,
  assay = "SCT",
  features = VariableFeatures(brain)[1:1000],
  selection.method = "moransi"
)
top_svg <- SpatiallyVariableFeatures(brain, selection.method = "moransi") |> head(12)
SpatialFeaturePlot(brain, features = top_svg[1:6], ncol = 3)
```

## 方法 2：SPARK-X

SPARK-X 在大样本下显著快于 SPARK / Moran's I。

```{r}
counts <- GetAssayData(brain, assay = "Spatial", slot = "counts")
coords <- GetTissueCoordinates(brain@images[[1]])
coords <- as.matrix(coords[, c("imagerow", "imagecol")])
sparkx_res <- sparkx(counts, coords, numCores = 2, option = "mixture")
sparkx_top <- rownames(sparkx_res$res_mtest)[order(sparkx_res$res_mtest$adjustedPval)][1:20]
SpatialFeaturePlot(brain, features = sparkx_top[1:6], ncol = 3)
```

## 方法 3：Voyager / Moran's I

```{r}
#| eval: false
sfe <- SpatialFeatureExperiment::toSpatialFeatureExperiment(brain)
sfe <- Voyager::runMoransI(sfe, features = rownames(sfe)[1:500])
```

## SVG 富集分析（含 Reactome）

```{r}
svg_genes <- top_svg
ego <- enrichGO(gene = svg_genes,
                OrgDb = "org.Mm.eg.db",
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH")
ekegg <- enrichKEGG(gene = clusterProfiler::bitr(svg_genes,
                                                 fromType = "SYMBOL",
                                                 toType = "ENTREZID",
                                                 OrgDb = "org.Mm.eg.db")$ENTREZID,
                    organism = "mmu")
ereact <- ReactomePA::enrichPathway(
  gene = clusterProfiler::bitr(svg_genes,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = "org.Mm.eg.db")$ENTREZID,
  organism = "mouse"
)
dotplot(ereact, showCategory = 15) + ggtitle("Reactome enrichment of SVGs")
```

## 三种方法对比

```{r}
seurat_svg <- top_svg
sparkx_svg <- sparkx_top[1:12]
overlap <- intersect(seurat_svg, sparkx_svg)
cat("Overlap:", length(overlap), "/", length(union(seurat_svg, sparkx_svg)), "\n")
```

## 注意事项

- SPARK-X 比 Seurat Moran's I 快很多，大样本（>5k spots）首选
- SVG 不等于 marker：marker 是 cluster 间高表达，SVG 是空间结构非随机
- 富集中必须包含 Reactome（repo 约定）

## 下游

- 找到 SVG 后可作为 BayesSpace / SpaGCN 的输入特征： `spatial/spatial_domain.qmd`
- SVG 模块共表达分析可对接 `hdWGCNA.qmd`（根目录）
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_SVG.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_SVG.qmd
git commit -m "add spatial_SVG.qmd: SVG detection with Seurat, SPARK-X, Voyager + Reactome enrichment"
```

---

## Task 8: Build `spatial_deconvolution.qmd`

**Files:**
- Create: `spatial/spatial_deconvolution.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_deconvolution.qmd`:

````markdown
---
title: "Spot 反卷积"
format: html
---

## 简介

Visium 每个 spot 直径约 55µm，包含多个细胞，因此需要用 scRNA 参考做反卷积，估计每个 spot 的细胞类型组成。

本文件涵盖三种工具：
- **RCTD（spacexr）**：泊松似然 + 双类型混合，鲁棒，目前最常用
- **CARD**：考虑空间相关性的非负回归
- **SPOTlight**：非负矩阵分解（NMF）

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(spacexr)     # RCTD
library(CARD)
library(SPOTlight)
library(qs2)
library(dplyr)
library(ggplot2)
source(here::here("singlecell_utiles.R"))
```

## 数据加载（ST + 参考 scRNA）

```{r}
# 空转
brain <- LoadData("stxBrain", type = "anterior1")
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")
brain <- subset(brain, subset = nFeature_Spatial > 200 & percent.mt < 20)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

# 参考 scRNA（Allen 脑图谱）
allen <- LoadData("allen_reference")
allen <- UpdateSeuratObject(allen)
# Seurat v5：JoinLayers 再取 counts
allen <- JoinLayers(allen)
```

## 方法 1：RCTD (spacexr)

```{r}
ref_counts <- GetAssayData(allen, assay = "RNA", slot = "counts")
ref_meta <- as.factor(allen$subclass)
names(ref_meta) <- colnames(allen)
reference <- Reference(ref_counts, ref_meta)

st_counts <- GetAssayData(brain, assay = "Spatial", slot = "counts")
st_coords <- GetTissueCoordinates(brain@images[[1]])[, c("imagerow", "imagecol")]
puck <- SpatialRNA(st_coords, st_counts)

rctd <- create.RCTD(puck, reference, max_cores = 2)
rctd <- run.RCTD(rctd, doublet_mode = "doublet")

# 提取权重并写入 Seurat
weights <- as.matrix(rctd@results$weights)
norm_w <- sweep(weights, 1, rowSums(weights), "/")
brain[["RCTD"]] <- CreateAssayObject(data = t(norm_w))

# 可视化
DefaultAssay(brain) <- "RCTD"
SpatialFeaturePlot(brain, features = head(rownames(brain[["RCTD"]]), 6),
                   ncol = 3, alpha = c(0.1, 1))
```

## 方法 2：CARD

```{r}
#| eval: false
card_obj <- createCARDObject(
  sc_count = ref_counts,
  sc_meta = data.frame(cellType = ref_meta, sampleInfo = "all"),
  spatial_count = st_counts,
  spatial_location = data.frame(x = st_coords[, 1], y = st_coords[, 2]),
  ct.varname = "cellType",
  ct.select = unique(ref_meta),
  sample.varname = "sampleInfo"
)
card_obj <- CARD_deconvolution(CARD_object = card_obj)
CARD.visualize.pie(proportion = card_obj@Proportion_CARD,
                   spatial_location = card_obj@spatial_location)
```

## 方法 3：SPOTlight

```{r}
#| eval: false
hvg <- VariableFeatures(allen)[1:3000]
mgs <- FindAllMarkers(allen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mgs_top <- mgs |> group_by(cluster) |> slice_max(avg_log2FC, n = 100)
spotlight_res <- SPOTlight(
  x = GetAssayData(allen, slot = "counts"),
  y = GetAssayData(brain, assay = "Spatial", slot = "counts"),
  groups = allen$subclass,
  mgs = mgs_top,
  hvg = hvg,
  weight_id = "avg_log2FC",
  group_id = "cluster",
  gene_id = "gene"
)
```

## 三种方法的可视化对比

```{r}
res <- PlotDeconvProportion(brain, decon_assay = "RCTD",
                            group_col = "seurat_clusters", top_n = 10)
res$bar
res$spatial
```

## 与原 `spacexr_STalign` 文件的关系

`spacexr` 即 RCTD 包；本文件已涵盖。原 `spatial/spacexr_STalign.qmd` 占位文件已删除。

## 注意事项

- Seurat v5 参考对象在 `Reference()` 之前必须 `JoinLayers()`
- `doublet_mode = "doublet"` 适合标准 Visium；高分辨率（Slide-seq、HD）用 `"full"`
- 序列化反卷积结果用 `qs2::qs_save`

## 下游

- niche / 邻域分析（基于反卷积比例）： `spatial/spatial_niche.qmd`
- 空间通讯（按反卷积细胞类型）： `spatial/spatial_celltalk.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_deconvolution.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_deconvolution.qmd
git commit -m "add spatial_deconvolution.qmd: RCTD + CARD + SPOTlight templates"
```

---

## Task 9: Build `spatial_domain.qmd` (skeleton for SpaGCN)

**Files:**
- Create: `spatial/spatial_domain.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_domain.qmd`:

````markdown
---
title: "空间 Domain 检测"
format: html
---

## 简介

"空间 domain" 是指组织上空间连续、转录组相似的区域。和普通聚类的区别：domain 检测**显式利用空间坐标**，结果更连续、更接近解剖区域。

主流工具：
- **BayesSpace**（R，可跑）— 在 PC 空间做 Markov 随机场聚类，利用 hex 邻接
- **SpaGCN**（Python，本文件仅注释，不执行）— 图卷积网络融合表达 + 坐标 + 组织图像

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
source(here::here("singlecell_utiles.R"))
```

## 数据加载

```{r}
brain <- LoadData("stxBrain", type = "anterior1")
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")
brain <- subset(brain, subset = nFeature_Spatial > 200 & percent.mt < 20)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
```

## 方法 1：BayesSpace（R，可跑）

转换为 SingleCellExperiment 并加上 spot 坐标：

```{r}
sce <- as.SingleCellExperiment(brain, assay = "Spatial")
coords <- GetTissueCoordinates(brain@images[[1]])
sce$row <- coords$imagerow
sce$col <- coords$imagecol
metadata(sce)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)
```

预处理 + 选择簇数：

```{r}
sce <- spatialPreprocess(sce, platform = "Visium",
                         n.PCs = 15, n.HVGs = 2000, log.normalize = FALSE)
sce <- qTune(sce, qs = seq(2, 12), platform = "Visium", d = 15)
qPlot(sce)
```

跑聚类：

```{r}
q <- 8  # 选 elbow / 业务知识
sce <- spatialCluster(sce, q = q, platform = "Visium", d = 15,
                      init.method = "mclust", model = "t",
                      gamma = 3, nrep = 1000, burn.in = 100)
clusterPlot(sce)
```

将 domain 写回 Seurat：

```{r}
brain$bayes_domain <- factor(sce$spatial.cluster)
SpatialDimPlot(brain, group.by = "bayes_domain", label = TRUE)
```

## 方法 2：SpaGCN（Python，仅参考）

SpaGCN 是 Python 工具，本仓库不标准化 Python 环境，因此本节**仅给出调用模板**，不在 quarto 渲染时执行。

```python
# python 环境（参考）
# pip install SpaGCN scanpy
import scanpy as sc
import SpaGCN as spg

adata = sc.read_visium("/path/to/spaceranger_outs")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 构建邻接图：表达 + 坐标 + 组织图像
adj = spg.calculate_adj_matrix(
    x=adata.obs["array_row"], y=adata.obs["array_col"],
    x_pixel=adata.obs["pxl_col_in_fullres"],
    y_pixel=adata.obs["pxl_row_in_fullres"],
    image=adata.uns["spatial"]["...your_sample..."]["images"]["hires"],
    beta=49, alpha=1, histology=True
)

clf = spg.SpaGCN()
clf.set_l(0.5)
clf.train(adata, adj, init_spa=True, init="louvain", res=0.4)
y_pred, _ = clf.predict()
adata.obs["spagcn_domain"] = y_pred
```

输出的 `spagcn_domain` 可通过 `Seurat_2_h5ad.qmd` 流程往返到 Seurat。

## BayesSpace vs SpaGCN

- BayesSpace：纯统计模型，结果稳定、可解释
- SpaGCN：能利用 H&E 图像、对噪声 spot 更鲁棒，但调参更敏感
- 实践上常先 BayesSpace，结果不理想或想用图像信息时换 SpaGCN

## 注意事项

- BayesSpace 用 `nrep = 50000` 才是 paper 推荐设置；本模板为加速用 1000
- SpaGCN 章节本仓库不执行；如需运行，自建 Python venv

## 下游

- domain 间差异分析等同于 cluster 间差异（见 `spatial_workflow.qmd`）
- domain 内的 niche 分析： `spatial_niche.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_domain.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_domain.qmd
git commit -m "add spatial_domain.qmd: BayesSpace runnable + SpaGCN reference"
```

---

## Task 10: Build `spatial_niche.qmd`

**Files:**
- Create: `spatial/spatial_niche.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_niche.qmd`:

````markdown
---
title: "Niche / 邻域分析"
format: html
---

## 简介

Niche / 邻域分析关心的是"哪些细胞类型 / domain 在空间上**共定位**"。本文件三种角度：
- **邻接共定位矩阵**（自定义 + helper）
- **hoodscanR**：基于 KNN 的邻域统计 + 富集 vs 随机
- **MISTy**：用 ML 模型把局部 / 中度 / 远距离邻域的影响拆开

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(hoodscanR)
library(SpatialExperiment)
library(mistyR)
library(future)
library(dplyr)
library(ggplot2)
source(here::here("singlecell_utiles.R"))
```

## 数据加载（已带反卷积结果）

假设 `spatial_deconvolution.qmd` 已跑完，加载结果：

```{r}
#| eval: false
brain <- qs2::qs_read("Results/brain_with_RCTD.qs")
DefaultAssay(brain) <- "RCTD"
```

为了独立可跑，这里 mock 一个 dominant cell type 标签：

```{r}
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(resolution = 0.4, verbose = FALSE)
brain$celltype <- paste0("CT", brain$seurat_clusters)
```

## 方法 1：邻接共定位（helper）

```{r}
nb <- NeighborComposition(brain, group_col = "celltype", k = 6)
agg <- nb |>
  dplyr::group_by(group, neighbor_group) |>
  dplyr::summarise(mean_frac = mean(fraction), .groups = "drop")
ggplot(agg, aes(x = neighbor_group, y = group, fill = mean_frac)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Neighbor cell type", y = "Center cell type", fill = "Mean fraction")
```

## 方法 2：hoodscanR

```{r}
spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = GetAssayData(brain, assay = "Spatial", slot = "counts")),
  colData = brain@meta.data,
  spatialCoords = as.matrix(GetTissueCoordinates(brain@images[[1]])[, c("imagerow", "imagecol")])
)
spe <- hoodscanR::readHoodData(spe, anno_col = "celltype")
fnc  <- hoodscanR::findNearCells(spe, k = 100)
pm   <- hoodscanR::scanHoods(fnc$distance)
hoods <- hoodscanR::mergeByGroup(pm, fnc$cells)
spe   <- hoodscanR::mergeHoodSpe(spe, hoods)
hoodscanR::plotHoodMat(spe, n = 20)
```

## 方法 3：MISTy（多视角邻域回归）

```{r}
#| eval: false
plan(multisession, workers = 2)
expr <- as.matrix(GetAssayData(brain, assay = "SCT", slot = "data"))
markers <- VariableFeatures(brain)[1:50]
expr_sub <- t(expr[markers, ])

coords_df <- as.data.frame(GetTissueCoordinates(brain@images[[1]])[, c("imagerow", "imagecol")])

views <- create_initial_view(as.data.frame(expr_sub)) |>
  add_juxtaview(positions = coords_df, neighbor.thr = 200) |>
  add_paraview(positions = coords_df, l = 1000)

run_misty(views, results.folder = "Results/misty_brain")
misty_results <- collect_results("Results/misty_brain")
plot_improvement_stats(misty_results)
plot_view_contributions(misty_results)
```

## 比较方向约定

本文件中每个矩阵的"row = center, column = neighbor"，所有标题与文字保持一致。

## 注意事项

- hoodscanR `k = 100` 是默认，组织小时调小
- MISTy 训练耗时长，示例只跑 50 个特征
- 反卷积比例输入 niche 分析需先 `JoinLayers()` 后再取 spot × cell-type 矩阵

## 下游

- niche × 通讯： `spatial_celltalk.qmd`
- niche × 区域 DE：参考 `spatial_workflow.qmd` Region 那节，把 Region 替换成 niche 标签
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_niche.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_niche.qmd
git commit -m "add spatial_niche.qmd: adjacency + hoodscanR + MISTy"
```

---

## Task 11: Build `spatial_celltalk.qmd`

**Files:**
- Create: `spatial/spatial_celltalk.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_celltalk.qmd`:

````markdown
---
title: "空间增强的细胞通讯"
format: html
---

## 简介

普通的细胞通讯分析（CellChat / NicheNet）不考虑空间，但实际上 ligand-receptor 信号通常发生在邻近 spot 间。本文件覆盖：
- **CellChat v2 spatial**（R，可跑）— 在 CellChat 中加入 spatial 距离约束
- **COMMOT**（Python，仅参考）— 用最优传输直接在空间上估计信号流

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(CellChat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
source(here::here("singlecell_utiles.R"))
```

## 数据准备

```{r}
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(resolution = 0.4, verbose = FALSE)
brain$celltype <- paste0("CT", brain$seurat_clusters)
```

## CellChat v2 spatial（可跑）

```{r}
data_input <- GetAssayData(brain, assay = "SCT", slot = "data")
meta <- data.frame(labels = brain$celltype, row.names = colnames(brain))

# 空间坐标
coords <- GetTissueCoordinates(brain@images[[1]])
spatial_locs <- as.matrix(coords[, c("imagerow", "imagecol")])

# Visium 的 spot 物理大小（µm）
scale_factors <- list(
  spot.diameter = 65,
  spot = 65,
  fiducial = 0,
  hires = brain@images[[1]]@scale.factors$hires,
  lowres = brain@images[[1]]@scale.factors$lowres
)

cellchat <- createCellChat(
  object = data_input,
  meta = meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = spatial_locs,
  scale.factors = scale_factors
)
cellchat@DB <- CellChatDB.mouse  # 小鼠数据
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,
                              type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE,
                              interaction.length = 200,
                              scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
```

可视化：

```{r}
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction weights")
```

## L-R pair 富集（含 Reactome）

```{r}
top_lr <- subsetCommunication(cellchat) |>
  arrange(desc(prob)) |>
  head(50)
genes_in_top <- unique(c(top_lr$ligand, top_lr$receptor))
entrez <- clusterProfiler::bitr(genes_in_top, fromType = "SYMBOL",
                                toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID
ereact <- ReactomePA::enrichPathway(gene = entrez, organism = "mouse")
dotplot(ereact, showCategory = 15) + ggtitle("Reactome enrichment of top L-R pairs")
```

## COMMOT（Python，仅参考）

```python
# 参考调用，本文件不执行
import scanpy as sc
import commot as ct

adata = sc.read_visium("/path/to/spaceranger_outs")
df_lr = ct.tl.ligand_receptor_database(species="mouse", signaling_type="Secreted Signaling")
ct.tl.spatial_communication(
    adata, database_name="cellchat", df_ligrec=df_lr,
    dis_thr=500, heteromeric=True, pathway_sum=True
)
ct.pl.plot_cell_communication(adata, database_name="cellchat",
                              pathway_name="WNT", plot_method="grid")
```

## 注意事项

- `interaction.length` 单位是像素（图像坐标）；按组织调整
- `distance.use = TRUE` 才启用空间约束，否则等价于普通 CellChat
- L-R 富集必须含 Reactome（repo 约定）

## 下游

- niche-by-niche 通讯比较：用 `spatial_niche.qmd` 的标签替换 `meta$labels`
- 通讯结果与 NicheNet 模型整合：参考根目录 `nichenet.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_celltalk.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_celltalk.qmd
git commit -m "add spatial_celltalk.qmd: CellChat v2 spatial runnable + COMMOT reference"
```

---

## Task 12: Build `spatial_scRNA_integration.qmd`

**Files:**
- Create: `spatial/spatial_scRNA_integration.qmd`

- [ ] **Step 1: Create file**

Create `spatial/spatial_scRNA_integration.qmd`:

````markdown
---
title: "ST × scRNA 整合 / 标签迁移"
format: html
---

## 简介

把 scRNA 的细胞类型注释迁移到 spot，是反卷积外的另一条路径：
- **反卷积**（`spatial_deconvolution.qmd`）：每个 spot 的混合比例
- **标签迁移**（本文件）：每个 spot 的"主要类型 + 置信度"，更接近单细胞分析的语义

主要工具：Seurat `FindTransferAnchors` + `TransferData` + `MapQuery`。

```{r}
#| include: false
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
source(here::here("singlecell_utiles.R"))
```

## 数据加载

```{r}
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE)

allen <- LoadData("allen_reference") |> UpdateSeuratObject()
allen <- JoinLayers(allen)  # Seurat v5 必需
allen <- SCTransform(allen, assay = "RNA", verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  RunUMAP(dims = 1:30, return.model = TRUE)
```

## 锚点 + 标签迁移

```{r}
anchors <- FindTransferAnchors(
  reference = allen,
  query = brain,
  normalization.method = "SCT",
  reduction = "pcaproject"
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = allen$subclass,
  weight.reduction = brain[["pca"]],
  dims = 1:30
)
brain <- AddMetaData(brain, metadata = predictions)
SpatialDimPlot(brain, group.by = "predicted.id", label = TRUE)
```

## 用 MapQuery 把 spot 投到参考 UMAP

```{r}
brain <- MapQuery(
  anchorset = anchors,
  reference = allen,
  query = brain,
  refdata = list(subclass = "subclass"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
DimPlot(brain, reduction = "ref.umap", group.by = "predicted.subclass", label = TRUE)
```

## 置信度过滤

每个 spot 一个 `predicted.id.score`（最大类别概率），低置信度 spot 单独标注：

```{r}
brain$confident_id <- ifelse(brain$predicted.id.score >= 0.5,
                             brain$predicted.id, "Uncertain")
SpatialDimPlot(brain, group.by = "confident_id", label = TRUE)
```

## 与反卷积结果交叉验证

```{r}
#| eval: false
# 假设 spatial_deconvolution.qmd 的结果已存
brain_decon <- qs2::qs_read("Results/brain_with_RCTD.qs")
# 比较：迁移得到的 predicted.id 与 RCTD argmax 的一致率
rctd_top <- apply(GetAssayData(brain_decon, assay = "RCTD"), 2, function(x) names(which.max(x)))
table(brain$predicted.id, rctd_top[colnames(brain)])
```

## 注意事项

- Seurat v5：参考对象必须 `JoinLayers()` 才能传给 `FindTransferAnchors`
- `weight.reduction` 用 query 的 `pca` 比 `pcaproject` 更稳
- 不要在多样本整合后忘记重做 SCTransform，否则 anchor 找不到

## 下游

- 与反卷积结果交叉： `spatial/spatial_deconvolution.qmd`
- 多切片标签传递：先 Harmony 整合 ST，然后用此流程批量迁移（参考 `spatial_workflow.qmd` 的整合段）
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/spatial_scRNA_integration.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/spatial_scRNA_integration.qmd
git commit -m "add spatial_scRNA_integration.qmd: label transfer with FindTransferAnchors + MapQuery"
```

---

## Task 13: Build `xenium_workflow.qmd`

**Files:**
- Create: `spatial/xenium_workflow.qmd`

- [ ] **Step 1: Create file**

Create `spatial/xenium_workflow.qmd`:

````markdown
---
title: "10x Xenium 工作流"
format: html
---

## 简介

Xenium 是 10x 的亚细胞分辨率（< 1µm）成像空转平台。和 Visium 相比：
- 每个 cell 一行（已有 segmentation），不是 spot
- 基因 panel 有限（~300–5000 基因）
- 数据量大（百万级细胞）

Seurat v5 原生支持 Xenium：`LoadXenium()`。

```{r}
#| include: false
library(Seurat)
library(future)
library(ggplot2)
library(dplyr)
library(qs2)
source(here::here("singlecell_utiles.R"))
plan("multisession", workers = 4)
options(future.globals.maxSize = 8 * 1024^3)
```

## 数据加载

```{r}
#| eval: false
# 10x 公开 Xenium 小鼠脑样本（Tiny Subset）
# 下载： https://www.10xgenomics.com/datasets （Xenium → Mouse Brain）
xen_dir <- "/path/to/xenium_output_bundle"
xen <- LoadXenium(xen_dir, fov = "fov")
```

替换为可演示的 mock（实际使用替换为上面的 `LoadXenium`）：

```{r}
xen <- readRDS(system.file("extdata", "xenium_demo.rds", package = "Seurat"))
# 若包内无 demo，跳过此 chunk
```

## QC

```{r}
xen[["nCount"]] <- xen$nCount_Xenium
xen[["nFeature"]] <- xen$nFeature_Xenium
SpatialQCPlot(xen, features = c("nCount", "nFeature"))
xen <- subset(xen, subset = nCount_Xenium > 10 & nFeature_Xenium > 5)
```

## 标准化与聚类

Xenium 基因数少，常用 `SCTransform` 或简单的 log-normalize：

```{r}
xen <- NormalizeData(xen) |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData() |>
  RunPCA(npcs = 30)
xen <- FindNeighbors(xen, dims = 1:20) |>
  FindClusters(resolution = 0.5)
xen <- RunUMAP(xen, dims = 1:20)
```

## 可视化

```{r}
ImageDimPlot(xen, fov = "fov", axes = TRUE, cols = "polychrome")
ImageFeaturePlot(xen, fov = "fov", features = head(VariableFeatures(xen), 4))
```

## 与 scRNA 参考做标签迁移

```{r}
#| eval: false
ref <- qs2::qs_read("/path/to/scRNA_reference.qs") |> JoinLayers()
ref <- SCTransform(ref) |> RunPCA() |> RunUMAP(dims = 1:30, return.model = TRUE)
anchors <- FindTransferAnchors(reference = ref, query = xen,
                               normalization.method = "LogNormalize",
                               reduction = "pcaproject")
preds <- TransferData(anchors, refdata = ref$celltype,
                      weight.reduction = xen[["pca"]], dims = 1:30)
xen <- AddMetaData(xen, preds)
ImageDimPlot(xen, group.by = "predicted.id", fov = "fov")
```

## 与 Visium 的关键差异

- **不需要反卷积**（已是单细胞分辨率）
- 邻域分析直接基于 cell 中心坐标，邻接定义不同（Delaunay / 半径 vs hex）
- 基因 panel 有限，HVG 选择基本无意义

## 注意事项

- `LoadXenium` 默认读 `cell_feature_matrix.h5` + `cells.csv.gz` + `morphology.ome.tif`
- 多 FOV 数据需用 `SubsetByImageBox`
- 大数据集用 BPCells 后端：`LoadXenium(..., load_bpcells = TRUE)`

## 下游

- niche / 邻域分析（基于 cell 邻接）： `spatial_niche.qmd` 中把 `k = 6` 改成更合适的 cell 邻接 k
- 空间通讯： `spatial_celltalk.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/xenium_workflow.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/xenium_workflow.qmd
git commit -m "add xenium_workflow.qmd: 10x Xenium imaging-based ST template"
```

---

## Task 14: Build `merfish_cosmx.qmd`

**Files:**
- Create: `spatial/merfish_cosmx.qmd`

- [ ] **Step 1: Create file**

Create `spatial/merfish_cosmx.qmd`:

````markdown
---
title: "MERFISH / CosMx 工作流"
format: html
---

## 简介

MERFISH（Vizgen）和 CosMx（Nanostring）是另两个亚细胞分辨率成像平台。和 Xenium 类似（每个 cell 一行），但 panel 大小、segmentation 策略、文件格式不同。

Seurat v5：
- MERFISH → `ReadVizgen()`
- CosMx → `ReadNanostring()` 或 `LoadNanostring()`

```{r}
#| include: false
library(Seurat)
library(future)
library(ggplot2)
library(dplyr)
source(here::here("singlecell_utiles.R"))
plan("multisession", workers = 4)
options(future.globals.maxSize = 8 * 1024^3)
```

## 数据加载：Vizgen MERFISH

```{r}
#| eval: false
mer_dir <- "/path/to/vizgen_output"
# Vizgen 输出：cell_by_gene.csv、cell_metadata.csv、cell_boundaries/
mer_data <- ReadVizgen(data.dir = mer_dir, type = "centroids")
mer <- CreateSeuratObject(counts = mer_data$transcripts, assay = "Vizgen")
mer[["fov"]] <- CreateFOV(mer_data$centroids, type = "centroids")
```

## 数据加载：Nanostring CosMx

```{r}
#| eval: false
cosmx_dir <- "/path/to/cosmx_output"
cosmx_data <- ReadNanostring(data.dir = cosmx_dir,
                             mtx.file = "exprMat_file.csv",
                             metadata.file = "metadata_file.csv",
                             molecules.file = "tx_file.csv",
                             segmentations.file = "polygons.csv",
                             type = c("centroids", "segmentations"))
cosmx <- CreateSeuratObject(counts = cosmx_data$matrix, assay = "Nanostring")
cosmx[["fov"]] <- CreateFOV(cosmx_data$centroids,
                            segmentations = cosmx_data$segmentations,
                            type = c("centroids", "segmentations"))
```

## QC + 标准化（共用）

```{r}
#| eval: false
mer[["nCount"]] <- mer$nCount_Vizgen
SpatialQCPlot(mer, features = c("nCount"))
mer <- subset(mer, subset = nCount_Vizgen > 20)
mer <- NormalizeData(mer) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 30) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(resolution = 0.5) |>
  RunUMAP(dims = 1:20)
```

## 可视化

```{r}
#| eval: false
ImageDimPlot(mer, fov = "fov", axes = TRUE)
ImageFeaturePlot(mer, fov = "fov", features = head(VariableFeatures(mer), 4))
```

## 平台差异速查

| 维度 | MERFISH (Vizgen) | CosMx (Nanostring) | Xenium (10x) |
|---|---|---|---|
| 典型 panel | 140–500 | 1000–6000 | 300–5000 |
| segmentation | DAPI + 膜 staining | DAPI + 膜 staining | DAPI（默认）/ multi-modal |
| 加载函数 | `ReadVizgen` | `ReadNanostring` | `LoadXenium` |
| 子细胞 transcript | `transcripts` | `molecules` | `transcripts` |

## 注意事项

- 三平台 segmentation 失败 spot 比例不同；先看 `nCount` 分布
- Vizgen 早期版本 cell ID 重复跨 FOV，加载后用 `paste0(fov, "_", cell)` 唯一化
- CosMx 数据量极大，强烈建议 `load_bpcells = TRUE` 或 disk-backed assay

## 下游

- 标签迁移： `spatial_scRNA_integration.qmd`（同样适用于 imaging 平台）
- 邻域分析： `spatial_niche.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/merfish_cosmx.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/merfish_cosmx.qmd
git commit -m "add merfish_cosmx.qmd: Vizgen MERFISH + Nanostring CosMx templates"
```

---

## Task 15: Build `visium_hd.qmd`

**Files:**
- Create: `spatial/visium_hd.qmd`

- [ ] **Step 1: Create file**

Create `spatial/visium_hd.qmd`:

````markdown
---
title: "Visium HD 工作流"
format: html
---

## 简介

Visium HD 是 10x 的高分辨率版（2µm × 2µm bins，分析常用 8µm 或 16µm bins）。和标准 Visium 的关键差异：
- 每张切片 ~10⁶ bin，spot 数从 5k 跳到 1M+
- 默认 bin 是 8µm
- 必须用 BPCells 等 disk-backed 后端

Seurat v5 原生支持： `Load10X_Spatial(..., bin.size = c(8, 16))`。

```{r}
#| include: false
library(Seurat)
library(BPCells)
library(future)
library(ggplot2)
library(dplyr)
source(here::here("singlecell_utiles.R"))
plan("multisession", workers = 4)
options(future.globals.maxSize = 16 * 1024^3)
```

## 数据加载

```{r}
#| eval: false
hd_dir <- "/path/to/visium_hd_output/binned_outputs"
hd <- Load10X_Spatial(
  data.dir = hd_dir,
  bin.size = c(8, 16),
  filter.matrix = TRUE
)
# bin 8 默认 active
DefaultAssay(hd) <- "Spatial.008um"
```

## QC

```{r}
#| eval: false
hd[["percent.mt"]] <- PercentageFeatureSet(hd, pattern = "^mt-",
                                           assay = "Spatial.008um")
SpatialQCPlot(hd, features = c("nCount_Spatial.008um",
                               "nFeature_Spatial.008um", "percent.mt"))
hd <- subset(hd, subset = nCount_Spatial.008um > 10 & percent.mt < 30)
```

## 标准化（log-norm 推荐，SCTransform 太慢）

```{r}
#| eval: false
hd <- NormalizeData(hd, assay = "Spatial.008um")
hd <- FindVariableFeatures(hd, assay = "Spatial.008um", nfeatures = 2000)
hd <- ScaleData(hd, assay = "Spatial.008um")
hd <- RunPCA(hd, assay = "Spatial.008um", reduction.name = "pca.008um", npcs = 30)
```

## 聚类（用 sketch / 子采样加速）

```{r}
#| eval: false
hd <- SketchData(hd, assay = "Spatial.008um", ncells = 50000,
                 method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(hd) <- "sketch"
hd <- FindVariableFeatures(hd) |> ScaleData() |> RunPCA(npcs = 30)
hd <- FindNeighbors(hd, dims = 1:20) |>
  FindClusters(resolution = 0.5)
hd <- RunUMAP(hd, dims = 1:20, return.model = TRUE)

# 投回完整数据
hd <- ProjectData(hd,
                  assay = "Spatial.008um",
                  full.reduction = "pca.008um",
                  sketched.assay = "sketch",
                  sketched.reduction = "pca",
                  umap.model = "umap",
                  dims = 1:20,
                  refdata = list(cluster_full = "seurat_clusters"))
```

## 可视化

```{r}
#| eval: false
DefaultAssay(hd) <- "Spatial.008um"
SpatialDimPlot(hd, group.by = "cluster_full", label = FALSE) +
  theme(legend.position = "right")
```

## 16µm bin（用于跨样本聚类比较）

```{r}
#| eval: false
DefaultAssay(hd) <- "Spatial.016um"
hd <- NormalizeData(hd, assay = "Spatial.016um")
# ... 重复上面的流程
```

## 注意事项

- Visium HD 必须 `library(BPCells)` 配合 `Load10X_Spatial`，否则内存爆掉
- 全量聚类不可行，用 `SketchData` + `ProjectData`
- `SCTransform` 在 1M+ bins 上太慢，用 `NormalizeData`
- 多个 bin size assay 同时挂在对象上，操作前显式 `DefaultAssay()`

## 下游

- 反卷积： Visium HD 的 8µm bin 已接近单细胞，可不必反卷积；如需反卷积用 RCTD `doublet_mode = "full"`
- 与 scRNA 整合： `spatial_scRNA_integration.qmd`
````

- [ ] **Step 2: Syntax check**

```bash
Rscript -e '
tmp <- tempfile(fileext = ".R")
knitr::purl("spatial/visium_hd.qmd", output = tmp, quiet = TRUE)
parse(tmp)
cat("OK\n")
'
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add spatial/visium_hd.qmd
git commit -m "add visium_hd.qmd: Visium HD 8/16um bin template with BPCells + sketch"
```

---

## Task 16: Update `README.Rmd` and regenerate `README.md`

**Files:**
- Modify: `README.Rmd`
- Regenerate: `README.md`

- [ ] **Step 1: Update `README.Rmd` spatial section**

Replace the current spatial bullet (around line 79–80 of `README.md`, find the matching place in `README.Rmd`):

```
- **`spatial/`** — `Spatial_Transcripts.qmd`, `10x_Genomics_Visium.R`,
  `spacexr_STalign.qmd`, `ST_test.qmd`, `spatial_scRNA_workflow.ipynb`.
```

with this expanded block:

```
- **`spatial/`** — 10 主题模板覆盖 Visium / Visium HD / Xenium / MERFISH /
  CosMx 全平台。按依赖关系分组：
  - 端到端流程： `spatial_workflow.qmd`
  - 空转专属分析： `spatial_SVG.qmd`, `spatial_deconvolution.qmd`,
    `spatial_domain.qmd`, `spatial_niche.qmd`, `spatial_celltalk.qmd`,
    `spatial_scRNA_integration.qmd`
  - 亚细胞分辨率平台： `xenium_workflow.qmd`, `merfish_cosmx.qmd`,
    `visium_hd.qmd`
  - Python 对照： `spatial_scRNA_workflow.ipynb`
```

- [ ] **Step 2: Regenerate `README.md`**

```bash
Rscript -e 'rmarkdown::render("README.Rmd")'
```

Expected: README.md updates without error.

- [ ] **Step 3: Verify**

```bash
grep -c "spatial_workflow.qmd" README.md
grep -c "xenium_workflow.qmd" README.md
```

Expected: both `>= 1`.

- [ ] **Step 4: Commit**

```bash
git add README.Rmd README.md
git commit -m "update README: expanded spatial/ section with 10-file template layout"
```

---

## Final integration check (manual, optional)

After all tasks complete, the user can do a full render of one or two `.qmd` files end-to-end as a sanity check (this is not part of the per-task verification because it requires `SeuratData` install + ~10 min compute):

```bash
quarto render spatial/spatial_workflow.qmd
quarto render spatial/spatial_SVG.qmd
```

If a render fails, fix the issue in that file and re-commit.

---

## Self-review notes

After writing this plan, I checked the following against the spec:

1. **Spec coverage:** 10 `.qmd` files (Tasks 5, 7–15) + 4 helpers (Tasks 1–4) + archive cleanup (Task 6) + README update (Task 16) = all spec deliverables covered.
2. **Placeholders:** No `TBD`, no "implement later", no `// fill in here`. Each `.qmd` task contains the full intended content of the file.
3. **Consistency:**
   - Helper signatures match across tasks (`LoadVisiumMulti(dirs, sample_names, filename, mt_pattern)` is called consistently).
   - All `source(here::here("singlecell_utiles.R"))` references match the actual repo path.
   - Reactome enrichment is present in `spatial_SVG.qmd` and `spatial_celltalk.qmd` (the two files where enrichment naturally appears).
   - `JoinLayers()` appears in every Seurat-v5 reference-loading step (`spatial_deconvolution.qmd`, `spatial_scRNA_integration.qmd`).
   - `qs2::` (not `qs::`) used everywhere serialization appears.
   - No `ggpubr` anywhere.
4. **Comparison direction:** `spatial_workflow.qmd` Region DE is "Cortex vs Striatum" with consistent variable name `de_cortex_vs_striatum`.

---

**End of plan.**
