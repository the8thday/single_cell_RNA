# Spatial Transcriptomics Workflow Templates — Design

**Date:** 2026-04-26
**Scope:** Reorganize and expand `spatial/` into a coherent set of template `.qmd` files covering the main spatial transcriptomics analyses, mirroring the style of the root-level scRNA workflows.

## Goal

Bring `spatial/` to the same template-quality bar as the root scRNA workflows (`cellchat.qmd`, `monocle3_trajectory.qmd`, etc.): topic-per-file, runnable on small public data, reusable as personal reference. Cover Visium (incl. Visium HD) and the imaging-based platforms (Xenium / MERFISH / CosMx), and add the spatial-specific analyses missing today (deconvolution, SVG, domain detection, niche, spatial cell-cell talk, ST × scRNA integration).

## Current state

`spatial/` today contains:

- `Spatial_Transcripts.qmd` (822 lines) — GBM Visium single-sample walkthrough
- `10x_Genomics_Visium.R` (562 lines) — multi-sample Visium pipeline (script form)
- `ST_test.qmd` (169 lines) — small Visium test
- `spacexr_STalign.qmd` (19 lines) — placeholder only
- `spatial_scRNA_workflow.ipynb` (1261 lines) — Python/Scanpy counterpart

Issues:
- The three R files overlap heavily (QC / SCT / cluster / markers) but share no helpers.
- Spatial-specific analyses are largely absent: spot deconvolution, spatially variable genes, spatial domains, niche analysis, spatial-aware cell–cell communication, ST × scRNA integration, imaging-platform handling.

## Decisions

- **Platform scope:** Visium + Visium HD + Xenium + MERFISH + CosMx. Stereo-seq / Slide-seqV2 / DBiT-seq are out of scope.
- **Existing files:** consolidate (option B). Merge the three overlapping R files into one end-to-end workflow; archive originals; expand from there.
- **Completeness:** option B for all files (runnable on a small public dataset). Two exceptions where Python-only tools fall back to option A (skeleton, annotated but not executed): SpaGCN in `spatial_domain.qmd`, and COMMOT in `spatial_celltalk.qmd`. The R-side sections of those two files (BayesSpace, CellChat v2 spatial) remain runnable.
- **Python:** the existing `.ipynb` is preserved as-is. No new Python files in this scope.

## File deliverables

### New / refactored `.qmd` (10 files, all in `spatial/`)

| File | Type | Completeness | Primary tools |
|---|---|---|---|
| `spatial_workflow.qmd` | refactored merge | B (runnable) | Seurat v5: QC / SCT / clustering / SpatialFeaturePlot |
| `spatial_deconvolution.qmd` | new | B | RCTD (spacexr) + CARD + SPOTlight |
| `spatial_SVG.qmd` | new | B | Moran's I / SPARK-X / `FindSpatiallyVariableFeatures` |
| `spatial_domain.qmd` | new | **A (skeleton)** | BayesSpace (R, runnable) + SpaGCN (annotated only, not executed) |
| `spatial_niche.qmd` | new | B | hoodscanR + MISTy + adjacency co-occurrence |
| `spatial_celltalk.qmd` | new | B (CellChat) + A (COMMOT) | CellChat v2 spatial runnable; COMMOT annotated only (Python) |
| `xenium_workflow.qmd` | new | B | Seurat v5 `LoadXenium` |
| `merfish_cosmx.qmd` | new | B | Seurat v5 `ReadVizgen` / `ReadNanostring` |
| `visium_hd.qmd` | new | B | Seurat v5 Visium HD (8µm bin) |
| `spatial_scRNA_integration.qmd` | new | B | `FindTransferAnchors` + label transfer |

### File adjustments

- `Spatial_Transcripts.qmd`, `ST_test.qmd`, `10x_Genomics_Visium.R` → contents merged into `spatial_workflow.qmd`; originals **moved to `archive/superseded_rmd/`** (preserve history; do not delete).
- `spacexr_STalign.qmd` (19-line placeholder) → relevant lines folded into `spatial_deconvolution.qmd`; original **deleted**.
- `spatial_scRNA_workflow.ipynb` → **untouched**.

## File internal structure

Every new `.qmd` follows the same template (consistent with root-level workflows for cross-reference):

1. Intro / when to use (1–2 short Chinese paragraphs)
2. `library()` + global options
3. Data loading (public sample + path variables)
4. Main analysis (topic-specific sections)
5. Visualization (`SpatialFeaturePlot` / `SpatialDimPlot` first-class)
6. Hand-off to downstream files (cross-link to other `spatial/` files or root scripts)
7. References / caveats

### Example data conventions

Choosing public datasets that load reliably and are small enough to render in reasonable time:

- Visium-class (`spatial_workflow`, `SVG`, `domain`, `niche`, `celltalk`, `scRNA_integration`) → `SeuratData::stxBrain` (mouse brain coronal, well-trodden in Seurat tutorials).
- `spatial_deconvolution.qmd` → `stxBrain` + `SeuratData::allen_reference` (mouse brain scRNA reference).
- `xenium_workflow.qmd` → 10x public Xenium mouse brain sample (path placeholder + download URL in comments).
- `merfish_cosmx.qmd` → Vizgen MERFISH MOp + Nanostring CosMx small sample (path placeholders).
- `visium_hd.qmd` → 10x public Visium HD mouse intestine sample (path placeholder).

For files using `SeuratData`, include the `InstallData()` hint in a commented setup block so the reader can install on first run.

## Repo-wide conventions (enforced per file where applicable)

These come from `CLAUDE.md` and saved feedback memories. Each is non-negotiable for files that touch the relevant code path:

- **Seurat v5:** call `JoinLayers()` before any tool that consumes raw counts directly (deconvolution inputs, CytoTRACE2-style tools, scVelo bridges).
- **Serialization:** `qs2::qs_save()` / `qs2::qs_read()`. Never `qs::`.
- **Statistical annotations:** base R (`wilcox.test`, `t.test`) + `ggplot2::annotate()`. **No `ggpubr`.**
- **Enrichment:** anywhere GO / KEGG enrichment appears (e.g. ligand–receptor pair enrichment in `spatial_celltalk.qmd`, SVG enrichment in `spatial_SVG.qmd`), Reactome (`ReactomePA::enrichPathway` or `clusterProfiler::gsePathway`) must also be run.
- **Comparison direction:** when a file computes DE / log2FC between conditions or niches, variable names, plot titles, and comments must agree on the same numerator/denominator throughout.
- **AUC-as-Assay (if it shows up):** clean feature names (no parens/spaces), no `NormalizeData` on the AUC assay, `FindMarkers(..., logfc.threshold = 0, min.pct = 0)`.
- **Multi-condition comparisons:** add exploratory analyses beyond the standard pipeline output (compositional shifts, condition-specific markers per cluster, pathway differences per cluster).

## Shared helpers (additions to `singlecell_utiles.R`)

Per the repo rule "new reusable helpers go into `singlecell_utiles.R`, not standalone `.R` scripts." All helpers use namespaced calls (`pkg::fn`) internally so they don't depend on `library()` order.

| Helper | Purpose | Used by |
|---|---|---|
| `LoadVisiumMulti(dirs, sample_names, ...)` | Batch-load multi-sample Visium + set `orig.ident` + add basic QC columns; returns `list<Seurat>` | `spatial_workflow`, `spatial_deconvolution`, `spatial_scRNA_integration` |
| `SpatialQCPlot(obj, features, ...)` | Combined `VlnPlot` + `SpatialFeaturePlot` panel (standard QC layout) | `spatial_workflow`, `xenium_workflow`, `visium_hd` |
| `PlotDeconvProportion(obj, decon_assay, ...)` | Stacked bar + spatial overlay of deconvolution proportions | `spatial_deconvolution`, `spatial_niche` |
| `NeighborComposition(obj, group_col, k = 6)` | k-neighborhood cell-type composition per spot; returns long-format `data.frame` | `spatial_niche`, `spatial_celltalk` |

### Explicitly NOT writing

- Thin wrappers around `SCTransform` / `RunPCA` — call Seurat directly.
- Platform-specific loaders for Xenium / MERFISH / CosMx — Seurat v5's `LoadXenium` / `ReadVizgen` / `ReadNanostring` are already sufficient.
- A unified SVG wrapper — different packages have substantially different APIs and assumptions; a wrapper would lose information.

## Dependencies and build order

```
spatial_workflow ──┬─→ spatial_SVG
                   ├─→ spatial_deconvolution ──┬─→ spatial_niche ──→ spatial_celltalk
                   │                           └─→ spatial_scRNA_integration
                   └─→ spatial_domain (independent)

xenium_workflow ──┐
merfish_cosmx ────┼─→ (imaging-based; cross-reference each other only)
visium_hd ────────┘
```

Implementation order follows this graph: `spatial_workflow` first (it depends on the helpers), then helpers, then the rest.

Concretely, the order of work is:

1. Add the four helpers to `singlecell_utiles.R`.
2. Build `spatial_workflow.qmd` (consumes 3 archived files' content).
3. Move `Spatial_Transcripts.qmd`, `ST_test.qmd`, `10x_Genomics_Visium.R` to `archive/superseded_rmd/`. Delete `spacexr_STalign.qmd`.
4. `spatial_SVG.qmd`, `spatial_deconvolution.qmd`, `spatial_domain.qmd` (in any order — independent).
5. `spatial_niche.qmd` (depends on deconvolution output).
6. `spatial_celltalk.qmd`, `spatial_scRNA_integration.qmd` (depend on niche / deconvolution respectively).
7. `xenium_workflow.qmd`, `merfish_cosmx.qmd`, `visium_hd.qmd` (independent imaging cluster).
8. README update.

## README update

In `README.md` (and `README.Rmd`, since `README.md` is regenerated from it), replace the current `spatial/` bullet under "Topic-specific subdirectories" with the grouped 10-file layout matching the dependency graph. Add a one-line note that `spatial/` has been expanded into a 10-topic template set.

## Out of scope (intentional non-goals)

- Multi-slice alignment / 3D reconstruction (`PASTE`, `PASTE2`, `STalign`, `STUtility` 3D) — toolchain still moving fast; would become stale templates.
- Stereo-seq, Slide-seqV2, DBiT-seq, GeoMx workflows.
- Any new Python notebook. The existing `.ipynb` is preserved as-is.
- Any change to `archive/legacy_2021_2023/` or other `archive/` content beyond moving the three superseded files in.
- Refactoring the helpers already in `singlecell_utiles.R`.

## Risks and mitigations

- **`SeuratData` install friction:** `stxBrain` and `allen_reference` require `InstallData()`. Mitigate by including the install command commented in setup, plus a fallback note pointing to the 10x download URL.
- **Imaging-platform datasets are large:** Xenium / MERFISH / CosMx full samples are multi-GB. Use the smallest official demo subsets and put paths behind variables so the user substitutes their own.
- **Python-only tools in two files:** SpaGCN (`spatial_domain.qmd`) and COMMOT (`spatial_celltalk.qmd`) are documented but not executed. The R-side analyses (BayesSpace; CellChat v2 spatial) run end-to-end so both files still produce output on render.
- **`.gitignore` pressure:** rendered outputs (`Results/`, `figures/`, `*.qs`, etc.) must not be committed. The existing `.gitignore` already covers this; verify before committing rendered HTMLs.
