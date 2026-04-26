<!-- README.md is generated from README.Rmd. Re-render with rmarkdown::render("README.Rmd"). -->

# single_cell_RNA

Personal collection of single-cell RNA-seq analysis scripts, workflows, and
notes — accumulated over years of practice. Organised by topic so the active
pipelines are easy to find and historical material stays out of the way.

## Directory layout

```
.
├── *.qmd / *.R / *.Rmd        # Active scRNA-seq workflows (see below)
├── projects/                  # Project- and dataset-specific analyses
├── spatial/                   # Spatial transcriptomics workflows
├── scATAC/                    # Single-cell ATAC-seq scripts
├── immune_repertoire/         # TCR/BCR (scRepertoire) scripts
├── archive/
│   ├── legacy_2021_2023/      # Early scripts, mostly superseded
│   └── superseded_rmd/        # Old .Rmd whose .qmd version is now in root
├── datasets/                  # Small demo data only (large intermediates gitignored)
└── README.Rmd / README.md
```

## Active workflows (root)

### End-to-end pipelines
- `2025scRNA_workflow.qmd` — current end-to-end scRNA-seq workflow
- `Seurat5.qmd` — Seurat v5 specifics (layers, BPCells, etc.)
- `scRNAseq_workflow.ipynb` — Python/Scanpy counterpart

### QC, normalisation, integration
- `scRNA_Analysis.qmd` — general analysis recipes
- `cell_anotation.qmd` — automated annotation (SingleR, etc.)
- `check-all-markers.R`, `tcell_markers_list.R` — marker gene utilities
- `inmmune_cells.qmd` — immune cell sub-clustering

### Differential expression & enrichment
- `SC_DESeq2.qmd` — pseudobulk DE with DESeq2
- `GSVA.qmd` — pathway activity scoring
- `AUCell.qmd` — gene-set scoring per cell
- `scMetabolism.qmd` — metabolic pathway scoring
- `progeny.qmd` — pathway activity (PROGENy)
  (the marker-module + enrichment heatmap helper now lives in
  `singlecell_utiles.R::plotMarkerEnrichmentHeatmap`)

### Trajectory / dynamics
- `monocle3_trajectory.qmd`, `monocle_lingang.qmd` — Monocle3 trajectories
- `cellrank.qmd` — CellRank fate mapping
- `CytoTRACE2.qmd` — differentiation potential

### Cell–cell communication
- `cellchat.qmd`, `cellphonedb.qmd`, `nichenet.qmd`

### Gene regulatory networks
- `SCENIC.qmd`, `pyscenic.qmd` — SCENIC / pySCENIC

### Co-expression & factor models
- `hdWGCNA.qmd`, `plot_hdWGCNA.qmd` — hdWGCNA modules
- `scNMF.qmd` — non-negative matrix factorisation
- `mofa2.qmd` — MOFA2 multi-omics factors

### Tumour-specific
- `inferCNV.qmd` — CNV inference from scRNA-seq

### Cell composition
- `cell_proportion_analysis.qmd` — proportion / compositional analysis

### Interop & utilities
- `Seurat_2_h5ad.qmd`, `anndataR.qmd`, `SC_ID_Convert.qmd` — format/ID conversions
- `singlecell_utiles.R` — shared helper functions
- `scRNA_Notes.Rmd` — personal notes (kept in root for quick reference)

## Topic-specific subdirectories

- **`projects/`** — analyses tied to specific datasets/projects:
  `BLM_six.qmd`, `BLM_mice.qmd`, `liver_cancer.qmd`, `GSE138794.qmd`,
  `GSE242889.qmd`, `临港实验室.qmd`.
- **`spatial/`** — `Spatial_Transcripts.qmd`, `10x_Genomics_Visium.R`,
  `spacexr_STalign.qmd`, `ST_test.qmd`, `spatial_scRNA_workflow.ipynb`.
- **`scATAC/`** — placeholders (`ArchR.R`, `Signnac.R`) to be expanded.
- **`immune_repertoire/`** — `scRepertoire.R`, to be expanded.

## Archive

`archive/legacy_2021_2023/` holds early Seurat v3/v4 era scripts, stub files,
and one-off helpers that have been replaced by the active workflows above.
`archive/superseded_rmd/` holds old `.Rmd` versions whose `.qmd` rewrites now
live in the root (Monocle3, inferCNV, GSVA, velocyto). Kept for reference; not
maintained.

## Usage

Open any `.qmd` / `.Rmd` in RStudio or VS Code with the Quarto extension and
render. Common dependencies: `Seurat (>=5)`, `tidyverse`,
`SingleCellExperiment`, `monocle3`, `harmony`, `CellChat`, `SCENIC`/`AUCell`,
`GSVA`, `infercnv`, `SingleR`, `hdWGCNA`. See the `library()` calls at the top
of each file for the exact list.
