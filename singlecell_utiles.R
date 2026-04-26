

SafeAddMetaData <- function(seurat_obj, metadata, col.name = NULL) {
  # 获取 Seurat 对象的细胞名
  cell_names <- colnames(seurat_obj)
  
  # 判断 metadata 类型
  if (is.vector(metadata)) {
    # 向量：必须有名字
    if (is.null(names(metadata))) {
      stop("Vector metadata must be a named vector with cell names as names.")
    }
    common_cells <- intersect(cell_names, names(metadata))
    if (length(common_cells) == 0) {
      stop("No matching cell names found between Seurat object and metadata vector.")
    }
    metadata <- metadata[cell_names]  # 顺序对齐
    if (is.null(col.name)) {
      stop("Please specify `col.name` when adding a vector.")
    }
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata, col.name = col.name)
    
  } else if (is.data.frame(metadata)) {
    # 数据框：行名必须是细胞名
    if (is.null(rownames(metadata))) {
      stop("Data frame metadata must have cell names as rownames.")
    }
    common_cells <- intersect(cell_names, rownames(metadata))
    if (length(common_cells) == 0) {
      stop("No matching cell names found between Seurat object and metadata data frame.")
    }
    metadata <- metadata[cell_names, , drop = FALSE]  # 顺序对齐
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
    
  } else {
    stop("Unsupported metadata type. Only vector and data.frame are supported.")
  }
  
  return(seurat_obj)
}


SafeAddMetaData <- function(seurat_obj, metadata, col.name = NULL, verbose = TRUE, 
                            drop.extra = TRUE) {
  cell_names <- colnames(seurat_obj)
  
  # 记录原来的 metadata 列名
  original_cols <- colnames(seurat_obj@meta.data)
  
  if (is.vector(metadata)) {
    if (is.null(names(metadata))) {
      stop("Vector metadata must have names (cell names).")
    }
    
    if (!all(cell_names %in% names(metadata))) {
      missing <- setdiff(cell_names, names(metadata))
      stop(paste0("Metadata vector missing ", length(missing), " cell(s): ", 
                  paste(head(missing), collapse = ", "), "..."))
    }
    
    metadata <- metadata[cell_names]
    if (is.null(col.name)) {
      stop("Please provide `col.name` when adding a vector.")
    }
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata, col.name = col.name)
    
    if (verbose) {
      message("✅ Added column '", col.name, "' to metadata.")
    }
    
  } else if (is.data.frame(metadata)) {
    if (is.null(rownames(metadata))) {
      stop("Data frame metadata must have rownames (cell names).")
    }
    
    # Warn about unmatched cell names
    missing_cells <- setdiff(cell_names, rownames(metadata))
    extra_cells <- setdiff(rownames(metadata), cell_names)
    
    if (length(missing_cells) > 0) {
      stop(paste0("Metadata missing ", length(missing_cells), " cell(s): ", 
                  paste(head(missing_cells), collapse = ", "), "..."))
    }
    
    if (length(extra_cells) > 0 && !drop.extra) {
      stop(paste0("Metadata contains ", length(extra_cells), " extra cell(s): ", 
                  paste(head(extra_cells), collapse = ", "), "..."))
    }
    
    # Align rows by Seurat cell order
    metadata <- metadata[cell_names, , drop = FALSE]
    
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
    
    if (verbose) {
      new_cols <- setdiff(colnames(seurat_obj@meta.data), original_cols)
      message("✅ Added metadata columns: ", paste(new_cols, collapse = ", "))
    }
    
  } else {
    stop("Unsupported metadata type. Must be a named vector or data.frame with rownames.")
  }
  
  return(seurat_obj)
}




#' Grouped Dot Plot for Seurat (x=genes, y=cell types, grouped by condition)
#'
#' @param obj Seurat object.
#' @param genes Character vector of genes for x-axis.
#' @param cell_var meta.data column for cell types (y-axis). Default "celltype".
#' @param condition_var meta.data column for grouping within each gene. Default "condition".
#' @param assay Assay to use; default DefaultAssay(obj).
#' @param layer Seurat v5 layer ("data","counts","scale.data"); auto-mapped to slot on v4.
#' @param map_size Which metric controls point size: "avg" (average expression) or "pct" (percent expressed). Default "avg".
#' @param detection_threshold Expression > threshold counts as expressed. Default 0.
#' @param standardize Whether to z-score per gene then rescale to `scale_range` for fill. Default TRUE.
#' @param scale_range Range to rescale standardized values to. Default c(-2.5, 2.5).
#' @param size_range Range of point sizes. Default c(1.5, 7).
#' @param palette Colors for fill (expression). Default c("#2166AC","#F7F7F7","#B2182B").
#' @param group_mode "dodge" (side-by-side within each gene) or "facet" (facet by condition). Default "dodge".
#' @param dodge_width Horizontal dodge width if group_mode="dodge". Default 0.6.
#' @param base_size ggplot base font size. Default 11.
#' @param legend_position "right","bottom","top","left","none". Default "right".
#' @param order_genes Optional explicit gene order on x-axis.
#' @param order_celltypes Optional explicit cell type order on y-axis.
#' @param order_conditions Optional explicit condition order for grouping.
#' @param show_zero_placeholders If TRUE, draw tiny dots when a combo has zero cells. Default FALSE.
#' @param zero_placeholder_size Size used for zero placeholders. Default 0.5.
#'
#' @return ggplot object
#' @export
scDotGrouped <- function(
  obj,
  genes,
  cell_var = "celltype",
  condition_var = "condition",
  assay = NULL,
  layer = "data",
  map_size = c("avg","pct"),
  detection_threshold = 0,
  standardize = TRUE,
  scale_range = c(-2.5, 2.5),
  size_range = c(1.5, 7),
  palette = c("#2166AC","#F7F7F7","#B2182B"),
  group_mode = c("dodge","facet"),
  dodge_width = 0.6,
  base_size = 11,
  legend_position = "right",
  order_genes = NULL,
  order_celltypes = NULL,
  order_conditions = NULL,
  show_zero_placeholders = FALSE,
  zero_placeholder_size = 0.5
) {
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("scales", quietly = TRUE)
  requireNamespace("tidyselect", quietly = TRUE)

  map_size <- match.arg(map_size)
  group_mode <- match.arg(group_mode)

  meta <- obj@meta.data
  if (!all(c(cell_var, condition_var) %in% colnames(meta))) {
    stop("meta.data 中缺少 ", paste(setdiff(c(cell_var, condition_var), colnames(meta)), collapse = ", "))
  }
  if (length(genes) == 0) stop("`genes` 不能为空。")

  # 去除空格，避免幽灵水平
  meta[[cell_var]]     <- trimws(meta[[cell_var]])
  meta[[condition_var]]<- trimws(meta[[condition_var]])

  # 选择 assay（不修改原对象）
  o <- obj
  if (!is.null(assay)) {
    if (!assay %in% names(o@assays)) stop(sprintf("Assay '%s' 不存在。", assay))
    Seurat::DefaultAssay(o) <- assay
  }

  # 仅保留存在的基因
  feats <- genes[genes %in% rownames(Seurat::GetAssay(o))]
  if (length(feats) == 0) stop("所有指定的基因都不存在于对象中。")
  if (length(feats) < length(genes)) {
    warning(sprintf("以下基因未找到并被忽略：%s", paste(setdiff(genes, feats), collapse = ", ")))
  }

  # 设定显示顺序
  cell_levels <- if (!is.null(order_celltypes)) order_celltypes else unique(meta[[cell_var]])
  cond_levels <- if (!is.null(order_conditions)) order_conditions else unique(meta[[condition_var]])
  gene_levels <- if (!is.null(order_genes)) order_genes else feats

  # Seurat v5 vs v4: FetchData 是否支持 layer
  fetch_formals <- names(formals(Seurat::FetchData))
  use_layer <- "layer" %in% fetch_formals
  vars_to_get <- unique(c(feats, cell_var, condition_var))

  df <- tryCatch(
    {
      if (use_layer) {
        Seurat::FetchData(o, vars = vars_to_get, layer = layer)
      } else {
        slot_fallback <- switch(layer,
                                "data" = "data",
                                "counts" = "counts",
                                "scale.data" = "scale.data",
                                layer)
        Seurat::FetchData(o, vars = vars_to_get, slot = slot_fallback)
      }
    },
    error = function(e) {
      stop("FetchData 失败，请检查 `layer/slot` 与 `assay`：", e$message)
    }
  )

  # 长表：每行一个 (cell × gene)
  long_df <- df |>
    tibble::as_tibble(rownames = "cell") |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(feats),
      names_to = "gene", values_to = "expr"
    ) |>
    dplyr::mutate(
      cell_type = .data[[cell_var]],
      condition = .data[[condition_var]]
    )

  # 聚合：每个 (gene × cell_type × condition)
  dot_data <- long_df |>
    dplyr::group_by(gene, cell_type, condition) |>
    dplyr::summarise(
      avg_exp = mean(expr, na.rm = TRUE),
      pct_exp = mean(expr > detection_threshold, na.rm = TRUE) * 100,
      .groups = "drop"
    )

  # 完整笛卡尔积，避免“只剩一个组”
  full_grid <- tidyr::expand_grid(
    gene       = gene_levels,
    cell_type  = cell_levels,
    condition  = cond_levels
  )
  dot_data <- full_grid |>
    dplyr::left_join(dot_data, by = c("gene","cell_type","condition")) |>
    dplyr::mutate(
      avg_exp = dplyr::coalesce(avg_exp, 0),
      pct_exp = dplyr::coalesce(pct_exp, 0)
    )

  # 颜色：按基因标准化（更利于跨细胞类型比较）
  if (isTRUE(standardize)) {
    dot_data <- dot_data |>
      dplyr::group_by(gene) |>
      dplyr::mutate(
        z = {
          m <- mean(avg_exp, na.rm = TRUE)
          s <- stats::sd(avg_exp, na.rm = TRUE)
          if (is.na(s) || s == 0) 0 else (avg_exp - m) / s
        },
        fill_val = {
          rng <- range(z, na.rm = TRUE)
          if (!all(is.finite(rng)) || diff(rng) == 0) mean(scale_range)
          else scales::rescale(z, to = scale_range, from = rng)
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-z)
  } else {
    dot_data$fill_val <- dot_data$avg_exp
  }

  # 尺寸：由 map_size 决定
  if (map_size == "avg") {
    dot_data$size_val <- dot_data$avg_exp
    size_name <- "Avg expr"
    size_limits <- NULL
  } else {
    dot_data$size_val <- dot_data$pct_exp
    size_name <- "Pct cells (%)"
    size_limits <- c(0, 100)
  }

  # 可选的零占位
  if (isTRUE(show_zero_placeholders)) {
    dot_data$size_val <- ifelse(dot_data$size_val == 0, zero_placeholder_size, dot_data$size_val)
  }

  # 因子顺序
  dot_data <- dot_data |>
    dplyr::mutate(
      gene = factor(gene, levels = gene_levels),
      cell_type = factor(cell_type, levels = cell_levels),
      condition = factor(condition, levels = cond_levels)
    )

  # 绘图
  pd <- ggplot2::position_dodge(width = dodge_width)

  base <- ggplot2::ggplot(
    dot_data,
    ggplot2::aes(x = gene, y = cell_type)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = size_val,
        fill = fill_val,
        group = condition,
        color = condition
      ),
      shape = 21, stroke = 0.25,
      position = if (group_mode == "dodge") pd else "identity",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_gradientn(colors = palette, name = if (standardize) "Scaled expr" else "Expr") +
    ggplot2::scale_size_continuous(range = size_range, limits = size_limits, name = size_name) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey80"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  p <- if (group_mode == "facet") {
    base + ggplot2::facet_grid(cols = ggplot2::vars(condition), scales = "fixed", space = "free_x")
  } else {
    base
  }

  return(p)
}




#' Summarize gene expression by cell type and condition from a Seurat object
#'
#' @param obj Seurat object.
#' @param genes Character vector of genes to summarize.
#' @param cell_var meta.data column for cell types (y-axis). Default "celltype".
#' @param condition_var meta.data column for groups/conditions (e.g., "disease"). Default "condition".
#' @param assay Assay to use; default DefaultAssay(obj).
#' @param layer Seurat v5 layer ("data","counts","scale.data"); auto-mapped to slot on v4.
#' @param detection_threshold Expression > threshold counts as expressed for pct. Default 0.
#' @param complete_grid Logical; if TRUE, include all gene × cell × condition combos
#'        (fill with n_cells=0, avg/median/pct=NA when no cells). Default TRUE.
#'
#' @return A tibble with columns: gene, cell_type, condition, n_cells, avg_expr, pct_expr, median_expr
#' @export
scSummarizeExpr <- function(
  obj,
  genes,
  cell_var = "celltype",
  condition_var = "condition",
  assay = NULL,
  layer = "data",
  detection_threshold = 0,
  complete_grid = TRUE
){
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("tidyselect", quietly = TRUE)

  # --- checks
  meta <- obj@meta.data
  if (!all(c(cell_var, condition_var) %in% colnames(meta))) {
    stop("meta.data 缺少列：", paste(setdiff(c(cell_var, condition_var), colnames(meta)), collapse = ", "))
  }
  if (length(genes) == 0) stop("`genes` 不能为空。")

  # 去除空格，避免幽灵水平
  meta[[cell_var]]      <- trimws(meta[[cell_var]])
  meta[[condition_var]] <- trimws(meta[[condition_var]])

  # --- choose assay without mutating input
  o <- obj
  if (!is.null(assay)) {
    if (!assay %in% names(o@assays)) stop(sprintf("Assay '%s' 不存在。", assay))
    Seurat::DefaultAssay(o) <- assay
  }

  # --- keep available genes
  feats <- genes[genes %in% rownames(Seurat::GetAssay(o))]
  if (!length(feats)) stop("指定的基因均不在对象中。")
  if (length(feats) < length(genes)) {
    warning(sprintf("以下基因未找到并被忽略：%s", paste(setdiff(genes, feats), collapse = ", ")))
  }

  # --- Seurat v5 vs v4: FetchData with layer/slot
  fetch_formals <- names(formals(Seurat::FetchData))
  use_layer <- "layer" %in% fetch_formals
  vars_to_get <- unique(c(feats, cell_var, condition_var))

  df <- tryCatch(
    {
      if (use_layer) {
        Seurat::FetchData(o, vars = vars_to_get, layer = layer)
      } else {
        slot_fallback <- switch(layer,
                                "data"       = "data",
                                "counts"     = "counts",
                                "scale.data" = "scale.data",
                                layer)
        Seurat::FetchData(o, vars = vars_to_get, slot = slot_fallback)
      }
    },
    error = function(e) stop("FetchData 失败：", e$message)
  )

  # --- long format: per cell × gene
  long_df <- df |>
    tibble::as_tibble(rownames = "cell") |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(feats),
      names_to = "gene", values_to = "expr"
    ) |>
    dplyr::mutate(
      cell_type = .data[[cell_var]],
      condition = .data[[condition_var]]
    )

  # --- aggregate
  out <- long_df |>
    dplyr::group_by(gene, cell_type, condition) |>
    dplyr::summarise(
      n_cells     = dplyr::n(),
      avg_expr    = mean(expr, na.rm = TRUE),
      pct_expr    = mean(expr > detection_threshold, na.rm = TRUE) * 100,
      median_expr = stats::median(expr, na.rm = TRUE),
      .groups = "drop"
    )

  # --- complete full grid if requested
  if (isTRUE(complete_grid)) {
    full_grid <- tidyr::expand_grid(
      gene      = feats,
      cell_type = unique(meta[[cell_var]]),
      condition = unique(meta[[condition_var]])
    )
    out <- full_grid |>
      dplyr::left_join(out, by = c("gene","cell_type","condition")) |>
      dplyr::mutate(
        n_cells     = dplyr::coalesce(n_cells, 0L),
        avg_expr    = ifelse(n_cells == 0, NA_real_, avg_expr),
        pct_expr    = ifelse(n_cells == 0, NA_real_, pct_expr),
        median_expr = ifelse(n_cells == 0, NA_real_, median_expr)
      )
  }

  dplyr::arrange(out, gene, cell_type, condition)
}




#' Create a Custom, Grouped Dot Plot for Seurat Objects
#'
#' This function generates a highly customizable dot plot, similar to `Seurat::DotPlot`,
#' but with enhanced and robust handling for splitting groups using `split.by`.
#' It solves common alignment issues with dodged points by ensuring a complete
#' data grid before plotting. The function uses efficient `dplyr` and `tidyr`
#' verbs for data processing.
#'
#' @param seurat_obj A Seurat object.
#' @param features A character vector of features (e.g., genes) to plot on the x-axis.
#' @param group.by A character string specifying the metadata column to group cells by
#'   on the y-axis (e.g., "seurat_clusters", "cell_type").
#' @param split.by An optional character string specifying a metadata column to split
#'   points by within each group. This will create dodged points with different shapes
#'   for each level of the `split.by` variable.
#' @param assay The assay to pull data from (e.g., "RNA", "SCT").
#' @param slot The data layer to use from the specified assay (e.g., "data", "counts", "scale.data").
#'   Corresponds to the `layer` argument in Seurat v5.
#' @param dot.scale A numeric value that controls the maximum size of the dots.
#' @param cols A character vector of 3 colors for the gradient scale (low, mid, high).
#' @param col.min The minimum value for the color scale after scaling. Values below this
#'   will be set to this minimum.
#' @param col.max The maximum value for the color scale after scaling. Values above this
#'   will be set to this maximum.
#' @param dot.min The minimum percentage of cells expressing a feature for a dot to be
#'   plotted. Dots representing a lower percentage will have a size of 0.
#' @param scale A logical value indicating whether to scale the average expression
#'   values per feature (z-score transformation). Defaults to `TRUE`.
#' @param cluster_rows A logical value indicating whether to hierarchically cluster the
#'   rows (y-axis groups).
#' @param cluster_cols A logical value indicating whether to hierarchically cluster the
#'   columns (features on x-axis).
#' @param angle The angle for the x-axis text labels (e.g., 45, 90).
#' @param shapes An optional numeric vector to specify custom shapes for the points when
#'   `split.by` is used. If `NULL` (default), a default set of shapes is used.
#'   If a single value is provided (e.g., 16), all groups will use that shape.
#'
#' @return A `ggplot` object representing the dot plot.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (require("Seurat") && require("SeuratObject")) {
#'   pbmc_small$groups <- sample(c("GroupA", "GroupB"), size = ncol(pbmc_small), replace = TRUE)
#'   features <- c("CD3D", "MS4A1", "CST3", "NKG7", "PPBP")
#'
#'   # Using the new 'shapes' parameter to make all points solid circles
#'   CustomGroupedDotPlot(
#'     seurat_obj = pbmc_small,
#'     features = features,
#'     group.by = "seurat_clusters",
#'     split.by = "groups",
#'     shapes = 16  # All groups will now be solid circles
#'   )
#'
#'   # Using default shapes (by not specifying the 'shapes' parameter)
#'   CustomGroupedDotPlot(
#'     seurat_obj = pbmc_small,
#'     features = features,
#'     group.by = "seurat_clusters",
#'     split.by = "groups"
#'   )
#' }
#' }
CustomGroupedDotPlot <- function(
    seurat_obj, 
    features, 
    group.by = "seurat_clusters", 
    split.by = NULL,
    assay = "RNA",
    slot = "data",
    dot.scale = 6,
    cols = c("blue", "white", "red"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    scale = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    angle = 45,
    shapes = NULL
) {
  
  # 1. Load necessary packages
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
  
  # 2. Input Validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input 'seurat_obj' must be a Seurat object.")
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }
  
  meta_cols <- c(group.by, split.by)
  missing_cols <- meta_cols[!meta_cols %in% colnames(seurat_obj@meta.data)]
  if (length(missing_cols) > 0) {
    stop(paste("Metadata columns not found:", paste(missing_cols, collapse = ", ")))
  }
  
  # 3. Data Extraction and Preparation
  exp_matrix <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  
  features <- intersect(features, rownames(exp_matrix))
  if (length(features) == 0) {
    stop("None of the specified features were found in the selected assay.")
  }
  
  grouping_vars <- c(group.by, split.by)
  meta_data_subset <- seurat_obj@meta.data %>%
    select(all_of(grouping_vars)) %>%
    rownames_to_column(var = "cell_id")
  
  plot_data <- as.matrix(exp_matrix[features, , drop = FALSE]) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell_id") %>%
    pivot_longer(
      cols = all_of(features),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    left_join(meta_data_subset, by = "cell_id")

  # 4. Calculate Summaries
  summary_df <- plot_data %>%
    group_by(across(all_of(c(grouping_vars, "Gene")))) %>%
    summarise(
      AvgExpression = mean(Expression, na.rm = TRUE),
      PctExpressed = sum(Expression > 0) / n() * 100,
      .groups = "drop"
    )

  # 5. Ensure complete grid for proper dodging
  if (!is.null(split.by)) {
    grid_list <- list(
      unique(meta_data_subset[[group.by]]),
      unique(meta_data_subset[[split.by]]),
      features
    )
    names(grid_list) <- c(group.by, split.by, "Gene")
    
    all_combos <- expand.grid(grid_list, stringsAsFactors = FALSE)
    
    summary_df <- left_join(all_combos, summary_df, by = c(group.by, split.by, "Gene")) %>%
      mutate(
        AvgExpression = replace_na(AvgExpression, 0),
        PctExpressed = replace_na(PctExpressed, 0)
      )
  }
  
  # 6. Data Scaling and Capping
  if (scale) {
    summary_df <- summary_df %>%
      group_by(Gene) %>%
      mutate(AvgExpressionScaled = scale(AvgExpression)[,1]) %>%
      ungroup()
  } else {
    summary_df$AvgExpressionScaled <- summary_df$AvgExpression
  }

  summary_df <- summary_df %>%
    mutate(
      AvgExpressionScaled = if_else(is.nan(AvgExpressionScaled), 0, AvgExpressionScaled),
      AvgExpressionScaled = pmax(pmin(AvgExpressionScaled, col.max), col.min),
      PctExpressed = pmax(PctExpressed, dot.min)
    )

  # 7. Clustering and Factor Ordering
  row_levels <- unique(as.character(summary_df[[group.by]]))
  col_levels <- features

  if (cluster_rows || cluster_cols) {
    value_var <- if (scale) "AvgExpressionScaled" else "AvgExpression"
    
    cluster_matrix_data <- summary_df %>%
        group_by(.data[[group.by]], Gene) %>%
        summarise(value = mean(.data[[value_var]], na.rm = TRUE), .groups = "drop")
    
    cluster_matrix <- cluster_matrix_data %>%
      pivot_wider(
        names_from = "Gene",
        id_cols = all_of(group.by),
        values_from = "value",
        values_fill = 0
      ) %>%
      column_to_rownames(var = group.by)

    if (cluster_rows && nrow(cluster_matrix) > 1) {
      row_order <- hclust(dist(cluster_matrix))$order
      row_levels <- rownames(cluster_matrix)[row_order]
    }
    
    if (cluster_cols && ncol(cluster_matrix) > 1) {
      col_order <- hclust(dist(t(cluster_matrix)))$order
      col_levels <- colnames(cluster_matrix)[col_order]
    }
  } 
  
  summary_df[[group.by]] <- factor(summary_df[[group.by]], levels = rev(row_levels))
  summary_df$Gene <- factor(summary_df$Gene, levels = col_levels)

  # 8. Create the ggplot object
  p <- ggplot(summary_df, aes(x = Gene, y = .data[[group.by]]))
  
  if (!is.null(split.by)) {
    geom_args <- aes(
      size = PctExpressed, 
      color = AvgExpressionScaled,
      shape = .data[[split.by]],
      group = .data[[split.by]]
    )
    p <- p + geom_point(geom_args, position = position_dodge(width = 0.8))
    
    # --- MODIFICATION FOR CUSTOM SHAPES ---
    shape_values <- if (is.null(shapes)) {
      # Use default shapes if user did not provide any
      c(16, 17, 15, 18, 8, 4, 3, 2, 0, 1) 
    } else {
      # Use user-provided shapes
      shapes
    }
    
    p <- p + scale_shape_manual(
      name = split.by,
      values = shape_values
    )
    # --- END MODIFICATION ---

  } else {
    geom_args <- aes(size = PctExpressed, color = AvgExpressionScaled)
    p <- p + geom_point(geom_args)
  }
  
  # 9. Apply Scales, Theme, and Labels
  color_lab <- if(scale) "Scaled Avg.\nExpression" else "Avg.\nExpression"
  
  p <- p +
    scale_color_gradientn(
      colors = cols, 
      name = color_lab,
      limits = c(col.min, col.max),
      oob = scales::squish
    ) +
    scale_size_continuous(
      range = c(0, dot.scale), 
      name = "Percent\nExpressed",
      limits = c(dot.min, 100)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  
  return(p)
}







#' Create a Faceted Dot Plot for Seurat Objects
#'
#' This function generates a series of dot plots, faceted by a specified metadata
#' column (`split.by`). It provides a clear side-by-side comparison of expression
#' patterns across different conditions or groups. This function is an alternative
#' to `CustomGroupedDotPlot` when dodged points become too crowded.
#' It uses the same robust data processing and clustering logic to ensure consistency.
#'
#' @param seurat_obj A Seurat object.
#' @param features A character vector of features (e.g., genes) to plot on the x-axis.
#' @param group.by A character string specifying the metadata column to group cells by
#'   on the y-axis (e.g., "seurat_clusters", "cell_type").
#' @param split.by A character string specifying a metadata column to create facets by.
#'   This parameter is required for this function.
#' @param assay The assay to pull data from (e.g., "RNA", "SCT").
#' @param slot The data layer to use from the specified assay (e.g., "data", "counts").
#' @param dot.scale A numeric value that controls the maximum size of the dots.
#' @param cols A character vector of 3 colors for the gradient scale (low, mid, high).
#' @param col.min The minimum value for the color scale after scaling.
#' @param col.max The maximum value for the color scale after scaling.
#' @param dot.min The minimum percentage of cells expressing a feature for a dot to be
#'   plotted.
#' @param scale A logical value indicating whether to scale the average expression
#'   values per feature. Defaults to `TRUE`.
#' @param cluster_rows A logical value indicating whether to hierarchically cluster the
#'   rows. The clustering is performed on data aggregated across all splits to ensure
#'   a consistent y-axis order across all facets.
#' @param cluster_cols A logical value indicating whether to hierarchically cluster the
#'   columns (features). The order is consistent across all facets.
#' @param facet_ncol An integer specifying the number of columns for the facet grid.
#' @param angle The angle for the x-axis text labels.
#'
#' @return A `ggplot` object with faceted dot plots.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (require("Seurat") && require("SeuratObject")) {
#'   pbmc_small$groups <- sample(c("GroupA", "GroupB"), size = ncol(pbmc_small), replace = TRUE)
#'   features <- c("CD3D", "MS4A1", "CST3", "NKG7", "PPBP")
#'
#'   # Create a faceted dot plot
#'   CustomFacetDotPlot(
#'     seurat_obj = pbmc_small,
#'     features = features,
#'     group.by = "seurat_clusters",
#'     split.by = "groups",
#'     cluster_cols = TRUE,
#'     facet_ncol = 2
#'   )
#' }
#' }
CustomFacetDotPlot <- function(
    seurat_obj, 
    features, 
    group.by = "seurat_clusters", 
    split.by,
    assay = "RNA",
    slot = "data",
    dot.scale = 6,
    cols = c("blue", "white", "red"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    scale = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    facet_ncol = NULL,
    angle = 45
) {
  
  # 1. Load necessary packages
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
  
  # 2. Input Validation
  if (missing(split.by) || is.null(split.by)) {
    stop("'split.by' parameter is required for CustomFacetDotPlot.")
  }
  # ... (Validation code is identical to the main function)
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input 'seurat_obj' must be a Seurat object.")
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }
  meta_cols <- c(group.by, split.by)
  missing_cols <- meta_cols[!meta_cols %in% colnames(seurat_obj@meta.data)]
  if (length(missing_cols) > 0) {
    stop(paste("Metadata columns not found:", paste(missing_cols, collapse = ", ")))
  }
  
  # 3. Data Extraction and Preparation (Identical to CustomGroupedDotPlot)
  exp_matrix <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  features <- intersect(features, rownames(exp_matrix))
  if (length(features) == 0) {
    stop("None of the specified features were found in the selected assay.")
  }
  grouping_vars <- c(group.by, split.by)
  meta_data_subset <- seurat_obj@meta.data %>%
    select(all_of(grouping_vars)) %>%
    rownames_to_column(var = "cell_id")
  plot_data <- as.matrix(exp_matrix[features, , drop = FALSE]) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell_id") %>%
    pivot_longer(cols = all_of(features), names_to = "Gene", values_to = "Expression") %>%
    left_join(meta_data_subset, by = "cell_id")

  # 4. Calculate Summaries (Identical)
  summary_df <- plot_data %>%
    group_by(across(all_of(c(grouping_vars, "Gene")))) %>%
    summarise(
      AvgExpression = mean(Expression, na.rm = TRUE),
      PctExpressed = sum(Expression > 0) / n() * 100,
      .groups = "drop"
    )

  # 5. Ensure complete grid (Good practice for consistent axes across facets)
  grid_list <- list(
    unique(meta_data_subset[[group.by]]),
    unique(meta_data_subset[[split.by]]),
    features
  )
  names(grid_list) <- c(group.by, split.by, "Gene")
  all_combos <- expand.grid(grid_list, stringsAsFactors = FALSE)
  summary_df <- left_join(all_combos, summary_df, by = c(group.by, split.by, "Gene")) %>%
    mutate(
      AvgExpression = replace_na(AvgExpression, 0),
      PctExpressed = replace_na(PctExpressed, 0)
    )
  
  # 6. Data Scaling and Capping (Identical)
  if (scale) {
    summary_df <- summary_df %>%
      group_by(Gene) %>%
      mutate(AvgExpressionScaled = scale(AvgExpression)[,1]) %>%
      ungroup()
  } else {
    summary_df$AvgExpressionScaled <- summary_df$AvgExpression
  }
  summary_df <- summary_df %>%
    mutate(
      AvgExpressionScaled = if_else(is.nan(AvgExpressionScaled), 0, AvgExpressionScaled),
      AvgExpressionScaled = pmax(pmin(AvgExpressionScaled, col.max), col.min),
      PctExpressed = pmax(PctExpressed, dot.min)
    )

  # 7. Clustering and Factor Ordering (Identical logic for consistent axes)
  row_levels <- unique(as.character(summary_df[[group.by]]))
  col_levels <- features
  if (cluster_rows || cluster_cols) {
    value_var <- if (scale) "AvgExpressionScaled" else "AvgExpression"
    cluster_matrix_data <- summary_df %>%
        group_by(.data[[group.by]], Gene) %>%
        summarise(value = mean(.data[[value_var]], na.rm = TRUE), .groups = "drop")
    cluster_matrix <- cluster_matrix_data %>%
      pivot_wider(names_from = "Gene", id_cols = all_of(group.by), values_from = "value", values_fill = 0) %>%
      column_to_rownames(var = group.by)
    if (cluster_rows && nrow(cluster_matrix) > 1) {
      row_order <- hclust(dist(cluster_matrix))$order
      row_levels <- rownames(cluster_matrix)[row_order]
    }
    if (cluster_cols && ncol(cluster_matrix) > 1) {
      col_order <- hclust(dist(t(cluster_matrix)))$order
      col_levels <- colnames(cluster_matrix)[col_order]
    }
  } 
  summary_df[[group.by]] <- factor(summary_df[[group.by]], levels = rev(row_levels))
  summary_df$Gene <- factor(summary_df$Gene, levels = col_levels)

  # 8. Create the ggplot object with Faceting
  p <- ggplot(summary_df, aes(x = Gene, y = .data[[group.by]])) +
    geom_point(aes(size = PctExpressed, color = AvgExpressionScaled)) +
    facet_wrap(vars(.data[[split.by]]), ncol = facet_ncol)

  # 9. Apply Scales, Theme, and Labels (Identical)
  color_lab <- if(scale) "Scaled Avg.\nExpression" else "Avg.\nExpression"
  p <- p +
    scale_color_gradientn(
      colors = cols, 
      name = color_lab,
      limits = c(col.min, col.max),
      oob = scales::squish
    ) +
    scale_size_continuous(
      range = c(0, dot.scale), 
      name = "Percent\nExpressed",
      limits = c(dot.min, 100)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "grey85", colour = "black"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}







#' Create Stacked Violin or Box Plots for Multiple Genes
#'
#' This function generates violin or box plots to visualize the expression distribution
#' of multiple features across different cell groups. Plots for each feature are
#' faceted, creating a "stacked" grid for easy comparison. The function allows for
#' further splitting of distributions by a second metadata variable (e.g., treatment condition).
#'
#' @param seurat_obj A Seurat object.
#' @param features A character vector of features (e.g., genes) to plot. Each feature
#'   will be a separate facet in the plot grid.
#' @param group.by A character string specifying the metadata column to group cells by
#'   on the x-axis (e.g., "seurat_clusters", "cell_type").
#' @param split.by An optional character string specifying a metadata column to split
#'   the violins/boxes by. This will create dodged plots colored by this variable.
#' @param plot.type A character string specifying the plot type: either `"violin"` (default)
#'   or `"boxplot"`.
#' @param assay The assay to pull data from (e.g., "RNA", "SCT").
#' @param slot The data layer to use from the specified assay (e.g., "data", "counts").
#' @param cols An optional character vector of colors to use for the `split.by` groups.
#'   If `NULL`, ggplot's default colors are used.
#' @param hide.outliers A logical value. If `TRUE`, for violin plots, the tails are
#'   trimmed to the data range (`trim=TRUE`). For box plots, outlier points are hidden
#'   (`outlier.shape=NA`). Defaults to `FALSE`.
#' @param pt.size A numeric value for the size of the jitter points to overlay on the plot.
#'   If `0` (default), no points are added.
#' @param facet.ncol An integer specifying the number of columns for the facet grid.
#' @param free.y.axis A logical value. If `TRUE` (default), each facet (gene) will have its
#'   own independent y-axis scale, which is recommended for genes with different
#'   expression ranges.
#' @param angle The angle for the x-axis text labels.
#'
#' @return A `ggplot` object.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if (require("Seurat") && require("SeuratObject")) {
#'   pbmc_small$groups <- sample(c("Control", "Stim"), size = ncol(pbmc_small), replace = TRUE)
#'   features <- c("MS4A1", "NKG7", "PPBP")
#'
#'   # Stacked violin plot, split by condition, with jitter points
#'   StackedVlnPlot(
#'     seurat_obj = pbmc_small,
#'     features = features,
#'     group.by = "seurat_clusters",
#'     split.by = "groups",
#'     pt.size = 0.1
#'   )
#'
#'   # Stacked box plot, without splitting, hiding outliers
#'   StackedVlnPlot(
#'     seurat_obj = pbmc_small,
#'     features = features,
#'     group.by = "seurat_clusters",
#'     plot.type = "boxplot",
#'     hide.outliers = TRUE,
#'     facet.ncol = 1 # Arrange in a single column
#'   )
#' }
#' }
StackedVlnPlot <- function(
    seurat_obj,
    features,
    group.by = "seurat_clusters",
    split.by = NULL,
    plot.type = "violin",
    assay = "RNA",
    assay.layer = "data",
    cols = NULL,
    hide.outliers = FALSE,
    pt.size = 0,
    facet.ncol = NULL,
    free.y.axis = TRUE,
    angle = 45
) {

  # 1. Input Validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input 'seurat_obj' must be a Seurat object.")
  }
  if (!plot.type %in% c("violin", "boxplot")) {
    stop("'plot.type' must be either 'violin' or 'boxplot'.")
  }
  
  meta_cols <- c(group.by, split.by)
  missing_cols <- meta_cols[!meta_cols %in% colnames(seurat_obj@meta.data)]
  if (length(missing_cols) > 0) {
    stop(paste("Metadata columns not found:", paste(missing_cols, collapse = ", ")))
  }
  
  # 2. Data Extraction and Preparation
  # exp_matrix <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  exp_matrix <- LayerData(seurat_obj, assay = assay, layer = assay.layer)
  
  features <- intersect(features, rownames(exp_matrix))
  if (length(features) == 0) {
    stop("None of the specified features were found in the selected assay.")
  }
  
  grouping_vars <- c(group.by, split.by)
  meta_data_subset <- seurat_obj@meta.data %>%
    select(all_of(grouping_vars)) %>%
    rownames_to_column(var = "cell_id")
  
  plot_data <- as.matrix(exp_matrix[features, , drop = FALSE]) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell_id") %>%
    pivot_longer(
      cols = all_of(features),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    left_join(meta_data_subset, by = "cell_id") %>%
    mutate(Gene = factor(Gene, levels = features)) # Preserve feature order

  # 3. Create the base ggplot object
  p <- ggplot(
    plot_data, 
    aes(
      x = .data[[group.by]], 
      y = Expression, 
      fill = if (!is.null(split.by)) .data[[split.by]] else .data[[group.by]]
    )
  )

  # 4. Add jitter points if requested
  if (pt.size > 0) {
    # Jitterdodge is used when split.by is active, otherwise simple jitter
    jitter_pos <- if (!is.null(split.by)) {
      position_jitter_dodge(jitter.width = 0.2, dodge.width = 0.9)
    } else {
      position_jitter(width = 0.2)
    }
    p <- p + geom_jitter(
      size = pt.size, 
      alpha = 0.5, 
      position = jitter_pos,
      aes(color = if (!is.null(split.by)) .data[[split.by]] else .data[[group.by]])
    )
  }

  # 5. Add main plot layer (violin or boxplot)
  dodge_pos <- position_dodge(width = 0.9)
  
  if (plot.type == "violin") {
    p <- p + geom_violin(
      trim = hide.outliers, 
      scale = "width", 
      position = dodge_pos,
      alpha = 0.8
    )
  } else { # boxplot
    outlier_shape <- if (hide.outliers) NA else 19
    p <- p + geom_boxplot(
      outlier.shape = outlier_shape, 
      position = dodge_pos,
      alpha = 0.8,
      width = 0.5
    )
  }

  # 6. Apply Faceting
  facet_scales <- if (free.y.axis) "free_y" else "fixed"
  p <- p + facet_wrap(~ Gene, ncol = facet.ncol, scales = facet_scales)

  # 7. Apply Aesthetics and Theming
  p <- p + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1, size = 10),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(face = "bold", size = 12),
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position = if (!is.null(split.by)) "right" else "none"
    ) +
    labs(y = paste(assay.layer, "Expression"), fill = split.by, color = split.by)

  # Apply custom colors if provided
  if (!is.null(cols)) {
    p <- p + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
  }
  
  return(p)
}




# a modified dot plot
scDotGrouped <- function(
  obj,
  genes,
  cell_var = "celltype",
  condition_var = "condition",
  assay = NULL,
  layer = "data",
  map_size = c("avg", "pct"),
  detection_threshold = 0,
  standardize = TRUE,
  scale_range = c(-2.5, 2.5), 
  size_range = c(0, 8),                         # Start at 0 to make 0% invisible
  palette = c("#2166AC", "#F7F7F7", "#B2182B"), # Blue-White-Red
  condition_palette = NULL,                     # Custom colors for group borders
  group_mode = c("dodge", "facet"),
  dodge_width = 0.8, 
  base_size = 11,
  legend_position = "right",
  order_genes = NULL,
  order_celltypes = NULL,
  order_conditions = NULL
) {
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("scales", quietly = TRUE)
  
  map_size <- match.arg(map_size)
  group_mode <- match.arg(group_mode)
  
  # --- 1. Data Prep ---
  o <- obj
  if (!is.null(assay)) Seurat::DefaultAssay(o) <- assay else assay <- Seurat::DefaultAssay(o)
  
  avail_feats <- rownames(o)
  valid_genes <- genes[genes %in% avail_feats]
  if (length(valid_genes) == 0) stop("No genes found.")
  
  meta <- o@meta.data
  meta[[cell_var]] <- trimws(as.character(meta[[cell_var]]))
  meta[[condition_var]] <- trimws(as.character(meta[[condition_var]]))
  
  if (is.null(order_celltypes)) order_celltypes <- sort(unique(meta[[cell_var]]))
  if (is.null(order_conditions)) order_conditions <- sort(unique(meta[[condition_var]]))
  if (is.null(order_genes)) order_genes <- valid_genes
  
  vars_to_fetch <- c(valid_genes, cell_var, condition_var)
  
  df <- tryCatch({
    Seurat::FetchData(o, vars = vars_to_fetch, layer = layer)
  }, error = function(e) {
    slot_use <- if(layer == "data") "data" else if(layer == "counts") "counts" else "scale.data"
    Seurat::FetchData(o, vars = vars_to_fetch, slot = slot_use)
  })
  
  # --- 2. Calculation ---
  long_df <- df %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(valid_genes), names_to = "gene", values_to = "expr") %>%
    dplyr::rename(cell_type = !!cell_var, condition = !!condition_var)
  
  stat_df <- long_df %>%
    dplyr::group_by(gene, cell_type, condition) %>%
    dplyr::summarise(
      avg_exp = mean(expr, na.rm = TRUE),
      pct_exp = mean(expr > detection_threshold, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Complete grid (needed for dodging to work correctly)
  full_grid <- tidyr::expand_grid(gene = order_genes, cell_type = order_celltypes, condition = order_conditions)
  plot_df <- full_grid %>%
    dplyr::left_join(stat_df, by = c("gene", "cell_type", "condition")) %>%
    dplyr::mutate(avg_exp = tidyr::replace_na(avg_exp, 0), pct_exp = tidyr::replace_na(pct_exp, 0))
  
  # Standardize
  if (standardize) {
    plot_df <- plot_df %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(
        z = (avg_exp - mean(avg_exp, na.rm=TRUE)) / sd(avg_exp, na.rm=TRUE),
        z = ifelse(is.na(z), 0, z),
        fill_val = scales::oob_squish(z, range = scale_range)
      ) %>%
      dplyr::ungroup()
    fill_lab <- "Scaled\nExpr"
  } else {
    plot_df$fill_val <- plot_df$avg_exp
    fill_lab <- "Avg\nExpr"
  }
  
  if (map_size == "avg") {
    plot_df$size_val <- plot_df$avg_exp; size_lab <- "Avg Expr"
  } else {
    plot_df$size_val <- plot_df$pct_exp; size_lab <- "% Expr"
  }
  
  # Set Factors
  plot_df$gene <- factor(plot_df$gene, levels = order_genes)
  plot_df$cell_type <- factor(plot_df$cell_type, levels = rev(order_celltypes))
  plot_df$condition <- factor(plot_df$condition, levels = order_conditions)
  
  # --- 3. Plotting ---
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = gene, y = cell_type))
  
  if (group_mode == "dodge") {
    p <- p + 
      ggplot2::geom_point(
        ggplot2::aes(
          size = size_val, 
          fill = fill_val, 
          color = condition, # Map color to condition
          group = condition
        ),
        shape = 21, 
        stroke = 0.2,        # <--- KEY FIX: Make border hairline thin on the plot
        position = ggplot2::position_dodge(width = dodge_width),
        na.rm = TRUE
      )
  } else {
    p <- p + 
      ggplot2::geom_point(
        ggplot2::aes(size = size_val, fill = fill_val, color = condition),
        shape = 21, stroke = 0.2
      ) +
      ggplot2::facet_grid(cols = ggplot2::vars(condition))
  }
  
  # --- 4. Aesthetics & Legends ---
  
  # Apply Custom Border Colors if provided
  if (!is.null(condition_palette)) {
    p <- p + ggplot2::scale_color_manual(values = condition_palette, name = condition_var)
  } else {
    p <- p + ggplot2::scale_color_discrete(name = condition_var)
  }
  
  p <- p +
    # Expression Fill
    ggplot2::scale_fill_gradient2(
      low = palette[1], mid = palette[2], high = palette[3], 
      midpoint = if(standardize) 0 else mean(range(plot_df$fill_val)), 
      limits = if(standardize) scale_range else NULL,
      name = fill_lab,
      oob = scales::squish
    ) +
    # Size Scale (Start range at 0 so 0% is invisible)
    ggplot2::scale_size_continuous(range = size_range, name = size_lab) +
    
    # Theme
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_blank(),
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_line(color = "grey95"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white")
    ) +
    
    # <--- KEY FIX: Force the LEGEND to be big and visible, even though plot borders are thin
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4, stroke = 1.5, fill = "transparent") 
      ),
      size = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_colorbar(order = 2)
    )
  
  return(p)
}




# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/
# group1 is the cell type/cluster annotation 
# group2 is any condition you want to further group, in this case, the stim 
GetMatrixFromSeuratByGroupSingle<- function(obj, feature, group1, group2){
  if (length(feature) != 1){
          stop("please only provide only one gene name")
  }
  exp_mat<- obj@assays$RNA@data[feature, ,drop=FALSE]
  count_mat<- obj@assays$RNA@counts[feature,,drop=FALSE ]
  
  meta<- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell")
        
  # get the average expression matrix
  exp_df<- as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene,{{group1}}, {{group2}}) %>%
    summarise(average_expression = mean(expression)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= average_expression) 
  
  exp_mat<- exp_df[, -c(1,2)] %>% as.matrix()
  rownames(exp_mat)<- exp_df %>% pull({{group2}})
  
  # get the percentage positive cell matrix
  count_df<- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, {{group1}}, {{group2}}) %>%
    summarise(percentage = mean(count >0)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= percentage) 

  percent_mat<- count_df[, -c(1,2)] %>% as.matrix()
  rownames(percent_mat)<- count_df %>% pull({{group2}})
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  
  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}


library(Seurat)
library(dplyr)
library(Matrix)

GetMatrixFromSeuratByGroupMulti <- function(obj, features, group1, group2, assay = NULL){
  
  # 1. 确定 Assay (兼容 V5)
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  
  # 检查 features 是否都在对象里
  available_features <- features[features %in% rownames(obj[[assay]])]
  if(length(available_features) < length(features)){
    warning(paste("部分基因未找到，已忽略:", paste(setdiff(features, available_features), collapse = ",")))
  }
  features <- available_features
  
  # 2. 安全获取矩阵 (Seurat V5 核心修改)
  # 使用 GetAssayData 会自动合并 V5 中的 layers (如 data.1, data.2)
  # 即使数据在磁盘上 (BPCells)，提取少量基因也不会太慢
  exp_mat <- GetAssayData(obj, assay = assay, layer = "data")[features, , drop=FALSE]
  count_mat <- GetAssayData(obj, assay = assay, layer = "counts")[features, , drop=FALSE]
  
  # 确保转换为内存中的稀疏矩阵或普通矩阵 (防止 BPCells 计算问题)
  exp_mat <- as(exp_mat, "CsparseMatrix")
  count_mat <- as(count_mat, "CsparseMatrix")
  
  # 3. 准备分组信息
  # 使用 rlang 处理非标准评估 ({{group}})
  meta <- obj@meta.data %>%
    dplyr::select({{group1}}, {{group2}}) %>%
    dplyr::mutate(
      g1 = as.character({{group1}}),
      g2 = as.character({{group2}}),
      # 拼接分组名，即原来的 pivot_wider 列名
      combined_group = paste(g1, g2, sep = "|") 
    )
  
  # 获取所有唯一的组合
  unique_groups <- sort(unique(meta$combined_group))
  
  # 4. 初始化结果矩阵
  res_exp_mat <- matrix(0, nrow = length(features), ncol = length(unique_groups))
  res_pct_mat <- matrix(0, nrow = length(features), ncol = length(unique_groups))
  
  rownames(res_exp_mat) <- features
  colnames(res_exp_mat) <- unique_groups
  rownames(res_pct_mat) <- features
  colnames(res_pct_mat) <- unique_groups
  
  # 5. 循环计算 (比 pivot_longer 快得多)
  for (g in unique_groups) {
    # 找到属于该组的细胞 ID
    cells_in_group <- rownames(meta)[meta$combined_group == g]
    
    if (length(cells_in_group) > 0) {
      # 提取子矩阵
      sub_exp <- exp_mat[, cells_in_group, drop=FALSE]
      sub_count <- count_mat[, cells_in_group, drop=FALSE]
      
      # 计算平均表达量 (mean of data slot)
      if (length(cells_in_group) == 1) {
        res_exp_mat[, g] <- as.vector(sub_exp)
        res_pct_mat[, g] <- as.vector(ifelse(sub_count > 0, 1, 0))
      } else {
        # Matrix::rowMeans 对稀疏矩阵非常快
        res_exp_mat[, g] <- Matrix::rowMeans(sub_exp)
        res_pct_mat[, g] <- Matrix::rowMeans(sub_count > 0)
      }
    }
  }
  
  # 6. 最后的检查 (逻辑上应该总是满足，但保留原有检查风格)
  if (!identical(dim(res_exp_mat), dim(res_pct_mat))) {
    stop("Dimension mismatch!")
  }
  
  return(list(exp_mat = res_exp_mat, percent_mat = res_pct_mat))
}





#' Grouped Dot Plot using ComplexHeatmap (Fixed Label Overlap)
#' 
#' @param obj Seurat object.
#' @param genes Character vector of genes.
#' @param cell_var Meta column for cell types.
#' @param condition_var Meta column for conditions.
#' @param assay Assay to use.
#' @param layer "data", "counts", or "scale.data".
#' @param scale_data Logical. If TRUE, z-score expression.
#' @param scale_range Vector c(min, max).
#' @param cluster_genes Logical.
#' @param exp_colors Vector of 3 colors (Low, Mid, High).
#' @param condition_colors Named vector for condition annotation.
#' @param dot_size_max Max dot size (mm).
#' @param order_celltypes Explicit order for cell types.
#' @param order_conditions Explicit order for conditions.
#' @param celltype_rot Rotation for top cell type labels (e.g., 0, 45, 90).
#' @param celltype_size Font size for top cell type labels.
#' 
#' @return ComplexHeatmap object.
#' @export
scDotComplex <- function(
  obj,
  genes,
  cell_var = "celltype",
  condition_var = "condition",
  assay = NULL,
  layer = "data",
  scale_data = TRUE,
  scale_range = c(-2.5, 2.5),
  cluster_genes = FALSE,
  exp_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
  condition_colors = NULL,
  dot_size_max = 4,
  order_celltypes = NULL,
  order_conditions = NULL,
  celltype_rot = 45,  # <--- 新增：默认旋转45度，避免重叠
  celltype_size = 10  # <--- 新增：控制字体大小
) {
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  requireNamespace("circlize", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)

  # --- 1. Data Prep ---
  if (!is.null(assay)) Seurat::DefaultAssay(obj) <- assay
  avail <- rownames(obj)
  genes <- genes[genes %in% avail]
  if(length(genes) == 0) stop("No valid genes found.")
  
  meta <- obj@meta.data
  meta[[cell_var]] <- trimws(as.character(meta[[cell_var]]))
  meta[[condition_var]] <- trimws(as.character(meta[[condition_var]]))
  
  vars <- c(genes, cell_var, condition_var)
  df <- tryCatch({
    Seurat::FetchData(obj, vars = vars, layer = layer)
  }, error = function(e) {
    slot_use <- if(layer == "data") "data" else if(layer == "counts") "counts" else "scale.data"
    Seurat::FetchData(obj, vars = vars, slot = slot_use)
  })
  
  # --- 2. Calculation ---
  long_df <- df %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(genes), names_to = "gene", values_to = "expr") %>%
    dplyr::rename(ct = !!cell_var, cond = !!condition_var) %>%
    dplyr::group_by(gene, ct, cond) %>%
    dplyr::summarise(
      avg = mean(expr, na.rm=TRUE),
      pct = mean(expr > 0, na.rm=TRUE),
      .groups = "drop"
    )

  if (is.null(order_celltypes)) order_celltypes <- sort(unique(long_df$ct))
  if (is.null(order_conditions)) order_conditions <- sort(unique(long_df$cond))
  
  full_grid <- tidyr::expand_grid(gene = genes, ct = order_celltypes, cond = order_conditions)
  
  plot_df <- full_grid %>%
    dplyr::left_join(long_df, by = c("gene", "ct", "cond")) %>%
    dplyr::mutate(
      avg = tidyr::replace_na(avg, 0),
      pct = tidyr::replace_na(pct, 0),
      col_id = paste(ct, cond, sep = "|||")
    )

  # --- 3. Matrix ---
  mat_exp <- plot_df %>% 
    dplyr::select(gene, col_id, avg) %>% 
    tidyr::pivot_wider(names_from = col_id, values_from = avg) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  mat_pct <- plot_df %>% 
    dplyr::select(gene, col_id, pct) %>% 
    tidyr::pivot_wider(names_from = col_id, values_from = pct) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  mat_exp <- mat_exp[genes, , drop=FALSE]
  mat_pct <- mat_pct[genes, , drop=FALSE]
  
  if (scale_data) {
    mat_scaled <- t(scale(t(mat_exp)))
    mat_scaled[mat_scaled > scale_range[2]] <- scale_range[2]
    mat_scaled[mat_scaled < scale_range[1]] <- scale_range[1]
    mat_scaled[is.na(mat_scaled)] <- 0
    mat_exp <- mat_scaled
    legend_title <- "Z-score"
  } else {
    legend_title <- "Avg Expr"
  }

  # --- 4. Annotation ---
  col_meta <- data.frame(id = colnames(mat_exp)) %>%
    tidyr::separate(id, into = c("ct", "cond"), sep = "\\|\\|\\|", remove = FALSE)
  col_meta$ct <- factor(col_meta$ct, levels = order_celltypes)
  col_meta$cond <- factor(col_meta$cond, levels = order_conditions)
  
  ord_idx <- order(col_meta$ct, col_meta$cond)
  mat_exp <- mat_exp[, ord_idx, drop=FALSE]
  mat_pct <- mat_pct[, ord_idx, drop=FALSE]
  col_meta <- col_meta[ord_idx, ]
  
  if (is.null(condition_colors)) {
    conds <- levels(col_meta$cond)
    condition_colors <- stats::setNames(scales::hue_pal()(length(conds)), conds)
  }
  
  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition = col_meta$cond,
    col = list(Condition = condition_colors),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(3, "mm")
  )
  
  # --- 5. Drawing ---
  col_fun <- circlize::colorRamp2(breaks = c(min(mat_exp), 0, max(mat_exp)), colors = exp_colors)
  
  cell_fun_logic <- function(j, i, x, y, width, height, fill) {
    pct_val <- mat_pct[i, j]
    if (pct_val > 0) {
      grid::grid.circle(
        x = x, y = y, 
        r = grid::unit(sqrt(pct_val) * (dot_size_max/2), "mm"), 
        gp = grid::gpar(fill = col_fun(mat_exp[i, j]), col = NA)
      )
    }
  }

  hm <- ComplexHeatmap::Heatmap(
    mat_exp,
    name = legend_title,
    
    # 分块设置
    column_split = col_meta$ct,
    cluster_column_slices = FALSE,
    
    # --- 这里修复标题重叠问题 ---
    column_title_rot = celltype_rot,   # 旋转角度 (0, 45, 90)
    column_title_gp = grid::gpar(fontface = "bold", fontsize = celltype_size),
    # ---------------------------
    
    cluster_rows = cluster_genes,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(type = "none"),
    cell_fun = cell_fun_logic,
    col = col_fun,
    top_annotation = top_anno,
    column_labels = col_meta$cond,
    row_names_side = "left",
    column_names_rot = 45, # 底部 Condition 的旋转角度
    column_names_gp = gpar(fontsize = 8),
    border = TRUE
  )
  
  lgd_size <- ComplexHeatmap::Legend(
    labels = c("25%", "50%", "75%", "100%"), title = "Pct Expr", type = "points", pch = 16,
    legend_gp = grid::gpar(col = "grey50"), size = grid::unit(sqrt(c(0.25, 0.5, 0.75, 1)) * dot_size_max, "mm"),
    background = "white"
  )

  ComplexHeatmap::draw(hm, annotation_legend_list = list(lgd_size))
}


#' Marker-module heatmap with side-by-side enrichment bar plot
#'
#' Builds the common single-cell figure of: rows = marker genes (split into
#' modules), columns = cell types, with a coloured block annotation per module
#' on the left and a horizontal bar plot of pathway / GO enrichment per module
#' on the right.
#'
#' @param expr_matrix Numeric matrix, rows = genes, cols = cell types. Typically
#'   row-scaled (z-scores) average expression.
#' @param gene_groups Named list of character vectors. Each element is one
#'   module; order of the list defines row block order.
#' @param enrichment_df data.frame with columns `gene_group` (matching
#'   names(gene_groups)), `term` (GO/pathway label), `neg_log10_padj`.
#' @param group_colors Named vector of fill colours for the modules. Default is
#'   a 7-colour palette recycled as needed.
#' @param expr_col ComplexHeatmap colour function for the expression scale.
#' @param column_colors Optional named vector of label colours for the column
#'   (cell-type) names.
#' @param row_height_mm Per-row height of the heatmap in mm. Default 4.5.
#' @param rel_widths Relative widths of (heatmap, bar plot) panels. Default
#'   c(4.2, 3).
#' @param outfile Optional path; if given, the combined figure is saved as PDF.
#'
#' @return A `cowplot` plot_grid object (also drawn to the active device).
plotMarkerEnrichmentHeatmap <- function(expr_matrix,
                                        gene_groups,
                                        enrichment_df,
                                        group_colors = NULL,
                                        expr_col = circlize::colorRamp2(c(-2, 0, 2),
                                                                        c("blue", "white", "red")),
                                        column_colors = NULL,
                                        row_height_mm = 4.5,
                                        rel_widths = c(4.2, 3),
                                        outfile = NULL) {
  stopifnot(is.matrix(expr_matrix),
            is.list(gene_groups), !is.null(names(gene_groups)),
            all(c("gene_group", "term", "neg_log10_padj") %in% colnames(enrichment_df)))

  if (is.null(group_colors)) {
    pal <- c("#F9A23D", "#F47C40", "#E7A56F", "#8CC37A",
             "#79BC5C", "#F7968C", "#ED6564")
    group_colors <- setNames(pal[((seq_along(gene_groups) - 1) %% length(pal)) + 1],
                             names(gene_groups))
  }

  row_groups <- rep(names(gene_groups), lengths(gene_groups))
  row_label_colors <- group_colors[row_groups]
  expr_matrix <- expr_matrix[unlist(gene_groups), , drop = FALSE]

  left_anno <- ComplexHeatmap::rowAnnotation(
    module = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = group_colors),
      labels = names(gene_groups),
      labels_gp = grid::gpar(col = "black", fontsize = 8)
    )
  )

  hm <- ComplexHeatmap::Heatmap(
    expr_matrix, name = "Expression", col = expr_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = factor(row_groups, levels = names(gene_groups)),
    row_gap = grid::unit(0, "mm"),
    show_row_names = TRUE, row_names_side = "left",
    column_names_side = "bottom",
    row_names_gp = grid::gpar(col = row_label_colors, fontsize = 8),
    column_names_gp = if (!is.null(column_colors))
      grid::gpar(col = column_colors[colnames(expr_matrix)], fontsize = 8)
    else grid::gpar(fontsize = 8),
    height = grid::unit(nrow(expr_matrix) * row_height_mm, "mm"),
    left_annotation = left_anno
  )

  enrichment_df$gene_group <- factor(enrichment_df$gene_group, levels = names(gene_groups))
  enrichment_df$term <- factor(enrichment_df$term, levels = rev(enrichment_df$term))

  bar <- ggplot2::ggplot(enrichment_df,
                         ggplot2::aes(x = neg_log10_padj, y = term, fill = gene_group)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  combined <- cowplot::plot_grid(
    grid::grid.grabExpr(ComplexHeatmap::draw(hm)),
    bar, ncol = 2, rel_widths = rel_widths, align = "h", axis = "tb"
  )

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, combined,
                    width = sum(rel_widths) + 1.5,
                    height = nrow(expr_matrix) * 0.25 + 1.5)
  }
  combined
}


#' Convert genes between human and mouse via Ensembl biomaRt
#'
#' @param genes character vector of gene IDs to convert
#' @param geneid one of "symbol", "ensembl", "entrez"
#' @param invert if TRUE, human -> mouse; if FALSE, mouse -> human
#' @param host biomaRt host. Use an archive (e.g. "https://aug2020.archive.ensembl.org")
#'   for a frozen GRCh38/GRCm38 mapping; default is the current Ensembl.
#' @return data.frame with the source and target gene IDs
convertHumanMouse <- function(genes,
                              geneid = c("symbol", "ensembl", "entrez"),
                              invert = FALSE,
                              host = "https://www.ensembl.org") {
  geneid <- match.arg(geneid)
  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
  mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)

  ids <- switch(geneid,
                symbol  = list(human = "hgnc_symbol",     mouse = "mgi_symbol"),
                ensembl = list(human = "ensembl_gene_id", mouse = "ensembl_gene_id"),
                entrez  = list(human = "entrezgene_id",   mouse = "entrezgene_id"))

  if (invert) {
    biomaRt::getLDS(attributes = ids$human, filters = ids$human, values = genes,
                    mart = human, attributesL = ids$mouse, martL = mouse,
                    uniqueRows = TRUE)
  } else {
    biomaRt::getLDS(attributes = ids$mouse, filters = ids$mouse, values = genes,
                    mart = mouse, attributesL = ids$human, martL = human,
                    uniqueRows = TRUE)
  }
}

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

