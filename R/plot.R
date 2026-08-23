#' Plot Quality Control Metrics
#'
#' @description Generates a grouped violin plot for quality control metrics (e.g., nFeature_RNA, nCount_RNA, percent_mt)
#' across specified cell groups. Automatically handles color interpolation for large numbers of groups.
#'
#' @param obj A \code{RNA} object containing calculated QC metrics in \code{meta.data}.
#' @param group_col Character. The column name in \code{meta.data} to group cells by. If NULL or not found, defaults to "orig.ident".
#' @param features Character vector. The QC metrics to plot. Defaults to \code{c("nFeature_RNA", "nCount_RNA", "percent_mt")}.
#'
#' @return A \code{ggplot} object representing the QC metrics violin plots.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter facet_wrap theme_bw labs theme element_text element_rect element_blank margin scale_fill_manual
#'
#' @examples
#' # 1. Create a minimal mock RNA object with required meta.data
#' mock_meta <- data.frame(
#'   orig.ident = rep(c("Sample_A", "Sample_B"), each = 20),
#'   nFeature_RNA = runif(40, min = 200, max = 2500),
#'   nCount_RNA = runif(40, min = 500, max = 8000),
#'   percent_mt = runif(40, min = 0, max = 15)
#' )
#' rownames(mock_meta) <- paste0("Cell_", 1:40)
#'
#' mock_obj <- list(meta.data = mock_meta)
#' class(mock_obj) <- "RNA"
#'
#' # 2. Generate QC plot
#' p_qc <- PlotQC_RNA(mock_obj, group_col = "orig.ident")
#' p_qc
#' @export
PlotQC_RNA <- function(obj,
                       group_col = NULL,
                       features = c("nFeature_RNA", "nCount_RNA", "percent_mt")) {

  if (!inherits(obj, "RNA")) stop("The input must be a RNA object!")

  meta <- obj$meta.data
  if (is.null(meta)) stop(" meta.data not found. Please check that the object is complete.")

  # --- 1. Check and extract the feature columns to be plotted ---
  missing_features <- setdiff(features, colnames(meta))
  if (length(missing_features) > 0) {
    stop(sprintf("Error: The feature column %s was not found in meta.data. Please run the QC calculation step first.",
                 paste(missing_features, collapse = ", ")))
  }

  # --- 2. Safely retrieve the grouping column name for plotting ---
  if (is.null(group_col) || !(group_col %in% colnames(meta))) {
    if ("orig.ident" %in% colnames(meta)) {
      group_col <- "orig.ident"
      message("No valid grouping column specified or column name does not exist; 'orig.ident' is used by default.")
    } else {
      # If even `orig.ident` is missing, create a virtual group for all cells
      meta$All_Cells <- "All_Cells"
      group_col <- "All_Cells"
      message("No grouping information found; all cells will be displayed as a single group.")
    }
  }

  # --- 3. Data Cleaning and Factorization ---
  # Extract the required columns to avoid loading the entire large metadata
  plot_data <- meta[, c(group_col, features), drop = FALSE]
  colnames(plot_data)[1] <- "Plot_Group" # Standardize naming for ggplot mapping

  # Prevent ggplot2 from crashing due to NA or empty strings
  plot_data$Plot_Group <- as.character(plot_data$Plot_Group)
  plot_data$Plot_Group[is.na(plot_data$Plot_Group) | plot_data$Plot_Group == ""] <- "Unknown"
  plot_data$Plot_Group <- as.factor(plot_data$Plot_Group)

  # --- 4. Converting a wide table to a long table (specifically designed for ggplot faceting) ---
  qc_long <- tidyr::pivot_longer(
    data = plot_data,
    cols = all_of(features),
    names_to = "Metric",
    values_to = "Value"
  )

  # Fix the panel order to match the order of user-input features
  qc_long$Metric <- factor(qc_long$Metric, levels = features)

  # --- 5. Dynamically expand the color palette ---
  sci_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                   "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                   "#7E6148FF", "#B09C85FF", "#FF7F00", "#6A3D9A")

  num_groups <- length(levels(qc_long$Plot_Group))

  if (num_groups <= length(sci_palette)) {
    # If the number of groups is <= 12, extract the primary colors strictly in order
    dynamic_colors <- sci_palette[1:num_groups]
  } else {
    # If the number of groups is > 12, generate a new set of colors with smooth transitions using an interpolation algorithm
    dynamic_colors <- grDevices::colorRampPalette(sci_palette)(num_groups)
  }

  # --- 6. Advanced ggplot2 native rendering output ---
  p_qc <- ggplot2::ggplot(qc_long, ggplot2::aes(x = Plot_Group, y = Value, fill = Plot_Group)) +
    ggplot2::geom_violin(trim = FALSE, color = "black", alpha = 0.9, scale = "width", linewidth = 0.6) +
    ggplot2::geom_jitter(width = 0.2, size = 0.5, alpha = 0.3, color = "grey30") +

    ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = length(features), strip.position = "top") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(x = group_col, y = "Metric Value") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(t = 10)),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),

      strip.background = ggplot2::element_rect(fill = "#f8f9fa", color = "black", linewidth = 1),
      strip.text = ggplot2::element_text(face = "bold", size = 14, color = "#2C3E50"),

      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1.2)
    ) +
    ggplot2::scale_fill_manual(values = dynamic_colors)

  return(p_qc)
}

#' Plot Dimensionality Reduction (PCA, t-SNE, UMAP)
#'
#' @description Generates scatter plots for dimensionality reduction coordinates.
#' Supports coloring by custom metadata and optionally juxtaposing with unsupervised clustering results.
#' Uses distinct color palettes for the metadata plot and the clustering plot to avoid visual confusion.
#'
#' @param obj A \code{RNA} object containing reductions and metadata.
#' @param reduction Character. The reduction to plot: "umap", "tsne", or "pca".
#' @param group_col Character. The column name in metadata to color the points by (e.g., "Condition", "Batch").
#' @param show_cluster Logical. If TRUE, plots a side-by-side comparison with 'Auto_Cluster' (requires patchwork package).
#'
#' @return A \code{ggplot} or \code{patchwork} object containing the generated plot(s).
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic labs theme element_text guides guide_legend
#'
#' @examples
#' # 1. Create minimal mock RNA object with reductions and metadata
#' set.seed(42)
#' mock_umap <- matrix(rnorm(80), ncol = 2)
#' rownames(mock_umap) <- paste0("Cell_", 1:40)
#' colnames(mock_umap) <- c("UMAP_1", "UMAP_2")
#'
#' mock_meta <- data.frame(
#'   orig.ident = rep(c("Control", "Treatment"), each = 20),
#'   Auto_Cluster = rep(paste0("Cluster_", 1:4), times = 10)
#' )
#' rownames(mock_meta) <- rownames(mock_umap)
#'
#' mock_obj <- list(
#'   reductions = list(umap = mock_umap),
#'   meta.data = mock_meta
#' )
#' class(mock_obj) <- "RNA"
#'
#' # 2. Generate plot
#' p_dim <- PlotDimRed_RNA(mock_obj, reduction = "umap", group_col = "orig.ident")
#' p_dim
#' @export
PlotDimRed_RNA <- function(obj,
                           reduction = "umap",
                           group_col = "orig.ident",
                           show_cluster = FALSE) {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # --- 1. Extract dimensionality reduction coordinates ---
  # Standardize naming format to prevent errors
  red_key <- tolower(reduction)
  if (red_key == "t-sne") red_key <- "tsne"

  coords <- obj$reductions[[red_key]]
  if (is.null(coords)) {
    stop(sprintf("Error: Cannot find coordinate data for '%s'. Please run RunDimReduction_RNA first.", reduction))
  }

  # --- 2. Extract Metadata and align with coordinates ---
  # Prioritize using filtered metadata
  meta <- if (!is.null(obj$filter_meta.data)) obj$filter_meta.data else obj$meta.data

  # Defensive programming: only take samples present in both to prevent cbind dimension mismatch errors
  common_cells <- intersect(rownames(coords), rownames(meta))
  if (length(common_cells) == 0) stop(" Sample names in the coordinate matrix and Metadata do not match at all!")

  coords_sub <- coords[common_cells, 1:2, drop = FALSE] # Force taking the first two dimensions for 2D plotting
  meta_sub <- meta[common_cells, , drop = FALSE]

  plot_df <- cbind(as.data.frame(coords_sub), meta_sub)
  dim_names <- colnames(coords_sub) # Get axis names, e.g., "UMAP_1", "UMAP_2"

  if (!(group_col %in% colnames(plot_df))) {
    stop(sprintf("Error: Column '%s' not found in Metadata.", group_col))
  }

  # --- 3. Dynamic dual independent color generation engine ---
  # Palette 1 (Classic NPG style, for plot 1: custom grouping)
  palette_1 <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                 "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                 "#7E6148FF", "#B09C85FF", "#FF7F00", "#6A3D9A")

  # Palette 2 (Classic D3 style, high contrast, for plot 2: Auto_Cluster)
  palette_2 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                 "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                 "#bcbd22", "#17becf", "#aec7e8", "#ffbb78")

  # Upgraded color extraction function, accepts specified palette as parameter
  get_dynamic_colors <- function(n, base_palette) {
    if (n <= length(base_palette)) {
      return(base_palette[1:n])
    } else {
      return(grDevices::colorRampPalette(base_palette)(n)) # Auto-interpolate when exceeding 12
    }
  }

  # --- 4. Build plot 1 (colored based on user-specified group, using palette_1) ---
  plot_df[[group_col]] <- as.factor(as.character(plot_df[[group_col]]))
  num_colors_1 <- length(levels(plot_df[[group_col]]))

  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[dim_names[1]]], y = .data[[dim_names[2]]], color = .data[[group_col]])) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8, stroke = 0) + # stroke = 0 to remove artifact edges around points
    ggplot2::scale_color_manual(values = get_dynamic_colors(num_colors_1, palette_1)) + # <- Apply palette 1
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title = paste(toupper(reduction), "by", group_col),
                  x = dim_names[1], y = dim_names[2]) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position = "right") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))

  # --- 5. Build plot 2 and combine (if cluster comparison is enabled, use palette_2) ---
  if (show_cluster) {
    if (!("Auto_Cluster" %in% colnames(plot_df))) {
      warning("Column 'Auto_Cluster' not found! Please run RunClustering_RNA first. Only a single plot will be output.")
      return(p1)
    }

    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Please install the 'patchwork' package (install.packages('patchwork')) to support side-by-side plotting. Only a single plot will be output.")
      return(p1)
    }

    plot_df$Auto_Cluster <- as.factor(as.character(plot_df$Auto_Cluster))
    num_colors_2 <- length(levels(plot_df$Auto_Cluster))

    p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[dim_names[1]]], y = .data[[dim_names[2]]], color = Auto_Cluster)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.8, stroke = 0) +
      ggplot2::scale_color_manual(values = get_dynamic_colors(num_colors_2, palette_2)) + # <- Apply palette 2 (Cluster-specific colors)
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::labs(title = "Unsupervised Clustering",
                    x = dim_names[1], y = dim_names[2]) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                     legend.position = "right") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))

    # Use patchwork for side-by-side rendering
    final_plot <- p1 + p2 + patchwork::plot_layout(ncol = 2, guides = "collect")
    return(final_plot)
  }

  # --- 6. Return single plot ---
  return(p1)
}

#' Plot Volcano Plot for Differential Expression Analysis
#'
#' @description Generates an academic-grade volcano plot based on differential expression analysis results.
#' Automatically highlights and labels top significant up- and down-regulated genes.
#'
#' @param res_df A data.frame containing DEA results. Must contain 'avg_log2FC' and 'p_val_adj' columns.
#' @param fc_cut Numeric. The absolute log2 fold change threshold. Default is 0.5.
#' @param p_cut Numeric. The adjusted p-value threshold. Default is 0.05.
#' @param top_n Integer. The number of top up/down regulated genes to label. Default is 10.
#' @param plot_title Character. Optional custom title. If NULL, auto-generates one.
#'
#' @return A \code{ggplot} object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline geom_point scale_fill_manual scale_size_manual scale_alpha_manual scale_color_manual theme_bw labs theme element_blank element_rect element_text element_line margin guide_legend arrow unit
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' # 1. Generate mock Differential Expression Analysis (DEA) results
#' set.seed(123)
#' mock_dea <- data.frame(
#'   gene = paste0("Gene_", 1:100),
#'   avg_log2FC = rnorm(100, mean = 0, sd = 1.5),
#'   p_val_adj = runif(100, min = 0, max = 0.1)
#' )
#' # Make a few genes highly significant for the plot
#' mock_dea$avg_log2FC[1:5] <- runif(5, 2, 4)
#' mock_dea$p_val_adj[1:5] <- runif(5, 1e-10, 1e-5)
#'
#' # 2. Generate Volcano Plot
#' p_volcano <- PlotVolcano_RNA(mock_dea, fc_cut = 1, p_cut = 0.05, top_n = 3)
#' p_volcano
#' @export
PlotVolcano_RNA <- function(res_df,
                            fc_cut = 1,
                            p_cut = 0.05,
                            top_n = 10,
                            plot_title = NULL) {

  # --- 1. Input data validation ---
  if (!is.data.frame(res_df)) stop("res_df must be a data.frame.")

  req_cols <- c("avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(req_cols, colnames(res_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Error: Data.frame is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Process gene names: if 'gene' column does not exist, try extracting rownames
  if (!"gene" %in% colnames(res_df)) {
    res_df$gene <- rownames(res_df)
  }

  # --- 2. Underlying data processing (specifically for visualization) ---
  # [Core Fix] Prevent log10 overflow causing plot crash when P-value is 0
  min_nonzero_p <- min(res_df$p_val_adj[res_df$p_val_adj > 0], na.rm = TRUE)
  if (is.infinite(min_nonzero_p) || is.na(min_nonzero_p)) min_nonzero_p <- 1e-300
  res_df$p_val_adj[res_df$p_val_adj == 0] <- min_nonzero_p

  # 3. Define color grouping logic
  res_df$Significance <- "NS"
  res_df$Significance[res_df$avg_log2FC > fc_cut & res_df$p_val_adj < p_cut] <- "UP"
  res_df$Significance[res_df$avg_log2FC < -fc_cut & res_df$p_val_adj < p_cut] <- "DOWN"
  res_df$Significance <- factor(res_df$Significance, levels = c("UP", "DOWN", "NS"))

  # 4. Smart extraction of Top genes for labeling (Base R approach, avoiding dplyr dependency)
  up_genes <- res_df[res_df$Significance == "UP", ]
  up_genes <- up_genes[order(up_genes$avg_log2FC, decreasing = TRUE), ]
  top_up <- head(up_genes, top_n)

  down_genes <- res_df[res_df$Significance == "DOWN", ]
  down_genes <- down_genes[order(down_genes$avg_log2FC, decreasing = FALSE), ]
  top_down <- head(down_genes, top_n)

  top_genes <- rbind(top_up, top_down)

  # 5. Build title
  if (is.null(plot_title)) {
    comparison_name <- if ("comparison" %in% colnames(res_df)) unique(res_df$comparison)[1] else "Differential Expression"
    plot_title <- paste("Volcano Plot:", comparison_name)
  }
  sub_title <- paste("Thresholds: |Log2FC| >", fc_cut, " &  adj.P <", p_cut)

  # 6. Define top-tier journal colors (similar to NPG - Nature Publishing Group)
  my_colors <- c("UP" = "#E64B35", "DOWN" = "#4DBBD5", "NS" = "#DFE6E9")

  # --- 3. Advanced visualization rendering ---
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj))) +

    # [Beautification] Draw threshold auxiliary lines at the bottom layer
    ggplot2::geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40", linewidth = 0.6, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey40", linewidth = 0.6, alpha = 0.8) +

    # [Beautification] Use shape=21 to give points a white outline, enhancing 3D effect
    ggplot2::geom_point(ggplot2::aes(fill = Significance, size = Significance, alpha = Significance),
                        shape = 21, color = "white", stroke = 0.3) +

    # Map colors, sizes, and transparencies
    ggplot2::scale_fill_manual(values = my_colors) +
    ggplot2::scale_size_manual(values = c("UP" = 2.5, "DOWN" = 2.5, "NS" = 1.2)) +
    ggplot2::scale_alpha_manual(values = c("UP" = 0.9, "DOWN" = 0.9, "NS" = 0.4)) +

    # [Beautification] Add non-overlapping gene name labels with halo (depends on ggrepel)
    ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene, color = Significance),
      size = 4.5,
      fontface = "bold.italic",    # Bold italic
      bg.color = "white",          # White halo edge
      bg.r = 0.15,
      box.padding = 0.8,
      point.padding = 0.3,
      segment.color = "grey30",
      segment.size = 0.6,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.015, "npc"), type = "closed"), # Leader line arrow
      min.segment.length = 0,      # Force display of leader lines
      show.legend = FALSE,
      max.overlaps = Inf
    ) +
    ggplot2::scale_color_manual(values = c("UP" = "#C0392B", "DOWN" = "#2980B9")) +

    # Academic theme beautification
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(title = plot_title,
                  subtitle = sub_title,
                  x = expression(bold("Log"[2]*" Fold Change")),
                  y = expression(bold("-Log"[10]*" (Adjusted P-value)"))) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.5),

      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, color = "#2C3E50", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "#7F8C8D", face = "italic", margin = ggplot2::margin(b = 15)),

      axis.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.8, color = "black"),

      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1.5, "cm")
    ) +
    # Override legend to make them solid circles of uniform size
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1, shape = 21, color = "white")))

  return(p)
}

#' Plot Trajectory and Pseudotime Inference
#'
#' @description Generates a dual-panel plot for pseudotime analysis.
#' The left panel shows the continuous pseudotime gradient, while the right panel shows the categorical clusters.
#' Overlays the inferred trajectory path and the root node (start cluster) if available.
#'
#' @param obj A \code{RNA} object containing pseudotime results in \code{reductions$pseudotime}.
#'
#' @return A \code{patchwork} object containing the combined trajectory plots.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_path scale_color_viridis_c scale_color_manual scale_fill_manual labs theme_minimal theme element_text element_blank element_rect arrow unit guide_legend guides
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @importFrom grDevices colorRampPalette
#' @importFrom tools toTitleCase
#'
#' @examples
#' # 1. Create mock RNA object with trajectory data
#' set.seed(123)
#' mock_umap <- matrix(rnorm(40), ncol = 2)
#' rownames(mock_umap) <- paste0("Cell_", 1:20)
#' colnames(mock_umap) <- c("UMAP_1", "UMAP_2")
#'
#' mock_meta <- data.frame(Auto_Cluster = rep(c("Cluster_1", "Cluster_2"), each = 10))
#' rownames(mock_meta) <- rownames(mock_umap)
#'
#' # Mock pseudotime structure
#' mock_pseudotime <- list(
#'   pseudotime = runif(20, 0, 1),
#'   dr_method = "umap",
#'   group_col = "Auto_Cluster",
#'   algorithm = "cluster",
#'   start_clus = "Cluster_1",
#'   curve_coords = matrix(rnorm(40), ncol = 2)
#' )
#'
#' mock_obj <- list(
#'   reductions = list(umap = mock_umap, pseudotime = mock_pseudotime),
#'   meta.data = mock_meta
#' )
#' class(mock_obj) <- "RNA"
#'
#' # 2. Plot trajectory (requires patchwork)
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   p_pseudo <- PlotPseudo_RNA(mock_obj)
#'   p_pseudo
#' }
#' @export
PlotPseudo_RNA <- function(obj) {

  # --- 1. Basic data validation ---
  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  if (is.null(obj$reductions$pseudotime)) {
    stop(" Pseudotime data not found. Please run the pseudotime inference function first.")
  }

  pt_res <- obj$reductions$pseudotime
  dr_method <- pt_res$dr_method
  group_col <- pt_res$group_col
  algo <- pt_res$algorithm
  start_clus <- pt_res$start_clus

  # Extract dimensionality reduction coordinates
  coords_raw <- obj$reductions[[dr_method]]
  if (is.null(coords_raw)) {
    stop(sprintf("Error: Cannot find the specified dimensionality reduction coordinates '%s'.", dr_method))
  }

  # --- 2. Data cleaning and merging ---
  coords <- as.data.frame(coords_raw[, 1:2])
  colnames(coords) <- c("Dim1", "Dim2")
  coords$Pseudotime <- pt_res$pseudotime

  # Prioritize using filter_meta.data
  meta <- if (!is.null(obj$filter_meta.data)) obj$filter_meta.data else obj$meta.data

  # Defensively match cell names
  common_cells <- intersect(rownames(coords), rownames(meta))
  if (length(common_cells) == 0) stop(" Coordinates and Metadata cell names do not match!")

  coords <- coords[common_cells, , drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]

  # Extract grouping information
  if (!is.null(group_col) && group_col %in% colnames(meta)) {
    groups <- as.character(meta[[group_col]])
    groups[is.na(groups) | groups == ""] <- "Unknown"
    coords$Cluster <- as.factor(groups)
  } else {
    coords$Cluster <- as.factor("All Cells")
    group_col <- "All Cells"
  }

  # --- 3. Dynamic palette and theme settings ---
  sci_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                   "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                   "#7E6148FF", "#B09C85FF", "#FF7F00", "#6A3D9A")

  n_clusters <- length(levels(coords$Cluster))
  dynamic_colors <- if (n_clusters <= length(sci_palette)) {
    sci_palette[1:n_clusters]
  } else {
    grDevices::colorRampPalette(sci_palette)(n_clusters)
  }

  my_theme <- ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16, color = "#2c3e50"),
      axis.title = ggplot2::element_text(face = "bold", size = 13, color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.2),
      legend.position = "right"
    )

  # Calculate the centroid of the starting cluster (used to place a prominent start marker)
  root_cells <- coords[coords$Cluster == start_clus, ]
  if (nrow(root_cells) > 0) {
    root_centroid <- data.frame(Dim1 = mean(root_cells$Dim1), Dim2 = mean(root_cells$Dim2))
  } else {
    root_centroid <- NULL
  }

  # --- 4. Build base scatter plots ---
  # [Left plot]: Continuous variable Pseudotime
  p_left <- ggplot2::ggplot(coords, ggplot2::aes(x = Dim1, y = Dim2)) +
    ggplot2::geom_point(ggplot2::aes(color = Pseudotime), size = 2.5, alpha = 0.8, shape = 16) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
    ggplot2::labs(title = paste0("Pseudotime Gradient (", tools::toTitleCase(algo), ")"),
                  x = paste0(toupper(dr_method), " 1"),
                  y = paste0(toupper(dr_method), " 2")) +
    my_theme

  # [Right plot]: Discrete variable Cluster
  p_right <- ggplot2::ggplot(coords, ggplot2::aes(x = Dim1, y = Dim2)) +
    ggplot2::geom_point(ggplot2::aes(color = Cluster), size = 2.5, alpha = 0.8, shape = 16) +
    ggplot2::scale_fill_manual(values = dynamic_colors) +
    ggplot2::scale_color_manual(values = dynamic_colors) +
    ggplot2::labs(title = paste0("Cluster (", group_col, ")"),
                  x = paste0(toupper(dr_method), " 1"),
                  y = "") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))) +
    my_theme

  # --- 5. Overlay trajectory lines and starting point ---
  # Overlay black main backbone curve
  if (algo == "cluster" && !is.null(pt_res$curve_coords)) {
    curve_coords <- as.data.frame(pt_res$curve_coords)
    colnames(curve_coords) <- c("Dim1", "Dim2")
    curve_coords$pt <- pt_res$pseudotime
    curve_coords <- curve_coords[order(curve_coords$pt), ] # Ensure lines are ordered by pseudotime

    p_left <- p_left + ggplot2::geom_path(data = curve_coords, ggplot2::aes(x = Dim1, y = Dim2),
                                          color = "black", linewidth = 1.5,
                                          arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.18, "inches"), ends = "last"))
    p_right <- p_right + ggplot2::geom_path(data = curve_coords, ggplot2::aes(x = Dim1, y = Dim2),
                                            color = "grey20", linewidth = 1, linetype = "dashed", alpha = 0.6)
  }

  # Uniformly place starting point marker (large yellow circle with black border)
  if (!is.null(root_centroid)) {
    p_left <- p_left + ggplot2::geom_point(data = root_centroid, ggplot2::aes(x = Dim1, y = Dim2),
                                           shape = 21, fill = "#f1c40f", color = "black", size = 6, stroke = 1.5)
    p_right <- p_right + ggplot2::geom_point(data = root_centroid, ggplot2::aes(x = Dim1, y = Dim2),
                                             shape = 21, fill = "#f1c40f", color = "black", size = 6, stroke = 1.5)
  }

  # --- 6. Combine plots using Patchwork ---
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("The 'patchwork' package is required to combine plots. Please install it using install.packages('patchwork').")
    return(list(Pseudotime_Plot = p_left, Cluster_Plot = p_right))
  }

  combined_plot <- patchwork::wrap_plots(p_left, p_right, ncol = 2) +
    patchwork::plot_annotation(
      title = paste("Trajectory Inference via", tools::toTitleCase(algo), "Algorithm"),
      subtitle = paste("Rooted at:", start_clus),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
                             plot.subtitle = ggplot2::element_text(size = 15, color = "#7f8c8d", hjust = 0.5))
    )

  return(combined_plot)
}

#' Plot Epigenetic Landscape (Ridge Plot) from Level Matrix
#'
#' @description Calculates the global mean across features for each sample from a
#' provided '.level' matrix and generates a high-quality ridge plot showing the
#' distribution across specified categorical groups.
#'
#' @param mat A numeric \code{matrix} or \code{data.frame} of extraction levels (Rows = Features/Regions, Columns = Samples).
#' @param meta A \code{data.frame} containing metadata/sample information.
#' @param sample_col Character. The column name in \code{meta} containing sample IDs matching the colnames of \code{mat}.
#' @param group_col Character. The column name in \code{meta} representing the groups/conditions for the Y-axis.
#' @param title Character. The main title of the plot.
#' @param subtitle Character. The subtitle of the plot.
#' @param xlab Character. The X-axis label.
#' @param ylab Character. The Y-axis label.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes scale_fill_manual labs theme_bw theme element_text element_line element_blank element_rect margin
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # 1. Create a mock epigenetic level matrix (Features x Samples)
#' set.seed(42)
#' mock_mat <- matrix(runif(150, min = 0, max = 1), nrow = 10, ncol = 15)
#' colnames(mock_mat) <- paste0("Sample_", 1:15)
#' rownames(mock_mat) <- paste0("Region_", 1:10)
#'
#' # 2. Create corresponding metadata
#' mock_meta <- data.frame(
#'   SampleID = paste0("Sample_", 1:15),
#'   Condition = rep(c("Tumor", "Normal", "Metastasis"), each = 5)
#' )
#'
#' # 3. Generate Ridge Plot (Requires ggridges)
#' if (requireNamespace("ggridges", quietly = TRUE)) {
#'   p_ridge <- PlotLandscape_Epi(
#'     mat = mock_mat, meta = mock_meta,
#'     sample_col = "SampleID", group_col = "Condition"
#'   )
#'   p_ridge
#' }
#' @export
PlotLandscape_Epi <- function(mat,
                              meta,
                              sample_col,
                              group_col,
                              title = "Global Epigenetic Landscape",
                              subtitle = "Distribution of Mean Epigenetic Levels per Single Cell",
                              xlab = "Global Mean Level",
                              ylab = "Group / Condition") {

  # --- 1. Basic validation and dependency checks (Enhanced defense mechanism) ---
  if (!requireNamespace("ggridges", quietly = TRUE)) stop("The 'ggridges' package is required")

  # If the input is a data.frame, check for non-numeric columns first
  if (is.data.frame(mat)) {
    # === Change: Replace `sapply` with `vapply` ===
    # Explicitly specify that the expected return value is a Boolean (TRUE/FALSE)
    is_num <- vapply(mat, is.numeric, FUN.VALUE = logical(1))

    if (!all(is_num)) {
      # Find non-numeric columns
      bad_cols <- names(mat)[!is_num]
      stop(sprintf(
        "Data format error: Your matrix contains non-numeric (text) columns, e.g., column: '%s'.\nPlease set this column as rownames and remove it (e.g., mat <- mat[, -1]) to ensure all inputs are strictly numeric!",
        bad_cols[1]
      ))
    }
    mat <- as.matrix(mat)
  }

  if (!is.matrix(mat) || !is.numeric(mat)) {
    # If still erroring, attempt forced coercion (for numeric matrices accidentally recognized as characters)
    storage.mode(mat) <- "numeric"
    if (!is.numeric(mat)) stop("Input 'mat' cannot be converted to a numeric matrix!")
  }

  if (!is.data.frame(meta)) stop("Input 'meta' must be a data.frame!")
  if (!(sample_col %in% colnames(meta))) stop(sprintf("Sample ID column '%s' not found in meta.", sample_col))
  if (!(group_col %in% colnames(meta))) stop(sprintf("Grouping column '%s' not found in meta.", group_col))

  # --- 2. Sample alignment (Strict and safe matching logic) ---
  mat_samples <- colnames(mat)
  meta_samples <- as.character(meta[[sample_col]])

  common_samples <- intersect(mat_samples, meta_samples)
  if (length(common_samples) < 3) {
    stop(sprintf("Matching error: The intersection of matrix column names and the '%s' column in Metadata is less than 3. Please check if sample names match!", sample_col))
  }

  mat_sub <- mat[, common_samples, drop = FALSE]
  meta_sub <- meta[match(common_samples, meta[[sample_col]]), , drop = FALSE]

  # --- 3. Core calculation: Calculate global mean ---
  cell_means <- colMeans(mat_sub, na.rm = TRUE)

  ridge_df <- data.frame(
    Global_Mean = cell_means,
    Group = as.factor(as.character(meta_sub[[group_col]]))
  )

  # --- 4. Color engine and theme configuration ---
  n_groups <- length(levels(ridge_df$Group))
  palette_1 <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                 "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                 "#7E6148FF", "#B09C85FF", "#FF7F00", "#6A3D9A")

  get_dynamic_colors <- function(n, base_palette) {
    if (n <= length(base_palette)) return(base_palette[1:n])
    return(grDevices::colorRampPalette(base_palette)(n))
  }
  dynamic_colors <- get_dynamic_colors(n_groups, palette_1)

  # --- 5. Render ridge plot ---
  p_ridge <- ggplot2::ggplot(ridge_df, ggplot2::aes(x = Global_Mean, y = Group, fill = Group)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01,
                                  color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = dynamic_colors) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, margin = ggplot2::margin(b = 10)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13, color = "grey30", margin = ggplot2::margin(b = 15)),
      axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text = ggplot2::element_text(color = "black", size = 13),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8)
    )

  return(p_ridge)
}

#' Plot Dimensionality Reduction Results
#'
#' @description Generates an elegant, publication-ready scatter plot from the output of `Reduce_epiDims`.
#' Features dynamic coloring, automatic confidence ellipses (if sample size allows), and clean themes.
#'
#' @param dr_data Data.frame. The output from `Reduce_epiDims` (must contain Dim1, Dim2, and the grouping column).
#' @param group_col Character. The name of the column in `dr_data` used for grouping/coloring.
#' @param method Character. The dimension reduction method used (e.g., "PCA", "UMAP"). Used for axis labels and title.
#' @param title Character. Optional custom title. If NULL, defaults to the `method` name.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline stat_ellipse geom_point labs scale_fill_manual scale_color_manual theme_bw theme element_text element_blank element_rect element_line unit margin
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # 1. Create mock dimensionality reduction output data
#' set.seed(123)
#' mock_dr <- data.frame(
#'   Dim1 = rnorm(15),
#'   Dim2 = rnorm(15),
#'   Group = rep(c("Tumor", "Normal", "Control"), each = 5)
#' )
#'
#' # 2. Plot PCA with confidence ellipses
#' p_pca <- PlotDimRed_Epi(
#'   dr_data = mock_dr,
#'   group_col = "Group",
#'   method = "PCA"
#' )
#' p_pca
#' @export
PlotDimRed_Epi <- function(dr_data,
                           group_col,
                           method = "PCA",
                           title = NULL) {

  # --- 1. Data validation ---
  if (!is.data.frame(dr_data)) stop("Input 'dr_data' must be a data.frame (usually the output of Reduce_epiDims)!")
  if (!all(c("Dim1", "Dim2") %in% colnames(dr_data))) stop("'dr_data' must contain 'Dim1' and 'Dim2' columns!")
  if (!(group_col %in% colnames(dr_data))) stop(sprintf("Grouping column '%s' not found in data!", group_col))

  # --- 2. Format plotting data ---
  # Convert to factor uniformly to ensure ggplot maps discrete colors
  dr_data$PlotGroup <- as.factor(dr_data[[group_col]])

  # --- 3. Label and title preparation ---
  # Note: Since Reduce_epiDims only returns coordinates, if it is PCA, display PCA 1 / PCA 2 uniformly here
  x_lab <- paste0(toupper(method), " 1")
  y_lab <- paste0(toupper(method), " 2")
  if (is.null(title)) title <- paste(toupper(method), "Plot")

  # If it is NMF, fine-tune the axis names
  if (toupper(method) == "NMF") {
    x_lab <- "Basis 1"
    y_lab <- "Basis 2"
  }

  # --- 4. Dynamic color scheme ---
  sci_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                   "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                   "#7E6148FF", "#B09C85FF", "#FF7F00", "#6A3D9A")

  n_groups <- length(levels(dr_data$PlotGroup))
  if (n_groups > length(sci_palette)) {
    dynamic_colors <- grDevices::colorRampPalette(sci_palette)(n_groups)
  } else {
    dynamic_colors <- sci_palette[1:n_groups]
  }

  # --- 5. Smart determination of whether to draw confidence ellipses ---
  # Draw ellipse only when sample size for all groups >= 3 and number of groups > 1 (otherwise stat_ellipse throws an error)
  group_counts <- table(dr_data$PlotGroup)
  can_draw_ellipse <- all(group_counts >= 3) && length(group_counts) > 1

  # --- 6. Build ggplot layers ---
  p <- ggplot2::ggplot(dr_data, ggplot2::aes(x = Dim1, y = Dim2, fill = PlotGroup, color = PlotGroup)) +
    # Bottom center dashed crosshairs
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5)

  # Confidence ellipses
  if (can_draw_ellipse) {
    p <- p +
      ggplot2::stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.1, show.legend = FALSE, color = NA) +
      ggplot2::stat_ellipse(level = 0.95, geom = "path", linewidth = 0.8, linetype = 2, show.legend = FALSE)
  }

  # Top layer scatter points with white outlines + beautified theme
  p <- p +
    ggplot2::geom_point(size = 4.5, alpha = 0.9, shape = 21, stroke = 0.6, color = "white") +
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    ggplot2::scale_fill_manual(values = dynamic_colors) +
    ggplot2::scale_color_manual(values = dynamic_colors) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, margin = ggplot2::margin(b = 15)),
      axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 10)),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      legend.position = "right",
      legend.title = ggplot2::element_blank(), # Hide legend title to make the plot cleaner
      legend.text = ggplot2::element_text(size = 13),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = ggplot2::unit(0.2, "cm")
    )

  return(p)
}

#' Plot Volcano Plot for Differential Epigenetic Analysis
#'
#' @description Generates a publication-quality volcano plot to visualize differential
#' methylation/epigenetic regions (DMRs). Supports dynamic effect size metrics
#' (e.g., logFC, Difference, Z-Scores).
#'
#' @param df A \code{data.frame} containing differential analysis results.
#' @param metric_col Character. The column name for the effect size (X-axis). Default is "Diff".
#' @param p_col Character. The column name for the P-values (Y-axis). Default is "P.Value".
#' @param feature_col Character. The column name for feature names (labels). Default is "chrdata".
#' @param th_effect Numeric. The absolute threshold for the effect size (e.g., logFC > 1). Default is 1.0.
#' @param th_p Numeric. The threshold for statistical significance (e.g., P-value < 0.05). Default is 0.05.
#' @param top_n Integer. Number of top significant features to label. Default is 10.
#' @param title Character. The main title of the plot. Default is "Volcano Plot".
#' @param custom_xlab Character or Expression. Custom label for the X-axis. If NULL, auto-generated based on \code{metric_col}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_point scale_fill_manual theme_bw labs theme element_text element_blank element_rect margin
#'
#' @examples
#' # 1. Create mock Differential Epigenetic Analysis (DEA) results
#' set.seed(42)
#' mock_epi_dea <- data.frame(
#'   chrdata = paste0("chr1:", 1000:(1000+99)),
#'   Diff = rnorm(100, mean = 0, sd = 0.8),
#'   P.Value = runif(100, min = 0, max = 0.5)
#' )
#' # Force some significant points
#' mock_epi_dea$Diff[1:3] <- c(1.5, 2.0, -1.8)
#' mock_epi_dea$P.Value[1:3] <- c(0.001, 0.0005, 0.002)
#'
#' # 2. Plot Epigenetic Volcano Plot
#' p_vol <- PlotVolcano_Epi(
#'   df = mock_epi_dea,
#'   metric_col = "Diff", p_col = "P.Value", feature_col = "chrdata",
#'   th_effect = 1.0, th_p = 0.05, top_n = 3
#' )
#' p_vol
#' @export
PlotVolcano_Epi <- function(df,
                            metric_col = "Diff",
                            p_col = "P.Value",
                            feature_col = "chrdata",
                            th_effect = 1.0,
                            th_p = 0.05,
                            top_n = 10,
                            title = "Volcano Plot",
                            custom_xlab = NULL) {

  # --- 1. Basic validation and dependency checks ---
  if (!is.data.frame(df)) stop("Input 'df' must be a data.frame.")
  required_cols <- c(metric_col, p_col, feature_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns not found in data.frame: %s", paste(missing_cols, collapse = ", ")))
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("The 'ggrepel' package is required to draw gene labels. Please run: install.packages('ggrepel')")
  }

  # --- 2. Dynamically mark significance (up/down regulated groups) ---
  df$Significance <- "Not Sig"
  effect_vals <- df[[metric_col]]
  p_vals <- df[[p_col]]

  df$Significance[effect_vals > th_effect & p_vals < th_p] <- "Up-regulated"
  df$Significance[effect_vals < -th_effect & p_vals < th_p] <- "Down-regulated"

  # Convert to factor, fix level order to ensure stable color mapping
  df$Significance <- factor(df$Significance, levels = c("Up-regulated", "Down-regulated", "Not Sig"))

  # --- 3. Sort and extract Top features for labeling ---
  # Ascending by P-value, descending by absolute Effect Size
  df <- df[order(df[[p_col]], -abs(df[[metric_col]])), ]
  top_features <- head(df[df$Significance != "Not Sig", ], top_n)

  # --- 4. Dynamically generate X-axis label ---
  if (is.null(custom_xlab)) {
    x_label <- switch(metric_col,
                      "logFC"  = expression(bold(Log[2]~"Fold Change")),
                      "Diff"   = expression(bold("Methylation Difference ("*Delta*")")),
                      "Scores" = expression(bold("Joint Z-Score")),
                      metric_col) # If not these three, use column name by default
  } else {
    x_label <- custom_xlab
  }

  # --- 5. Render publication-quality elegant volcano plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric_col]], y = -log10(.data[[p_col]]), fill = Significance)) +
    # Auxiliary lines
    ggplot2::geom_hline(yintercept = -log10(th_p), linetype = "dashed", color = "black", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = c(-th_effect, th_effect), linetype = "dashed", color = "black", linewidth = 0.6) +

    # Scatter plot layer (hollow dots with white outlines)
    ggplot2::geom_point(alpha = 0.85, size = 2.5, shape = 21, stroke = 0.3, color = "white") +

    # Custom color mapping (classic NPG red-blue palette, drop=FALSE ensures no error when groups are missing)
    ggplot2::scale_fill_manual(
      values = c("Up-regulated" = "#DC0000FF", "Down-regulated" = "#4DBBD5FF", "Not Sig" = "#CCCCCC"),
      drop = FALSE
    ) +

    # Axes and title
    ggplot2::labs(
      title = title,
      x = x_label,
      y = expression(bold(-Log[10]~italic("P")*"-value"))
    ) +

    # Advanced theme beautification
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18, margin = ggplot2::margin(b = 15)),
      legend.position = "top",
      legend.title = ggplot2::element_blank(), # Hide legend title
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.text = ggplot2::element_text(color = "black")
    )

  # --- 6. Add smart non-overlapping text labels ---
  if (nrow(top_features) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_features,
      ggplot2::aes(label = .data[[feature_col]]),
      size = 4.5,
      color = "black",
      box.padding = 0.6,
      max.overlaps = Inf,
      fontface = "italic",
      bg.color = "white",
      bg.r = 0.15 # Text white edge to prevent illegibility when overlapping with background dots
    )
  }

  return(p)
}

#' Plot Multi-Omics Correlation Matrix
#'
#' @description Generates a matrix plot showing the distribution (density) of each omics layer on the diagonal,
#' pairwise scatter plots on the lower triangle, and Spearman correlation coefficients on the upper triangle.
#' Automatically adapts to dual-omics (2x2) or tri-omics (3x3) data based on the input dataframe.
#'
#' @param df Data.frame. The integrated multi-omics dataframe output by \code{IntegrateMultiOmics}.
#' @param title Character. Global title for the plot.
#' @param cor_method Character. Correlation method to use ("spearman", "pearson", or "kendall").
#' @param point_color Character. Color of the scatter plot points.
#' @param density_fill Character. Fill color for the density plots.
#' @param alpha Numeric. Transparency for points and density fills (0 to 1).
#'
#' @return A \code{ggplot} object representing the correlation plots.
#'
#' @import ggplot2
#' @importFrom stats cor
#' @importFrom patchwork plot_annotation plot_layout
#'
#' @examples
#' # 1. Create mock multi-omics integrated dataframe
#' set.seed(123)
#' mock_omics <- data.frame(
#'   RNA_Exp = rnorm(50, mean = 10, sd = 2),
#'   CpG_level = runif(50, 0, 1),
#'   GpC_level = runif(50, 0, 1)
#' )
#'
#' # 2. Generate Omics Correlation Matrix Plot (requires patchwork)
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   p_matrix <- PlotOmicsMatrix(mock_omics, cor_method = "spearman")
#'   p_matrix
#' }
#' @export
PlotOmicsMatrix <- function(
    df,
    title = "Omics Correlation Matrix",
    cor_method = "spearman",
    point_color = "#3C5488",
    density_fill = "#B09C85",
    alpha = 0.5
) {
  # 1. Dynamically detect included omics columns and rename them for plotting
  cols_sel <- c()
  if ("RNA_Exp" %in% names(df)) cols_sel <- c(cols_sel, "RNA" = "RNA_Exp")
  if ("CpG_level" %in% names(df)) cols_sel <- c(cols_sel, "CpG" = "CpG_level")
  if ("GpC_level" %in% names(df)) cols_sel <- c(cols_sel, "GpC" = "GpC_level")

  if (length(cols_sel) < 2) {
    stop("Plotting failed: The input data.frame must contain at least two omics data columns (RNA_Exp, CpG_level, GpC_level)!")
  }

  # Extract needed columns and rename them to simpler names (RNA, CpG, GpC)
  plot_df <- df[, cols_sel, drop = FALSE]
  colnames(plot_df) <- names(cols_sel)

  # 2. Internal closure function: Draw diagonal density plots (upgraded to modern .data[[]] syntax)
  plot_density <- function(data, var, title_text) {
    ggplot2::ggplot(data, aes(x = .data[[var]])) +
      geom_density(fill = density_fill, alpha = alpha, color = "black") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      annotate("text", x = Inf, y = Inf, label = title_text,
               hjust = 1.2, vjust = 1.5, size = 5, fontface = "bold")
  }

  # 3. Internal closure function: Draw lower-left scatter plots
  plot_scatter_matrix <- function(data, x_var, y_var) {
    ggplot2::ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(alpha = alpha, size = 1, color = point_color) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  }

  # 4. Internal closure function: Draw upper-right correlation coefficient text
  plot_cor_text <- function(data, x_var, y_var) {
    r_val <- cor(data[[x_var]], data[[y_var]], method = cor_method, use = "complete.obs")
    lbl <- sprintf("R = %.2f", r_val)
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = lbl, size = 6, fontface = "bold", color = "black") +
      theme_void() +
      theme(panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1))
  }

  # 5. Dynamically build Patchwork layout based on the number of omics
  if (ncol(plot_df) == 3) {
    # Tri-omics (3x3 matrix)
    p11 <- plot_density(plot_df, "RNA", "RNA Exp")
    p12 <- plot_cor_text(plot_df, "RNA", "CpG")
    p13 <- plot_cor_text(plot_df, "RNA", "GpC")

    p21 <- plot_scatter_matrix(plot_df, "RNA", "CpG") + labs(x = "", y = "CpG")
    p22 <- plot_density(plot_df, "CpG", "CpG Level")
    p23 <- plot_cor_text(plot_df, "CpG", "GpC")

    p31 <- plot_scatter_matrix(plot_df, "RNA", "GpC") + labs(x = "RNA", y = "GpC")
    p32 <- plot_scatter_matrix(plot_df, "CpG", "GpC") + labs(x = "CpG", y = "")
    p33 <- plot_density(plot_df, "GpC", "GpC Level")

    # Assemble 3x3 layout
    p_mat <- (p11 | p12 | p13) /
      (p21 | p22 | p23) /
      (p31 | p32 | p33)

  } else if (ncol(plot_df) == 2) {
    # Dual-omics (2x2 matrix)
    v1 <- names(plot_df)[1]
    v2 <- names(plot_df)[2]

    p11 <- plot_density(plot_df, v1, paste(v1, "Level"))
    p12 <- plot_cor_text(plot_df, v1, v2)

    p21 <- plot_scatter_matrix(plot_df, v1, v2) + labs(x = v1, y = v2)
    p22 <- plot_density(plot_df, v2, paste(v2, "Level"))

    # Assemble 2x2 layout
    p_mat <- (p11 | p12) /
      (p21 | p22)
  }

  # 6. Add global main title
  p_mat <- p_mat + patchwork::plot_annotation(
    title = title,
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

  return(p_mat)
}

#' Plot Multi-Omics Pairwise Scatter Plots
#'
#' @description Generates pairwise scatter plots for integrated multi-omics data.
#' Automatically adds linear regression trend lines (LM) and calculates correlation
#' coefficients (Spearman by default) using \code{ggpubr::stat_cor}.
#' Dynamically generates 1 to 3 plots depending on the available omics layers in the dataframe.
#'
#' @param df Data.frame. The integrated multi-omics dataframe output by \code{IntegrateMultiOmics}.
#' @param cor_method Character. Correlation method: "spearman" (default) or "pearson".
#' @param point_size Numeric. Size of the scatter points.
#' @param point_alpha Numeric. Transparency of the scatter points (0 to 1).
#' @param color_cpg_rna Character. Color for CpG vs RNA plot.
#' @param color_gpc_rna Character. Color for GpC vs RNA plot.
#' @param color_cpg_gpc Character. Color for CpG vs GpC plot.
#'
#' @return A \code{ggplot} object representing the scatter plots.
#'
#' @import ggplot2
#' @importFrom ggpubr stat_cor
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' # 1. Create mock multi-omics integrated dataframe
#' set.seed(42)
#' mock_omics <- data.frame(
#'   RNA_Exp = rnorm(30, mean = 5, sd = 1),
#'   CpG_level = runif(30, 0, 1)
#'   # Omitted GpC_level to test the 2-omics dynamic layout
#' )
#'
#' # Add artificial correlation
#' mock_omics$RNA_Exp <- mock_omics$RNA_Exp - (mock_omics$CpG_level * 2)
#'
#' # 2. Generate pairwise scatter plots (requires ggpubr and patchwork)
#' if (requireNamespace("ggpubr", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   p_scatter <- PlotOmicsScatter(mock_omics)
#'   p_scatter
#' }
#' @export
PlotOmicsScatter <- function(
    df,
    cor_method = "spearman",
    point_size = 1.5,
    point_alpha = 0.6,
    color_cpg_rna = "#4DBBD5",
    color_gpc_rna = "#00A087",
    color_cpg_gpc = "#E64B35"
) {

  # Define globally unified theme
  my_theme <- theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank())

  plot_list <- list()

  # 1. Plot CpG vs RNA (if both columns exist)
  if ("CpG_level" %in% names(df) && "RNA_Exp" %in% names(df)) {
    p1 <- ggplot(df, aes(x = CpG_level, y = RNA_Exp)) +
      geom_point(alpha = point_alpha, size = point_size, color = color_cpg_rna) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      ggpubr::stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 5) +
      labs(x = "CpG Methylation", y = "RNA Expression") +
      my_theme

    plot_list <- append(plot_list, list(p1))
  }

  # 2. Plot GpC vs RNA (if both columns exist)
  if ("GpC_level" %in% names(df) && "RNA_Exp" %in% names(df)) {
    p2 <- ggplot(df, aes(x = GpC_level, y = RNA_Exp)) +
      geom_point(alpha = point_alpha, size = point_size, color = color_gpc_rna) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      ggpubr::stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 5) +
      labs(x = "GpC Accessibility", y = "RNA Expression") +
      my_theme

    plot_list <- append(plot_list, list(p2))
  }

  # 3. Plot CpG vs GpC (if both columns exist)
  if ("CpG_level" %in% names(df) && "GpC_level" %in% names(df)) {
    p3 <- ggplot(df, aes(x = CpG_level, y = GpC_level)) +
      geom_point(alpha = point_alpha, size = point_size, color = color_cpg_gpc) +
      geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
      ggpubr::stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 5) +
      labs(x = "CpG Methylation", y = "GpC Accessibility") +
      my_theme

    plot_list <- append(plot_list, list(p3))
  }

  # 4. Safety check: If no plots were generated (column names are incorrect)
  if (length(plot_list) == 0) {
    stop("Plotting failed: No matching omics columns found in the data.frame (must include at least two of RNA_Exp, CpG_level, GpC_level).")
  }

  # 5. Use patchwork for dynamic layout (arranged in a single row)
  final_p <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))

  return(final_p)
}
