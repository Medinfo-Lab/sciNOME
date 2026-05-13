#' Create a Lightweight RNA Object
#'
#' @description Internal core function to initialize a RNA object from a raw counts matrix.
#' Automatically calculates basic QC metrics including percent_mt.
#'
#' @param counts A sparse matrix or data.frame containing raw expression counts.
#' @param project_name Character string specifying the project name. Default is "RNA_Project".
#' @param mt_pattern Character string. Regex pattern to identify mitochondrial genes. Default is "^[Mm][Tt]-".
#'
#' @return A \code{RNA} object containing initialized assays and basic QC metrics.
#'
#' @importFrom Matrix colSums
#' @importFrom methods as
#' @export
CreateRNAObject <- function(counts, project_name = "RNA_Project", mt_pattern = "^[Mm][Tt]-") {

  # Force conversion to a sparse matrix
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(as.matrix(counts), "CsparseMatrix")
  }

  cell_names <- colnames(counts)

  # Basic QC calculation
  nCount_RNA <- as.numeric(Matrix::colSums(counts))
  nFeature_RNA <- as.numeric(Matrix::colSums(counts > 0))

  # New: Extremely fast calculation of mitochondrial gene proportion (percent_mt)
  mt_genes <- grep(mt_pattern, rownames(counts), value = TRUE)

  if (length(mt_genes) > 0) {
    # Extract the sub-matrix of mitochondrial genes and sum
    mt_counts <- as.numeric(Matrix::colSums(counts[mt_genes, , drop = FALSE]))
    percent_mt <- (mt_counts / nCount_RNA) * 100

    # Prevent 0/0=NaN caused by nCount_RNA=0 for some samples
    percent_mt[is.na(percent_mt)] <- 0
  } else {
    # If no mitochondrial genes are found (e.g., extracted matrix), directly assign 0
    percent_mt <- rep(0, length(cell_names))
    message(sprintf("Note: No mitochondrial genes matching the '%s' pattern were detected in the matrix. percent_mt is automatically set to 0.", mt_pattern))
  }

  # Build meta.data
  meta_data <- data.frame(
    row.names = cell_names,
    nCount_RNA = nCount_RNA,
    nFeature_RNA = nFeature_RNA,
    percent_mt = percent_mt  # Add percent_mt column
  )

  # Assemble the returned object
  object <- list(
    project = project_name,
    assays = list(
      RNA = list(
        counts = counts,
        filter_counts = NULL,
        data = NULL,
        scale.data = NULL)
    ),
    meta.data = meta_data,
    filter_meta.data = NULL,
    reductions = list(
      pca = NULL,
      umap = NULL,
      tsne = NULL
    ),
    var.genes = NULL
  )

  class(object) <- "RNA"
  return(object)
}

#' Add Metadata to a RNA Object
#'
#' @description Core engine for injecting metadata into the \code{meta.data} slot of a RNA object.
#' Supports both strict ID matching and ordered vector injection.
#'
#' @param object A \code{RNA} object.
#' @param meta_data A data.frame containing metadata to add.
#' @param id_col Character. Column name in \code{meta_data} to match with cell names. If NULL or "MATCH_BY_ORDER", matches by row order.
#' @param group_cols Character vector of column names to add. If NULL, adds all columns except \code{id_col}.
#'
#' @return An updated \code{RNA} object with integrated metadata.
#'
#' @export
AddMetaData_RNA <- function(object, meta_data, id_col = NULL, group_cols = NULL) {
  current_meta <- object$meta.data
  cell_names <- rownames(current_meta)

  # If group_cols is not specified, add all columns in metadata by default
  if (is.null(group_cols)) {
    group_cols <- colnames(meta_data)
    if (!is.null(id_col) && id_col != "MATCH_BY_ORDER") {
      group_cols <- setdiff(group_cols, id_col) # Exclude the ID column itself
    }
  }

  # Case A: Vector mode / hard match by order
  if (is.null(id_col) || id_col == "MATCH_BY_ORDER") {
    for (col in group_cols) {
      vec <- meta_data[[col]]
      if (length(vec) != length(cell_names)) {
        stop(paste("The length of column", col, "(", length(vec), ") does not match the number of samples (", length(cell_names), ")."))
      }
      current_meta[[col]] <- vec
    }
  }
  # Case B: Dataframe precise ID matching
  else {
    clean_meta <- meta_data[!duplicated(meta_data[[id_col]]), ]
    rownames(clean_meta) <- as.character(clean_meta[[id_col]])
    common_cells <- intersect(cell_names, rownames(clean_meta))

    if (length(common_cells) == 0) {
      stop(sprintf("ID matching failed! Sample names in expression matrix: %s... \nID names in Metadata: %s...",
                   paste(head(cell_names, 3), collapse=", "),
                   paste(head(rownames(clean_meta), 3), collapse=", ")))
    }
    for (col in group_cols) {
      current_meta[[col]] <- NA # Initialize as NA
      current_meta[common_cells, col] <- clean_meta[common_cells, col]
    }
  }
  object$meta.data <- current_meta
  return(object)
}

#' Build a Complete RNA Object
#'
#' @description One-step function to generate a RNA object directly from an expression matrix and a metadata data.frame.
#' Automatically handles data cleaning, NA replacement, pre-filtering, and QC metrics calculation.
#'
#' @param expr_mat Expression matrix (data.frame or matrix), rows are genes, columns are samples.
#' @param meta_data Metadata data.frame for sample grouping.
#' @param meta_id_col Character. Column in \code{meta_data} matching matrix column names. Defaults to ordered matching if NULL.
#' @param meta_group_cols Character vector. Specific metadata columns to add (Defaults to NULL, adding all columns).
#' @param min_cells Integer. Pre-filtering: keep genes expressed in at least this many cells.
#' @param min_features Integer. Pre-filtering: keep cells expressing at least this many genes.
#' @param mt_pattern Character string. Regex pattern for mitochondrial genes. Default is "^[Mm][Tt]-".
#' @param project_name Character. Project name.
#'
#' @return A built and pre-filtered \code{RNA} object.
#'
#' @importFrom Matrix colSums rowSums
#' @export
Build_RNAObject <- function(expr_mat,
                            meta_data = NULL,
                            meta_id_col = NULL,
                            meta_group_cols = NULL,
                            min_cells = 0,
                            min_features = 0,
                            mt_pattern = "^[Mm][Tt]-",
                            project_name = "RNA_Project") {

  # --- 1. Smart data cleaning and matrix transformation ---
  message("Checking and cleaning expression matrix format...")
  if (is.data.frame(expr_mat)) {
    if (!is.numeric(expr_mat[[1]])) {
      message("Detected non-numeric first column (suspected gene names), automatically converting it to row names...")
      rownames(expr_mat) <- make.unique(as.character(expr_mat[[1]]))
      expr_mat <- expr_mat[, -1, drop = FALSE]
    }

    # === Change 1: Replace `sapply` with `vapply` ===
    # The `is.numeric` function should return a logical value (TRUE/FALSE)
    numeric_cols <- vapply(expr_mat, is.numeric, FUN.VALUE = logical(1))

    if (sum(numeric_cols) == 0) stop("Error: No numeric columns found in the expression matrix!")
    expr_mat <- as.matrix(expr_mat[, numeric_cols, drop = FALSE])
  } else {
    expr_mat <- matrix(as.numeric(expr_mat), nrow = nrow(expr_mat), dimnames = dimnames(expr_mat))
  }

  # --- 2. Check missing values and row/col names ---
  if (any(is.na(expr_mat))) {
    message("Warning: The matrix contains NA missing values, automatically replacing NA with 0...")
    expr_mat[is.na(expr_mat)] <- 0
  }
  if (is.null(colnames(expr_mat))) colnames(expr_mat) <- paste0("Sample_", 1:ncol(expr_mat))
  if (is.null(rownames(expr_mat))) rownames(expr_mat) <- paste0("Gene_", 1:nrow(expr_mat))
  if (any(is.na(colnames(expr_mat)))) {
    colnames(expr_mat)[is.na(colnames(expr_mat))] <- paste0("UnknownSample_", 1:sum(is.na(colnames(expr_mat))))
  }

  # --- 3. Pre-filter ---
  message("Performing filtering: min_cells = ", min_cells, ", min_features = ", min_features)
  n_features_per_cell <- Matrix::colSums(expr_mat > 0)
  keep_cells <- n_features_per_cell >= min_features
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]

  n_cells_per_gene <- Matrix::rowSums(expr_mat > 0)
  keep_genes <- n_cells_per_gene >= min_cells
  expr_mat <- expr_mat[keep_genes, , drop = FALSE]

  if(ncol(expr_mat) == 0 || nrow(expr_mat) == 0) {
    stop("Filtering is too strict! 0 samples or 0 genes remaining, please decrease min_cells / min_features.")
  }

  # --- 4. Create base object (pass mt_pattern here) ---
  message("Creating underlying object to calculate QC (including mitochondrial proportion)...")
  sci_obj <- CreateRNAObject(expr_mat, project_name = project_name, mt_pattern = mt_pattern)

  # --- 5. Integrate Metadata ---
  if (!is.null(meta_data)) {
    message("Merging Metadata grouping information...")
    if (!is.data.frame(meta_data)) meta_data <- as.data.frame(meta_data)
    sci_obj <- AddMetaData_RNA(object = sci_obj,
                               meta_data = meta_data,
                               id_col = meta_id_col,
                               group_cols = meta_group_cols)
  }

  message(sprintf("Build complete! Retained %d genes and %d samples.", nrow(expr_mat), ncol(expr_mat)))
  return(sci_obj)
}

#' Process QC, Filtering, and Normalization for RNA Object
#'
#' @description Performs mitochondrial proportion calculation, rigorous cell filtering,
#' sparse matrix normalization (LogNormalize, LogCPM, or TPM), and optional data scaling.
#'
#' @param obj A \code{RNA} object.
#' @param mt_pattern Character. Regex pattern to identify mitochondrial genes (e.g., "^MT-" for human, "^mt-" for mouse).
#' @param min_nCount Numeric. Minimum UMI count threshold.
#' @param max_nCount Numeric. Maximum UMI count threshold.
#' @param min_nFeature Numeric. Minimum number of detected genes.
#' @param max_nFeature Numeric. Maximum number of detected genes.
#' @param max_mt Numeric. Maximum allowed mitochondrial percentage (e.g., 15 means 15\%).
#' @param norm_method Character. Normalization method: "LogNormalize", "LogCPM", or "TPM".
#' @param do_scale Logical. Whether to center and scale the data (highly memory-intensive).
#'
#' @return A \code{RNA} object after QC, filtering, and normalization.
#'
#' @importFrom Matrix colSums t
#' @export
ProcessQC_RNA <- function(obj,
                          mt_pattern = "^MT-",
                          min_nCount = 500,
                          max_nCount = 50000,
                          min_nFeature = 200,
                          max_nFeature = 8000,
                          max_mt = 15,
                          norm_method = "LogNormalize",
                          do_scale = TRUE) {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # --- Stage 1: Calculate mitochondrial proportion ---
  message(sprintf("1. Calculating mitochondrial proportion (Pattern: %s)...", mt_pattern))
  mt_genes <- grep(mt_pattern, rownames(obj$assays$RNA$counts), value = TRUE)

  if(length(mt_genes) > 0) {
    mt_counts <- Matrix::colSums(obj$assays$RNA$counts[mt_genes, , drop = FALSE])
    obj$meta.data$percent_mt <- (mt_counts / obj$meta.data$nCount_RNA) * 100
  } else {
    obj$meta.data$percent_mt <- 0
    message("Warning: No mitochondrial genes matching '", mt_pattern, "' found. percent_mt is set to 0.")
  }

  # --- Stage 2: Strict QC filtering ---
  message(sprintf("2. Performing filtering: nCount(%d-%d), nFeature(%d-%d), MT(<%d%%)...",
                  min_nCount, max_nCount, min_nFeature, max_nFeature, max_mt))
  meta <- obj$meta.data

  # Vectorized condition filtering
  keep_cells <- rownames(meta)[
    meta$nCount_RNA >= min_nCount &
      meta$nCount_RNA <= max_nCount &
      meta$nFeature_RNA >= min_nFeature &
      meta$nFeature_RNA <= max_nFeature &
      meta$percent_mt <= max_mt
  ]

  if (length(keep_cells) == 0) {
    stop("Error: Filtering conditions are too strict! 0 cells remaining. Please relax the filtering thresholds.")
  }

  # Save filtered data
  obj$filter_meta.data <- meta[keep_cells, , drop = FALSE]
  obj$assays$RNA$filter_counts <- obj$assays$RNA$counts[, keep_cells, drop = FALSE]
  message(sprintf("   -> Filtering complete. Retained %d / %d samples.", length(keep_cells), nrow(meta)))

  # --- Stage 3: Extremely fast sparse matrix normalization ---
  message(sprintf("3. Performing normalization (Method: %s)...", norm_method))
  counts_sparse <- obj$assays$RNA$filter_counts
  col_sums <- Matrix::colSums(counts_sparse)
  norm_mat <- NULL

  if (norm_method == "LogNormalize") {
    # ln( (count / total_count) * 10000 + 1 )
    norm_mat <- Matrix::t(Matrix::t(counts_sparse) / col_sums) * 10000
    norm_mat@x <- log1p(norm_mat@x)
  } else if (norm_method == "LogCPM") {
    # log2( (count / total_count) * 1,000,000 + 1 )
    norm_mat <- Matrix::t(Matrix::t(counts_sparse) / col_sums) * 1e6
    norm_mat@x <- log2(norm_mat@x + 1)
  } else if (norm_method == "TPM") {
    # (count / total_count) * 1,000,000
    norm_mat <- Matrix::t(Matrix::t(counts_sparse) / col_sums) * 1e6
  } else {
    stop("Unsupported norm_method, please select LogNormalize, LogCPM, or TPM.")
  }

  obj$assays$RNA$data <- norm_mat

  # --- Stage 4: Data centering and scaling ---
  if (do_scale) {
    message("4. Performing data centering and scaling - Note: This step consumes a large amount of memory...")
    tryCatch({
      # Sparse to dense, scale by row (gene)
      scale_mat <- t(scale(t(as.matrix(norm_mat))))
      # Remove NaN calculated due to entirely 0 rows
      scale_mat[is.na(scale_mat)] <- 0
      obj$assays$RNA$scale.data <- scale_mat
      message("   -> Scaling complete!")
    }, error = function(err) {
      message("Warning: Out of memory (OOM), Scaling step skipped to prevent program crash.")
      obj$assays$RNA$scale.data <- NULL
    })
  } else {
    message("4. User chose to skip Scaling step.")
    obj$assays$RNA$scale.data <- NULL
  }

  # Force garbage collection to release memory
  gc()

  message("All QC and data processing complete!")
  return(obj)
}

#' Run Dimensionality Reduction (PCA, t-SNE, UMAP)
#'
#' @description Computes Highly Variable Genes (HVG) and performs PCA, t-SNE, or UMAP on the specified data layer.
#'
#' @param obj A \code{RNA} object.
#' @param method Character. Dimensionality reduction method: "PCA", "t-SNE", or "UMAP".
#' @param layer_name Character. Data layer to use: "data" (recommended), "scale.data", "counts", or "filter_counts".
#' @param n_hvg Integer. Number of highly variable genes to calculate (for PCA only).
#' @param pca_rank Integer. Number of principal components to compute.
#' @param use_dims Integer vector. Which PCs to use for t-SNE/UMAP downstream (e.g., 1:20).
#' @param tsne_perp Numeric. Perplexity parameter for t-SNE.
#' @param tsne_iter Integer. Maximum iterations for t-SNE.
#' @param umap_neighbors Integer. Number of nearest neighbors for UMAP.
#' @param umap_mindist Numeric. Minimum distance parameter for UMAP.
#'
#' @return A \code{RNA} object with embedded coordinates.
#'
#' @importFrom Matrix rowMeans t
#' @export
RunDimReduction_RNA <- function(obj,
                                method = "PCA",
                                layer_name = "data",
                                n_hvg = 2000,
                                pca_rank = 50,
                                use_dims = 1:20,
                                tsne_perp = 30,
                                tsne_iter = 1000,
                                umap_neighbors = 30,
                                umap_mindist = 0.3) {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # Get specified data layer
  mat <- obj$assays$RNA[[layer_name]]
  if(is.null(mat)) stop(paste("Error: Cannot find data layer", layer_name, ", please check if normalization or Scale has been performed!"))

  # Algorithm A: PCA
  if (method == "PCA") {
    message(sprintf("Calculating highly variable genes (HVG) on %s layer...", layer_name))
    # Extremely fast calculation of HVG
    rM <- Matrix::rowMeans(mat)
    rVar <- Matrix::rowMeans(mat^2) - rM^2
    hvgs <- names(sort(rVar, decreasing = TRUE))[1:min(n_hvg, length(rVar))]
    obj$var.genes <- hvgs

    message(sprintf("Running truncated SVD (PCA) using %d feature genes...", length(hvgs)))
    mat_hvg <- mat[hvgs, , drop = FALSE]
    pca_res <- irlba::prcomp_irlba(Matrix::t(mat_hvg), n = pca_rank, center = TRUE, scale. = TRUE)

    pca_coords <- pca_res$x
    rownames(pca_coords) <- colnames(mat)
    colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))

    obj$reductions$pca <- pca_coords
    obj$reductions$current_plot <- pca_coords[, 1:2]
    colnames(obj$reductions$current_plot) <- c("Dim1", "Dim2")
    obj$reductions$current_method <- "PCA"
    message("PCA complete!")
  }

  # Algorithm B: t-SNE
  else if (method == "t-SNE") {
    if(is.null(obj$reductions$pca)) stop("PCA results not found! Please run PCA before running t-SNE.")
    message(sprintf("Running t-SNE... (Using PC dimensions: %d to %d)", min(use_dims), max(use_dims)))

    pca_input <- obj$reductions$pca[, use_dims, drop = FALSE]
    tsne_res <- Rtsne::Rtsne(pca_input, pca = FALSE,
                             perplexity = tsne_perp,
                             max_iter = tsne_iter,
                             check_duplicates = FALSE)

    tsne_coords <- tsne_res$Y
    rownames(tsne_coords) <- rownames(pca_input)
    colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")

    obj$reductions$tsne <- tsne_coords
    obj$reductions$current_plot <- tsne_coords
    colnames(obj$reductions$current_plot) <- c("Dim1", "Dim2")
    obj$reductions$current_method <- "t-SNE"
    message("t-SNE coordinates generated successfully!")
  }

  # Algorithm C: UMAP
  else if (method == "UMAP") {
    if(is.null(obj$reductions$pca)) stop("PCA results not found! Please run PCA before running UMAP.")
    message(sprintf("Running UMAP... (Using PC dimensions: %d to %d)", min(use_dims), max(use_dims)))

    pca_input <- obj$reductions$pca[, use_dims, drop = FALSE]
    umap_res <- uwot::umap(pca_input,
                           n_neighbors = umap_neighbors,
                           min_dist = umap_mindist,
                           n_components = 2,
                           fast_sgd = TRUE)

    umap_coords <- umap_res
    rownames(umap_coords) <- rownames(pca_input)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

    obj$reductions$umap <- umap_coords
    obj$reductions$current_plot <- umap_coords
    colnames(obj$reductions$current_plot) <- c("Dim1", "Dim2")
    obj$reductions$current_method <- "UMAP"
    message("UMAP coordinates generated successfully!")
  } else {
    stop("Unsupported dimensionality reduction method, please select 'PCA', 't-SNE', or 'UMAP'.")
  }

  return(obj)
}

#' Run Clustering for RNA Object
#'
#' @description Performs cell/sample clustering using either a Graph-based (Leiden) approach or Hierarchical clustering.
#'
#' @param obj A \code{RNA} object (Requires prior dimensionality reduction).
#' @param reduction Character. Which reduction to use: "pca" (recommended), "umap", or "tsne".
#' @param method Character. Clustering algorithm: "graph" (KNN + Leiden) or "hierarchical".
#' @param use_dims Integer vector. Dimensions to extract from the reduction (e.g., 1:20).
#' @param cluster_res Numeric. Resolution parameter for Leiden algorithm (larger values = more clusters).
#' @param k_neighbors Integer. Number of neighbors for K-NN graph construction.
#' @param cluster_k Integer. Fixed number of clusters (k) for hierarchical clustering.
#'
#' @return A \code{RNA} object with \code{Auto_Cluster} added to its \code{filter_meta.data}.
#'
#' @importFrom stats dist hclust cutree
#' @importFrom igraph graph_from_edgelist simplify cluster_leiden
#' @export
RunClustering_RNA <- function(obj,
                              reduction = "pca",
                              method = "graph",
                              use_dims = 1:20,
                              cluster_res = 0.8,k_neighbors = 20,
                              cluster_k = 5) {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # --- 1. Dynamically extract specified dimensionality reduction data ---
  # Standardize name format to prevent inconsistent casing or hyphens from user input
  reduction_key <- tolower(reduction)
  if (reduction_key == "t-sne") reduction_key <- "tsne"

  coord_mat <- obj$reductions[[reduction_key]]

  if (is.null(coord_mat)) {
    stop(sprintf("Error: Cannot find dimensionality reduction results for '%s'! Please run RunDimReduction_RNA first to generate this data.", reduction))
  }

  # --- 2. Dimension safety truncation ---
  # Prevent specified dimensions from exceeding actual existing dimensions (e.g., UMAP usually has only 2 dimensions, but use_dims might be 1:20)
  start_dim <- max(1, min(use_dims))
  end_dim <- min(max(use_dims), ncol(coord_mat))

  if (start_dim > end_dim) {
    stop(sprintf("Error: Specified use_dims is out of range. '%s' only has %d dimensions.", reduction, ncol(coord_mat)))
  }

  cluster_input <- coord_mat[, start_dim:end_dim, drop = FALSE]
  message(sprintf("Extracting %s (dimensions %d to %d) for clustering...", toupper(reduction), start_dim, end_dim))

  # --- 3. Algorithm A: Graph-based (KNN + Leiden) ---
  if (method == "graph") {
    message(sprintf("Constructing K-NN graph (FNN, k=%d)...", k_neighbors))
    # Prevent cases where sample size is less than k
    k_n <- min(k_neighbors, nrow(cluster_input) - 1)
    if (k_n < 1) stop("Sample size is too small to construct a K-NN graph!")

    knn_res <- FNN::get.knn(cluster_input, k = k_n)
    edges <- do.call(rbind, lapply(1:nrow(cluster_input), function(i) {
      cbind(rep(i, k_n), knn_res$nn.index[i, ])
    }))

    message(sprintf("Running Leiden community detection algorithm (Resolution=%.2f)...", cluster_res))
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    g <- igraph::simplify(g)
    leiden_res <- igraph::cluster_leiden(g, resolution_parameter = cluster_res)
    cluster_labels <- paste0("Cluster_", leiden_res$membership)

  }
  # --- 4. Algorithm B: Hierarchical Clustering ---
  else if (method == "hierarchical") {
    message(sprintf("Running hierarchical clustering (ward.D2, k=%d)...", cluster_k))
    d_mat <- stats::dist(cluster_input)
    hc_res <- stats::hclust(d_mat, method = "ward.D2")
    # Prevent k from being greater than sample size
    actual_k <- min(cluster_k, nrow(cluster_input))
    cluster_labels <- paste0("Cluster_", stats::cutree(hc_res, k = actual_k))

  } else {
    stop("Unsupported clustering method, please select 'graph' or 'hierarchical'.")
  }

  # --- 5. Factorization: Elegant sorting ---
  # Ensure Cluster_2 is ordered before Cluster_10
  cluster_nums <- as.numeric(gsub("Cluster_", "", cluster_labels))
  sorted_levels <- paste0("Cluster_", sort(unique(cluster_nums)))
  cluster_factor <- factor(cluster_labels, levels = sorted_levels)

  # Safely write labels into Metadata
  if (!is.null(obj$filter_meta.data)) {
    obj$filter_meta.data$Auto_Cluster <- cluster_factor
  }
  # Also write to the original meta.data to maintain global data structure integrity
  common_cells <- intersect(rownames(obj$meta.data), rownames(cluster_input))
  obj$meta.data[common_cells, "Auto_Cluster"] <- cluster_factor

  message(sprintf("Clustering complete! Found %d clusters. Written to filter_meta.data$Auto_Cluster.", length(sorted_levels)))
  return(obj)
}

#' Run Differential Expression Analysis (Wilcoxon Rank Sum Test)
#'
#' @description Extremely fast implementation of Wilcoxon Rank Sum Test for discovering differentially expressed genes between groups.
#' Automatically handles 1-vs-1 or 1-vs-Rest comparisons.
#'
#' @param obj A \code{RNA} object.
#' @param layer_name Character. Data layer to use: "data" (recommended for log2FC) or "scale.data".
#' @param group_col Character. Column name in \code{filter_meta.data} used for grouping cells.
#' @param ident_1 Character. The target cluster/group name (Required).
#' @param ident_2 Character. The control group name. If \code{NULL}, computes \code{ident_1} vs all other groups (Rest).
#' @param min_pct Numeric. Minimum fraction of cells expressing the gene in either group to retain for testing.
#' @param logfc_thresh Numeric. Minimum absolute log2 fold change (or mean difference) required.
#'
#' @return A sorted \code{data.frame} of DEA results containing p-values, adjusted p-values, FC/Diff, and expression percentages.
#'
#' @importFrom Matrix rowSums rowMeans
#' @importFrom stats wilcox.test p.adjust
#' @export
RunDEA_RNA <- function(obj,
                       layer_name = "data",
                       group_col = "Auto_Cluster",
                       ident_1 = NULL,ident_2 = NULL,
                       min_pct = 0.1,logfc_thresh = 0.25) {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # --- 1. Get Metadata and check grouping ---
  meta <- obj$filter_meta.data
  if (is.null(meta)) stop("filter_meta.data not found, please run QC filtering first!")
  if (!(group_col %in% colnames(meta))) stop(sprintf("Column name '%s' is not in filter_meta.data!", group_col))

  group_vec <- as.character(meta[[group_col]])

  # Mandatory requirement: ident_1 cannot be null
  if (is.null(ident_1)) {
    stop("Error: ident_1 cannot be null! A target group must be specified.")
  }
  if (!(ident_1 %in% group_vec)) {
    stop(sprintf("Error: ident_1 ('%s') does not exist in column %s!", ident_1, group_col))
  }

  cells_1 <- rownames(meta)[which(group_vec == ident_1)]

  # --- 2. Smartly divide cell groups (automatically determine comparison mode) ---
  if (is.null(ident_2)) {
    # If ident_2 is null, automatically perform ident_1 vs Rest (all other groups)
    cells_2 <- rownames(meta)[which(group_vec != ident_1 & !is.na(group_vec))]
    ident_2_name <- "Rest"
    message(sprintf("ident_2 not specified, automatically performing analysis of '%s' vs 'all other groups (Rest)'.", ident_1))

  } else {
    # If ident_2 has a value, automatically perform ident_1 vs ident_2
    if (ident_1 == ident_2) stop("Error: ident_1 and ident_2 cannot be the same group!")
    if (!(ident_2 %in% group_vec)) stop(sprintf("Error: ident_2 ('%s') does not exist in column %s!", ident_2, group_col))

    cells_2 <- rownames(meta)[which(group_vec == ident_2)]
    ident_2_name <- ident_2
    message(sprintf("ident_2 detected, performing precise comparison: '%s' vs '%s'.", ident_1, ident_2))
  }

  # Check cell count
  if (length(cells_1) < 3 || length(cells_2) < 3) {
    stop("The number of cells/samples in one group is less than 3, unable to perform reliable statistical tests!")
  }

  # --- 3. Get matrix layer ---
  if (!(layer_name %in% c("data", "scale.data"))) {
    stop("layer_name must be 'data' or 'scale.data'")
  }
  mat <- obj$assays$RNA[[layer_name]]
  if (is.null(mat)) stop(sprintf("Layer '%s' not found, please run normalization or Scale first.", layer_name))

  message(sprintf("Extracting data for Group 1 (%s, N=%d) and Group 2 (%s, N=%d)...", ident_1, length(cells_1), ident_2_name, length(cells_2)))

  # Split matrix
  mat_1 <- mat[, cells_1, drop = FALSE]
  mat_2 <- mat[, cells_2, drop = FALSE]

  # --- 4. Extremely fast calculation of expression percentages (pct.1, pct.2) ---
  message("Calculating expression percentages...")
  pct_1 <- round(Matrix::rowSums(mat_1 > 0) / length(cells_1), 3)
  pct_2 <- round(Matrix::rowSums(mat_2 > 0) / length(cells_2), 3)

  # --- 5. Calculate Fold Change / Difference ---
  message("Calculating fold change magnitude...")
  if (layer_name == "data") {
    # Assuming data is log1p transformed, calculating rigorous log2FC:
    # First use expm1 to restore original expression and calculate mean, then calculate log2.
    mat_1_exp <- mat_1; mat_1_exp@x <- expm1(mat_1_exp@x)
    mat_2_exp <- mat_2; mat_2_exp@x <- expm1(mat_2_exp@x)

    mean_1 <- Matrix::rowMeans(mat_1_exp)
    mean_2 <- Matrix::rowMeans(mat_2_exp)
    # Add pseudo-count to prevent log2(0)
    log2fc <- log2(mean_1 + 1e-9) - log2(mean_2 + 1e-9)
    fc_col_name <- "avg_log2FC"
  } else {
    # If scale.data (containing negative z-scores), calculating FC will produce NaN
    # In this case, directly calculate Mean Difference
    mean_1 <- Matrix::rowMeans(mat_1)
    mean_2 <- Matrix::rowMeans(mat_2)
    log2fc <- mean_1 - mean_2
    fc_col_name <- "mean_diff" # Change column name to prevent misleading
  }

  # --- 6. Pre-filtering (significantly accelerates subsequent p-value calculation) ---
  genes_pass <- names(which(
    (pct_1 >= min_pct | pct_2 >= min_pct) & abs(log2fc) >= logfc_thresh
  ))

  if (length(genes_pass) == 0) {
    stop("No genes passed the pre-filtering of min_pct and logfc_thresh! Please try lowering the thresholds.")
  }

  # --- 7. Wilcoxon Rank Sum Test ---
  message(sprintf("Performing Wilcoxon Rank Sum Test on %d candidate genes...", length(genes_pass)))

  # Extract sub-matrix to be tested and convert to ordinary matrix (improves speed of row-by-row extraction)
  test_mat_1 <- as.matrix(mat_1[genes_pass, , drop = FALSE])
  test_mat_2 <- as.matrix(mat_2[genes_pass, , drop = FALSE])

  # === Changes 2: Replace `sapply` with `vapply`, and replace `1:length()` with `seq_along()` ===
  # The goal is to return a `p.value` (numeric type) for each iteration
  p_vals <- vapply(seq_along(genes_pass), function(i) {
    x <- test_mat_1[i, ]
    y <- test_mat_2[i, ]
    # Disable exact and correct to accelerate large-scale calculation
    res <- stats::wilcox.test(x, y, exact = FALSE, correct = FALSE)
    return(as.numeric(res$p.value))
  }, FUN.VALUE = numeric(1))

  # --- 8. Format results and multiple testing correction ---
  message("Formatting output results...")
  res_df <- data.frame(
    gene = genes_pass,
    p_val = p_vals,
    FC_metric = log2fc[genes_pass],
    pct.1 = pct_1[genes_pass],
    pct.2 = pct_2[genes_pass],
    cluster = ident_1,
    comparison = paste0(ident_1, " vs ", ident_2_name),
    stringsAsFactors = FALSE
  )

  # Dynamically name FC column
  colnames(res_df)[3] <- fc_col_name

  # Calculate FDR (BH correction)
  res_df$p_val_adj <- stats::p.adjust(res_df$p_val, method = "BH")

  # Adjust column order: gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, comparison
  res_df <- res_df[, c("gene", "p_val", fc_col_name, "pct.1", "pct.2", "p_val_adj", "cluster", "comparison")]

  # Sort ascending by p_val
  res_df <- res_df[order(res_df$p_val), ]
  rownames(res_df) <- NULL

  message(sprintf("DEA complete! Found %d significantly differentially expressed genes.", nrow(res_df)))
  return(res_df)
}

#' Run Pseudotime and Trajectory Inference
#'
#' @description Computes cellular pseudotime using either a cluster-based Principal Curve method or a graph-based K-NN shortest path method.
#'
#' @param obj A \code{RNA} object.
#' @param reduction Character. Embedding to use: "umap" (recommended), "tsne", or "pca".
#' @param group_col Character. Column name in \code{filter_meta.data} defining clusters. Defaults to "Auto_Cluster".
#' @param start_clus Character. Name of the root/starting cluster (Required).
#' @param algorithm Character. Trajectory algorithm: "cluster" (Principal Curve) or "graph" (KNN Shortest Path).
#'
#' @return A \code{RNA} object with \code{Pseudotime} stored in \code{filter_meta.data} and detailed trajectory info in \code{reductions$pseudotime}.
#'
#' @importFrom igraph graph_from_edgelist simplify distances V
#' @export
RunPseudotime_RNA <- function(obj,
                              reduction = "umap",
                              group_col = "Auto_Cluster",
                              start_clus = NULL,algorithm = "cluster") {

  if (!inherits(obj, "RNA")) stop("Input must be a RNA object!")

  # --- 1. Strictly lock filter_meta.data ---
  meta <- obj$filter_meta.data
  if (is.null(meta)) stop("Error: filter_meta.data not found! Please run QC filtering first.")
  if (!(group_col %in% colnames(meta))) stop(sprintf("Error: Column name '%s' is not in filter_meta.data!", group_col))
  if (is.null(start_clus)) stop("Error: Developmental root (start_clus) must be specified! E.g., start_clus = 'Cluster_1'")

  # --- 2. Get root samples ---
  group_vec <- as.character(meta[[group_col]])
  start_cells <- rownames(meta)[which(group_vec == start_clus)]

  if (length(start_cells) == 0) {
    stop(sprintf("Error: In the '%s' column of filter_meta.data, the group named '%s' cannot be found!", group_col, start_clus))
  }
  message(sprintf("Identified %d root samples (belonging to %s).", length(start_cells), start_clus))

  # --- 3. Get dimensionality reduction coordinates and strictly align with filter_meta.data ---
  reduction_key <- tolower(reduction)
  if (reduction_key == "t-sne") reduction_key <- "tsne"

  raw_coords <- obj$reductions[[reduction_key]]
  if (is.null(raw_coords)) stop(sprintf("Error: Cannot find dimensionality reduction data '%s'! Please run dimensionality reduction analysis first.", reduction))

  # Core alignment: Only keep coordinates of samples present in filter_meta.data
  valid_cells <- intersect(rownames(meta), rownames(raw_coords))
  if (length(valid_cells) == 0) stop("Error: Sample names in filter_meta.data and dimensionality reduction coordinates do not match at all!")

  coords <- raw_coords[valid_cells, 1:2, drop = FALSE] # Force restriction to 2D space for pseudotime
  start_cells <- intersect(start_cells, valid_cells)   # Ensure root cells are valid

  # --- 4. Core algorithm branching ---
  if (algorithm == "cluster") {
    # Engine A: Cluster-based (Principal Curve fitting, similar to Slingshot)
    message("Running Cluster-based principal curve fitting (princurve)...")
    if (!requireNamespace("princurve", quietly = TRUE)) stop("Please install the princurve package (install.packages('princurve'))!")

    fit <- princurve::principal_curve(coords, smoother = "smooth.spline")
    pt <- fit$lambda
    names(pt) <- rownames(coords)
    curve_out <- fit$s # Extract smoothed principal curve backbone coordinates

    message("Calibrating the biological time arrow...")
    mean_pt_start <- mean(pt[names(pt) %in% start_cells], na.rm = TRUE)
    mean_pt_all <- mean(pt, na.rm = TRUE)

    # If the mean time of the root is greater than the global mean, the curve is reversed and needs to be flipped
    if (mean_pt_start > mean_pt_all) {
      pt <- max(pt, na.rm = TRUE) - pt
    }

  } else if (algorithm == "graph") {
    # Engine B: Graph-based (Graph theory shortest path, similar to Monocle3)
    message("Running Graph-based K-NN graph network inference...")

    k_neighbors <- min(15, nrow(coords) - 1)
    knn_res <- FNN::get.knn(coords, k = k_neighbors)
    edges <- do.call(rbind, lapply(1:nrow(coords), function(i) {
      cbind(rep(i, k_neighbors), knn_res$nn.index[i, ])
    }))

    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    g <- igraph::simplify(g)

    message("Calculating shortest path distance of the network graph (Pseudotime)...")
    start_indices <- which(rownames(coords) %in% start_cells)

    dist_matrix <- igraph::distances(g, v = start_indices, to = igraph::V(g))
    pt <- apply(dist_matrix, 2, min)
    names(pt) <- rownames(coords)

    # Handle isolated cells (unreachable nodes)
    finite_max <- max(pt[is.finite(pt)])
    pt[is.infinite(pt)] <- finite_max * 1.1

    curve_out <- NULL # Dendrograms do not output a single backbone coordinate

  } else {
    stop("Unsupported algorithm parameter! Please select 'cluster' or 'graph'.")
  }

  # --- 5. Result formatting and saving ---
  message("Normalizing pseudotime to between 0 and 1...")
  pt <- (pt - min(pt, na.rm = TRUE)) / (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))

  # Save to standard underlying location
  obj$reductions$pseudotime <- list(
    pseudotime = pt,
    curve_coords = curve_out,
    dr_method = reduction_key,
    start_clus = start_clus,
    group_col = group_col,
    algorithm = algorithm
  )

  # Key step: Write directly to filter_meta.data for fast downstream access
  # Since the intersection was done earlier with valid_cells, it is very safe to directly assign by name to filter_meta.data.
  obj$filter_meta.data$Pseudotime <- NA # Fill entirely with NA first as a fallback
  obj$filter_meta.data[names(pt), "Pseudotime"] <- pt

  message(sprintf("Trajectory inference complete (%s)! Pseudotime data has been saved to filter_meta.data$Pseudotime.", toupper(algorithm)))
  return(obj)
}
