#' Create an RNA storage object
#'
#' @param file_data RNA Expression Matrix
#' @param min_cells Include features detected in at least this many cells.
#' @param min_features 	Include cells where at least this many features are detected
#'
#' @returns RNA object
#' @export
#'
#' @examples
CreatRNAObject <- function(file_data,min_cells=50,min_features=0){
  RNA_pbmc <- CreateSeuratObject(counts = file_data,
                                 min.cells = min_cells, min.features = min_features)

  RNA_pbmc <- NormalizeData(RNA_pbmc, normalization.method = "LogNormalize")
  RNA_pbmc <- FindVariableFeatures(RNA_pbmc, selection.method = "vst", nfeatures = 2000)
  RNA_pbmc <- VariableFeatures(RNA_pbmc)
  RNA_pbmc <- ScaleData(RNA_pbmc,features = rownames(RNA_pbmc))
  return(RNA_pbmc)
}


#' Dimensionality reduction analysis
#'
#' @param file_data Transposed expression matrix
#' @param method Behavioral cells, classified as genes
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
#' @param min_dist The effective minimum distance between embedded points.
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#'
#' @returns Dimensionality-Reduced Dataset
#' @export
#'
#' @examples
DR_process <- function(file_data,method="RNA_UMAP",n_neighbors=100,min_dist=0.5,
                       center=T,scale=T,cor_method="pearson"){
  if(method=="RNA_UMAP"){
    umap_result <- uwot::umap(file_data,n_neighbors = n_neighbors, min_dist = min_dist, metric = "cosine")
    umap_df <- as.data.frame(umap_result)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    return(umap_df)
  }else if (method=="Epi_UMAP") {
    umap_result <- uwot::umap(file_data,n_neighbors=n_neighbors, min_dist = min_dist)
    umap_df <- as.data.frame(umap_result)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    return(umap_df)
  }else if (method=="PCA") {
    pca_res <- stats::prcomp(file_data, center = center, scale. = center, rank. = 30)
    var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
    pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
    pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
    pca_embedding <- as.data.frame(pca_res$x)

    pca_list <- list()
    pca_list$pca_embedding <- pca_embedding
    pca_list$pc1_lab <- pc1_lab
    pca_list$pc2_lab <- pc2_lab
    return(pca_list)
  }else if (method=="tSNE") {
    tsne_out <- Rtsne::Rtsne(file_data,
                             dims = 2,
                             perplexity = 30,
                             verbose = TRUE,
                             max_iter = 1000,
                             check_duplicates = F,
                             pca = TRUE)
    tsne_plot_data <- data.frame(
      tSNE_1 = tsne_out$Y[, 1],
      tSNE_2 = tsne_out$Y[, 2]
    )
    return(tsne_plot_data)
  }else if(method=="MDS"){
    cor_matrix <- cor(file_data, use = "pairwise.complete.obs", method = cor_method)
    dist_matrix <- as.dist(1 - cor_matrix)
    dist_matrix[is.na(dist_matrix)] <- 1
    mds_points <- cmdscale(dist_matrix, k = 2, eig = TRUE)
    return(mds_points)
  }
}

