#' UMAP dimensionality reduction data was used for quasi-temporal analysis
#'
#' @param rds_tmp RDS objection
#'
#' @return monocle3 object and Pseudotime data
#' @export
#'
#' @examples
Pseudotime_analysis <- function(rds_tmp){
  pbmc <- rds_tmp
  expression_matrix <- GetAssayData(pbmc,assay="RNA",slot="counts")
  cell_metadata <- pbmc@meta.data
  gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix),
                              row.names = rownames(expression_matrix))
  cds <- new_cell_data_set(expression_data = expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)
  #preprocess_cds==seurat NormalizeData+ScaleData+RunPCA
  # cds <- preprocess_cds(cds)
  # cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP")

  int.embed <- Embeddings(pbmc,reduction = "umap")
  # int.embed <- int.embed[rownames(cds.embed),]
  reducedDims(cds)$UMAP <- int.embed


  cds <- cluster_cells(cds,reduction_method = "UMAP")
  cds <- learn_graph(cds,
                     learn_graph_control = list(
                       minimal_branch_len = 1,
                       euclidean_distance_ratio = 1))
  cds <- order_cells(cds)

  pseudotime_values <- pseudotime(cds, reduction_method = "UMAP")
  pseudo_seurat <- as.data.frame(pseudotime_values)
  pseudo_seurat_UMAP <- cds@int_colData$reducedDims$UMAP
  pseudo_seurat_UMAP <- as.data.frame(pseudo_seurat_UMAP)
  pseudo_seurat$UMAP1 <- pseudo_seurat_UMAP$umap_1
  pseudo_seurat$UMAP2 <- pseudo_seurat_UMAP$umap_2
  pseudo_seurat$group <- pbmc$seurat_clusters

  pseudo_list <- list()
  pseudo_list$pseudo_object <- cds
  pseudo_list$pseudo_data <- pseudo_seurat
  return(pseudo_list)
}


#' RNA standard data
#'
#' @param file_tmp_length The length of each gene entered
#' @param method Methods for data standardization, TPM, FPKM, CPM
#' @param file_tmp Enter the RNA expression matrix
#'
#' @return Expression matrix after normalization
#' @export
#'
#' @examples
RNA_data_standard <- function(file_tmp,file_tmp_length,method="tpm"){
  file_data <- STAD_seurat_data
  file_gene_length <- STAD_expr_length_choose

  if (method=="tpm") {
    file_data_count <- file_data
    file_gene_length_kb <- file_gene_length[,1]/1000

    rpk <- file_data_count/file_gene_length_kb
    tpm_matrix <- t(t(rpk) / colSums(rpk) * 1e6)
    return(tpm_matrix)
  }else if(method=="fpkm"){
    file_data_count <- file_data
    file_gene_length_kb <- file_gene_length[,1]/1000

    rpk <- file_data_count/file_gene_length_kb
    fpkm_matrix <- t(t(rpk) / colSums(rpk) * 1e9)
    return(fpkm_matrix)
  }else if(method=="cpm"){
    file_data_count <- file_data
    total_counts <- colSums(file_data_count)

    cpm_matrix <- t(t(file_data_count) / total_counts) * 1e6
    return(cpm_matrix)
  }
}


#' LASSO gene
#'
#' @param file_tmp Enter the RNA expression matrix
#' @param file_group_tmp Grouping information for the sample
#' @param test_size Training-to-testing ratio
#'
#' @return Returns the expression matrix after selection
#' @export
#'
#' @examples
LASSO_feature <- function(file_tmp,file_group_tmp,test_size=0.3){
  source_python_file <- system.file(package = "scMethylCA")
  source_def <- paste(source_python_file,"python/lasso_feature.py",sep = "/")
  reticulate::source_python(source_def)
  # reticulate::source_python("inst/python/lasso_feature.py")

  result <- reticulate::py$lasso_feature(file_tmp,file_group_tmp,test_size)

  return(result)
}



#' Grouping model prediction
#'
#' @param file_tmp Enter the RNA expression matrix
#' @param test_size Training-to-testing ratio
#' @param model Model selection
#' @param file_group Grouping information for the sample
#'
#' @return Post-training data
#' @export
#'
#' @examples
Group_model_predict <- function(file_tmp,file_group,test_size=0.3,model="lgbm"){
  source_python_file <- system.file(package = "scMethylCA")
  source_def <- paste(source_python_file,"python/model_algorithm.py",sep = "/")
  reticulate::source_python(source_def)
  # reticulate::source_python("inst/python/model_algorithm.py")

  result <- reticulate::py$model_algorithm(file_tmp,file_group,test_size,model)

  return(result)
}
