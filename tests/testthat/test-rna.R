# tests/testthat/test-rna.R

test_that("RNA Object workflow (Build, QC, DimRed, Cluster, DEA, Pseudotime)", {
  set.seed(42)

  # 1. Generate a simulated matrix (20 genes x 20 cells)
  mock_counts <- matrix(rnbinom(400, mu = 10, size = 1), nrow = 20, ncol = 20)
  rownames(mock_counts) <- c("MT-ND1", "MT-ND2", paste0("Gene_", 3:20)) # Contains mitochondrial genes
  colnames(mock_counts) <- paste0("Cell_", 1:20)

  mock_meta <- data.frame(ID = colnames(mock_counts), Group = rep(c("A", "B"), each = 10))

  # 2. Test build object
  rna_obj <- Build_RNAObject(mock_counts, meta_data = mock_meta, meta_id_col = "ID")
  expect_s3_class(rna_obj, "RNA")
  expect_true("percent_mt" %in% colnames(rna_obj$meta.data))

  # 3. Test QC and standardization
  rna_obj <- ProcessQC_RNA(rna_obj, mt_pattern = "^MT-", min_nCount = 0, min_nFeature = 0, max_mt = 100, do_scale = FALSE)
  expect_false(is.null(rna_obj$assays$RNA$data)) # Check whether the normalization layer is generated

  # 4. Test dimensionality reduction (adaptive small sample parameters)
  rna_obj <- RunDimReduction_RNA(rna_obj, method = "PCA", n_hvg = 10, pca_rank = 3)
  expect_false(is.null(rna_obj$reductions$pca))

  # 5. Test clustering
  rna_obj <- RunClustering_RNA(rna_obj, reduction = "pca", method = "hierarchical", cluster_k = 2)
  expect_true("Auto_Cluster" %in% colnames(rna_obj$filter_meta.data))

  # 6. Test differential expression analysis (DEA)
  dea_res <- RunDEA_RNA(rna_obj, group_col = "Group", ident_1 = "A", ident_2 = "B", min_pct = 0, logfc_thresh = 0)
  expect_true(is.data.frame(dea_res))
  expect_true("p_val_adj" %in% colnames(dea_res))

  # 7. Test pseudotime analysis
  # Force the use of PCA space instead of UMAP for quick testing
  rna_obj <- RunPseudotime_RNA(rna_obj, reduction = "pca", group_col = "Group", start_clus = "A", algorithm = "cluster")
  expect_true("Pseudotime" %in% colnames(rna_obj$filter_meta.data))
})
