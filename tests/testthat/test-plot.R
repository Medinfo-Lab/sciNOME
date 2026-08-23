# tests/testthat/test-plot.R

test_that("Plotting functions return ggplot/patchwork objects", {

  # --- Prepare basic simulation objects ---
  mock_meta <- data.frame(
    orig.ident = c("A", "B"), Auto_Cluster = c("C1", "C2"),
    nFeature_RNA = c(1000, 2000), nCount_RNA = c(3000, 4000), percent_mt = c(5, 10),
    row.names = c("Cell_1", "Cell_2")
  )
  mock_umap <- matrix(runif(4), ncol = 2, dimnames = list(c("Cell_1", "Cell_2"), c("UMAP_1", "UMAP_2")))

  mock_rna <- list(
    meta.data = mock_meta,
    filter_meta.data = mock_meta,
    reductions = list(umap = mock_umap)
  )
  class(mock_rna) <- "RNA"

  # 1. PlotQC_RNA
  p_qc <- PlotQC_RNA(mock_rna)
  expect_s3_class(p_qc, "ggplot")

  # 2. PlotDimRed_RNA
  p_dim <- PlotDimRed_RNA(mock_rna, reduction = "umap")
  expect_s3_class(p_dim, "ggplot")

  # 3. PlotVolcano_RNA
  mock_dea <- data.frame(gene = "GeneA", avg_log2FC = 2, p_val_adj = 0.01)
  p_volc <- PlotVolcano_RNA(mock_dea)
  expect_s3_class(p_volc, "ggplot")

  # 4. PlotDimRed_Epi
  mock_dr_epi <- data.frame(Dim1 = c(1,2), Dim2 = c(3,4), Group = c("A","B"))
  p_epi_pca <- PlotDimRed_Epi(mock_dr_epi, group_col = "Group", method = "PCA")
  expect_s3_class(p_epi_pca, "ggplot")

  # 5. PlotVolcano_Epi
  mock_epi_dea <- data.frame(chrdata = "chr1:1-2", Diff = 1.5, P.Value = 0.01)
  p_epi_volc <- PlotVolcano_Epi(mock_epi_dea, feature_col = "chrdata")
  expect_s3_class(p_epi_volc, "ggplot")

  # 6. PlotOmicsMatrix
  skip_if_not_installed("patchwork")
  mock_omics <- data.frame(RNA_Exp = c(1,2), CpG_level = c(0.1,0.2), GpC_level = c(0.8,0.9))
  p_matrix <- PlotOmicsMatrix(mock_omics)
  expect_true(inherits(p_matrix, "ggplot") || inherits(p_matrix, "patchwork"))
})
