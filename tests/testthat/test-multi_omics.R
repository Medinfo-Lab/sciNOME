# tests/testthat/test-multi_omics.R

test_that("Multi-Omics Integration Logic works", {
  # 1. Generate simulated data
  meta_df <- data.frame(
    group1 = rep("Tumor", 3),
    RNA_ID = paste0("cell", 1:3),
    CpG_ID = paste0("samp", 1:3),
    GpC_ID = paste0("samp", 1:3)
  )

  region_df <- data.frame(
    chr = c("chr1", "chr2"),
    start = c(100, 200),
    end = c(150, 250),
    gene_id = c("GeneA", "GeneB")
  )

  rna_mat <- matrix(c(10, 2, 5, 12, 6, 1), nrow = 2, byrow = TRUE,
                    dimnames = list(c("GeneA", "GeneB"), paste0("cell", 1:3)))
  rna_obj <- list(assays = list(RNA = list(data = rna_mat)))

  # rownames must be completely consistent with chrdata formed by combining region_df
  epi_rownames <- c("chr1:100-150", "chr2:200-250")
  cpg_mat <- matrix(runif(6), nrow = 2, dimnames = list(epi_rownames, paste0("samp", 1:3)))
  gpc_mat <- matrix(runif(6), nrow = 2, dimnames = list(epi_rownames, paste0("samp", 1:3)))

  # 2. Test Integration Function (Tri-omics Model)
  res_tri <- Integrate_MultiOmics(
    mode = "tri", target_group = "Tumor", meta_df = meta_df, group_col = "group1",
    region_df = region_df,
    rna_obj = rna_obj, rna_id_col = "RNA_ID",
    cpg_mat = cpg_mat, cpg_id_col = "CpG_ID",
    gpc_mat = gpc_mat, gpc_id_col = "GpC_ID"
  )

  expect_true(is.data.frame(res_tri))
  expect_true(all(c("RNA_Exp", "CpG_level", "GpC_level") %in% colnames(res_tri)))
  expect_equal(nrow(res_tri), 2) # Two genes should be integrated
})
