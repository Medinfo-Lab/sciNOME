# tests/testthat/test-epi.R

test_that("Epigenetic Matrix Operations (Aggregate, Extract, QC, Reduce, DEA)", {
  set.seed(42)

  # --- 1. Test Aggregate_epiRegions (simulated file read/write) ---
  tmp_dir <- tempdir()

  bed_df <- data.frame(chr = "chr1", start = 100, end = 200)
  bed_path <- file.path(tmp_dir, "test.bed")
  write.table(bed_df, bed_path, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

  cov_df <- data.frame(chr = "chr1", pos = 150, dummy = 0, dummy2 = 0, meth = 8, unmeth = 2)
  cov_path <- file.path(tmp_dir, "Sample1.cov")
  write.table(cov_df, cov_path, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

  agg_res <- Aggregate_epiRegions(cov_dir = tmp_dir, bed_file = bed_path, n_cores = 1)
  expect_true(is.data.frame(agg_res))
  expect_equal(agg_res$Sample1.meth[1], 8)

  # --- 2. Test Extract and QC ---
  mock_agg <- data.frame(
    region_id = c("reg1", "reg2"),
    S1.meth = c(8, NA), S1.nonmeth = c(2, NA), S1.level = c(0.8, NA),
    S2.meth = c(9, 1), S2.nonmeth = c(1, 9), S2.level = c(0.9, 0.1)
  )

  level_df <- Extract_epiData(mock_agg, suffix = ".level")
  expect_equal(ncol(level_df), 3) # region, S1.level, S2.level

  # Test QC (allowing a maximum NA proportion of 0.6, S1 has half NA, will be retained)
  qc_df <- QC_epiData(mock_agg, top_n_rows = 2, max_col_na_ratio = 0.6)
  expect_true("S1.level" %in% colnames(qc_df))

  # --- 3. Test Reduce_epiDims (dimensionality reduction) ---
  mock_lvl <- data.frame(
    region_id = paste0("reg", 1:5),
    S1 = runif(5), S2 = runif(5), S3 = runif(5), S4 = runif(5)
  )
  mock_meta <- data.frame(ID = paste0("S", 1:4), Group = rep(c("T", "N"), each = 2))

  pca_res <- Reduce_epiDims(mock_lvl, mock_meta, sample_col = "ID", group_col = "Group", dr_method = "PCA")
  expect_equal(nrow(pca_res), 4) # four samples
  expect_true("Dim1" %in% colnames(pca_res))

  # --- 4. Test Epigenetic Differential Analysis Run_Diffanalysis ---
  mock_raw <- data.frame(
    region_id = c("reg1", "reg2", "reg3"),
    S1.lvl=c(0.9,0.9,0.1), S1.m=c(9,9,1), S1.nm=c(1,1,9),
    S2.lvl=c(0.9,0.9,0.1), S2.m=c(9,9,1), S2.nm=c(1,1,9),
    C1.lvl=c(0.1,0.1,0.9), C1.m=c(1,1,9), C1.nm=c(9,9,1),
    C2.lvl=c(0.1,0.1,0.9), C2.m=c(1,1,9), C2.nm=c(9,9,1)
  )
  mock_meta2 <- data.frame(
    Sample = c("S1", "S2", "C1", "C2"), Group = c("Tumor", "Tumor", "Normal", "Normal"),
    L = c("S1.lvl", "S2.lvl", "C1.lvl", "C2.lvl"),
    M = c("S1.m", "S2.m", "C1.m", "C2.m"),
    N = c("S1.nm", "S2.nm", "C1.nm", "C2.nm")
  )

  diff_res <- Run_Diffanalysis(mock_raw, mock_meta2, group_col = "Group", target_group = "Tumor",
                               col_level = "L", col_meth = "M", col_nonmeth = "N", verbose = FALSE)
  expect_true("FDR" %in% colnames(diff_res))
  expect_true("Significance" %in% colnames(diff_res))
})
