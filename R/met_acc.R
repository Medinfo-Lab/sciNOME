#' Aggregate Bismark Coverage Files into Defined Genomic Regions
#'
#' @description This function takes a directory of Bismark `.cov.gz` files and a BED file containing
#' genomic regions. It aggregates the methylation data across these regions in parallel.
#'
#' @param cov_dir Character. Path to the directory containing `.cov.gz` files.
#' @param bed_file Character. Path to the BED file. Must contain at least 3 columns: chr, start, end.
#' @param n_cores Integer. Number of cores to use for parallel processing. Default is 1.
#'
#' @return A \code{data.table} containing the aggregated results. The first column is \code{region_id},
#' followed by \code{.meth}, \code{.nonmeth}, and \code{.level} columns for each sample.
#'
#' @importFrom parallel makePSOCKcluster stopCluster clusterEvalQ
#' @importFrom pbapply pblapply
#' @importFrom dplyr filter mutate select
#' @importFrom data.table fread setDT
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @export
Aggregate_epiRegions <- function(cov_dir, bed_file, n_cores = 1) {
  # 1. Get the list of all cov files
  cov_files <- list.files(path = cov_dir, pattern = "\\.cov(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(cov_files) == 0) {
    stop("No .cov or .cov.gz files found in the specified directory.")
  }
  # 2. Read and process BED file (build master GRanges)
  message(sprintf("Loading BED file: %s", bed_file))
  bed_df <- data.table::fread(bed_file, select = 1:3, nThread = 1)
  # bed_df <- bed_data[,c(1:3)]
  data.table::setnames(bed_df, 1:3, c("chr", "start", "end"))
  bed_df[, region_id := paste0(chr, ":", start, "-", end)]
  master_bed_gr <- GenomicRanges::GRanges(
    seqnames = bed_df$chr,
    ranges = IRanges::IRanges(start = bed_df$start, end = bed_df$end)
  )
  master_region_ids <- bed_df$region_id
  n_regions <- length(master_bed_gr)
  message(sprintf("Found %d regions in BED file.", n_regions))

  # 3. Define single file processing logic (Worker Function)
  process_single_file <- function(f_path, local_master_gr, local_n_regions) {
    # Extract and clean sample names
    f_name <- basename(f_path)
    cell_name <- gsub("\\.cov(\\.gz)?$|\\.bed(\\.gz)?$|\\.txt(\\.gz)?$|\\.csv$", "", f_name, ignore.case = TRUE)
    cell_name <- gsub("[-_\\. ]+", ".", cell_name)
    cell_name <- gsub("^\\.|\\.$", "", cell_name)

    # Read Coverage data (assuming columns 1, 2, 5, 6 are chr, pos, meth, unmeth respectively)
    cov_df <- data.table::fread(f_path, select = c(1, 2, 5, 6), nThread = 1)
    data.table::setnames(cov_df, c("chr", "pos", "meth", "unmeth"))

    # Convert to GRanges and find overlapping regions
    cov_gr <- GenomicRanges::GRanges(
      seqnames = cov_df$chr,
      ranges = IRanges::IRanges(start = cov_df$pos, width = 1)
    )
    hits <- GenomicRanges::findOverlaps(cov_gr, local_master_gr)

    # Aggregate data
    hit_dt <- data.table::data.table(
      region_idx = S4Vectors::subjectHits(hits),
      meth = cov_df$meth[S4Vectors::queryHits(hits)],
      unmeth = cov_df$unmeth[S4Vectors::queryHits(hits)]
    )
    agg <- hit_dt[, .(meth = sum(meth, na.rm = TRUE), unmeth = sum(unmeth, na.rm = TRUE)), by = region_idx]

    # Initialize empty vectors (fill uncovered regions with NA)
    res_meth <- rep(NA_real_, local_n_regions)
    res_unmeth <- rep(NA_real_, local_n_regions)
    res_meth[agg$region_idx] <- agg$meth
    res_unmeth[agg$region_idx] <- agg$unmeth

    # Calculate Level
    res_level <- res_meth / (res_meth + res_unmeth)
    res_level[is.nan(res_level)] <- NA_real_

    # Pack into a list
    cell_vecs <- list(res_meth, res_unmeth, res_level)
    names(cell_vecs) <- c(paste0(cell_name, ".meth"),
                          paste0(cell_name, ".nonmeth"),
                          paste0(cell_name, ".level"))
    return(cell_vecs)
  }

  # 4. Setup parallel computing
  message(sprintf("Starting aggregation on %d files using %d cores...", length(cov_files), n_cores))
  if (n_cores > 1) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE) # Ensure cluster is stopped even if exiting with an error
    # Load necessary packages into child processes
    parallel::clusterEvalQ(cl, {
      library(data.table)
      library(GenomicRanges)
      library(IRanges)
    })
    # Use pbapply to show progress bar
    results_list <- pbapply::pblapply(
      cov_files,
      FUN = process_single_file,
      local_master_gr = master_bed_gr,
      local_n_regions = n_regions,
      cl = cl
    )
  } else {
    # Single-thread mode
    results_list <- pbapply::pblapply(
      cov_files,
      FUN = process_single_file,
      local_master_gr = master_bed_gr,
      local_n_regions = n_regions
    )
  }

  # 5. Assemble final matrix
  message("Finalizing matrix...")
  # Flatten the list of lists
  flat_list <- unlist(results_list, recursive = FALSE)
  # Merge region_id
  final_list <- c(list(region_id = master_region_ids), flat_list)
  final_dt <- data.table::setDT(final_list)
  message("Done!")
  return(final_dt)
}

#' Extract Specific Metric Matrix from Aggregated Epigenetic Data
#'
#' @description Extracts columns ending with a specific suffix (e.g., ".meth", ".nonmeth", ".level")
#' from the aggregated results, while always retaining the 'region_id' column.
#'
#' @param data A \code{data.frame} or \code{data.table} generated by \code{aggregate_epi_regions}.
#' @param suffix Character. The suffix to extract. Examples: ".meth", ".nonmeth", ".level".
#' @param clean_names Logical. If \code{TRUE} (default), removes the suffix from the column names
#' in the returned data, leaving only the clean sample names.
#'
#' @return A \code{data.table} containing 'region_id' and the extracted metric columns.
#'
#' @importFrom data.table as.data.table
#' @export
Extract_epiData <- function(data, suffix, clean_names = TRUE) {
  # 1. Input data check and conversion
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame or data.table.")
  }
  dt <- data.table::as.data.table(data)

  # 2. Check if the required region_id column exists
  if (!"region_id" %in% names(dt)) {
    stop("The input data must contain a 'region_id' column.")
  }

  # 3. Smart suffix handling (automatically add '.' if missing to prevent unexpected string matching)
  if (!startsWith(suffix, ".")) {
    suffix <- paste0(".", suffix)
  }

  # 4. Identify target columns
  all_cols <- names(dt)
  # Use endsWith for exact suffix matching (faster than regex and no need to escape '.')
  target_cols <- all_cols[endsWith(all_cols, suffix)]
  if (length(target_cols) == 0) {
    stop(sprintf("No columns found ending with the suffix '%s'.", suffix))
  }

  # 5. Extract data
  keep_cols <- c("region_id", target_cols)
  # Native data.table subset extraction
  res_dt <- dt[, ..keep_cols]

  # 6. (Optional) Clean column names, restore 'Sample1.level' to 'Sample1'
  # if (clean_names) {
  #   # Only modify non-region_id column names
  #   new_names <- gsub(paste0(suffix, "$"), "", target_cols)
  #   data.table::setnames(res_dt, old = target_cols, new = new_names)
  # }
  return(res_dt)
}

#' Two-Step Quality Control for Epigenetic Matrix
#'
#' @description Performs a two-step QC on a methylation level matrix.
#' Step 1: Sorts rows (regions) by the number of NA values in ascending order based ONLY on 'level' columns,
#'         and selects the top N rows (highest quality regions).
#' Step 2: Using ONLY the selected top rows, calculates the NA ratio for each sample's 'level' column.
#'         Retains samples with an NA ratio strictly less than the specified threshold.
#' Finally, it keeps the 'meth', 'nonmeth', and 'level' columns for all samples that passed the QC.
#'
#' @param data A \code{data.frame} or \code{data.table} containing the extracted
#'             metric (e.g., from \code{extract_epi_metric}). Must contain a 'region_id' column.
#' @param top_n_rows Integer. The number of top rows (regions with the fewest NAs) to retain.
#'                   Default is 5000.
#' @param max_col_na_ratio Numeric. The maximum allowed NA ratio for columns (samples).
#'                         Should be between 0 and 1. Default is 0.8 (i.e., < 80\% NAs).
#' @param level_suffix Character. The suffix identifying the 'level' columns. Default is "_level".
#' @param meth_suffix Character. The suffix identifying the 'meth' columns. Default is "_meth".
#' @param nonmeth_suffix Character. The suffix identifying the 'nonmeth' columns. Default is "_nonmeth".
#'
#' @return A filtered \code{data.table} retaining meth, nonmeth, and level columns for passed samples.
#'
#' @importFrom data.table as.data.table
#' @export
QC_epiData <- function(data,
                       top_n_rows = 5000,
                       max_col_na_ratio = 0.8,
                       level_suffix = ".level",
                       meth_suffix = ".meth",
                       nonmeth_suffix = ".nonmeth") {
  # 1. Input validation
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame or data.table.")
  }
  if (!"region_id" %in% names(data)) {
    stop("Input data must contain a 'region_id' column.")
  }
  if (max_col_na_ratio > 1 || max_col_na_ratio < 0) {
    if (max_col_na_ratio > 1 && max_col_na_ratio <= 100) {
      max_col_na_ratio <- max_col_na_ratio / 100 # If user inputs 80, automatically convert to 0.8
    } else {
      stop("'max_col_na_ratio' must be between 0 and 1 (e.g., 0.8 for 80%).")
    }
  }

  dt <- data.table::as.data.table(data)

  # 2. Identify all level columns and extract base sample names
  all_cols <- setdiff(names(dt), "region_id")
  level_cols <- grep(paste0(level_suffix, "$"), all_cols, value = TRUE)

  if (length(level_cols) == 0) {
    stop(sprintf("No level columns found! Please check if your columns end with '%s'.", level_suffix))
  }

  # Extract base sample names (e.g.: "Sample1_level" -> "Sample1")
  sample_names <- gsub(paste0(level_suffix, "$"), "", level_cols)

  message(sprintf("Original data: %d rows (regions) and %d samples (Total %d columns).",
                  nrow(dt), length(sample_names), ncol(dt)))

  # Step 1: Row QC - [Based ONLY on level columns] Sort by NA count and subset
  message(sprintf("Step 1: Sorting rows by NA count in 'level' columns and selecting top %d...", top_n_rows))
  # Calculate NA count per row (only for level columns)
  row_na_counts <- rowSums(is.na(as.matrix(dt[, ..level_cols])))
  # Get indices ordered by NA count in ascending order (fewer NAs, higher rank)
  ordered_indices <- order(row_na_counts, decreasing = FALSE)
  # Extract the top N rows
  keep_n <- min(top_n_rows, nrow(dt))
  top_indices <- ordered_indices[1:keep_n]
  dt_row_filtered <- dt[top_indices]
  message(sprintf(" -> Kept %d highest-quality regions.", nrow(dt_row_filtered)))

  # Step 2: Column QC - [Based ONLY on level columns] Calculate missing rate
  message(sprintf("Step 2: Filtering samples with NA ratio < %.1f%% based on selected regions...",
                  max_col_na_ratio * 100))

  # === 修改点 1：将 sapply 替换为 vapply ===
  # 明确声明：每次循环都期望返回一个长度为 1 的数值型数据 numeric(1)
  col_na_ratios <- vapply(dt_row_filtered[, ..level_cols], function(x) {
    sum(is.na(x)) / nrow(dt_row_filtered)
  }, FUN.VALUE = numeric(1))

  # Find passing level column names
  passed_level_cols <- names(col_na_ratios)[col_na_ratios < max_col_na_ratio]
  failed_level_cols <- names(col_na_ratios)[col_na_ratios >= max_col_na_ratio]

  # Map back to base sample names
  passed_samples <- gsub(paste0(level_suffix, "$"), "", passed_level_cols)
  failed_samples <- gsub(paste0(level_suffix, "$"), "", failed_level_cols)

  if (length(passed_samples) == 0) {
    stop("All samples failed the column QC! Try increasing 'max_col_na_ratio' or 'top_n_rows'.")
  }

  # Print filtering information
  if (length(failed_samples) > 0) {
    message(sprintf(" -> Removed %d sample(s) due to high NA ratio: %s",
                    length(failed_samples), paste(failed_samples, collapse = ", ")))
  } else {
    message(" -> All samples passed column QC.")
  }

  # Step 3: Assemble final data - [Bundle and keep meth, nonmeth, level columns]

  # === Change 2: Replace `unlist(lapply)` with the safer `as.vector(vapply)` ===
  # Explicitly specify: Each iteration is expected to return a vector containing 3 strings (`character(3)`)
  passed_sample_cols <- as.vector(vapply(passed_samples, function(s) {
    c(paste0(s, meth_suffix),
      paste0(s, nonmeth_suffix),
      paste0(s, level_suffix))
  }, FUN.VALUE = character(3)))

  # Intersect with original data column names (prevent errors due to missing columns in original data)
  valid_sample_cols <- intersect(passed_sample_cols, names(dt))

  final_cols <- c("region_id", valid_sample_cols)
  final_dt <- dt_row_filtered[, ..final_cols]

  message(sprintf("Final QC'd data: %d rows (regions) and %d samples (%d target columns retained).",
                  nrow(final_dt), length(passed_samples), length(valid_sample_cols)))

  return(final_dt)
}

#' Run Dimensionality Reduction on Epigenetic Level Data
#'
#' @description Preprocesses methylation/accessibility level matrices, performs smart sample
#' name matching with metadata, imputes missing values, and calculates dimensionality
#' reduction coordinates (PCA, UMAP, MDS, or NMF).
#'
#' @param mat Data.frame. The level data matrix. The first column MUST be the genomic region (e.g., Chromosome regions). Rows are regions, columns are samples.
#' @param meta Data.frame. Metadata containing sample IDs and grouping information.
#' @param sample_col Character. Column name in `meta` representing the sample IDs.
#' @param group_col Character. Column name in `meta` representing the sample groups/conditions.
#' @param dr_method Character. Dimensionality reduction method. Options: "PCA", "UMAP", "MDS", "NMF". (Default: "PCA").
#' @param impute_method Character. Method for handling NA values. Options: "knn", "mean". (Default: "mean").
#' @param knn_k Integer. Number of neighbors to use if \code{impute_method = "knn"}. (Default: 10).
#' @param umap_neighbors Integer. Number of neighbors for UMAP. (Default: 15).
#' @param nmf_rank Integer. Rank for NMF factorization. (Default: 2).
#'
#' @return A data.frame containing sample IDs, Dim1, Dim2, and original metadata columns, ready for plotting.
#'
#' @importFrom stats prcomp cor as.dist cmdscale
#' @export
Reduce_epiDims <- function(
    mat,
    meta,
    sample_col,
    group_col,
    dr_method = c("PCA", "UMAP", "MDS", "NMF"),
    impute_method = c("mean", "knn"),
    knn_k = 10,
    umap_neighbors = 15,
    nmf_rank = 2
) {

  # Parameter matching validation
  dr_method <- match.arg(toupper(dr_method), choices = c("PCA", "UMAP", "MDS", "NMF"))
  impute_method <- match.arg(tolower(impute_method), choices = c("mean", "knn"))

  # Step A: Extract matrix and process column names
  # Set the first column (chromosome regions) as row names, extract pure numeric matrix
  region_names <- as.character(mat[[1]])
  meth_mat <- as.matrix(mat[, -1, drop = FALSE])
  rownames(meth_mat) <- region_names

  # Ensure matrix is numeric (prevent conversion to character due to formatting issues during file reading)
  if(!is.numeric(meth_mat)) {
    mode(meth_mat) <- "numeric"
  }

  # Step B: Ultimate Smart Sample Matching (Smart Fuzzy Matching)
  meta_data <- as.data.frame(meta)

  if(!sample_col %in% colnames(meta_data)) stop(sprintf("Sample column not found in Metadata: %s", sample_col))
  if(!group_col %in% colnames(meta_data)) stop(sprintf("Group column not found in Metadata: %s", group_col))

  # Core cleaning function: strip all possible suffixes and special characters
  smart_clean <- function(x) {
    x <- gsub("\\.level$", "", x, ignore.case = TRUE)
    x <- gsub("\\.meth$", "", x, ignore.case = TRUE)
    x <- gsub("\\.nonmeth$", "", x, ignore.case = TRUE)
    x <- gsub("\\.cov$", "", x, ignore.case = TRUE)
    x <- gsub("\\.gz$", "", x, ignore.case = TRUE)
    x <- gsub("[-_ ]", ".", x)
    x <- gsub("\\.+", ".", x)
    x <- gsub("^\\.|\\.$", "", x)
    return(x)
  }

  mat_clean_names <- smart_clean(colnames(meth_mat))
  meta_clean_names <- smart_clean(as.character(meta_data[[sample_col]]))

  colnames(meth_mat) <- mat_clean_names
  meta_data$CleanID <- meta_clean_names

  # Get intersecting samples
  common_samples <- intersect(meta_data$CleanID, colnames(meth_mat))

  if(length(common_samples) < 3) {
    stop(paste0(
      "Sample matching failed! Less than 3 intersecting samples.\n",
      "Top 3 samples after matrix cleaning: ", paste(head(mat_clean_names, 3), collapse = ", "), "\n",
      "Top 3 samples after Meta cleaning: ", paste(head(meta_clean_names, 3), collapse = ", ")
    ))
  }

  # Align matrix and Metadata
  meth_mat <- meth_mat[, common_samples, drop = FALSE]
  meta_data <- meta_data[match(common_samples, meta_data$CleanID), , drop = FALSE]
  meta_data$SampleID <- meta_data$CleanID # Create standard ID column for subsequent Merge

  # Step C: Missing Value Imputation
  if (any(is.na(meth_mat))) {
    message("NA values detected, performing imputation...")

    if (impute_method == "knn") {
      if (!requireNamespace("impute", quietly = TRUE)) {
        stop("Using KNN imputation requires the 'impute' package. Please run BiocManager::install('impute') to install it.")
      }
      # Ensure k value is not greater than the number of samples
      safe_k <- min(knn_k, ncol(meth_mat) - 1)
      meth_mat <- impute::impute.knn(meth_mat, k = safe_k)$data

    } else if (impute_method == "mean") {
      row_means <- rowMeans(meth_mat, na.rm = TRUE)
      idx_na <- which(is.na(meth_mat), arr.ind = TRUE)
      meth_mat[idx_na] <- row_means[idx_na[, 1]]
    }
  }

  # Step D: Run Dimensionality Reduction Algorithm
  message(sprintf("Running %s dimensionality reduction...", dr_method))
  coords <- NULL

  if (dr_method == "PCA") {
    # PCA: Needs features as columns (transpose)
    pca_res <- prcomp(t(meth_mat), center = TRUE, scale. = FALSE)
    coords <- data.frame(
      SampleID = rownames(pca_res$x),
      Dim1 = pca_res$x[, 1],
      Dim2 = pca_res$x[, 2]
    )

  } else if (dr_method == "UMAP") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Using UMAP requires the 'uwot' package. Please run install.packages('uwot') to install it.")
    }
    n_samples <- ncol(meth_mat)
    safe_n_neighbors <- min(umap_neighbors, n_samples - 1)
    if (safe_n_neighbors < 2) safe_n_neighbors <- 2 # Fallback protection

    umap_res <- uwot::umap(
      X = t(meth_mat),
      n_neighbors = safe_n_neighbors,
      min_dist = 0.1,
      n_threads = 1,
      pca = NULL
    )
    coords <- data.frame(
      SampleID = colnames(meth_mat),
      Dim1 = umap_res[, 1],
      Dim2 = umap_res[, 2]
    )

  } else if (dr_method == "MDS") {
    # MDS: Based on correlation distance between samples
    dist_mat <- stats::as.dist(1 - stats::cor(meth_mat, method = "pearson"))
    mds_res <- stats::cmdscale(dist_mat, k = 2)
    coords <- data.frame(
      SampleID = rownames(mds_res),
      Dim1 = mds_res[, 1],
      Dim2 = mds_res[, 2]
    )

  } else if (dr_method == "NMF") {
    if (!requireNamespace("NMF", quietly = TRUE)) {
      stop("Using NMF requires the 'NMF' package. Please run install.packages('NMF') to install it.")
    }
    if (min(meth_mat, na.rm = TRUE) < 0) {
      stop("NMF algorithm requires the input matrix to be non-negative. Your data contains negative values!")
    }

    # NMF calculation (meth_mat keeps Features x Samples format)
    nmf_res <- NMF::nmf(meth_mat, rank = nmf_rank, seed = "random")
    # Extract coefficient matrix H (Rank x Samples), and transpose to get (Samples x Rank)
    h_mat <- t(NMF::coef(nmf_res))
    coords <- data.frame(
      SampleID = rownames(h_mat),
      Dim1 = h_mat[, 1],
      Dim2 = h_mat[, 2]
    )
  }

  # Step E: Integrate coordinates with grouping information and return
  # Merge coordinates with original metadata
  meta_subset <- meta_data[, c("SampleID", group_col), drop = FALSE]
  res_df <- merge(coords, meta_subset, by = "SampleID", all.x = TRUE)

  message("Dimensionality reduction analysis complete!")
  return(res_df)
}

#' Core Algorithm for Epigenetic Differential Region Analysis
#'
#' @description Performs differential methylation/epigenetic analysis between a target group and a control group (or the rest of the samples).
#' It calculates level variances, fold changes, and utilizes Fisher's Exact Test on count data to determine statistical significance.
#' Finally, it classifies regions as "Up-regulated", "Down-regulated", or "Not Sig" based on user-defined thresholds.
#'
#' @param raw_mat A data.frame or matrix containing the raw epigenetic data. If the first column is of type character, it is automatically converted to row names.
#' @param meta A data.frame containing metadata. This maps sample groups to their corresponding data columns.
#' @param group_col Character. The column name in \code{meta} representing sample grouping.
#' @param target_group Character. The name of the target (experimental) group.
#' @param control_group Character or \code{NULL}. The name of the control group. If \code{NULL}, a "Target vs Rest" (1 vs All) analysis is performed.
#' @param col_level Character. The column name in \code{meta} indicating the data columns for methylation 'level' (0-1 values).
#' @param col_meth Character. The column name in \code{meta} indicating the data columns for 'meth' sites (counts).
#' @param col_nonmeth Character. The column name in \code{meta} indicating the data columns for 'nonmeth' sites (counts).
#' @param effect_metric Character. The metric used to determine up/down regulation. Options are \code{"logFC"}, \code{"Diff"}, or \code{"Scores"}. Default is \code{"Diff"}.
#' @param effect_th Numeric. The absolute threshold for the selected \code{effect_metric}. Default is 1.
#' @param p_th Numeric. The P-value threshold for significance. Default is 0.05.
#' @param verbose Logical. Whether to print progress information to the console. Default is \code{TRUE}.
#'
#' @return A \code{data.frame} containing the comprehensive differential analysis results, sorted dynamically by P-value (ascending) and effect size (descending).
#'
#' @importFrom stats var fisher.test p.adjust
#' @export
Run_Diffanalysis <- function(raw_mat,
                             meta,
                             group_col,target_group,
                             control_group = NULL,
                             col_level,col_meth,col_nonmeth,
                             effect_metric = "Diff",
                             effect_th = 1,
                             p_th = 0.05,
                             verbose = TRUE) {

  start_time <- Sys.time()

  # --- 0. Preprocessing ---
  raw_mat <- as.data.frame(raw_mat)
  # Set the first column as row names
  if(is.character(raw_mat[[1]])) {
    rownames(raw_mat) <- raw_mat[[1]]
    raw_mat <- raw_mat[, -1, drop = FALSE] # Remove the first column which is already used as row names
  }

  # --- Step A: Extract Group Indices ---
  if(verbose) message("Step 1/4: Extracting group indices and mapping columns...")
  group_vec <- as.character(meta[[group_col]])
  target_idx <- which(group_vec == target_group)

  # Core logic: Determine if it's 1vs1 or 1vsRest
  if (!is.null(control_group)) {
    # 1 vs 1 mode
    control_idx <- which(group_vec == control_group)
  } else {
    # 1 vs Rest mode
    control_group <- "Rest"
    control_idx <- which(group_vec != target_group & !is.na(group_vec) & trimws(group_vec) != "")
  }

  if (length(target_idx) < 2 || length(control_idx) < 2) {
    stop("Need at least 2 samples in both Target and Control groups to calculate variance!")
  }

  # Clean and match column names
  clean_vector <- function(x) { trimws(as.character(x[!is.na(x)])) }

  t_level <- clean_vector(meta[[col_level]][target_idx]); t_level <- t_level[t_level %in% colnames(raw_mat)]
  t_site  <- clean_vector(meta[[col_meth]][target_idx]);  t_site  <- t_site[t_site %in% colnames(raw_mat)]
  t_non   <- clean_vector(meta[[col_nonmeth]][target_idx]); t_non <- t_non[t_non %in% colnames(raw_mat)]

  c_level <- clean_vector(meta[[col_level]][control_idx]); c_level <- c_level[c_level %in% colnames(raw_mat)]
  c_site  <- clean_vector(meta[[col_meth]][control_idx]);  c_site  <- c_site[c_site %in% colnames(raw_mat)]
  c_non   <- clean_vector(meta[[col_nonmeth]][control_idx]); c_non <- c_non[c_non %in% colnames(raw_mat)]

  if(length(t_level) < 2 || length(c_level) < 2) {
    stop("Sample names from metadata do not match matrix column names.")
  }

  # --- Step B: Level Data Processing (Highly Vectorized) ---
  if(verbose) message("Step 2/4: Calculating Level Variances & Fold Changes...")
  t_mat_level <- as.matrix(raw_mat[, t_level, drop=FALSE])
  c_mat_level <- as.matrix(raw_mat[, c_level, drop=FALSE])

  # 1. Strict filtering: The number of non-NAs in both groups must be >= 2
  keep_idx <- (rowSums(!is.na(t_mat_level)) >= 2) & (rowSums(!is.na(c_mat_level)) >= 2)
  t_mat_level <- t_mat_level[keep_idx, , drop=FALSE]
  c_mat_level <- c_mat_level[keep_idx, , drop=FALSE]
  valid_regions <- rownames(t_mat_level)

  # 2. Vectorized calculation
  t_mean <- rowMeans(t_mat_level, na.rm = TRUE)
  c_mean <- rowMeans(c_mat_level, na.rm = TRUE)
  t_var <- apply(t_mat_level, 1, var, na.rm = TRUE)
  c_var <- apply(c_mat_level, 1, var, na.rm = TRUE)

  # 3. Calculate difference formulas
  diff_val <- t_mean - c_mean
  scores <- diff_val / (1 + sqrt(t_var^2 + c_var^2))
  fc <- t_mean / c_mean
  log_fc <- log2((t_mean + 1) / (c_mean + 1))

  # --- Step C: Fisher's Exact Test on Count Data ---
  if(verbose) message(sprintf("Step 3/4: Running Fisher's Exact Test on %d valid sites...", length(valid_regions)))

  t_mat_site <- as.matrix(raw_mat[valid_regions, t_site, drop=FALSE])
  t_mat_non  <- as.matrix(raw_mat[valid_regions, t_non, drop=FALSE])
  c_mat_site <- as.matrix(raw_mat[valid_regions, c_site, drop=FALSE])
  c_mat_non  <- as.matrix(raw_mat[valid_regions, c_non, drop=FALSE])

  # Summarize counts
  t_site_sum <- rowSums(t_mat_site, na.rm=TRUE)
  t_non_sum  <- rowSums(t_mat_non, na.rm=TRUE)
  c_site_sum <- rowSums(c_mat_site, na.rm=TRUE)
  c_non_sum  <- rowSums(c_mat_non, na.rm=TRUE)

  p_values <- rep(1, length(valid_regions))
  n_test <- length(valid_regions)

  for (i in 1:n_test) {
    mat <- matrix(c(t_site_sum[i], t_non_sum[i], c_site_sum[i], c_non_sum[i]), nrow=2)
    if(sum(mat) > 0) {
      try({ p_values[i] <- fisher.test(mat)$p.value }, silent = TRUE)
    }
  }

  # --- Step D: Merge Calculation Results and Return ---
  if(verbose) message("Step 4/4: Generating Final Table...")

  final_df <- data.frame(
    chrdata = valid_regions,
    P.Value = p_values,
    FDR = p.adjust(p_values, method = "fdr"),
    var_group = t_var,
    var_nogroup = c_var,
    Diff = diff_val,
    Scores = scores,
    FC = fc,
    logFC = log_fc,
    stringsAsFactors = FALSE
  )

  # Dynamically mark significance
  final_df$Significance <- "Not Sig"
  effect_vals <- final_df[[effect_metric]]

  final_df$Significance[effect_vals > effect_th & final_df$P.Value < p_th] <- "Up-regulated"
  final_df$Significance[effect_vals < -effect_th & final_df$P.Value < p_th] <- "Down-regulated"

  # Dynamically sort (ascending by P-value, descending by absolute Effect value)
  final_df <- final_df[order(final_df$P.Value, -abs(final_df[[effect_metric]])), ]

  end_time <- Sys.time()
  if(verbose) message(sprintf("Analysis Complete! Time elapsed: %.2f seconds.", as.numeric(difftime(end_time, start_time, units = "secs"))))

  # Return dataframe only
  return(final_df)
}
