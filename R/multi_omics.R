#' Integrate RNA, DNA Methylation (CpG), and Chromatin Accessibility (GpC) Data
#'
#' @description Merges multiple omics layers based on genomic region annotations and Gene IDs.
#' Uses a single global metadata table to safely map samples across different omics matrices.
#' Optionally filters features based on differential analysis results.
#'
#' @param mode Character. Integration mode: "tri" (all 3), "rna_cpg", "rna_gpc", or "cpg_gpc".
#' @param target_group Character. The biological group/condition to calculate averages for.
#' @param meta_df Data.frame. The global metadata linking all omics samples (e.g., "ALL" sheet).
#' @param group_col Character. The column name in \code{meta_df} containing group info (e.g., "group1").
#' @param region_df Data.frame. Genomic region annotations (must contain chr, start, end, and gene id).
#'
#' @param rna_obj Object. The scRNA-seq object.
#' @param rna_id_col Character. Column in \code{meta_df} matching RNA cell/sample names.
#' @param rna_diff_df Data.frame. Optional. RNA differential analysis results.
#' @param rna_pval_col Character. P-value column in \code{rna_diff_df}.
#' @param rna_logfc_col Character. logFC column in \code{rna_diff_df}.
#' @param rna_pval_th Numeric. P-value threshold for RNA.
#' @param rna_logfc_th Numeric. Absolute logFC threshold for RNA.
#'
#' @param cpg_mat Matrix/Data.frame. CpG level matrix (Rows = Regions, Cols = Samples).
#' @param cpg_id_col Character. Column in \code{meta_df} matching CpG matrix columns.
#' @param cpg_diff_df Data.frame. Optional. CpG DMR analysis results.
#' @param cpg_pval_col Character. P-value column in \code{cpg_diff_df}.
#' @param cpg_diff_col Character. Difference column in \code{cpg_diff_df}.
#' @param cpg_pval_th Numeric. P-value threshold for CpG.
#' @param cpg_diff_th Numeric. Absolute difference threshold for CpG.
#'
#' @param gpc_mat Matrix/Data.frame. GpC level matrix.
#' @param gpc_id_col Character. Column in \code{meta_df} matching GpC matrix columns.
#' @param gpc_diff_df Data.frame. Optional. GpC DMR analysis results.
#' @param gpc_pval_col Character. P-value column in \code{gpc_diff_df}.
#' @param gpc_diff_col Character. Difference column in \code{gpc_diff_df}.
#' @param gpc_pval_th Numeric. P-value threshold for GpC.
#' @param gpc_diff_th Numeric. Absolute difference threshold for GpC.
#'
#' @return A merged \code{data.frame} containing integrated omics levels.
#'
#' @examples
#' # 1. Create mock metadata
#' meta_df <- data.frame(
#'   group1 = c("Tumor", "Tumor", "Normal"),
#'   sample = c("cell1", "cell2", "cell3"),
#'   CpG_level = c("samp1", "samp2", "samp3")
#' )
#'
#' # 2. Create mock region annotation
#' region_df <- data.frame(
#'   chr = c("chr1", "chr2"),
#'   start = c(100, 200),
#'   end = c(150, 250),
#'   gene_id = c("GeneA", "GeneB")
#' )
#'
#' # 3. Create mock RNA object
#' rna_mat <- matrix(
#'   c(10, 12, 2, 5, 6, 1), nrow = 2, byrow = TRUE,
#'   dimnames = list(c("GeneA", "GeneB"), c("cell1", "cell2", "cell3"))
#' )
#' rna_obj <- list(assays = list(RNA = list(data = rna_mat)))
#'
#' # 4. Create mock CpG matrix (rownames match chr:start-end from region_df)
#' cpg_mat <- matrix(
#'   c(0.8, 0.9, 0.1, 0.7, 0.8, 0.2), nrow = 2, byrow = TRUE,
#'   dimnames = list(c("chr1:100-150", "chr2:200-250"), c("samp1", "samp2", "samp3"))
#' )
#'
#' # 5. Run Integration
#' merged_res <- Integrate_MultiOmics(
#'   mode = "rna_cpg",
#'   target_group = "Tumor",
#'   meta_df = meta_df,
#'   group_col = "group1",
#'   region_df = region_df,
#'   rna_obj = rna_obj, rna_id_col = "sample",
#'   cpg_mat = cpg_mat, cpg_id_col = "CpG_level"
#' )
#'
#' head(merged_res)
#'
#' @importFrom dplyr select filter mutate distinct group_by summarise ungroup inner_join left_join rename first pull
#' @importFrom Matrix rowMeans
#' @importFrom utils head
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom tidyselect all_of
#' @export
Integrate_MultiOmics <- function(
    mode = c("tri", "rna_cpg", "rna_gpc", "cpg_gpc"),
    target_group,
    meta_df,
    group_col = "group1",
    region_df,

    # --- RNA parameters ---
    rna_obj = NULL, rna_id_col = "sample",
    rna_diff_df = NULL, rna_pval_col = "p_val", rna_logfc_col = "avg_log2FC",
    rna_pval_th = 0.05, rna_logfc_th = 0.5,

    # --- CpG parameters ---
    cpg_mat = NULL, cpg_id_col = "CpG_level",
    cpg_diff_df = NULL, cpg_pval_col = "P.Value", cpg_diff_col = "Diff",
    cpg_pval_th = 0.05, cpg_diff_th = 0.1,

    # --- GpC parameters ---
    gpc_mat = NULL, gpc_id_col = "GpC_level",
    gpc_diff_df = NULL, gpc_pval_col = "P.Value", gpc_diff_col = "Diff",
    gpc_pval_th = 0.05, gpc_diff_th = 0.1
) {

  mode <- match.arg(mode)
  target_group <- trimws(as.character(target_group))
  if (target_group == "") stop("Please provide a valid target_group.")
  if (!is.data.frame(meta_df)) stop("meta_df must be a data.frame.")

  # --- New: Support logic for "ALL" global grouping ---
  if (toupper(target_group) == "ALL") {
    meta_sub <- meta_df
    if (nrow(meta_sub) == 0) stop("The provided meta_df is empty.")
  } else {
    if (!group_col %in% colnames(meta_df)) stop(sprintf("Group column not found in meta_df: %s", group_col))
    meta_sub <- meta_df[trimws(as.character(meta_df[[group_col]])) == target_group, , drop = FALSE]
    if (nrow(meta_sub) == 0) stop(sprintf("Target group '%s' not found in column %s of meta_df.", target_group, group_col))
  }

  # STEP A: Region annotation processing
  col_names <- colnames(region_df)
  id_col <- grep("^gene_id$|^ensembl_id$|^id$", col_names, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(id_col)) stop("gene_id column not found in region_df.")
  name_col <- grep("^gene_name$|^symbol$|^genename$", col_names, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(name_col)) name_col <- id_col

  region_df$final_gene_id <- region_df[[id_col]]
  region_df$final_gene_name <- region_df[[name_col]]

  if (!"chrdata" %in% colnames(region_df)) {
    region_df$chrdata <- paste0(trimws(region_df$chr), ":", trimws(region_df$start), "-", trimws(region_df$end))
  }
  region_df$chrdata <- trimws(region_df$chrdata)
  region_df$chrdata_nochr <- gsub("^chr", "", region_df$chrdata, ignore.case = TRUE)

  region_map <- region_df %>%
    dplyr::select(chrdata, chrdata_nochr, GeneID = final_gene_id, GeneName = final_gene_name) %>%
    dplyr::mutate(GeneID = toupper(gsub("\\..*|_.*", "", trimws(as.character(GeneID)))),
                  GeneName = trimws(as.character(GeneName))) %>%
    dplyr::filter(GeneID != "", !is.na(GeneID))

  gene_id_to_name <- region_map %>%
    dplyr::select(GeneID, GeneName) %>%
    dplyr::distinct(GeneID, .keep_all = TRUE)

  # STEP B: RNA Data Processing
  rna_res <- NULL
  if (mode %in% c("tri", "rna_cpg", "rna_gpc")) {
    if (is.null(rna_obj)) stop("The selected mode requires rna_obj input.")

    target_rna_cells <- trimws(as.character(meta_sub[[rna_id_col]]))
    target_rna_cells <- target_rna_cells[!is.na(target_rna_cells) & target_rna_cells != "nan"]

    rna_assay <- rna_obj$assays$RNA
    expr_mat <- if (!is.null(rna_assay$data)) rna_assay$data else if (!is.null(rna_assay$filter_counts)) rna_assay$filter_counts else rna_assay$counts

    valid_cells <- intersect(target_rna_cells, colnames(expr_mat))
    if (length(valid_cells) == 0) stop("No intersection between extracted RNA sample names from the master table and column names of the rna_obj matrix!")

    raw_rna_genes <- rownames(expr_mat)
    clean_mat_genes <- toupper(gsub("\\..*|_.*", "", trimws(raw_rna_genes)))

    if (!is.null(rna_diff_df)) {
      rna_diff_df[[rna_pval_col]] <- as.numeric(as.character(rna_diff_df[[rna_pval_col]]))
      rna_diff_df[[rna_logfc_col]] <- as.numeric(as.character(rna_diff_df[[rna_logfc_col]]))
      rna_diff_df <- rna_diff_df[which(rna_diff_df[[rna_pval_col]] < rna_pval_th & abs(rna_diff_df[[rna_logfc_col]]) > rna_logfc_th), ]

      gene_col <- grep("gene|id|ensembl", colnames(rna_diff_df), ignore.case = TRUE, value = TRUE)[1]
      if (is.na(gene_col)) gene_col <- colnames(rna_diff_df)[1]

      diff_genes_raw <- unique(trimws(as.character(rna_diff_df[[gene_col]])))
      diff_genes_clean <- toupper(gsub("\\..*|_.*", "", diff_genes_raw))

      common_clean_genes <- intersect(clean_mat_genes, diff_genes_clean)

      if (length(common_clean_genes) == 0) {
        stop(sprintf("No gene match! Please check the ID format.\nTop 3 genes in expression matrix: %s\nTop 3 genes in differential table: %s",
                     paste(head(raw_rna_genes, 3), collapse = ", "),
                     paste(head(diff_genes_raw, 3), collapse = ", ")))
      }
      valid_rna_genes <- raw_rna_genes[clean_mat_genes %in% common_clean_genes]
    } else {
      valid_rna_genes <- raw_rna_genes
    }

    sub_mat <- expr_mat[valid_rna_genes, valid_cells, drop = FALSE]
    # Take the average (when it's ALL, this is the overall average of all samples)
    rna_vals <- Matrix::rowMeans(sub_mat, na.rm = TRUE)

    final_rna_ids <- toupper(gsub("\\..*|_.*", "", trimws(rownames(sub_mat))))

    rna_res <- data.frame(GeneID = final_rna_ids, RNA_Exp = as.numeric(rna_vals), stringsAsFactors = FALSE) %>%
      dplyr::filter(RNA_Exp > 0) %>%
      dplyr::group_by(GeneID) %>%
      dplyr::summarise(RNA_Exp = sum(RNA_Exp, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  # STEP C: Methylation/Accessibility common processing closure
  process_meth_layer <- function(mat, id_col, type_name, diff_df, pval_col, diff_col, p_th, d_th) {
    target_mat_cols <- trimws(as.character(meta_sub[[id_col]]))
    target_mat_cols <- target_mat_cols[!is.na(target_mat_cols) & target_mat_cols != "nan"]

    if (is.data.frame(mat)) {
      if (is.character(mat[[1]]) || is.factor(mat[[1]])) {
        rnames <- as.character(mat[[1]])
        mat <- as.matrix(mat[, -1, drop = FALSE])
        rownames(mat) <- rnames
      } else {
        mat <- as.matrix(mat)
      }
    }

    valid_cols <- intersect(target_mat_cols, colnames(mat))
    if (length(valid_cols) == 0) stop(sprintf("Matching sample columns not found for %s!", type_name))

    mat_feats <- rownames(mat)
    match_std <- sum(mat_feats %in% region_map$chrdata)
    match_nochr <- sum(mat_feats %in% region_map$chrdata_nochr)
    use_col <- if (match_nochr > match_std) "chrdata_nochr" else "chrdata"

    if (!is.null(diff_df)) {
      diff_df[[pval_col]] <- as.numeric(as.character(diff_df[[pval_col]]))
      diff_df[[diff_col]] <- as.numeric(as.character(diff_df[[diff_col]]))
      diff_df <- diff_df[which(diff_df[[pval_col]] < p_th & abs(diff_df[[diff_col]]) > d_th), ]
      common_feats <- intersect(trimws(as.character(diff_df[[1]])), rownames(mat))
      if (length(common_feats) == 0) stop(sprintf("%s differential regions have no intersection with the matrix.", type_name))
      mat <- mat[common_feats, , drop = FALSE]
    }

    sub_mat <- mat[, valid_cols, drop = FALSE]
    # Take the average (when it's ALL, this is the overall average of all samples)
    avg_vals <- rowMeans(sub_mat, na.rm = TRUE)

    res_df <- data.frame(RegionID = rownames(sub_mat), value = as.numeric(avg_vals), stringsAsFactors = FALSE)
    res_df[[paste0(type_name, "_level")]] <- res_df$value; res_df$value <- NULL

    join_vec <- stats::setNames(use_col, "RegionID")
    res_df <- res_df %>% dplyr::inner_join(region_map, by = join_vec)

    if (use_col == "chrdata") res_df <- res_df %>% dplyr::rename(chrdata = RegionID)
    res_df <- res_df %>% dplyr::select(chrdata, GeneID, !!dplyr::sym(paste0(type_name, "_level")))
    return(res_df)
  }

  # STEP D & E: Execute processing and merging
  cpg_res <- NULL; gpc_res <- NULL
  if (mode %in% c("tri", "rna_cpg", "cpg_gpc")) cpg_res <- process_meth_layer(cpg_mat, cpg_id_col, "CpG", cpg_diff_df, cpg_pval_col, cpg_diff_col, cpg_pval_th, cpg_diff_th)
  if (mode %in% c("tri", "cpg_gpc", "rna_gpc")) gpc_res <- process_meth_layer(gpc_mat, gpc_id_col, "GpC", gpc_diff_df, gpc_pval_col, gpc_diff_col, gpc_pval_th, gpc_diff_th)

  final_region_dict <- region_map %>% dplyr::group_by(GeneID) %>% dplyr::summarise(Associated_Regions = paste(unique(chrdata), collapse = ";")) %>% dplyr::ungroup()

  if (mode == "tri") {
    cpg_agg <- cpg_res %>% dplyr::group_by(GeneID) %>% dplyr::summarise(CpG_level = mean(CpG_level, na.rm = TRUE))
    gpc_agg <- gpc_res %>% dplyr::group_by(GeneID) %>% dplyr::summarise(GpC_level = mean(GpC_level, na.rm = TRUE))
    merged_df <- rna_res %>% dplyr::inner_join(cpg_agg, by = "GeneID") %>% dplyr::inner_join(gpc_agg, by = "GeneID") %>% dplyr::left_join(final_region_dict, by = "GeneID") %>% dplyr::left_join(gene_id_to_name, by = "GeneID") %>% dplyr::select(Associated_Regions, GeneID, GeneName, RNA_Exp, CpG_level, GpC_level)
  } else if (mode == "rna_cpg") {
    cpg_agg <- cpg_res %>% dplyr::group_by(GeneID) %>% dplyr::summarise(CpG_level = mean(CpG_level, na.rm = TRUE))
    merged_df <- rna_res %>% dplyr::inner_join(cpg_agg, by = "GeneID") %>% dplyr::left_join(final_region_dict, by = "GeneID") %>% dplyr::left_join(gene_id_to_name, by = "GeneID") %>% dplyr::select(Associated_Regions, GeneID, GeneName, RNA_Exp, CpG_level)
  } else if (mode == "rna_gpc") {
    gpc_agg <- gpc_res %>% dplyr::group_by(GeneID) %>% dplyr::summarise(GpC_level = mean(GpC_level, na.rm = TRUE))
    merged_df <- rna_res %>% dplyr::inner_join(gpc_agg, by = "GeneID") %>% dplyr::left_join(final_region_dict, by = "GeneID") %>% dplyr::left_join(gene_id_to_name, by = "GeneID") %>% dplyr::select(Associated_Regions, GeneID, GeneName, RNA_Exp, GpC_level)
  } else if (mode == "cpg_gpc") {
    cpg_reg <- cpg_res %>% dplyr::group_by(chrdata) %>% dplyr::summarise(CpG_level = mean(CpG_level, na.rm=TRUE), GeneID = dplyr::first(GeneID))
    gpc_reg <- gpc_res %>% dplyr::group_by(chrdata) %>% dplyr::summarise(GpC_level = mean(GpC_level, na.rm=TRUE), GeneID = dplyr::first(GeneID))
    merged_df <- dplyr::inner_join(cpg_reg, gpc_reg, by = "chrdata") %>% dplyr::rename(Associated_Regions = chrdata, GeneID = GeneID.x) %>% dplyr::left_join(gene_id_to_name, by = "GeneID") %>% dplyr::select(Associated_Regions, GeneID, GeneName, CpG_level, GpC_level)
  }

  if (nrow(merged_df) == 0) stop("Integration failed: No common genes/regions found that successfully match across all three omics!")
  return(merged_df)
}
