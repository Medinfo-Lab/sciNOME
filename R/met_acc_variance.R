#' Distinguish chrXX:xxx-xxx as chrXX, xxx, xxx
#'
#' @param file_tmp chromosome physical segment data, chrXX:xxx-xxx format or chrXX, xxx, xxx format
#' @param method chromosome fragment processing method
#'
#' @return Processed physical fragment data, chrXX, xxx, xxx format or chrXX:xxx-xxx format
#' @export
#'
#' @examples
Chr_region_process <- function(file_tmp,method){
  if(method=="split"){
    chr_data <- file_tmp
    chr_data_frame <- data.frame(region=chr_data)
    chr_split <- chr_data_frame %>%
      tidyr::separate(region, into = c("chr", "start", "end"), sep = ":|-")
    chr_split$start <- as.numeric(chr_split$start)
    chr_split$end <- as.numeric(chr_split$end)
    return(chr_split)
  }else if(method=="paste"){
    chr_data <- file_tmp
    chr_data_paste <- sprintf("%s:%d-%d",chr_data$chr,chr_data$start,chr_data$end)
    chr_data_paste_frame <- as.data.frame(chr_data_paste)
    colnames(chr_data_paste_frame) <- "chrdata"
    return(chr_data_paste_frame)
  }
}


#' Methlevel group variance analysis
#'
#' @param methlevel_data DNA methylation or chromatin accessibility methlevel data
#' @param group_data DNA methylation or chromatin accessibility group data
#' @param group_suff grouping field
#' @param methlevel_group_suff methlevel grouping field
#' @param suff grouping suffix name
#'
#' @return methlevel Difference analysis results
#' @export
#'
#' @examples
Methlevel_group_variance_analysis <- function(methlevel_data,group_data,group_suff,methlevel_group_suff,suff){
  group_methlevel <- group_data %>%
    filter(.data[[group_suff]] %in% suff)
  nogroup_methlevel <- group_data %>%
    filter(!.data[[group_suff]] %in% suff)

  group_methlevel_data <- methlevel_data[,group_methlevel[[methlevel_group_suff]]]
  nogroup_methlevel_data <- methlevel_data[,nogroup_methlevel[[methlevel_group_suff]]]

  group_methlevel_data_clean <- group_methlevel_data[!apply(is.na(group_methlevel_data), 1, all), ]
  nogroup_methlevel_data_clean <- nogroup_methlevel_data[!apply(is.na(nogroup_methlevel_data), 1, all), ]

  rowname_ids <- intersect(rownames(group_methlevel_data_clean), rownames(nogroup_methlevel_data_clean))

  group_methlevel_data_clean_choose <- group_methlevel_data_clean[rowname_ids,]
  group_methlevel_data_clean_choose <- as.matrix(group_methlevel_data_clean_choose)

  nogroup_methlevel_data_clean_choose <- nogroup_methlevel_data_clean[rowname_ids,]
  nogroup_methlevel_data_clean_choose <- as.matrix(nogroup_methlevel_data_clean_choose)

  group_methlevel_data_clean_choose_rowmean <- rowMeans(group_methlevel_data_clean_choose,na.rm = T)
  nogroup_methlevel_data_clean_choose_rowmean <- rowMeans(nogroup_methlevel_data_clean_choose,na.rm = T)


  DEG_data <- data.frame(chr=rowname_ids)

  for (i in 1:length(rowname_ids)) {
    result <- wilcox.test(group_methlevel_data_clean_choose[i,],
                          nogroup_methlevel_data_clean_choose[i,]
    )
    DEG_data$p.value[i] <- result$p.value
    DEG_data$logFC[i] <- log2(group_methlevel_data_clean_choose_rowmean[i]/nogroup_methlevel_data_clean_choose_rowmean[i])
  }
  DEG_data$adj_p.value_fdr <- p.adjust(DEG_data$p.value, method = "fdr")
  return(DEG_data)
}


#' Meth group variance analysis
#'
#' @param meth_data DNA methylation or chromatin accessibility meth data
#' @param UNmeth_data DNA methylation or chromatin accessibility UNmeth data
#' @param group_data DNA methylation or chromatin accessibility group data
#' @param group_suff grouping field
#' @param title_suff meth or UNmeth grouping field
#' @param suff group suffix name
#'
#' @return meth differential data
#' @export
#'
#' @examples
Meth_group_variance_analysis <- function(meth_data,UNmeth_data,group_data,group_suff,title_suff,suff){
  group_sample <- group_data %>%
    filter(.data[[group_suff]] %in% suff)
  no_group_sample <- group_data %>%
    filter(!.data[[group_suff]] %in% suff)


  meth_data_meth_target <- meth_data[,group_sample[[title_suff]]]
  UNmeth_data_UNmeth_target <- UNmeth_data[,group_sample[[title_suff]]]

  meth_data_meth_control <- meth_data[,no_group_sample[[title_suff]]]
  UNmeth_data_UNmeth_control <- UNmeth_data[,no_group_sample[[title_suff]]]


  meth_data_meth_target_sum <- rowSums(meth_data_meth_target)
  meth_data_UNmeth_target_sum <- rowSums(UNmeth_data_UNmeth_target)

  UNmeth_data_meth_control_sum <- rowSums(meth_data_meth_control)
  UNmeth_data_UNmeth_control_sum <- rowSums(UNmeth_data_UNmeth_control)


  fisher_data <- data.frame(chr=rownames(meth_data_meth_target))

  for (i in 1:nrow(fisher_data)) {
    x = matrix(c(meth_data_meth_target_sum[i],
                 meth_data_UNmeth_target_sum[i],
                 UNmeth_data_meth_control_sum[i],
                 UNmeth_data_UNmeth_control_sum[i]), nrow=2, ncol=2)
    p <- fisher.test(x)
    fisher_data[i,1] <- p$p.value
  }
  fisher_data$adj_p.value_fdr <- p.adjust(fisher_data$p.value,method = "fdr")
  fisher_data$log10_adj_p.value_fdr <- -log10(fisher_data$adj_p.value_fdr)
  return(fisher_data)
}

