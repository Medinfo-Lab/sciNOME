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
#' @param methlevel_group_suff level grouping field
#' @param suff grouping suffix name
#'
#' @return methlevel Difference analysis results
#' @export
#'
#' @examples
Level_group_variance_analysis <- function(methlevel_data,group_data,
                                          group_suff="group",
                                          methlevel_group_suff="level",
                                          suff){
  # methlevel_data <- GSE121690_CpG_genetss2k_methlevel_data
  # group_data <- GSE121690_CpG_Epis_sample
  # group_suff="Developmental_stage"
  # methlevel_group_suff="level"
  # suff="E4.5"

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


  DEG_data <- data.frame(chrdata=rowname_ids)

  for (i in 1:length(rowname_ids)) {
    # result <- wilcox.test(group_methlevel_data_clean_choose[i,],
    #                       nogroup_methlevel_data_clean_choose[i,]
    # )
    # DEG_data$p.value[i] <- result$p.value
    DEG_data$var_group[i] <- var(group_methlevel_data_clean_choose[i,],na.rm = T)
    DEG_data$var_nogroup[i] <- var(nogroup_methlevel_data_clean_choose[i,],na.rm = T)
    DEG_data$Diff[i] <- group_methlevel_data_clean_choose_rowmean[i]-nogroup_methlevel_data_clean_choose_rowmean[i]
    DEG_data$Scores[i] <- DEG_data$Diff[i]/(1+sqrt((DEG_data$var_group[i])^2+(DEG_data$var_nogroup[i])^2))
    DEG_data$FC[i] <- group_methlevel_data_clean_choose_rowmean[i]/nogroup_methlevel_data_clean_choose_rowmean[i]
    DEG_data$logFC[i] <- log2((group_methlevel_data_clean_choose_rowmean[i]+1)/
                                (nogroup_methlevel_data_clean_choose_rowmean[i]+1))
  }
  # DEG_data$adj_p.value_fdr <- p.adjust(DEG_data$p.value, method = "fdr")
  return(DEG_data)
}


#' Meth group variance analysis
#'
#' @param meth_data DNA methylation or chromatin accessibility meth data
#' @param UNmeth_data DNA methylation or chromatin accessibility UNmeth data
#' @param group_data DNA methylation or chromatin accessibility group data
#' @param group_suff grouping field
#' @param meth_title_suff site grouping field
#' @param UNmeth_title_suff UNsite grouping field
#' @param suff group suffix name
#'
#' @return meth differential data
#' @export
#'
#' @examples
Site_group_variance_analysis <- function(meth_data,UNmeth_data,group_data,
                                         group_suff="group",
                                         meth_title_suff="site",
                                         UNmeth_title_suff="UNsite",
                                         suff){
  # meth_data <- GSE121690_CpG_genetss2k_meth_data
  # UNmeth_data <- GSE121690_CpG_genetss2k_UNmeth_data
  # group_data <- GSE121690_CpG_Epis_sample
  # group_suff="Developmental_stage"
  # meth_title_suff="site"
  # UNmeth_title_suff="UNsite"
  # suff="E4.5"

  group_sample <- group_data %>%
    filter(.data[[group_suff]] %in% suff)
  no_group_sample <- group_data %>%
    filter(!.data[[group_suff]] %in% suff)


  meth_data_meth_target <- meth_data[,group_sample[[meth_title_suff]]]
  UNmeth_data_UNmeth_target <- UNmeth_data[,group_sample[[UNmeth_title_suff]]]

  meth_data_meth_control <- meth_data[,no_group_sample[[meth_title_suff]]]
  UNmeth_data_UNmeth_control <- UNmeth_data[,no_group_sample[[UNmeth_title_suff]]]

  meth_data_meth_target_sum <- rowSums(meth_data_meth_target)
  meth_data_UNmeth_target_sum <- rowSums(UNmeth_data_UNmeth_target)

  UNmeth_data_meth_control_sum <- rowSums(meth_data_meth_control)
  UNmeth_data_UNmeth_control_sum <- rowSums(UNmeth_data_UNmeth_control)

  fisher_data <- data.frame(chrdata=rownames(meth_data_meth_target))

  for (i in 1:nrow(fisher_data)) {
    x = matrix(c(meth_data_meth_target_sum[i],
                 meth_data_UNmeth_target_sum[i],
                 UNmeth_data_meth_control_sum[i],
                 UNmeth_data_UNmeth_control_sum[i]), nrow=2, ncol=2)
    p <- fisher.test(x)
    fisher_data$p.value[i] <- p$p.value
  }
  fisher_data$adj_p.value_fdr <- p.adjust(fisher_data$p.value,method = "fdr")
  return(fisher_data)
}

