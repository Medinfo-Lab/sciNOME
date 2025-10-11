#' Read meth, UNmeth, and meth_level data, the corresponding column name is \\.name$
#'
#' @param file_tmp DNA methylation or chromatin accessibility matrix
#' @param string_data sample type information, such as meth, UNmeth, methlevel
#'
#' @return sample data after selecting columns
#' @export
#'
#' @examples
Read_file_meth_colname <- function(file_tmp,string_data){
  file_data <- file_tmp
  string_tmp <- string_data
  string_col <- paste0("\\.",string_tmp)
  string_col2 <- paste0(string_col,"$")

  column_indices <- grep(string_col2, names(file_data), value = TRUE)
  file_data_meth <- file_data[,column_indices]
  # file_data_meth <- cbind(rownames(file_data),file_data_meth)
  rownames(file_data_meth) <- rownames(file_tmp)
  return(file_data_meth)
}

#' Data quality control
#'
#' @param file_tmp DNA methylation or chromatin accessibility matrix
#' @param flag choose to perform quality control by column or row
#'
#' @return data after quality control
#' @export
#'
#' @examples
QC_meth_unmeth <- function(file_tmp,flag="col"){
  file_data <- file_tmp

  file_meth_name <- grep("\\.meth$", colnames(file_data), value = TRUE)
  file_unmeth_name <- grep("\\.UNmeth$", colnames(file_data), value = TRUE)

  file_meth_data <- file_data[,file_meth_name]
  file_unmeth_data <- file_data[,file_unmeth_name]

  file_meth_unmeth <- file_meth_data+file_unmeth_data
  rownames(file_meth_unmeth) <- rownames(file_data)
  file_name <- sub("\\.meth$","",file_meth_name)
  colnames(file_meth_unmeth) <- paste0(file_name,".meth_UNmeth")

  if (flag=="row") {
    file_data_num_count <- data.frame()
    file_data_num_count <- rownames(file_meth_unmeth)
    file_data_num_count <- as.data.frame(file_data_num_count)
    colnames(file_data_num_count)[1] <- "chr"

    file_meth_data_sum <- apply(file_meth_data, MARGIN = 1, sum)
    file_data_num_count$chr_meth_UNmeth_SUM <- apply(file_meth_unmeth, MARGIN = 1, sum)
    file_data_num_count$chr_methlevel <- file_meth_data_sum/file_data_num_count$chr_meth_UNmeth_SUM
    file_data_num_count_sorted <- arrange(file_data_num_count, -chr_meth_UNmeth_SUM)
    return(file_data_num_count_sorted)
  }else if (flag=="col") {
    file_data_num_count <- data.frame()
    file_data_num_count <- colnames(file_meth_unmeth)
    file_data_num_count <- as.data.frame(file_data_num_count)
    colnames(file_data_num_count)[1] <- "sample"

    file_meth_data_sum <- apply(file_meth_data, MARGIN = 2, sum)
    file_data_num_count$sample_meth_UNmeth_SUM <- apply(file_meth_unmeth, MARGIN = 2, sum)
    file_data_num_count$sample_methlevel <- file_meth_data_sum/file_data_num_count$sample_meth_UNmeth_SUM
    file_data_num_count_sorted <- arrange(file_data_num_count, -sample_meth_UNmeth_SUM)
    return(file_data_num_count_sorted)
  }
}

#' Cov files to text files
#'
#' @param cov_file site file name
#' @param cov_file_data site data
#' @param chr_tmp chromosome data
#' @param region chromosomal physical fragment
#' @param suffixname file extension
#' @param datatmp select the meth or meth-level method
#'
#' @return DNA methylation or chromatin accessibility matrix
#' @export
#'
#' @examples
Cov_to_data <- function(cov_file, cov_file_data, chr_tmp, region, suffixname, datatmp) {
  region_chr <- region %>% filter(chr %in% chr_tmp)
  cov_data <- cov_file_data
  colnames(cov_data) <- c("chr", "start", "end", "methlevel", "meth", "UNmeth")
  cov_data_chr <- cov_data %>% filter(chr %in% chr_tmp)
  cov_data_chr_order <- cov_data_chr %>%
    group_by(chr) %>%
    arrange(start,.by_group = T)


  if(datatmp=="meth"){
    # 生成列名
    cov_file_name <- basename(cov_file)
    split_result <- strsplit(cov_file_name, suffixname)[[1]][1]
    col_names <- paste0(split_result, c(".meth", ".UNmeth"))

    result_data <- data.frame()

    # 使用向量化操作代替嵌套循环
    for (i in 1:length(chr_tmp)) {
      cov_data_chr_order_chr <- cov_data_chr_order %>%
        filter(chr %in% chr_tmp[i])
      region_chr_chr <- region_chr %>%
        filter(chr %in% chr_tmp[i])

      region_chr_paste <- sprintf("%s:%s-%s", region_chr_chr$chr, region_chr_chr$start, region_chr_chr$end)
      region_chr_paste <- as.data.frame(region_chr_paste)
      colnames(region_chr_paste) <- "chr"

      df_meth <- data.frame(matrix(nrow = nrow(region_chr_paste), ncol = 2))
      rownames(df_meth) <- region_chr_paste$chr
      colnames(df_meth) <- col_names

      for (j in 1:nrow(region_chr_chr)) {
        # 使用向量化操作替代双重循环
        # chr_match <- (cov_data_chr_order_chr$chr == region_chr_chr$chr[j])
        start_match <- (cov_data_chr_order_chr$start >= region_chr_chr$start[j])
        end_match <- (cov_data_chr_order_chr$start <= region_chr_chr$end[j])
        valid <- start_match & end_match

        df_meth[j, 1] <- sum(cov_data_chr_order_chr$meth[valid])
        df_meth[j, 2] <- sum(cov_data_chr_order_chr$UNmeth[valid])
      }
      result_data <- rbind(result_data,df_meth)
    }
    return(result_data)
  }

  if(datatmp=="methlevel"){
    # 生成列名
    cov_file_name <- basename(cov_file)
    split_result <- strsplit(cov_file_name, suffixname)[[1]][1]
    col_names <- paste0(split_result, c(".site", ".methlevel"))

    result_data <- data.frame()

    # 使用向量化操作代替嵌套循环
    for (i in 1:length(chr_tmp)) {
      cov_data_chr_order_chr <- cov_data_chr_order %>%
        filter(chr %in% chr_tmp[i])
      region_chr_chr <- region_chr %>%
        filter(chr %in% chr_tmp[i])

      region_chr_paste <- sprintf("%s:%s-%s", region_chr_chr$chr, region_chr_chr$start, region_chr_chr$end)
      region_chr_paste <- as.data.frame(region_chr_paste)
      colnames(region_chr_paste) <- "chr"

      df_methlevel <- data.frame(matrix(nrow = nrow(region_chr_paste), ncol = 2))
      rownames(df_methlevel) <- region_chr_paste$chr
      colnames(df_methlevel) <- col_names

      for (j in 1:nrow(region_chr_chr)) {
        # 过滤出当前 region_chr 相关的 cov_data_chr
        filtered_cov_data <- cov_data_chr_order_chr[cov_data_chr_order_chr$start >= region_chr_chr$start[j] & cov_data_chr_order_chr$start <= region_chr_chr$end[j], ]

        # 计算 site_count 和 methlevel_sum
        site_count <- nrow(filtered_cov_data)
        methlevel_sum <- ifelse(site_count > 0, mean(filtered_cov_data$methlevel), NA)

        # 赋值到 df_meth
        df_methlevel[j, ] <- c(site_count, methlevel_sum)
      }
      result_data <- rbind(result_data,df_methlevel)
    }
    return(result_data)
  }
}

#' Filtering of site information
#'
#' @param files Site information file path
#' @param save_path Save folder path
#'
#' @return nothing
#' @export
#'
#' @examples
Cov_QC <- function(files,save_path){
  for (i in 1:length(files)) {
    name <- basename(files[i])
    result <- sub("\\.cov$", "", name)
    result2 <- paste0(result, ".csv")

    file_data <- read.table(files[i])
    colnames(file_data) <- c("chr","start","end","methlevel","meth","UNmeth")

    file_data_filter <- file_data %>%
      filter(methlevel>=70 | methlevel<=30)
    file_data_save_path <- paste0(save_path, result2)

    # write.csv(file_data_filter,file_data_save_path,row.names = F)
    cat(files[i],"Success\n")
    fwrite(file_data_filter,file_data_save_path,sep = ',')
  }
}

#' Differential level data processing
#'
#' @param methlevel_data Chromosome-level locus data
#' @param group_data group data
#' @param suff group information
#' @param group_suff grouping field
#' @param methlevel_group_suff methlevel grouping field
#'
#' @return methlevel results
#' @export
#'
#' @examples
Methlevel_diff <- function(methlevel_data,group_data,group_suff,methlevel_group_suff,suff){
  group_sample <- group_data %>%
    filter(.data[[group_suff]] == suff)
  no_group_sample <- group_data %>%
    filter(.data[[group_suff]] != suff)

  methlevel_data_meth_target <- methlevel_data[,group_sample[[methlevel_group_suff]]]
  methlevel_data_meth_control <- methlevel_data[,no_group_sample[[methlevel_group_suff]]]

  methlevel_diff_data <- data.table()

  methlevel_diff_data_target <- rowMeans(methlevel_data_meth_target,na.rm = T)
  methlevel_diff_data_control <- rowMeans(methlevel_data_meth_control,na.rm = T)

  methlevel_diff_data$rate_diff <- methlevel_diff_data_target-methlevel_diff_data_control
  methlevel_diff_data_frame <- as.data.frame(methlevel_diff_data)
  methlevel_diff_data_frame$logFC <- log2(methlevel_diff_data_target/methlevel_diff_data_control)
  return(methlevel_diff_data_frame)
}

#' Difference level chart
#'
#' @param CpG_DEG_file_data Methylation differences
#' @param GpC_DEG_file_data Open difference information
#'
#' @return CpG and GpC differences data
#' @export
#'
#' @examples
DEG_diff_plot <- function(CpG_DEG_file_data,GpC_DEG_file_data){
  # CpG_DEG_file_data <- read.csv(CpG_DEG_files)
  # GpC_DEG_file_data <- read.csv(GpC_DEG_files)

  CpG_DEG_file_data_clean_df_nan_na <- CpG_DEG_file_data %>%
    filter(!is.nan(rate_diff) & !is.na(rate_diff))


  GpC_DEG_file_data_clean_df_nan_na <- GpC_DEG_file_data %>%
    filter(!is.nan(rate_diff) & !is.na(rate_diff))

  CpG_GpC_DEG_union <- intersect(GpC_DEG_file_data_clean_df_nan_na$chr,
                                 CpG_DEG_file_data_clean_df_nan_na$chr)
  # CpG_GpC_DEG_geneid <- intersect(GpC_DEG_file_data_clean_df_nan_na$genename,
  #                                 CpG_DEG_file_data_clean_df_nan_na$genename)
  CpG_GpC_DEG_genename <- GpC_DEG_file_data_clean_df_nan_na %>%
    filter(chr %in% CpG_GpC_DEG_union)

  CpG_DEG_file_data_choose <- subset(CpG_DEG_file_data_clean_df_nan_na, chr %in% CpG_GpC_DEG_union)
  GpC_DEG_file_data_choose <- subset(GpC_DEG_file_data_clean_df_nan_na, chr %in% CpG_GpC_DEG_union)


  gene_diff <- data.frame(chr=CpG_DEG_file_data_choose$chr,
                          CpG_p.value=CpG_DEG_file_data_choose$P.value,
                          CpG_adj_p.value=CpG_DEG_file_data_choose$adj_P.value_fdr,
                          GpC_p.value=GpC_DEG_file_data_choose$P.value,
                          GpC_adj_p.value=GpC_DEG_file_data_choose$adj_P.value_fdr,
                          CpG_rate_diff=CpG_DEG_file_data_choose$rate_diff,
                          GpC_rate_diff=GpC_DEG_file_data_choose$rate_diff,
                          geneid=CpG_GpC_DEG_genename$geneid,
                          genename=CpG_GpC_DEG_genename$genename
  )
  return(gene_diff)
}

