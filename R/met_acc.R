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


#' Coverage files to text files
#'
#' @param cov_file site file name
#' @param cov_file_data site data
#' @param chr_data chromosome data
#' @param region chromosomal physical fragment
#' @param suffixname file extension
#' @param datatmp select the meth or meth-level method
#'
#' @return DNA methylation or chromatin accessibility matrix
#' @export
#'
#' @examples
Coverage_to_data <- function(cov_file, cov_file_data, region_data, chr_data, suffixname_data=".cov.gz", type) {
  region_chr <- region_data %>% filter(chr %in% chr_data)
  cov_data <- cov_file_data
  colnames(cov_data) <- c("chr", "start", "end", "methlevel", "meth", "UNmeth")
  cov_data_chr <- cov_data %>% filter(chr %in% chr_data)
  cov_data_chr_order <- cov_data_chr %>%
    group_by(chr) %>%
    arrange(start,.by_group = T)
  suffixname <- paste0("\\",suffixname_data)

  if(type=="meth"){
    # 生成列名
    cov_file_name <- basename(cov_file)
    split_result <- strsplit(cov_file_name, suffixname)[[1]][1]
    col_names <- paste0(split_result, c(".meth", ".UNmeth"))

    result_data <- data.frame()

    # 使用向量化操作代替嵌套循环
    for (i in 1:length(chr_data)) {
      cov_data_chr_order_chr <- cov_data_chr_order %>%
        filter(chr %in% chr_data[i])
      region_chr_chr <- region_chr %>%
        filter(chr %in% chr_data[i])

      region_chr_paste <- sprintf("%s:%s-%s", region_chr_chr$chr, region_chr_chr$start, region_chr_chr$end)
      region_chr_paste <- as.data.frame(region_chr_paste)
      colnames(region_chr_paste) <- "chrdata"

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

  if(type=="methlevel"){
    # 生成列名
    cov_file_name <- basename(cov_file)
    split_result <- strsplit(cov_file_name, suffixname)[[1]][1]
    col_names <- paste0(split_result, c(".site", ".methlevel"))

    result_data <- data.frame()

    # 使用向量化操作代替嵌套循环
    for (i in 1:length(chr_data)) {
      cov_data_chr_order_chr <- cov_data_chr_order %>%
        filter(chr %in% chr_data[i])
      region_chr_chr <- region_chr %>%
        filter(chr %in% chr_data[i])

      region_chr_paste <- sprintf("%s:%s-%s", region_chr_chr$chr, region_chr_chr$start, region_chr_chr$end)
      region_chr_paste <- as.data.frame(region_chr_paste)
      colnames(region_chr_paste) <- "chrdata"

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

#' Methlevel quality control
#'
#' @param file_data
#' @param sample_cutoff
#'
#' @returns
#' @export
#'
#' @examples
QC_methlevel_data <- function(file_data,sample_cutoff=0.8){
  samples_keep <- colMeans(is.na(file_data)) < sample_cutoff
  CpG_data_clean <- file_data[, samples_keep]

  cat("After filtering the samples, the matrix dimensions:", dim(CpG_data_clean), "\n")
  cat("Number of excluded samples:", ncol(file_data) - ncol(CpG_data_clean), "\n")

  return(CpG_data_clean)
}


#' Difference level chart
#'
#' @param CpG_DEG_file_data Methylation differences information
#' @param GpC_DEG_file_data Chromitan difference information
#'
#' @return CpG and GpC differences data
#' @export
#'
#' @examples
DEG_diff_process <- function(CpG_DEG_file_data,GpC_DEG_file_data){
  # CpG_DEG_file_data <- read.csv(CpG_DEG_files)
  # GpC_DEG_file_data <- read.csv(GpC_DEG_files)

  CpG_GpC_chrdata <- intersect(CpG_DEG_file_data$chrdata,GpC_DEG_file_data$chrdata)

  CpG_DEG_data_choose_union <- CpG_DEG_file_data %>%
    filter(chrdata %in% CpG_GpC_chrdata)

  GpC_DEG_data_choose_union <- GpC_DEG_file_data %>%
    filter(chrdata %in% CpG_GpC_chrdata)

  colnames(CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(CpG_DEG_data_choose_union))
  colnames(GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(GpC_DEG_data_choose_union))

  DEG_data_union <- cbind(CpG_DEG_data_choose_union,GpC_DEG_data_choose_union)
  return(DEG_data_union)
}


#' Dimensionality reduction processing
#'
#' @param file_data Enter level data
#' @param method Dimension reduction methods
#'
#' @returns Dimension reduction results
#' @export
#'
#' @examples
DR_process <- function(file_data,method="MDS"){
  if (method=="MDS") {
    cor_matrix <- cor(file_data, use = "pairwise.complete.obs", method = "spearman")
    # 将相关性转换为距离 (1 - correlation)
    dist_matrix <- as.dist(1 - cor_matrix)
    # 如果 dist_matrix 里有 NA，需要把 NA 替换为最大距离
    dist_matrix[is.na(dist_matrix)] <- 1
    # MDS 降维 (cmdscale)
    mds_points <- cmdscale(dist_matrix, k = 2, eig = TRUE)
    return(mds_points)
  }else if (method=="PCA") {
    pca_res <- prcomp(t(file_data),
                      center = T,
                      scale. = F,
                      rank. = 30)
    return(pca_res)
  }else if (method=="UMAP") {
    data_umap_res <- uwot::umap(file_data,n_neighbors=100,min_dist=0.5)
    return(data_umap_res)
  }
}


#' Calculate the silhouette coefficient
#'
#' @param file_data Dimensionality reduction data
#' @param file_type Dimensionality-Reduced data types
#'
#' @returns Mean Silhouette coefficient value
#' @export
#'
#' @examples
SC_process <- function(file_data,file_type="UMAP"){
  if (file_type=="UMAP") {
    umap_coords <- file_data[,c(1,2)]
    umap_groups <- file_data$Stage
    umap_dist <- dist(umap_coords)
    sil_umap <- silhouette(as.numeric(factor(umap_groups)), umap_dist)
    mean_sil_umap <- mean(sil_umap[, 3])
    cat("UMAP Mean Silhouette Coefficient Value:", round(mean_sil_umap, 3), "\n")
    return(mean_sil_umap)
  }else if(file_type=="PCA"){
    pca_coords <- file_data[,c(1,2)]
    pca_groups <- file_data$Stage
    pca_dist <- dist(pca_coords)
    sil_pca <- silhouette(as.numeric(factor(pca_groups)), pca_dist)
    mean_sil_pca <- mean(sil_pca[, 3])
    cat("PCA Mean Silhouette Coefficient Value:", round(mean_sil_pca, 3), "\n")
    return(mean_sil_pca)
  }else if(file_type=="MDS"){
    mds_coords <- file_data[,c(1,2)]
    mds_groups <- file_data$Stage
    mds_dist <- dist(mds_coords)
    sil_mds <- silhouette(as.numeric(factor(mds_groups)), mds_dist)
    mean_sil_mds <- mean(sil_mds[, 3])
    cat("MDS Mean Silhouette Coefficient Value:", round(mean_sil_mds, 3), "\n")
    return(mean_sil_mds)
  }
}
