#' count the number of NAs in the ranks
#'
#' @param file_tmp Methlevel data of DNA methylation or chromatin accessibility
#' @param flag Choose to count by column or row
#'
#' @return statistical information
#' @export
#'
#' @examples
NA_count_order <- function(file_tmp,flag="row"){
  file_data_meth_level <- file_tmp

  if(flag=="row"){
    file_data_NA_count <- rowMeans(is.na(file_data_meth_level))
    file_data_NA_count <- as.data.frame(file_data_NA_count)
    colnames(file_data_NA_count)[1] <- "NA_percentage"

    file_data_NA_count$NA_count <- rowSums(is.na(file_data_meth_level))

    file_data_NA_count_sorted <- arrange(file_data_NA_count, NA_percentage)
    return(file_data_NA_count_sorted)
  }else if(flag=="col"){
    file_data_NA_count <- colMeans(is.na(file_data_meth_level))
    file_data_NA_count <- as.data.frame(file_data_NA_count)
    colnames(file_data_NA_count)[1] <- "NA_percentage"

    file_data_NA_count$NA_count <- colSums(is.na(file_data_meth_level))

    file_data_NA_count_sorted <- arrange(file_data_NA_count, NA_percentage)
    return(file_data_NA_count_sorted)
  }
}


#' The NA value is populated by the specified value and by the average value
#'
#' @param file_tmp Methlevel data of DNA methylation or chromatin accessibility
#' @param flag Choose to count by column or row
#' @param num Enter the number to fill
#'
#' @return After filling the matrix
#' @export
#'
#' @examples
NA_padding_num_mean <- function(file_tmp,num=0.5,flag="row"){
  if(flag=="row"){
    file_data_meth <- file_tmp
    file_data_meth2 <- file_tmp

    #The selected data is filled in with the average NA value according to the corresponding chr
    file_data_meth[is.na(file_data_meth)] <- num

    # rownames(file_data_meth) <- file_data_meth$chr
    # file_data_meth <- file_data_meth[,-1]

    file_data_mean_meth <- data.frame()
    file_data_mean_meth <- rownames(file_data_meth)
    file_data_mean_meth <- as.data.frame(file_data_mean_meth)
    colnames(file_data_mean_meth)[1] <- "chr"

    file_data_mean_meth$mean <- apply(file_data_meth,1,mean)

    # rownames(file_data_meth2) <- file_data_meth2$chr
    # file_data_meth2 <- file_data_meth2[,-1]

    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    cat("\n正在处理数据...\n")

    for (i in 1:nrow(file_data_meth2)) {
      for (j in 1:ncol(file_data_meth2)) {
        if (is.na(file_data_meth2[i,j])) {
          file_data_meth2[i,j] <- file_data_mean_meth$mean[i]
        }
      }
      setTxtProgressBar(pb, i / nrow(file_data_meth2) * 100)
    }
    return(file_data_meth2)
    close(pb)

  }else if(flag=="col"){
    file_data_meth <- file_tmp
    file_data_meth2 <- file_tmp

    #The selected data is filled in with the average NA value according to the corresponding chr
    file_data_meth[is.na(file_data_meth)] <- num

    # rownames(file_data_meth) <- file_data_meth$chr
    # file_data_meth <- file_data_meth[,-1]

    file_data_mean_meth <- data.frame()
    file_data_mean_meth <- colnames(file_data_meth)
    file_data_mean_meth <- as.data.frame(file_data_mean_meth)
    colnames(file_data_mean_meth)[1] <- "sample"

    file_data_mean_meth$mean <- apply(file_data_meth,2,mean)

    # rownames(file_data_meth2) <- file_data_meth2$chr
    # file_data_meth2 <- file_data_meth2[,-1]

    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    cat("\n正在处理数据...\n")

    for (i in 1:ncol(file_data_meth2)) {
      for (j in 1:nrow(file_data_meth2)) {
        if (is.na(file_data_meth2[j,i])) {
          file_data_meth2[j,i] <- file_data_mean_meth$mean[i]
        }
      }
      setTxtProgressBar(pb, i / ncol(file_data_meth2) * 100)
    }
    return(file_data_meth2)
    close(pb)
  }
}


#' The NA value is populated by the average
#'
#' @param file_tmp Methlevel data of DNA methylation or chromatin accessibility
#'
#' @return After filling the matrix
#' @export
#'
#' @examples
NA_padding_mean <- function(file_tmp,flag="row"){
  if(flag=="row"){

    file_data_meth_tmp <- file_tmp
    # rownames(file_data_meth_tmp) <- file_data_meth_tmp[,1]
    # file_data_meth_tmp <- file_data_meth_tmp[,-1]
    file_data_meth_mean <- rowMeans(file_data_meth_tmp,na.rm=T)

    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    cat("\n正在处理数据...\n")

    for (i in 1:nrow(file_data_meth_tmp)) {
      for (j in 1:ncol(file_data_meth_tmp)) {
        file_data_meth_tmp[i,j] <- ifelse(is.na(file_data_meth_tmp[i,j]), file_data_meth_mean[i], file_data_meth_tmp[i,j])
      }
      setTxtProgressBar(pb, i / nrow(file_data_meth_tmp) * 100)
    }
    return(file_data_meth_tmp)
    close(pb)

  }else if(flag=="col"){
    file_data_meth_tmp <- file_tmp
    # rownames(file_data_meth_tmp) <- file_data_meth_tmp[,1]
    # file_data_meth_tmp <- file_data_meth_tmp[,-1]
    file_data_meth_mean <- colMeans(file_data_meth_tmp,na.rm=T)

    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    cat("\n正在处理数据...\n")

    for (i in 1:ncol(file_data_meth_tmp)) {
      for (j in 1:nrow(file_data_meth_tmp)) {
        file_data_meth_tmp[j,i] <- ifelse(is.na(file_data_meth_tmp[j,i]), file_data_meth_mean[i], file_data_meth_tmp[j,i])
      }
      setTxtProgressBar(pb, i / ncol(file_data_meth_tmp) * 100)
    }
    return(file_data_meth_tmp)
    close(pb)
  }
}
