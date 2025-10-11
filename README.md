# sciNOME: Integrative Analysis of Nucleosome Occupancy, Genome Methylation and Expression across Single-Cells

# Installation instructions

*sciNOME: A Region-Centric Framework for Integrative Analysis of Single-Cell Epigenomic and Transcriptomic Data.* 

***sciNOME*** is an R package for jointly analyzing transcriptomic, DNA methylation, and chromatin accessibility data. The package is designed to process sequencing data from scNOMe-seq or scNMT-seq experiments. Additionally, ***sciNOME*** converts sequencing data into region-based formats for efficient storage and subsequent analysis. The R package incorporates features such as differential analysis, dimensionality reduction analysis, and differential analysis of horizontal sites.

```R
#install devtools if you don't have it already for easy installation
install.packages("devtools")
library(devtools)
install_github("Medinfo-Lab/sciNOME")
```

If you prefer to build the package by hand, follow these steps:

- Make sure that you have the dependencies from CRAN ("dplyr","reticulate","utils","MOFA","data.table")

- Download and build from source:

```R
git clone git@github.com:Medinfo-Lab/sciNOME.git
R CMD build sciNOME
R CMD INSTALL sciNOME-0.1.4.tar.gz
```

# Usage

**The Workflow:**

![](https://imgur.com/HzTMLit.png)

**The Epigenomic Processing Flow:**

![](https://imgur.com/SQdvvC8.png)

**The Transcriptomic Processing Flow:**

![](https://imgur.com/gVY0VJ0.png)

```R
library(sciNOME)
library(dplyr)
library(data.table)

#First provide a coverage data, bed data and chromosome data
merge_coverage <- list.files(
  coverage_path,
  full.names = TRUE,
  pattern = "\\.cov.gz$"
)

load("data/list.RData")
CPU_cores <- 10

#site methlevel
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  lines <- readLines(merge_coverage[i], warn = FALSE)
  if (length(lines)<1) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }
  
  cov_data <- fread(merge_coverage[i])
  b <- data.frame()
  # Parallel processing of each chromosome using mclapply
  merge_list <- mclapply(chromosome_data, function(chr_tmp) {
    merge_chr <- Cov_to_data(merge_coverage[i], cov_data, chr_tmp, bed_data, suffixname, "methlevel")
    merge_name <- paste0("merge_", chr_tmp)
    return(list(merge_name = merge_chr))
  }, mc.cores = CPU_cores)  # Number of CPU cores used
  
  b <- do.call(rbind, lapply(merge_list, function(x) x$merge_name))
  sample_name <- colnames(b)
  bed_data_paste_methlevel[, sample_name] <- b[, sample_name]
  cat(i,"bed data ",merge_coverage[i],'\n')
  cat("File processing time:", round(Sys.time() - start_time, 1), "秒\n")
}
                             
#meth UNmeth
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  lines <- readLines(merge_coverage[i], warn = FALSE)
  if (length(lines)<1) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }
  
  cov_data <- fread(merge_coverage[i])
  b <- data.frame()
  # Parallel processing of each chromosome using mclapply
  merge_list <- mclapply(chromosome_data, function(chr_tmp) {
    merge_chr <- Cov_to_data(merge_coverage[i], cov_data, chr_tmp, bed_data, suffixname, "meth")
    merge_name <- paste0("merge_", chr_tmp)
    return(list(merge_name = merge_chr))
  }, mc.cores = CPU_cores)  # Number of CPU cores used
  
  b <- do.call(rbind, lapply(merge_list, function(x) x$merge_name))
  sample_name <- colnames(b)
  bed_data_paste_meth[, sample_name] <- b[, sample_name]
  cat(i,"bed data ",merge_coverage[i],'\n')
  cat("File processing time:", round(Sys.time() - start_time, 1), "秒\n")
}
                             
bed_data_paste_methlevel_methleveldata <- Read_file_meth_colname(bed_data_paste_methlevel,"methlevel")
bed_data_paste_meth_methdata <- Read_file_meth_colname(bed_data_paste_meth,"meth")
bed_data_paste_meth_UNmethdata <- Read_file_meth_colname(bed_data_paste_meth,"UNmeth")
```
