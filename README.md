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
R CMD INSTALL sciNOME-0.2.4.tar.gz
```

# Usage

**The Workflow:**

![](https://imgur.com/tQcn3RL.png)

**The Epigenomic Processing Flow:**

![](https://imgur.com/9B3kq6u.png)

**The Transcriptomic Processing Flow:**

![](https://imgur.com/gVY0VJ0.png)

```R

library(sciNOME)
library(dplyr)
library(data.table)
library(parallel)

#First provide a coverage data, bed data, chromosome data and group data
merge_coverage <- list.files(
  coverage_path,
  full.names = TRUE,
  pattern = "\\.cov.gz$"
)

load("data/List_Data.RData")
load("data/Epi_Group_Data.RData")

bed_data_paste_Site <- sprintf("%s:%s-%s", bed_data$chr, bed_data$start, bed_data$end)
bed_data_paste_Level <- sprintf("%s:%s-%s", bed_data$chr, bed_data$start, bed_data$end)

bed_data_paste_Site <- as.data.frame(bed_data_paste_Site)
colnames(bed_data_paste_Site) <- "chrdata"
bed_data_paste_Level <- as.data.frame(bed_data_paste_Level)
colnames(bed_data_paste_Level) <- "chrdata"

# example chromosome data
chr_data <- c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
              "chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4",
              "chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")

#Level
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  if (file.size(merge_CpG[i]) == 0) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }

  cov_data <- fread(merge_CpG[i]) 
  result_df <- tryCatch({
    Coverage_to_data(
      cov_file = merge_coverage[i],     # File path used for extracting sample names
      cov_file_data = cov_data,         # Data content after reading
      region_data = bed_data,    		# Region data
      chr_data = chr_data,    			# All chromosome vectors requiring processing
      suffixname_data = "suffixname",   # Suffix
      method_type = "Level"             # The first letter must be capitalized, consistent with the function definition.
    )
  }, error = function(e) {
    message("Error in file processing: ", merge_CpG[i], " - ", e$message)
    return(NULL)
  })

  if (is.null(result_df)) next
  new_cols <- colnames(result_df)
  if (nrow(result_df) == nrow(bed_data_paste_Level)) {
    bed_data_paste_Level[new_cols] <- result_df
  } else {
    warning(paste("Warning: Document", merge_coverage[i], "The number of processed rows does not match the target table; skip assignment."))
  }

  cat(i, "Bed Data Level Processed:", merge_coverage[i], '\n')
  cat("Time:", round(difftime(Sys.time(), start_time, units = "secs"), 1), "second\n")
}

#Site
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  if (file.size(merge_CpG[i]) == 0) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }

  cov_data <- fread(merge_CpG[i]) 
  result_df <- tryCatch({
    Coverage_to_data(
      cov_file = merge_coverage[i],     # File path used for extracting sample names
      cov_file_data = cov_data,         # Data content after reading
      region_data = bed_data,    		# Region data
      chr_data = chr_data,    			# All chromosome vectors requiring processing
      suffixname_data = "suffixname",   # Suffix
      method_type = "Site"              # The first letter must be capitalized, consistent with the function definition.
    )
  }, error = function(e) {
    message("Error in file processing: ", merge_coverage[i], " - ", e$message)
    return(NULL)
  })

  if (is.null(result_df)) next
  new_cols <- colnames(result_df)
  if (nrow(result_df) == nrow(bed_data_paste_Site)) {
    bed_data_paste_Site[new_cols] <- result_df
  } else {
    warning(paste("Warning: Document", merge_coverage[i], "The number of processed rows does not match the target table; skip assignment."))
  }

  cat(i, "Bed Data Site Processed:", merge_coverage[i], '\n')
  cat("Time:", round(difftime(Sys.time(), start_time, units = "secs"), 1), "second\n")
}

bed_data_paste_leveldata <- Read_file_colname(bed_data_paste_Level,"level")
bed_data_paste_sitedata <- Read_file_colname(bed_data_paste_Site,"site")
bed_data_paste_nonsitedata <- Read_file_colname(bed_data_paste_Site,"nonsite")

group_levels <- levels(factor(group$Developmental_stage))

for(i in 1:length(group_levels)){
  group_suffix <- group_levels[i]
  cat(group_suffix,"\n")
  level_DEG_data <- Level_group_variance_analysis(bed_data_paste_leveldata, Group_data,
                                                  "Developmental_stage", "level", group_suffix)
  level_DEG_data <- merge(level_DEG_data, list_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)
  site_DEG_data <- Site_group_variance_analysis(bed_data_paste_sitedata, bed_data_paste_nonsitedata, Group_data,
                                                "Developmental_stage", "site","UNsite", group_suffix)
  site_DEG_data <- merge(site_DEG_data, list_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)
}
```
