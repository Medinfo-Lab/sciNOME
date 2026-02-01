library(dplyr)
library(readxl)
library(sciNOME)


load("data/mouse/Level_Data.RData")
load("data/mouse/Site_Data.RData")
load("data/mouse/List_Data.RData")
load("data/mouse/Epi_Group_Data.RData")



#site----
#CpG
GSE121690_CpG_genetss2k_meth_data <- Read_file_colname(GSE121690_CpG_genetss2k_meth,"site")
GSE121690_CpG_encode_meth_data <- Read_file_colname(GSE121690_CpG_enhancer_meth,"site")

GSE121690_CpG_genetss2k_UNmeth_data <- Read_file_colname(GSE121690_CpG_genetss2k_meth,"nonsite")
GSE121690_CpG_encode_UNmeth_data <- Read_file_colname(GSE121690_CpG_enhancer_meth,"nonsite")

#GpC
GSE121690_GpC_genetss2k_meth_data <- Read_file_colname(GSE121690_GpC_genetss2k_meth,"site")
GSE121690_GpC_encode_meth_data <- Read_file_colname(GSE121690_GpC_enhancer_meth,"site")

GSE121690_GpC_genetss2k_UNmeth_data <- Read_file_colname(GSE121690_GpC_genetss2k_meth,"nonsite")
GSE121690_GpC_encode_UNmeth_data <- Read_file_colname(GSE121690_GpC_enhancer_meth,"nonsite")


group_levels <- levels(factor(GSE121690_CpG_Epis_sample$Developmental_stage))


#CpG----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Site_group_variance_analysis(
    GSE121690_CpG_genetss2k_meth_data,
    GSE121690_CpG_genetss2k_UNmeth_data,
    GSE121690_CpG_Epis_sample,
    "Developmental_stage","site","UNsite",
    group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/meth_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_CpG.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  # write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Site_group_variance_analysis(
    GSE121690_CpG_encode_meth_data,
    GSE121690_CpG_encode_UNmeth_data,
    GSE121690_CpG_Epis_sample,
    "Developmental_stage","site","UNsite",
    group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/meth_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_CpG.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  # write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}

group_levels <- levels(factor(GSE121690_GpC_Epis_sample$Developmental_stage))

#GpC----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Site_group_variance_analysis(
    GSE121690_GpC_genetss2k_meth_data,
    GSE121690_GpC_genetss2k_UNmeth_data,
    GSE121690_GpC_Epis_sample,
    "Developmental_stage","site","UNsite",
    group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/meth_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_GpC.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  # write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Site_group_variance_analysis(
    GSE121690_GpC_encode_meth_data,
    GSE121690_GpC_encode_UNmeth_data,
    GSE121690_GpC_Epis_sample,
    "Developmental_stage","site","UNsite",
    group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/meth_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_GpC.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  # write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}




#level----
#CpG
GSE121690_CpG_genetss2k_methlevel_data <- Read_file_colname(GSE121690_CpG_genetss2k_methlevel,"level")
GSE121690_CpG_encode_methlevel_data <- Read_file_colname(GSE121690_CpG_enhancer_methlevel,"level")

GSE121690_CpG_genetss2k_methlevel_data <- GSE121690_CpG_genetss2k_methlevel_data/100
GSE121690_CpG_encode_methlevel_data <- GSE121690_CpG_encode_methlevel_data/100

#GpC
GSE121690_GpC_genetss2k_methlevel_data <- Read_file_colname(GSE121690_GpC_genetss2k_methlevel,"level")
GSE121690_GpC_encode_methlevel_data <- Read_file_colname(GSE121690_GpC_enhancer_methlevel,"level")

GSE121690_GpC_genetss2k_methlevel_data <- GSE121690_GpC_genetss2k_methlevel_data/100
GSE121690_GpC_encode_methlevel_data <- GSE121690_GpC_encode_methlevel_data/100


group_levels <- levels(factor(GSE121690_CpG_Epis_sample$Developmental_stage))

#CpG----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Level_group_variance_analysis(GSE121690_CpG_genetss2k_methlevel_data,GSE121690_CpG_Epis_sample,
                                                        "Developmental_stage","level",group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/methlevel_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_CpG.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  # write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Level_group_variance_analysis(GSE121690_CpG_encode_methlevel_data,GSE121690_CpG_Epis_sample,
                                                     "Developmental_stage","level",group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/methlevel_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_CpG.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  # write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}

group_levels <- levels(factor(GSE121690_GpC_Epis_sample$Developmental_stage))

#GpC----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Level_group_variance_analysis(GSE121690_GpC_genetss2k_methlevel_data,GSE121690_GpC_Epis_sample,
                                                        "Developmental_stage","level",group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/methlevel_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_GpC.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  # write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Level_group_variance_analysis(GSE121690_GpC_encode_methlevel_data,GSE121690_GpC_Epis_sample,
                                                     "Developmental_stage","level",group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/methlevel_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_GpC.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  # write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}
