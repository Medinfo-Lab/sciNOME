library(dplyr)
library(readxl)
library(sciNOME)


load("data/List_Data.RData")
load("data/Epi_Group.RData")
load("data/Methlevel.RData")
load("data/Meth.RData")



#meth----
#CpG----
GSE121690_CpG_genetss2k_meth_data <- Read_file_meth_colname(GSE121690_CpG_genetss2k_meth,"meth")
GSE121690_CpG_genebody_meth_data <- Read_file_meth_colname(GSE121690_CpG_genebody_meth,"meth")
GSE121690_CpG_encode_meth_data <- Read_file_meth_colname(GSE121690_CpG_encode_meth,"meth")

GSE121690_CpG_genetss2k_UNmeth_data <- Read_file_meth_colname(GSE121690_CpG_genetss2k_meth,"UNmeth")
GSE121690_CpG_genebody_UNmeth_data <- Read_file_meth_colname(GSE121690_CpG_genebody_meth,"UNmeth")
GSE121690_CpG_encode_UNmeth_data <- Read_file_meth_colname(GSE121690_CpG_encode_meth,"UNmeth")



#GpC----
GSE121690_GpC_genetss2k_meth_data <- Read_file_meth_colname(GSE121690_GpC_genetss2k_meth,"meth")
GSE121690_GpC_genebody_meth_data <- Read_file_meth_colname(GSE121690_GpC_genebody_meth,"meth")
GSE121690_GpC_encode_meth_data <- Read_file_meth_colname(GSE121690_GpC_encode_meth,"meth")

GSE121690_GpC_genetss2k_UNmeth_data <- Read_file_meth_colname(GSE121690_GpC_genetss2k_meth,"UNmeth")
GSE121690_GpC_genebody_UNmeth_data <- Read_file_meth_colname(GSE121690_GpC_genebody_meth,"UNmeth")
GSE121690_GpC_encode_UNmeth_data <- Read_file_meth_colname(GSE121690_GpC_encode_meth,"UNmeth")


group_levels <- levels(factor(GSE121690_CpGgroup$Developmental_stage))


#CpG----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Meth_group_variance_analysis(
    GSE121690_CpG_genetss2k_meth_data,
    GSE121690_CpG_genetss2k_UNmeth_data,
    GSE121690_CpGgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/meth_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_CpG.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  # write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("genebody ")
  genebody_DEG_result <- Meth_group_variance_analysis(
    GSE121690_CpG_genebody_meth_data,
    GSE121690_CpG_genebody_UNmeth_data,
    GSE121690_CpGgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  genebody_DEG_data <- merge(genebody_DEG_result, genebody_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genebody_file_suffix <- "data/meth_DEG/genebody/"
  genebody_file_name <- paste0(group_suffix,"_CpG.csv")
  genebody_file_suffix_name <- paste0(genebody_file_suffix,genebody_file_name)

  # write.csv(genebody_DEG_data,genebody_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Meth_group_variance_analysis(
    GSE121690_CpG_encode_meth_data,
    GSE121690_CpG_encode_UNmeth_data,
    GSE121690_CpGgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/meth_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_CpG.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  # write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}



#GpC----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Meth_group_variance_analysis(
    GSE121690_GpC_genetss2k_meth_data,
    GSE121690_GpC_genetss2k_UNmeth_data,
    GSE121690_GpCgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/meth_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_GpC.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("genebody ")
  genebody_DEG_result <- Meth_group_variance_analysis(
    GSE121690_GpC_genebody_meth_data,
    GSE121690_GpC_genebody_UNmeth_data,
    GSE121690_GpCgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  genebody_DEG_data <- merge(genebody_DEG_result, genebody_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genebody_file_suffix <- "data/meth_DEG/genebody/"
  genebody_file_name <- paste0(group_suffix,"_GpC.csv")
  genebody_file_suffix_name <- paste0(genebody_file_suffix,genebody_file_name)

  write.csv(genebody_DEG_data,genebody_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Meth_group_variance_analysis(
    GSE121690_GpC_encode_meth_data,
    GSE121690_GpC_encode_UNmeth_data,
    GSE121690_GpCgroup,
    "Developmental_stage","meth","UNmeth",
    group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/meth_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_GpC.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}




#methlevel----
#CpG----
GSE121690_CpG_genetss2k_methlevel_data <- Read_file_meth_colname(GSE121690_CpG_genetss2k_methlevel,"methlevel")
GSE121690_CpG_genebody_methlevel_data <- Read_file_meth_colname(GSE121690_CpG_genebody_methlevel,"methlevel")
GSE121690_CpG_encode_methlevel_data <- Read_file_meth_colname(GSE121690_CpG_encode_methlevel,"methlevel")

GSE121690_CpG_genetss2k_methlevel_data <- GSE121690_CpG_genetss2k_methlevel_data/100
GSE121690_CpG_genebody_methlevel_data <- GSE121690_CpG_genebody_methlevel_data/100
GSE121690_CpG_encode_methlevel_data <- GSE121690_CpG_encode_methlevel_data/100


#GpC----
GSE121690_GpC_genetss2k_methlevel_data <- Read_file_meth_colname(GSE121690_GpC_genetss2k_methlevel,"methlevel")
GSE121690_GpC_genebody_methlevel_data <- Read_file_meth_colname(GSE121690_GpC_genebody_methlevel,"methlevel")
GSE121690_GpC_encode_methlevel_data <- Read_file_meth_colname(GSE121690_GpC_encode_methlevel,"methlevel")

GSE121690_GpC_genetss2k_methlevel_data <- GSE121690_GpC_genetss2k_methlevel_data/100
GSE121690_GpC_genebody_methlevel_data <- GSE121690_GpC_genebody_methlevel_data/100
GSE121690_GpC_encode_methlevel_data <- GSE121690_GpC_encode_methlevel_data/100


group_levels <- levels(factor(GSE121690_CpGgroup$Developmental_stage))

#CpG----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Methlevel_group_variance_analysis(GSE121690_CpG_genetss2k_methlevel_data,GSE121690_CpGgroup,
                                                            "Developmental_stage","methlevel",group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/methlevel_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_CpG.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("genebody ")
  genebody_DEG_result <- Methlevel_group_variance_analysis(GSE121690_CpG_genebody_methlevel_data,GSE121690_CpGgroup,
                                                           "Developmental_stage","methlevel",group_suffix)
  genebody_DEG_data <- merge(genebody_DEG_result, genebody_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genebody_file_suffix <- "data/methlevel_DEG/genebody/"
  genebody_file_name <- paste0(group_suffix,"_CpG.csv")
  genebody_file_suffix_name <- paste0(genebody_file_suffix,genebody_file_name)

  write.csv(genebody_DEG_data,genebody_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Methlevel_group_variance_analysis(GSE121690_CpG_encode_methlevel_data,GSE121690_CpGgroup,
                                                         "Developmental_stage","methlevel",group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/methlevel_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_CpG.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}


#GpC----
for (i in 1:length(group_levels)) {
  start_time <- Sys.time()
  group_suffix <- group_levels[i]

  cat(group_suffix,"\n")
  cat("genetss2k ")
  genetss2k_DEG_result <- Methlevel_group_variance_analysis(GSE121690_GpC_genetss2k_methlevel_data,GSE121690_GpCgroup,
                                                            "Developmental_stage","methlevel",group_suffix)
  genetss2k_DEG_data <- merge(genetss2k_DEG_result, genetss2k_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genetss2k_file_suffix <- "data/methlevel_DEG/genetss2k/"
  genetss2k_file_name <- paste0(group_suffix,"_GpC.csv")
  genetss2k_file_suffix_name <- paste0(genetss2k_file_suffix,genetss2k_file_name)

  write.csv(genetss2k_DEG_data,genetss2k_file_suffix_name,row.names = F)


  cat("genebody ")
  genebody_DEG_result <- Methlevel_group_variance_analysis(GSE121690_GpC_genebody_methlevel_data,GSE121690_GpCgroup,
                                                           "Developmental_stage","methlevel",group_suffix)
  genebody_DEG_data <- merge(genebody_DEG_result, genebody_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  genebody_file_suffix <- "data/methlevel_DEG/genebody/"
  genebody_file_name <- paste0(group_suffix,"_GpC.csv")
  genebody_file_suffix_name <- paste0(genebody_file_suffix,genebody_file_name)

  write.csv(genebody_DEG_data,genebody_file_suffix_name,row.names = F)


  cat("enhancer  ")
  encode_DEG_result <- Methlevel_group_variance_analysis(GSE121690_GpC_encode_methlevel_data,GSE121690_GpCgroup,
                                                         "Developmental_stage","methlevel",group_suffix)
  encode_DEG_data <- merge(encode_DEG_result, encode_data, by.x = "chrdata", by.y = "chrdata", all.x = TRUE)

  encode_file_suffix <- "data/methlevel_DEG/enhancer/"
  encode_file_name <- paste0(group_suffix,"_GpC.csv")
  encode_file_suffix_name <- paste0(encode_file_suffix,encode_file_name)

  write.csv(encode_DEG_data,encode_file_suffix_name,row.names = F)
  cat("The ",group_suffix," ", as.numeric((Sys.time() - start_time),units = "mins"), "mins\n")
}
