library(ggtext)
library(ggrepel)
library(fpc)
library(cluster)
library(clusterProfiler)
library(readxl)
library(dplyr)
library(ggplot2)
library(sciNOME)
library(data.table)
library(org.Mm.eg.db)
library(ggVennDiagram)
library(VennDiagram)
library(enrichplot)
library(tidydr)
library(uwot)
library(impute)
library(MOFA2)
library(mixOmics)

load("data/GSE121690_Epi_GeneTSS2k_Dimensionality.RData")


#NA process----
GSE121690_CpG_methlevel_na_count <- NA_count_order(GSE121690_CpG_methlevel)
GSE121690_GpC_methlevel_na_count <- NA_count_order(GSE121690_GpC_methlevel)

GSE121690_CpG_methlevel_na_count_nona <- GSE121690_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.37)

GSE121690_CpG_methlevel_choose <- GSE121690_CpG_methlevel[rownames(GSE121690_CpG_methlevel_na_count_nona),]

GSE121690_GpC_methlevel_na_count_nona <- GSE121690_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.33)

GSE121690_GpC_methlevel_choose <- GSE121690_GpC_methlevel[rownames(GSE121690_GpC_methlevel_na_count_nona),]



#CpG umap----
sample_cutoff <- 0.8
samples_keep <- colMeans(is.na(GSE121690_CpG_methlevel_choose)) < sample_cutoff
CpG_data_clean <- GSE121690_CpG_methlevel_choose[, samples_keep]

cat("After filtering the samples, the matrix dimensions:", dim(CpG_data_clean), "\n")
cat("Number of excluded samples:", ncol(GSE121690_CpG_methlevel_choose) - ncol(CpG_data_clean), "\n")

CpG_data_clean <- as.matrix(CpG_data_clean)
GSE121690_CpG_methlevel_choose_imputed_res <- impute.knn(CpG_data_clean, k = 10)
GSE121690_CpG_methlevel_choose_imputed_res_data <- GSE121690_CpG_methlevel_choose_imputed_res$data
GSE121650_CpG_methlevel_sample_group <- GSE121690_CpG_Epis_sample %>%
  filter(Sample_title %in% colnames(CpG_data_clean))


#不插补 + 相关性距离 + MDS (PCoA) —— 最适合高稀疏数据
# 计算细胞间的相关性矩阵
CpG_cor_matrix <- cor(CpG_data_clean, use = "pairwise.complete.obs", method = "spearman")
# 将相关性转换为距离 (1 - correlation)
CpG_dist_matrix <- as.dist(1 - CpG_cor_matrix)
# 如果 dist_matrix 里有 NA，需要把 NA 替换为最大距离
CpG_dist_matrix[is.na(CpG_dist_matrix)] <- 1
# MDS 降维 (cmdscale)
CpG_mds_points <- cmdscale(CpG_dist_matrix, k = 2, eig = TRUE)

CpG_mds_data <- data.frame(
  Dim1 = CpG_mds_points$points[,1],
  Dim2 = CpG_mds_points$points[,2],
  Stage = GSE121650_CpG_methlevel_sample_group$Developmental_stage
)
# write.csv(CpG_mds_data,"plot/genetss2k/UMAP_DATA/genetss2k_CpG_mds.csv")
# CpG_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_mds.csv")

ggplot(CpG_mds_data, aes(x=Dim1, y=Dim2, color=Stage)) +
  geom_point(CpG_mds_data,mapping=aes(Dim1,Dim2,color=Stage),size = 3.5,alpha = 0.85)+
  theme_classic()+
  labs(x = "Dim1",y = "Dim2",title = "GeneTSS2k CpG MDS")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# "#E69F00", "#56B4E9", "#009E73", "#CC79A7"

#pca
GSE121690_CpG_methlevel_choose_mean <- NA_padding_mean(CpG_data_clean)

pca_res <- prcomp(t(GSE121690_CpG_methlevel_choose_mean),
                  center = T,
                  scale. = T,
                  rank. = 30)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
pca_df <- as.data.frame(pca_res$x)
pca_df$Stage <- GSE121650_CpG_methlevel_sample_group$Developmental_stage # 加上分组信息
# write.csv(pca_df,"plot/genetss2k/UMAP_DATA/genetss2k_CpG_pca.csv")

ggplot(pca_df)+
  geom_point(pca_df,mapping=aes(PC1,PC2,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = pc1_lab,y = pc2_lab,title = "GeneTSS2k CpG PCA")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# theme_dr()


#umap
# GSE121690_CpG_methlevel_choose_imputed_res_data_umap_res <- umap(t(GSE121690_CpG_methlevel_choose_imputed_res_data),
#                                                                  n_neighbors = 100, min_dist = 0.5)
GSE121690_CpG_methlevel_choose_imputed_res_data_umap_res <- umap(pca_df,
                                                                 n_neighbors = 100, min_dist = 1)

plot_df <- data.frame(
  UMAP1 = GSE121690_CpG_methlevel_choose_imputed_res_data_umap_res[,1],
  UMAP2 = GSE121690_CpG_methlevel_choose_imputed_res_data_umap_res[,2],
  Stage = GSE121650_CpG_methlevel_sample_group$Developmental_stage
)
# write.csv(plot_df,"plot/genetss2k/UMAP_DATA/genetss2k_CpG_umap.csv")

ggplot(plot_df)+
  geom_point(plot_df,mapping=aes(UMAP1,UMAP2,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "UMAP1",y = "UMAP2",title = "GeneTSS2k CpG UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
  # theme_dr()



#CpG mofa----
GSE121690_CpG_methlevel_choose_mean_matrix <- as.matrix(CpG_data_clean)
mofa_list <- list()
mofa_list$file1 <- GSE121690_CpG_methlevel_choose_mean_matrix

GSE121690_CpG_mofa <- create_mofa(mofa_list,group = GSE121650_CpG_methlevel_sample_group$Developmental_stage)

ModelOptions <- get_default_model_options(GSE121690_CpG_mofa)
TrainOptions <- get_default_training_options(GSE121690_CpG_mofa)
DataOptions <- get_default_data_options(GSE121690_CpG_mofa)

# TrainOptions$convergence_mode <- "medium"
# DataOptions$scale_views <- TRUE
# DataOptions$center_groups <- TRUE

GSE121690_CpG_mofa <- prepare_mofa(GSE121690_CpG_mofa,
                                   model_options = ModelOptions,
                                   training_options = TrainOptions,
                                   data_options = DataOptions)

GSE121690_CpG_mofa <- run_mofa(GSE121690_CpG_mofa)

# GSE121690_CpG_mofa <- run_umap(GSE121690_CpG_mofa)
GSE121690_CpG_mofa <- run_umap(GSE121690_CpG_mofa,
                               n_neighbors = 100, min_dist = 0.5)

plot_dimred(GSE121690_CpG_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())


GSE121690_CpG_umap_data <- GSE121690_CpG_mofa@dim_red$UMAP
GSE121690_CpG_umap_data <- cbind(GSE121690_CpG_umap_data,GSE121690_CpG_mofa@samples_metadata$group)
colnames(GSE121690_CpG_umap_data)[4] <- "group"
rownames(GSE121690_CpG_umap_data) <- 1:nrow(GSE121690_CpG_umap_data)

# write.csv(GSE121690_CpG_umap_data,"plot/genetss2k/UMAP_DATA/genetss2k_CpG_mofaumap.csv",row.names = F)
# GSE121690_CpGgene_umap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_mofaumap.csv")

ggplot(GSE121690_CpG_umap_data)+
  geom_point(GSE121690_CpG_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 4,alpha = 0.85)+
  theme_bw()+
  labs(x = "UMAP1",y = "UMAP2",title = "GeneTSS2k CpG MOFA UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# theme_dr()
# theme_dr(xlength = 0.2,
#          ylength = 0.2,
#          arrow = arrow(length = unit(0.2, "inches"),type = "closed"))
# "#E69F00", "#56B4E9", "#009E73", "#CC79A7"




#GpC umap----
sample_cutoff <- 0.8
samples_keep <- colMeans(is.na(GSE121690_GpC_methlevel_choose)) < sample_cutoff
GpC_data_clean <- GSE121690_GpC_methlevel_choose[, samples_keep]

cat("After filtering the samples, the matrix dimensions:", dim(GpC_data_clean), "\n")
cat("Number of excluded samples:", ncol(GSE121690_GpC_methlevel_choose) - ncol(GpC_data_clean), "\n")

GpC_data_clean <- as.matrix(GpC_data_clean)
GSE121690_GpC_methlevel_choose_imputed_res <- impute.knn(GpC_data_clean, k = 10)
GSE121690_GpC_methlevel_choose_imputed_res_data <- GSE121690_GpC_methlevel_choose_imputed_res$data
GSE121650_GpC_methlevel_sample_group <- GSE121690_GpC_Epis_sample %>%
  filter(Sample_title %in% colnames(GpC_data_clean))


#不插补 + 相关性距离 + MDS (PCoA) —— 最适合高稀疏数据
# 计算细胞间的相关性矩阵
GpC_cor_matrix <- cor(GpC_data_clean, use = "pairwise.complete.obs", method = "spearman")
# 将相关性转换为距离 (1 - correlation)
GpC_dist_matrix <- as.dist(1 - GpC_cor_matrix)
# 如果 dist_matrix 里有 NA，需要把 NA 替换为最大距离
GpC_dist_matrix[is.na(GpC_dist_matrix)] <- 1
# MDS 降维 (cmdscale)
GpC_mds_points <- cmdscale(GpC_dist_matrix, k = 2, eig = TRUE)

GpC_mds_data <- data.frame(
  Dim1 = GpC_mds_points$points[,1],
  Dim2 = GpC_mds_points$points[,2],
  Stage = GSE121650_GpC_methlevel_sample_group$Developmental_stage
)
# write.csv(GpC_mds_data,"plot/genetss2k/UMAP_DATA/genetss2k_GpC_mds.csv")
# GpC_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_mds.csv")

ggplot(GpC_mds_data, aes(x=Dim1, y=Dim2, color=Stage)) +
  geom_point(GpC_mds_data,mapping=aes(Dim1,Dim2,color=Stage),size = 3.5,alpha = 0.85)+
  theme_classic()+
  labs(x = "Dim1",y = "Dim2",title = "GeneTSS2k GpC MDS")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# "#E69F00", "#56B4E9", "#009E73", "#CC79A7"

#pca
GSE121690_GpC_methlevel_choose_mean <- NA_padding_mean(GpC_data_clean)

pca_res <- prcomp(t(GSE121690_GpC_methlevel_choose_mean),
                  center = T,
                  scale. = T,
                  rank. = 30)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
pca_df <- as.data.frame(pca_res$x)
pca_df$Stage <- GSE121650_GpC_methlevel_sample_group$Developmental_stage # 加上分组信息
# write.csv(pca_df,"plot/genetss2k/UMAP_DATA/genetss2k_GpC_pca.csv")

ggplot(pca_df)+
  geom_point(pca_df,mapping=aes(PC1,PC2,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = pc1_lab,y = pc2_lab,title = "GeneTSS2k GpC PCA")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# theme_dr()


#umap
# GSE121690_GpC_methlevel_choose_imputed_res_data_umap_res <- umap(t(GSE121690_GpC_methlevel_choose_imputed_res_data),
#                                                                  n_neighbors = 100, min_dist = 0.5,metric = "cosine")
GSE121690_GpC_methlevel_choose_imputed_res_data_umap_res <- umap(pca_df,n_neighbors = 100, min_dist = 0.5)

plot_df <- data.frame(
  UMAP1 = GSE121690_GpC_methlevel_choose_imputed_res_data_umap_res[,1],
  UMAP2 = GSE121690_GpC_methlevel_choose_imputed_res_data_umap_res[,2],
  Stage = GSE121650_GpC_methlevel_sample_group$Developmental_stage
)
# write.csv(plot_df,"plot/genetss2k/UMAP_DATA/genetss2k_GpC_umap.csv")


ggplot(plot_df)+
  geom_point(plot_df,mapping=aes(UMAP1,UMAP2,color=Stage),size = 4,alpha = 0.9)+
  theme_bw()+
  labs(x = "UMAP1",y = "UMAP2",title = "GeneTSS2k GpC UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# theme_dr()


#GpC mofa----
GSE121690_GpC_methlevel_choose_mean_matrix <- as.matrix(GpC_data_clean)
mofa_list <- list()
mofa_list$file1 <- GSE121690_GpC_methlevel_choose_mean_matrix

GSE121690_GpC_mofa <- create_mofa(mofa_list,group = GSE121650_GpC_methlevel_sample_group$Developmental_stage)

ModelOptions <- get_default_model_options(GSE121690_GpC_mofa)
TrainOptions <- get_default_training_options(GSE121690_GpC_mofa)
DataOptions <- get_default_data_options(GSE121690_GpC_mofa)

# TrainOptions$convergence_mode <- "slow"
# DataOptions$scale_views <- TRUE
# DataOptions$center_groups <- TRUE

GSE121690_GpC_mofa <- prepare_mofa(GSE121690_GpC_mofa,
                                   model_options = ModelOptions,
                                   training_options = TrainOptions,
                                   data_options = DataOptions)

GSE121690_GpC_mofa <- run_mofa(GSE121690_GpC_mofa)

# GSE121690_GpC_mofa <- run_umap(GSE121690_GpC_mofa,min_dist = 1,spread = 2)
GSE121690_GpC_mofa <- run_umap(GSE121690_GpC_mofa,
                               n_neighbors = 100, min_dist = 0.5)

plot_dimred(GSE121690_GpC_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

GSE121690_GpC_umap_data <- GSE121690_GpC_mofa@dim_red$UMAP
GSE121690_GpC_umap_data <- cbind(GSE121690_GpC_umap_data,GSE121690_GpC_mofa@samples_metadata$group)
colnames(GSE121690_GpC_umap_data)[4] <- "group"
rownames(GSE121690_GpC_umap_data) <- 1:nrow(GSE121690_GpC_umap_data)

# write.csv(GSE121690_GpC_umap_data,"plot/genetss2k/UMAP_DATA/genetss2k_GpC_mofaumap.csv",row.names = F)
# GSE121690_GpCgene_umap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_mofaumap.csv")

ggplot(GSE121690_GpC_umap_data)+
  geom_point(GSE121690_GpC_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 4,alpha = 0.85)+
  theme_bw()+
  labs(x = "UMAP1",y = "UMAP2",title = "GeneTSS2k GpC MOFA UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
        )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))

# "#E69F00", "#56B4E9", "#009E73", "#CC79A7"
# load("GSE121690_Epi_GeneTSS2k_Dimensionality.RData")


#Silhouette Score----
CpG_pca_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_pca.csv",row.names = 1)
CpG_umap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_umap.csv",row.names = 1)
CpG_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_mds.csv",row.names = 1)
CpG_mofaumap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_mofaumap.csv",row.names = 1)

GpC_pca_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_pca.csv",row.names = 1)
GpC_umap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_umap.csv",row.names = 1)
GpC_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_mds.csv",row.names = 1)
GpC_mofaumap_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_mofaumap.csv",row.names = 1)


pca_data <- GpC_pca_data
umap_data <- GpC_umap_data
mds_data <- GpC_mds_data
mofa_umap_data <- GpC_mofaumap_data

#pca
pca_coords <- pca_data[,c(1,2)]
pca_groups <- pca_data$Stage
pca_dist <- dist(pca_coords)
sil_pca <- silhouette(as.numeric(factor(pca_groups)), pca_dist)
mean_sil_pca <- mean(sil_pca[, 3])
cat("PCA平均轮廓系数:", round(mean_sil_pca, 3), "\n")


#umap
umap_coords <- umap_data[,c(1,2)]
umap_groups <- umap_data$Stage
umap_dist <- dist(umap_coords)
sil_umap <- silhouette(as.numeric(factor(umap_groups)), umap_dist)
mean_sil_umap <- mean(sil_umap[, 3])
cat("UMAP平均轮廓系数:", round(mean_sil_umap, 3), "\n")


#mds
mds_coords <- mds_data[,c(1,2)]
mds_groups <- mds_data$Stage
mds_dist <- dist(mds_coords)
sil_mds <- silhouette(as.numeric(factor(mds_groups)), mds_dist)
mean_sil_mds <- mean(sil_mds[, 3])
cat("MDS平均轮廓系数:", round(mean_sil_mds, 3), "\n")


#mofa umap
mofaumap_coords <- mofa_umap_data[,c(1,2)]
mofaumap_groups <- mofa_umap_data$group
mofaumap_dist <- dist(mofaumap_coords)
sil_mofaumap <- silhouette(as.numeric(factor(mofaumap_groups)), mofaumap_dist)
mean_sil_mofaumap <- mean(sil_mofaumap[, 3])
cat("UMAP平均轮廓系数:", round(mean_sil_mofaumap, 3), "\n")



PCA_SC_data <- read_xlsx("SC.xlsx",sheet = "PCA")
UMAP_SC_data <- read_xlsx("SC.xlsx",sheet = "UMAP")
MDS_SC_data <- read_xlsx("SC.xlsx",sheet = "MDS")
MOFAUMAP_SC_data <- read_xlsx("SC.xlsx",sheet = "MOFA")


# 绘制柱状图并美化
ggplot(PCA_SC_data, aes(x = type, y = SC.value, fill = type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8) +  # 柱状图基础设置
  geom_text(aes(label = round(SC.value, 3)), vjust = -0.5, size = 4, color = "black") +  # 添加数值标签
  scale_fill_brewer(palette = "Set2") +  # 使用柔和配色
  labs(
    title = "Silhouette Coefficient PCA",  # 标题
    x = "",                                    # x = "Type",                                    # x轴标签
    y = "Silhouette Coefficient Value"            # y轴标签
    # caption = "Data Source: Silhouette_Coefficient.xlsx"  # 数据来源
  ) +
  theme_minimal() +  # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题 = 16),  # 标题居中加粗
    axis.title = element_text(size = 15),             # 坐标轴标题大小
    axis.text = element_text(size = 12,angle = 15),              # 坐标轴文字大小
    legend.position = "none",                         # 隐藏图例（类别已在x轴显示）
    panel.grid.major = element_blank(),               # 隐藏主要网格线
    # panel.grid.minor = element_blank(),                # 隐藏次要网格线
    axis.ticks.y = element_line(color = "black")
  ) +
  ylim(min(PCA_SC_data$SC.value) - 0.1, max(PCA_SC_data$SC.value) + 0.1)  # 扩展y轴范围以显示数值标签

ggplot(UMAP_SC_data, aes(x = type, y = SC.value, fill = type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8) +  # 柱状图基础设置
  geom_text(aes(label = round(SC.value, 3)), vjust = -0.5, size = 4, color = "black") +  # 添加数值标签
  scale_fill_brewer(palette = "Set2") +  # 使用柔和配色
  labs(
    title = "Silhouette Coefficient UMAP",  # 标题
    x = "",                                    # x = "Type",                                    # x轴标签
    y = "Silhouette Coefficient Value"            # y轴标签
    # caption = "Data Source: Silhouette_Coefficient.xlsx"  # 数据来源
  ) +
  theme_minimal() +  # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题 = 16),  # 标题居中加粗
    axis.title = element_text(size = 15),             # 坐标轴标题大小
    axis.text = element_text(size = 12,angle = 15),              # 坐标轴文字大小
    legend.position = "none",                         # 隐藏图例（类别已在x轴显示）
    panel.grid.major = element_blank(),               # 隐藏主要网格线
    # panel.grid.minor = element_blank(),                # 隐藏次要网格线
    axis.ticks.y = element_line(color = "black")
  ) +
  ylim(min(UMAP_SC_data$SC.value) - 0.1, max(UMAP_SC_data$SC.value) + 0.1)  # 扩展y轴范围以显示数值标签

MDS_SC_data <- MDS_SC_data[c(-2,-1),]
ggplot(MDS_SC_data, aes(x = type, y = SC.value, fill = type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8) +  # 柱状图基础设置
  geom_text(aes(label = round(SC.value, 3)), vjust = -0.5, size = 4, color = "black") +  # 添加数值标签
  scale_fill_brewer(palette = "Set2") +  # 使用柔和配色
  labs(
    title = "Silhouette Coefficient MDS",  # 标题
    x = "",                                    # x = "Type",                                    # x轴标签
    y = "Silhouette Coefficient Value"            # y轴标签
    # caption = "Data Source: Silhouette_Coefficient.xlsx"  # 数据来源
  ) +
  theme_minimal() +  # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题 = 16),  # 标题居中加粗
    axis.title = element_text(size = 15),             # 坐标轴标题大小
    axis.text = element_text(size = 12,angle = 15),              # 坐标轴文字大小
    legend.position = "none",                         # 隐藏图例（类别已在x轴显示）
    panel.grid.major = element_blank(),               # 隐藏主要网格线
    # panel.grid.minor = element_blank(),                # 隐藏次要网格线
    axis.ticks.y = element_line(color = "black")
  ) +
  ylim(min(MDS_SC_data$SC.value) - 0.1, max(MDS_SC_data$SC.value) + 0.1)  # 扩展y轴范围以显示数值标签

ggplot(MOFAUMAP_SC_data, aes(x = type, y = SC.value, fill = type)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8) +  # 柱状图基础设置
  geom_text(aes(label = round(SC.value, 3)), vjust = -0.5, size = 4, color = "black") +  # 添加数值标签
  scale_fill_brewer(palette = "Set2") +  # 使用柔和配色
  labs(
    title = "Silhouette Coefficient MOFA",  # 标题
    x = "",                                    # x = "Type",                                    # x轴标签
    y = "Silhouette Coefficient Value"            # y轴标签
    # caption = "Data Source: Silhouette_Coefficient.xlsx"  # 数据来源
  ) +
  theme_minimal() +  # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 标题 = 16),  # 标题居中加粗
    axis.title = element_text(size = 15),             # 坐标轴标题大小
    axis.text = element_text(size = 12,angle = 15),              # 坐标轴文字大小
    legend.position = "none",                         # 隐藏图例（类别已在x轴显示）
    panel.grid.major = element_blank(),               # 隐藏主要网格线
    # panel.grid.minor = element_blank(),                # 隐藏次要网格线
    axis.ticks.y = element_line(color = "black")
  ) +
  ylim(min(MOFAUMAP_SC_data$SC.value) - 0.1, max(MOFAUMAP_SC_data$SC.value) + 0.1)  # 扩展y轴范围以显示数值标签



#CpG GpC mean meth group----
#CpG
GSE121690_CpG_methlevel_choose_colmean <- colMeans(CpG_data_clean,na.rm = T)
GSE121690_CpG_methlevel_choose_colmean <- as.data.frame(GSE121690_CpG_methlevel_choose_colmean)
GSE121690_CpG_methlevel_choose_colmean$group <- GSE121650_CpG_methlevel_sample_group$Developmental_stage
GSE121690_CpG_methlevel_choose_colmean$group <- factor(GSE121690_CpG_methlevel_choose_colmean$group)

GSE121690_CpG_methlevel_choose_colmean_mean_data <- GSE121690_CpG_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(GSE121690_CpG_methlevel_choose_colmean, na.rm = TRUE))

ggplot(GSE121690_CpG_methlevel_choose_colmean,
       aes(x = group, y = GSE121690_CpG_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.9,scale = "width",aes(linetype = NA))+
  geom_segment(data = GSE121690_CpG_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.05,
                   xend = as.numeric(group) + 0.05,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_classic()+
  labs(y="Percentage (%)",
       x="",
       title = "GeneTSS2k CpG LevelData Means")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")+
  scale_fill_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                    values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))


#GpC
GSE121690_GpC_methlevel_choose_colmean <- colMeans(GpC_data_clean,na.rm = T)
GSE121690_GpC_methlevel_choose_colmean <- as.data.frame(GSE121690_GpC_methlevel_choose_colmean)
GSE121690_GpC_methlevel_choose_colmean$group <- GSE121650_GpC_methlevel_sample_group$Developmental_stage
GSE121690_GpC_methlevel_choose_colmean$group <- factor(GSE121690_GpC_methlevel_choose_colmean$group)

GSE121690_GpC_methlevel_choose_colmean_mean_data <- GSE121690_GpC_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(GSE121690_GpC_methlevel_choose_colmean, na.rm = TRUE))

ggplot(GSE121690_GpC_methlevel_choose_colmean,
       aes(x = group, y = GSE121690_GpC_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.9,scale = "width",aes(linetype = NA))+
  geom_segment(data = GSE121690_GpC_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.05,
                   xend = as.numeric(group) + 0.05,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_classic()+
  labs(y="Percentage (%)",
       x="",
       title = "GeneTSS2k GpC LevelData Means")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")+
  scale_fill_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                    values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))



#DEG data process----
CpG_meth_DEG_files <- list.files("data/meth_DEG/genetss2k/",full.names = T,pattern = "\\CpG.csv*")
GpC_meth_DEG_files <- list.files("data/meth_DEG/genetss2k/",full.names = T,pattern = "\\GpC.csv*")

CpG_methlevel_DEG_files <- list.files("data/methlevel_DEG/genetss2k/",full.names = T,pattern = "\\CpG.csv*")
GpC_methlevel_DEG_files <- list.files("data/methlevel_DEG/genetss2k/",full.names = T,pattern = "\\GpC.csv*")

#4.5
CpG_meth_DEG_data <- read.csv(CpG_meth_DEG_files[1])
GpC_meth_DEG_data <- read.csv(GpC_meth_DEG_files[1])

CpG_methlevel_DEG_data <- read.csv(CpG_methlevel_DEG_files[1])
GpC_methlevel_DEG_data <- read.csv(GpC_methlevel_DEG_files[1])

E4.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E4.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)

#5.5
CpG_meth_DEG_data <- read.csv(CpG_meth_DEG_files[2])
GpC_meth_DEG_data <- read.csv(GpC_meth_DEG_files[2])

CpG_methlevel_DEG_data <- read.csv(CpG_methlevel_DEG_files[2])
GpC_methlevel_DEG_data <- read.csv(GpC_methlevel_DEG_files[2])

E5.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E5.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)

#6.5
CpG_meth_DEG_data <- read.csv(CpG_meth_DEG_files[3])
GpC_meth_DEG_data <- read.csv(GpC_meth_DEG_files[3])

CpG_methlevel_DEG_data <- read.csv(CpG_methlevel_DEG_files[3])
GpC_methlevel_DEG_data <- read.csv(GpC_methlevel_DEG_files[3])

E6.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E6.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)

#7.5
CpG_meth_DEG_data <- read.csv(CpG_meth_DEG_files[4])
GpC_meth_DEG_data <- read.csv(GpC_meth_DEG_files[4])

CpG_methlevel_DEG_data <- read.csv(CpG_methlevel_DEG_files[4])
GpC_methlevel_DEG_data <- read.csv(GpC_methlevel_DEG_files[4])

E7.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E7.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)


#CpG GpC volcano----
E4.5_CpG_DEG_data_choose <- E4.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05)

# E4.5_GpC_DEG_data <- E4.5_GpC_DEG_data %>%
#   mutate(numeric_logFC = parse_number(logFC)) %>%
#   filter(!is.na(numeric_logFC)) %>%
#   filter(adj_p.value_fdr < 0.05) %>%
#   slice_max(numeric_logFC, n = 25000)



volcano_data <- E4.5_CpG_DEG_data_choose %>%
  mutate(
    # chr = rownames(SRP151137_CpG_DEG_P$chr),  # 确保有基因名列
    log_pval = -log10(adj_p.value_fdr),  # 转换校正p值
    direction = case_when(         # 标记上下调
      numeric_logFC > 1 & adj_p.value_fdr < 1e-5 ~ "Up",
      numeric_logFC < -1 & adj_p.value_fdr < 1e-5 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

ggplot(volcano_data, aes(x = Scores, y = log_pval, colour = direction)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (Adj P-value)")+
  theme_bw()+
  # ggtitle("GpC ENCODE DEG Volcano")+
  # 图例
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))



#CpG GpC DEG diff methlevel----
#E4.5
E4.5_CpG_DEG_data_choose <- E4.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-10)

E4.5_GpC_DEG_data_choose <- E4.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-10)

E4.5_CpG_GpC_chrdata <- intersect(E4.5_CpG_DEG_data_choose$chrdata,E4.5_GpC_DEG_data_choose$chrdata)

E4.5_CpG_DEG_data_choose_union <- E4.5_CpG_DEG_data_choose %>%
  filter(chrdata %in% E4.5_CpG_GpC_chrdata)

E4.5_GpC_DEG_data_choose_union <- E4.5_GpC_DEG_data_choose %>%
  filter(chrdata %in% E4.5_CpG_GpC_chrdata)

colnames(E4.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E4.5_CpG_DEG_data_choose_union))
colnames(E4.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E4.5_GpC_DEG_data_choose_union))

E4.5_DEG_data_union <- cbind(E4.5_CpG_DEG_data_choose_union,E4.5_GpC_DEG_data_choose_union)
E4.5_DEG_data_union <- E4.5_DEG_data_union %>%
  arrange(CpG.adj_p.value_fdr,GpC.adj_p.value_fdr) %>%
  arrange(CpG.numeric_logFC,GpC.numeric_logFC)
E4.5_DEG_data_union$top10  <- FALSE
E4.5_DEG_data_union$top10[1:10] <- TRUE  # 前10行标记为TRUE

ggplot(data = E4.5_DEG_data_union,mapping = aes(x = CpG.Diff, y = GpC.Diff))+
  geom_point(size = 3,alpha = 0.9,color = "grey75")+
  geom_point(data = subset(E4.5_DEG_data_union, top10),
             color = "black", size = 3)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "#FF7F00", linetype = "dashed") +  # 垂直参考线（x=50）
  geom_hline(yintercept = 0, color = "#FF7F00", linetype = "dashed")+   # 水平参考线（y=50）
  labs(x = "Differential methylation",y = "Differential accessibility",
       title = "E4.5 GeneTSS2k Differential")+
  theme(legend.title=element_blank(),
        # plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        plot.title = element_textbox(
          # family = "serif",  # 字体
          face = "bold",
          size = 20,
          linewidth = 0,    # 边框粗细
          linetype = 0,       # 实线边框
          lineheight = 1.5,
          color = "black",    # 文字颜色
          fill = "#E69F00",   # 背景填充色(浅蓝)
          box.color = "#E69F00",  # 边框颜色(深蓝)
          hjust = 0.5,       # 水平居中
          # vjust = 0.5,      # 垂直居中
          padding = margin(8, 15, 8, 15),  # 内边距（上、右、下、左）
          margin = margin(b = 10),   # 标题下方外边距
        ),
        axis.text.x = element_text(face = "bold",size = 13),
        axis.text.y = element_text(face = "bold",size = 13),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  geom_label_repel(
    data = subset(E4.5_DEG_data_union, top10),
    aes(label = CpG.gene_name.x),
    color = "#E69F00",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 1,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


#E5.5
E5.5_CpG_DEG_data_choose <- E5.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-5)

E5.5_GpC_DEG_data_choose <- E5.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-5)

E5.5_CpG_GpC_chrdata <- intersect(E5.5_CpG_DEG_data_choose$chrdata,E5.5_GpC_DEG_data_choose$chrdata)

E5.5_CpG_DEG_data_choose_union <- E5.5_CpG_DEG_data_choose %>%
  filter(chrdata %in% E5.5_CpG_GpC_chrdata)

E5.5_GpC_DEG_data_choose_union <- E5.5_GpC_DEG_data_choose %>%
  filter(chrdata %in% E5.5_CpG_GpC_chrdata)

colnames(E5.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E5.5_CpG_DEG_data_choose_union))
colnames(E5.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E5.5_GpC_DEG_data_choose_union))

E5.5_DEG_data_union <- cbind(E5.5_CpG_DEG_data_choose_union,E5.5_GpC_DEG_data_choose_union)
E5.5_DEG_data_union <- E5.5_DEG_data_union %>%
  arrange(CpG.adj_p.value_fdr,GpC.adj_p.value_fdr) %>%
  arrange(CpG.numeric_logFC,GpC.numeric_logFC)
E5.5_DEG_data_union$top10  <- FALSE
E5.5_DEG_data_union$top10[1:10] <- TRUE  # 前10行标记为TRUE

ggplot(data = E5.5_DEG_data_union,mapping = aes(x = CpG.Diff, y = GpC.Diff))+
  geom_point(size = 3,alpha = 0.9,color = "grey75")+
  geom_point(data = subset(E5.5_DEG_data_union, top10),
             color = "black", size = 3)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "#FF7F00", linetype = "dashed") +  # 垂直参考线（x=50）
  geom_hline(yintercept = 0, color = "#FF7F00", linetype = "dashed")+   # 水平参考线（y=50）
  labs(x = "Differential methylation",y = "Differential accessibility",
       title = "E5.5 GeneTSS2k Differential")+
  theme(legend.title=element_blank(),
        # plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        plot.title = element_textbox(
          # family = "serif",  # 字体
          face = "bold",
          size = 20,
          linewidth = 0,    # 边框粗细
          linetype = 0,       # 实线边框
          lineheight = 1.5,
          color = "black",    # 文字颜色
          fill = "#56B4E9",   # 背景填充色(浅蓝)
          box.color = "#56B4E9",  # 边框颜色(深蓝)
          hjust = 0.5,       # 水平居中
          # vjust = 0.5,      # 垂直居中
          padding = margin(8, 15, 8, 15),  # 内边距（上、右、下、左）
          margin = margin(b = 10),   # 标题下方外边距
        ),
        axis.text.x = element_text(face = "bold",size = 13),
        axis.text.y = element_text(face = "bold",size = 13),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  geom_label_repel(
    data = subset(E5.5_DEG_data_union, top10),
    aes(label = CpG.gene_name.x),
    color = "#56B4E9",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


#6.5
E6.5_CpG_DEG_data_choose <- E6.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-6)

E6.5_GpC_DEG_data_choose <- E6.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-6)

E6.5_CpG_GpC_chrdata <- intersect(E6.5_CpG_DEG_data_choose$chrdata,E6.5_GpC_DEG_data_choose$chrdata)

E6.5_CpG_DEG_data_choose_union <- E6.5_CpG_DEG_data_choose %>%
  filter(chrdata %in% E6.5_CpG_GpC_chrdata)

E6.5_GpC_DEG_data_choose_union <- E6.5_GpC_DEG_data_choose %>%
  filter(chrdata %in% E6.5_CpG_GpC_chrdata)

colnames(E6.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E6.5_CpG_DEG_data_choose_union))
colnames(E6.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E6.5_GpC_DEG_data_choose_union))

E6.5_DEG_data_union <- cbind(E6.5_CpG_DEG_data_choose_union,E6.5_GpC_DEG_data_choose_union)
E6.5_DEG_data_union <- E6.5_DEG_data_union %>%
  arrange(CpG.adj_p.value_fdr,GpC.adj_p.value_fdr) %>%
  arrange(CpG.numeric_logFC,GpC.numeric_logFC)
E6.5_DEG_data_union$top10  <- FALSE
E6.5_DEG_data_union$top10[1:10] <- TRUE  # 前10行标记为TRUE

ggplot(data = E6.5_DEG_data_union,mapping = aes(x = CpG.Diff, y = GpC.Diff))+
  geom_point(size = 3,alpha = 0.9,color = "grey75")+
  geom_point(data = subset(E6.5_DEG_data_union, top10),
             color = "black", size = 3)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "#FF7F00", linetype = "dashed") +  # 垂直参考线（x=50）
  geom_hline(yintercept = 0, color = "#FF7F00", linetype = "dashed")+   # 水平参考线（y=50）
  labs(x = "Differential methylation",y = "Differential accessibility",
       title = "E6.5 GeneTSS2k Differential")+
  theme(legend.title=element_blank(),
        # plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        plot.title = element_textbox(
          # family = "serif",  # 字体
          face = "bold",
          size = 20,
          linewidth = 0,    # 边框粗细
          linetype = 0,       # 实线边框
          lineheight = 1.5,
          color = "black",    # 文字颜色
          fill = "#009E73",   # 背景填充色(浅蓝)
          box.color = "#009E73",  # 边框颜色(深蓝)
          hjust = 0.5,       # 水平居中
          # vjust = 0.5,      # 垂直居中
          padding = margin(8, 15, 8, 15),  # 内边距（上、右、下、左）
          margin = margin(b = 10),   # 标题下方外边距
        ),
        axis.text.x = element_text(face = "bold",size = 13),
        axis.text.y = element_text(face = "bold",size = 13),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  geom_label_repel(
    data = subset(E6.5_DEG_data_union, top10),
    aes(label = CpG.gene_name.x),
    color = "#009E73",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


#7.5
E7.5_CpG_DEG_data_choose <- E7.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-10)

E7.5_GpC_DEG_data_choose <- E7.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 1e-10)

E7.5_CpG_GpC_chrdata <- intersect(E7.5_CpG_DEG_data_choose$chrdata,E7.5_GpC_DEG_data_choose$chrdata)

E7.5_CpG_DEG_data_choose_union <- E7.5_CpG_DEG_data_choose %>%
  filter(chrdata %in% E7.5_CpG_GpC_chrdata)

E7.5_GpC_DEG_data_choose_union <- E7.5_GpC_DEG_data_choose %>%
  filter(chrdata %in% E7.5_CpG_GpC_chrdata)

colnames(E7.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E7.5_CpG_DEG_data_choose_union))
colnames(E7.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E7.5_GpC_DEG_data_choose_union))

E7.5_DEG_data_union <- cbind(E7.5_CpG_DEG_data_choose_union,E7.5_GpC_DEG_data_choose_union)
E7.5_DEG_data_union <- E7.5_DEG_data_union %>%
  arrange(CpG.adj_p.value_fdr,GpC.adj_p.value_fdr) %>%
  arrange(CpG.numeric_logFC,GpC.numeric_logFC)
E7.5_DEG_data_union$top10  <- FALSE
E7.5_DEG_data_union$top10[1:10] <- TRUE  # 前10行标记为TRUE

ggplot(data = E7.5_DEG_data_union,mapping = aes(x = CpG.Diff, y = GpC.Diff))+
  geom_point(size = 3,alpha = 0.9,color = "grey75")+
  geom_point(data = subset(E7.5_DEG_data_union, top10),
             color = "black", size = 3)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "#FF7F00", linetype = "dashed") +  # 垂直参考线（x=50）
  geom_hline(yintercept = 0, color = "#FF7F00", linetype = "dashed")+   # 水平参考线（y=50）
  labs(x = "Differential methylation",y = "Differential accessibility",
       title = "E7.5 GeneTSS2k Differential")+
  theme(legend.title=element_blank(),
        # plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        plot.title = element_textbox(
          # family = "serif",  # 字体
          face = "bold",
          size = 20,
          linewidth = 0,    # 边框粗细
          linetype = 0,       # 实线边框
          lineheight = 1.5,
          color = "black",    # 文字颜色
          fill = "#CC79A7",   # 背景填充色(浅蓝)
          box.color = "#CC79A7",  # 边框颜色(深蓝)
          hjust = 0.5,       # 水平居中
          # vjust = 0.5,      # 垂直居中
          padding = margin(8, 15, 8, 15),  # 内边距（上、右、下、左）
          margin = margin(b = 10),   # 标题下方外边距
        ),
        axis.text.x = element_text(face = "bold",size = 13),
        axis.text.y = element_text(face = "bold",size = 13),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  geom_label_repel(
    data = subset(E7.5_DEG_data_union, top10),
    aes(label = CpG.gene_name.x),
    color = "#CC79A7",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )



