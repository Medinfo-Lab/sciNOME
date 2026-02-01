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
library(mixOmics)
library(readr)

load("data/human/Level_Data.RData")
load("data/human/List_Data.RData")
load("data/human/Epi_Group_Data.RData")



#epis data----
GSE255339_CpG_meth_level_data <- GSE255339_CpG_enhancer_methlevel
GSE255339_GpC_meth_level_data <- GSE255339_GpC_enhancer_methlevel
GSE255339_CpG_meth_data <- GSE255339_CpG_enhancer_meth
GSE255339_GpC_meth_data <- GSE255339_GpC_enhancer_meth

GSE255339_CpG_methlevel <- Read_file_colname(GSE255339_CpG_meth_level_data,"level")
GSE255339_GpC_methlevel <- Read_file_colname(GSE255339_GpC_meth_level_data,"level")
GSE255339_CpG_meth <- Read_file_colname(GSE255339_CpG_meth_data,"site")
GSE255339_CpG_UNmeth <- Read_file_colname(GSE255339_CpG_meth_data,"nonsite")
GSE255339_GpC_meth <- Read_file_colname(GSE255339_GpC_meth_data,"site")
GSE255339_GpC_UNmeth <- Read_file_colname(GSE255339_GpC_meth_data,"nonsite")

CpG_sample_title <- GSE255339_CpG_Epis_sample %>%
  filter(level %in% colnames(GSE255339_CpG_methlevel))
GpC_sample_title <- GSE255339_GpC_Epis_sample %>%
  filter(level %in% colnames(GSE255339_GpC_methlevel))

colnames(GSE255339_CpG_methlevel) <- CpG_sample_title$sample
colnames(GSE255339_GpC_methlevel) <- GpC_sample_title$sample

colnames(GSE255339_CpG_meth) <- CpG_sample_title$sample
colnames(GSE255339_CpG_UNmeth) <- CpG_sample_title$sample
colnames(GSE255339_GpC_meth) <- GpC_sample_title$sample
colnames(GSE255339_GpC_UNmeth) <- GpC_sample_title$sample



#NA process----
GSE255339_CpG_methlevel_na_count <- NA_count_order(GSE255339_CpG_methlevel)
GSE255339_GpC_methlevel_na_count <- NA_count_order(GSE255339_GpC_methlevel)

GSE255339_CpG_methlevel_na_count_nona <- GSE255339_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.62)

GSE255339_CpG_methlevel_choose <- GSE255339_CpG_methlevel[rownames(GSE255339_CpG_methlevel_na_count_nona),]

GSE255339_GpC_methlevel_na_count_nona <- GSE255339_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.58)

GSE255339_GpC_methlevel_choose <- GSE255339_GpC_methlevel[rownames(GSE255339_GpC_methlevel_na_count_nona),]


#CpG umap----
sample_cutoff <- 0.8
samples_keep <- colMeans(is.na(GSE255339_CpG_methlevel_choose)) < sample_cutoff
CpG_data_clean <- GSE255339_CpG_methlevel_choose[, samples_keep]

cat("After filtering the samples, the matrix dimensions:", dim(CpG_data_clean), "\n")
cat("Number of excluded samples:", ncol(GSE255339_CpG_methlevel_choose) - ncol(CpG_data_clean), "\n")

CpG_data_clean <- as.matrix(CpG_data_clean)
GSE255339_CpG_methlevel_choose_imputed_res <- impute.knn(CpG_data_clean, k = 10)
GSE255339_CpG_methlevel_choose_imputed_res_data <- GSE255339_CpG_methlevel_choose_imputed_res$data
GSE121650_CpG_methlevel_sample_group <- GSE255339_CpG_Epis_sample %>%
  filter(sample %in% colnames(CpG_data_clean))


#GpC umap----
sample_cutoff <- 0.8
samples_keep <- colMeans(is.na(GSE255339_GpC_methlevel_choose)) < sample_cutoff
GpC_data_clean <- GSE255339_GpC_methlevel_choose[, samples_keep]

cat("After filtering the samples, the matrix dimensions:", dim(GpC_data_clean), "\n")
cat("Number of excluded samples:", ncol(GSE255339_GpC_methlevel_choose) - ncol(GpC_data_clean), "\n")

GpC_data_clean <- as.matrix(GpC_data_clean)
GSE255339_GpC_methlevel_choose_imputed_res <- impute.knn(GpC_data_clean, k = 10)
GSE255339_GpC_methlevel_choose_imputed_res_data <- GSE255339_GpC_methlevel_choose_imputed_res$data
GSE121650_GpC_methlevel_sample_group <- GSE255339_GpC_Epis_sample %>%
  filter(sample %in% colnames(GpC_data_clean))


#CpG GpC mean meth group----
#CpG
GSE255339_CpG_methlevel_choose_colmean <- colMeans(CpG_data_clean,na.rm = T)
GSE255339_CpG_methlevel_choose_colmean <- as.data.frame(GSE255339_CpG_methlevel_choose_colmean)
GSE255339_CpG_methlevel_choose_colmean$group <- GSE121650_CpG_methlevel_sample_group$group1
GSE255339_CpG_methlevel_choose_colmean$group <- factor(GSE255339_CpG_methlevel_choose_colmean$group)

GSE255339_CpG_methlevel_choose_colmean_mean_data <- GSE255339_CpG_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(GSE255339_CpG_methlevel_choose_colmean, na.rm = TRUE))

ggplot(GSE255339_CpG_methlevel_choose_colmean,
       aes(x = group, y = GSE255339_CpG_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.9,scale = "width",aes(linetype = NA))+
  geom_segment(data = GSE255339_CpG_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.05,
                   xend = as.numeric(group) + 0.05,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_classic()+
  labs(y="Percentage (%)",
       x="",
       title = "Enhancer CpG LevelData Means")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  scale_fill_manual(breaks = c("Unt","AZA", "DAC"),
                     values = c("#868686", "#E64B35", "#3C5488"))+
  guides(fill="none")



#GpC
GSE255339_GpC_methlevel_choose_colmean <- colMeans(GpC_data_clean,na.rm = T)
GSE255339_GpC_methlevel_choose_colmean <- as.data.frame(GSE255339_GpC_methlevel_choose_colmean)
GSE255339_GpC_methlevel_choose_colmean$group <- GSE121650_GpC_methlevel_sample_group$group1
GSE255339_GpC_methlevel_choose_colmean$group <- factor(GSE255339_GpC_methlevel_choose_colmean$group)

GSE255339_GpC_methlevel_choose_colmean_mean_data <- GSE255339_GpC_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(GSE255339_GpC_methlevel_choose_colmean, na.rm = TRUE))

ggplot(GSE255339_GpC_methlevel_choose_colmean,
       aes(x = group, y = GSE255339_GpC_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.9,scale = "width",aes(linetype = NA))+
  geom_segment(data = GSE255339_GpC_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.05,
                   xend = as.numeric(group) + 0.05,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_classic()+
  labs(y="Percentage (%)",
       x="",
       title = "Enhancer GpC LevelData Means")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")+
  scale_fill_manual(breaks = c("Unt","AZA", "DAC"),
                    values = c("#868686", "#E64B35", "#3C5488"))



#DEG data process----
CpG_enhancer_meth_DEG_files <- list.files("site_DEG/enhancer/",full.names = T,pattern = "\\CpG.csv*")
GpC_enhancer_meth_DEG_files <- list.files("site_DEG/enhancer/",full.names = T,pattern = "\\GpC.csv*")

CpG_enhancer_methlevel_DEG_files <- list.files("level_DEG/enhancer/",full.names = T,pattern = "\\CpG.csv*")
GpC_enhancer_methlevel_DEG_files <- list.files("level_DEG/enhancer/",full.names = T,pattern = "\\GpC.csv*")

#AZA DAC
CpG_enhancer_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[1])
GpC_enhancer_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[1])

CpG_enhancer_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[1])
GpC_enhancer_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[1])

E4.5_CpG_enhancer_DEG_data <- merge(CpG_enhancer_methlevel_DEG_data,CpG_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E4.5_GpC_enhancer_DEG_data <- merge(GpC_enhancer_methlevel_DEG_data,GpC_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)

#AZA Unt
CpG_enhancer_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[2])
GpC_enhancer_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[2])

CpG_enhancer_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[2])
GpC_enhancer_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[2])

E5.5_CpG_enhancer_DEG_data <- merge(CpG_enhancer_methlevel_DEG_data,CpG_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E5.5_GpC_enhancer_DEG_data <- merge(GpC_enhancer_methlevel_DEG_data,GpC_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)

#DAC Unt
CpG_enhancer_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[3])
GpC_enhancer_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[3])

CpG_enhancer_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[3])
GpC_enhancer_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[3])

E6.5_CpG_enhancer_DEG_data <- merge(CpG_enhancer_methlevel_DEG_data,CpG_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E6.5_GpC_enhancer_DEG_data <- merge(GpC_enhancer_methlevel_DEG_data,GpC_enhancer_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)

# save.image("Enhancer_Differential_Regional_Data.RData")


#CpG GpC DEG diff methlevel----
#AZA DAC
E4.5_CpG_DEG_data_choose <- E4.5_CpG_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)
# slice_max(numeric_logFC, n = 25000)

E4.5_GpC_DEG_data_choose <- E4.5_GpC_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)
# slice_max(numeric_logFC, n = 25000)

E4.5_CpG_GpC_chrdata <- intersect(E4.5_CpG_DEG_data_choose$chr,E4.5_GpC_DEG_data_choose$chr)

E4.5_CpG_DEG_data_choose_union <- E4.5_CpG_DEG_data_choose %>%
  filter(chr %in% E4.5_CpG_GpC_chrdata)

E4.5_GpC_DEG_data_choose_union <- E4.5_GpC_DEG_data_choose %>%
  filter(chr %in% E4.5_CpG_GpC_chrdata)

colnames(E4.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E4.5_CpG_DEG_data_choose_union))
colnames(E4.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E4.5_GpC_DEG_data_choose_union))

E4.5_DEG_data_union <- cbind(E4.5_CpG_DEG_data_choose_union,E4.5_GpC_DEG_data_choose_union)
E4.5_DEG_data_union <- E4.5_DEG_data_union %>%
  arrange(CpG.p.value,GpC.p.value) %>%
  arrange(CpG.logFC,GpC.logFC)
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
       title = "AZA vs DAC Enhancer Differential")+
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
    aes(label = CpG.geneone.x),
    color = "#E69F00",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


#E5.5
E5.5_CpG_DEG_data_choose <- E5.5_CpG_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)

E5.5_GpC_DEG_data_choose <- E5.5_GpC_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)

E5.5_CpG_GpC_chrdata <- intersect(E5.5_CpG_DEG_data_choose$chr,E5.5_GpC_DEG_data_choose$chr)

E5.5_CpG_DEG_data_choose_union <- E5.5_CpG_DEG_data_choose %>%
  filter(chr %in% E5.5_CpG_GpC_chrdata)

E5.5_GpC_DEG_data_choose_union <- E5.5_GpC_DEG_data_choose %>%
  filter(chr %in% E5.5_CpG_GpC_chrdata)

colnames(E5.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E5.5_CpG_DEG_data_choose_union))
colnames(E5.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E5.5_GpC_DEG_data_choose_union))

E5.5_DEG_data_union <- cbind(E5.5_CpG_DEG_data_choose_union,E5.5_GpC_DEG_data_choose_union)
E5.5_DEG_data_union <- E5.5_DEG_data_union %>%
  arrange(CpG.p.value,GpC.p.value) %>%
  arrange(CpG.logFC,GpC.logFC)
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
       title = "AZA vs Unt Enhancer Differential")+
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
    aes(label = CpG.geneone.x),
    color = "#56B4E9",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


#6.5
E6.5_CpG_DEG_data_choose <- E6.5_CpG_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)

E6.5_GpC_DEG_data_choose <- E6.5_GpC_enhancer_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 1e-5)

E6.5_CpG_GpC_chrdata <- intersect(E6.5_CpG_DEG_data_choose$chr,E6.5_GpC_DEG_data_choose$chr)

E6.5_CpG_DEG_data_choose_union <- E6.5_CpG_DEG_data_choose %>%
  filter(chr %in% E6.5_CpG_GpC_chrdata)

E6.5_GpC_DEG_data_choose_union <- E6.5_GpC_DEG_data_choose %>%
  filter(chr %in% E6.5_CpG_GpC_chrdata)

colnames(E6.5_CpG_DEG_data_choose_union) <- paste0("CpG.", colnames(E6.5_CpG_DEG_data_choose_union))
colnames(E6.5_GpC_DEG_data_choose_union) <- paste0("GpC.", colnames(E6.5_GpC_DEG_data_choose_union))

E6.5_DEG_data_union <- cbind(E6.5_CpG_DEG_data_choose_union,E6.5_GpC_DEG_data_choose_union)
E6.5_DEG_data_union <- E6.5_DEG_data_union %>%
  arrange(CpG.p.value,GpC.p.value) %>%
  arrange(CpG.logFC,GpC.logFC)
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
       title = "DAC vs Unt Enhancer Differential")+
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
    aes(label = CpG.geneone.x),
    color = "#009E73",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )


