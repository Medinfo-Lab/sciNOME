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
library(Seurat)
library(ggpubr)
library(patchwork)
library(ggalluvial)
library(ggridges)
library(GGally)
library(SNFtool)
library(readr)


load("data/multi_omics_sourcedata.RData")


#RNA----
#4.5
RNA_marker_data_E4.5 <- RNA_marker_data %>%
  filter(cluster == "E4.5")
RNA_marker_data_E4.5 <- RNA_marker_data_E4.5 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E4.5_counts <- GSE121650_pbmc_counts_frame[RNA_marker_data_E4.5$gene,]

RNA_cells_E4.5 <- GSE121650_RNA_sample$Sample_title[GSE121650_RNA_sample$Developmental_stage == "E4.5"]
valid_cells <- intersect(RNA_cells_E4.5, colnames(RNA_marker_data_E4.5_counts))
sub_mat_E4.5 <- RNA_marker_data_E4.5_counts[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "Expression"

E4.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E4.5), Expression = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE121650_RNA_data_gene, by = c("GeneID" = "Geneid"))



#E5.5
RNA_marker_data_E5.5 <- RNA_marker_data %>%
  filter(cluster == "E5.5")
RNA_marker_data_E5.5 <- RNA_marker_data_E5.5 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E5.5_counts <- GSE121650_pbmc_counts_frame[RNA_marker_data_E5.5$gene,]

RNA_cells_E5.5 <- GSE121650_RNA_sample$Sample_title[GSE121650_RNA_sample$Developmental_stage == "E5.5"]
valid_cells <- intersect(RNA_cells_E5.5, colnames(RNA_marker_data_E5.5_counts))
sub_mat_E5.5 <- RNA_marker_data_E5.5_counts[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "Expression"

E5.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E5.5), Expression = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE121650_RNA_data_gene, by = c("GeneID" = "Geneid"))



#E6.5
RNA_marker_data_E6.5 <- RNA_marker_data %>%
  filter(cluster == "E6.5")
RNA_marker_data_E6.5 <- RNA_marker_data_E6.5 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E6.5_counts <- GSE121650_pbmc_counts_frame[RNA_marker_data_E6.5$gene,]

RNA_cells_E6.5 <- GSE121650_RNA_sample$Sample_title[GSE121650_RNA_sample$Developmental_stage == "E6.5"]
valid_cells <- intersect(RNA_cells_E6.5, colnames(RNA_marker_data_E6.5_counts))
sub_mat_E6.5 <- RNA_marker_data_E6.5_counts[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "Expression"

E6.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E6.5), Expression = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE121650_RNA_data_gene, by = c("GeneID" = "Geneid"))



#E7.5
RNA_marker_data_E7.5 <- RNA_marker_data %>%
  filter(cluster == "E7.5")
RNA_marker_data_E7.5 <- RNA_marker_data_E7.5 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E7.5_counts <- GSE121650_pbmc_counts_frame[RNA_marker_data_E7.5$gene,]

RNA_cells_E7.5 <- GSE121650_RNA_sample$Sample_title[GSE121650_RNA_sample$Developmental_stage == "E7.5"]
valid_cells <- intersect(RNA_cells_E7.5, colnames(RNA_marker_data_E7.5_counts))
sub_mat_E7.5 <- RNA_marker_data_E7.5_counts[, valid_cells, drop = FALSE]
avg_E7.5 <- as.data.frame(rowMeans(sub_mat_E7.5, na.rm = TRUE))
colnames(avg_E7.5) <- "Expression"

E7.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E7.5), Expression = avg_E7.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE121650_RNA_data_gene, by = c("GeneID" = "Geneid"))


#genetss2k----
GSE121690_CpG_genetss2k_methlevel <- Read_file_meth_colname(GSE121690_CpG_genetss2k_meth_level_data,"methlevel")
GSE121690_GpC_genetss2k_methlevel <- Read_file_meth_colname(GSE121690_GpC_genetss2k_meth_level_data,"methlevel")

CpG_sample_title <- GSE121690_CpG_Epis_sample %>%
  filter(methlevel %in% colnames(GSE121690_CpG_genetss2k_methlevel))
GpC_sample_title <- GSE121690_GpC_Epis_sample %>%
  filter(methlevel %in% colnames(GSE121690_GpC_genetss2k_methlevel))

colnames(GSE121690_CpG_genetss2k_methlevel) <- CpG_sample_title$Sample_title
colnames(GSE121690_GpC_genetss2k_methlevel) <- GpC_sample_title$Sample_title



#4.5----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[1])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[1])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[1])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[1])

E4.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E4.5_CpG_DEG_data <- E4.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E4.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E4.5_GpC_DEG_data <- E4.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E4.5_common_data <- unique(intersect(intersect(E4.5_CpG_DEG_data$gene_name.x, E4.5_GpC_DEG_data$gene_name.x),
                                E4.5_RNA_df_result$gene_name))

E4.5_RNA_df_result_common <- E4.5_RNA_df_result %>%
  filter(gene_name %in% E4.5_common_data)
E4.5_RNA_df_result_common_unique <- E4.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E4.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E4.5"]
valid_cells <- intersect(CpG_cells_E4.5, colnames(GSE121690_CpG_genetss2k_methlevel))
sub_mat_E4.5 <- GSE121690_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "CpG_level"

E4.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E4.5), CpG_level = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E4.5_CpG_df_result_common <- E4.5_CpG_df_result %>%
  filter(genename %in% E4.5_common_data)
E4.5_CpG_df_result_common_unique <- E4.5_CpG_df_result_common %>%
  distinct(genename, .keep_all = TRUE)



#GpC
GpC_cells_E4.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E4.5"]
valid_cells <- intersect(GpC_cells_E4.5, colnames(GSE121690_GpC_genetss2k_methlevel))
sub_mat_E4.5 <- GSE121690_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "GpC_level"

E4.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E4.5), GpC_level = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E4.5_GpC_df_result_common <- E4.5_GpC_df_result %>%
  filter(genename %in% E4.5_common_data)
E4.5_GpC_df_result_common_unique <- E4.5_GpC_df_result_common %>%
  distinct(genename, .keep_all = TRUE)


#CpG GpC RNA
E4.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E4.5_common_data,
  RNA_Exp = E4.5_RNA_df_result_common_unique$Expression,
  CpG_level = E4.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E4.5_GpC_df_result_common_unique$GpC_level
)


#小提琴图 + 箱中箱
plot_data <- E4.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
               names_to = "Omics_Type",
               values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "GeneTSS2k E4.5 Violin Plot", y = "Value")


#相关性散点图
# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k E4.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E4.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E4.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3


#相关性矩阵图
ggpairs(E4.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k E4.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()


#5.5----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[2])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[2])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[2])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[2])

E5.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E5.5_CpG_DEG_data <- E5.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E5.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E5.5_GpC_DEG_data <- E5.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E5.5_common_data <- unique(intersect(intersect(E5.5_CpG_DEG_data$gene_name.x, E5.5_GpC_DEG_data$gene_name.x),
                                     E5.5_RNA_df_result$gene_name))

E5.5_RNA_df_result_common <- E5.5_RNA_df_result %>%
  filter(gene_name %in% E5.5_common_data)

E5.5_RNA_df_result_common_unique <- E5.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)



#CpG
CpG_cells_E5.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E5.5"]
valid_cells <- intersect(CpG_cells_E5.5, colnames(GSE121690_CpG_genetss2k_methlevel))
sub_mat_E5.5 <- GSE121690_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "CpG_level"

E5.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E5.5), CpG_level = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E5.5_CpG_df_result_common <- E5.5_CpG_df_result %>%
  filter(genename %in% E5.5_common_data)
E5.5_CpG_df_result_common_unique <- E5.5_CpG_df_result_common %>%
  distinct(genename, .keep_all = TRUE)



#GpC
GpC_cells_E5.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E5.5"]
valid_cells <- intersect(GpC_cells_E5.5, colnames(GSE121690_GpC_genetss2k_methlevel))
sub_mat_E5.5 <- GSE121690_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "GpC_level"

E5.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E5.5), GpC_level = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E5.5_GpC_df_result_common <- E5.5_GpC_df_result %>%
  filter(genename %in% E5.5_common_data)
E5.5_GpC_df_result_common_unique <- E5.5_GpC_df_result_common %>%
  distinct(genename, .keep_all = TRUE)


#CpG GpC RNA
E5.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E5.5_common_data,
  RNA_Exp = E5.5_RNA_df_result_common_unique$Expression,
  CpG_level = E5.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E5.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E5.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "GeneTSS2k E5.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k E5.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E5.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E5.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3


#相关性矩阵图
ggpairs(E5.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k E5.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()



#6.5----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[3])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[3])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[3])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[3])

E6.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E6.5_CpG_DEG_data <- E6.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E6.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E6.5_GpC_DEG_data <- E6.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E6.5_common_data <- unique(intersect(intersect(E6.5_CpG_DEG_data$gene_name.x, E6.5_GpC_DEG_data$gene_name.x),
                                     E6.5_RNA_df_result$gene_name))

E6.5_RNA_df_result_common <- E6.5_RNA_df_result %>%
  filter(gene_name %in% E6.5_common_data)
E6.5_RNA_df_result_common_unique <- E6.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E6.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E6.5"]
valid_cells <- intersect(CpG_cells_E6.5, colnames(GSE121690_CpG_genetss2k_methlevel))
sub_mat_E6.5 <- GSE121690_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "CpG_level"

E6.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E6.5), CpG_level = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E6.5_CpG_df_result_common <- E6.5_CpG_df_result %>%
  filter(genename %in% E6.5_common_data)
E6.5_CpG_df_result_common_unique <- E6.5_CpG_df_result_common %>%
  distinct(genename, .keep_all = TRUE)



#GpC
GpC_cells_E6.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E6.5"]
valid_cells <- intersect(GpC_cells_E6.5, colnames(GSE121690_GpC_genetss2k_methlevel))
sub_mat_E6.5 <- GSE121690_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "GpC_level"

E6.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E6.5), GpC_level = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E6.5_GpC_df_result_common <- E6.5_GpC_df_result %>%
  filter(genename %in% E6.5_common_data)
E6.5_GpC_df_result_common_unique <- E6.5_GpC_df_result_common %>%
  distinct(genename, .keep_all = TRUE)


#CpG GpC RNA
E6.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E6.5_common_data,
  RNA_Exp = E6.5_RNA_df_result_common_unique$Expression,
  CpG_level = E6.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E6.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E6.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "GeneTSS2k E6.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k E6.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E6.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E6.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E6.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k E6.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()






#7.5----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[4])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[4])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[4])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[4])

E7.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E7.5_CpG_DEG_data <- E7.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E7.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E7.5_GpC_DEG_data <- E7.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E7.5_common_data <- unique(intersect(intersect(E7.5_CpG_DEG_data$gene_name.x, E7.5_GpC_DEG_data$gene_name.x),
                                     E7.5_RNA_df_result$gene_name))

E7.5_RNA_df_result_common <- E7.5_RNA_df_result %>%
  filter(gene_name %in% E7.5_common_data)
E7.5_RNA_df_result_common_unique <- E7.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E7.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E7.5"]
valid_cells <- intersect(CpG_cells_E7.5, colnames(GSE121690_CpG_genetss2k_methlevel))
sub_mat_E7.5 <- GSE121690_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E7.5 <- as.data.frame(rowMeans(sub_mat_E7.5, na.rm = TRUE))
colnames(avg_E7.5) <- "CpG_level"

E7.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E7.5), CpG_level = avg_E7.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E7.5_CpG_df_result_common <- E7.5_CpG_df_result %>%
  filter(genename %in% E7.5_common_data)
E7.5_CpG_df_result_common_unique <- E7.5_CpG_df_result_common %>%
  distinct(genename, .keep_all = TRUE)



#GpC
GpC_cells_E7.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E7.5"]
valid_cells <- intersect(GpC_cells_E7.5, colnames(GSE121690_GpC_genetss2k_methlevel))
sub_mat_E7.5 <- GSE121690_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
avg_E7.5 <- as.data.frame(rowMeans(sub_mat_E7.5, na.rm = TRUE))
colnames(avg_E7.5) <- "GpC_level"

E7.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E7.5), GpC_level = avg_E7.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(genetss2k_data_paste, by = c("chrdata" = "chrdata"))

E7.5_GpC_df_result_common <- E7.5_GpC_df_result %>%
  filter(genename %in% E7.5_common_data)
E7.5_GpC_df_result_common_unique <- E7.5_GpC_df_result_common %>%
  distinct(genename, .keep_all = TRUE)


#CpG GpC RNA
E7.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E7.5_common_data,
  RNA_Exp = E7.5_RNA_df_result_common_unique$Expression,
  CpG_level = E7.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E7.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E7.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "GeneTSS2k E7.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k E7.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E7.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k E7.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E6.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k E7.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()







#enhancer----
GSE121690_CpG_enhancer_methlevel <- Read_file_meth_colname(GSE121690_CpG_enhancer_meth_level_data,"methlevel")
GSE121690_GpC_enhancer_methlevel <- Read_file_meth_colname(GSE121690_GpC_enhancer_meth_level_data,"methlevel")

CpG_sample_title <- GSE121690_CpG_Epis_sample %>%
  filter(methlevel %in% colnames(GSE121690_CpG_enhancer_methlevel))
GpC_sample_title <- GSE121690_GpC_Epis_sample %>%
  filter(methlevel %in% colnames(GSE121690_GpC_enhancer_methlevel))

colnames(GSE121690_CpG_enhancer_methlevel) <- CpG_sample_title$Sample_title
colnames(GSE121690_GpC_enhancer_methlevel) <- GpC_sample_title$Sample_title


#4.5----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[1])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[1])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[1])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[1])

E4.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E4.5_CpG_DEG_data <- E4.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E4.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E4.5_GpC_DEG_data <- E4.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E4.5_common_data <- unique(intersect(intersect(E4.5_CpG_DEG_data$geneone.x, E4.5_GpC_DEG_data$geneone.x),
                                     E4.5_RNA_df_result$gene_name))

E4.5_RNA_df_result_common <- E4.5_RNA_df_result %>%
  filter(gene_name %in% E4.5_common_data)
E4.5_RNA_df_result_common_unique <- E4.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E4.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E4.5"]
valid_cells <- intersect(CpG_cells_E4.5, colnames(GSE121690_CpG_enhancer_methlevel))
sub_mat_E4.5 <- GSE121690_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "CpG_level"

E4.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E4.5), CpG_level = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E4.5_CpG_df_result_common <- E4.5_CpG_df_result %>%
  filter(geneone %in% E4.5_common_data)
E4.5_CpG_df_result_common_unique <- E4.5_CpG_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)



#GpC
GpC_cells_E4.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E4.5"]
valid_cells <- intersect(GpC_cells_E4.5, colnames(GSE121690_GpC_enhancer_methlevel))
sub_mat_E4.5 <- GSE121690_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "GpC_level"

E4.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E4.5), GpC_level = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E4.5_GpC_df_result_common <- E4.5_GpC_df_result %>%
  filter(geneone %in% E4.5_common_data)
E4.5_GpC_df_result_common_unique <- E4.5_GpC_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)


#CpG GpC RNA
E4.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E4.5_common_data,
  RNA_Exp = E4.5_RNA_df_result_common_unique$Expression,
  CpG_level = E4.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E4.5_GpC_df_result_common_unique$GpC_level
)


#小提琴图 + 箱中箱
plot_data <- E4.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "Enhancer E4.5 Violin Plot", y = "Value")


#相关性散点图
# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer E4.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E4.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E4.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E4.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer E4.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()



#5.5----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[2])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[2])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[2])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[2])

E5.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E5.5_CpG_DEG_data <- E5.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E5.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E5.5_GpC_DEG_data <- E5.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E5.5_common_data <- unique(intersect(intersect(E5.5_CpG_DEG_data$geneone.x, E5.5_GpC_DEG_data$geneone.x),
                                     E5.5_RNA_df_result$gene_name))

E5.5_RNA_df_result_common <- E5.5_RNA_df_result %>%
  filter(gene_name %in% E5.5_common_data)
E5.5_RNA_df_result_common_unique <- E5.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E5.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E5.5"]
valid_cells <- intersect(CpG_cells_E5.5, colnames(GSE121690_CpG_enhancer_methlevel))
sub_mat_E5.5 <- GSE121690_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "CpG_level"

E5.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E5.5), CpG_level = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E5.5_CpG_df_result_common <- E5.5_CpG_df_result %>%
  filter(geneone %in% E5.5_common_data)
E5.5_CpG_df_result_common_unique <- E5.5_CpG_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)



#GpC
GpC_cells_E5.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E5.5"]
valid_cells <- intersect(GpC_cells_E5.5, colnames(GSE121690_GpC_enhancer_methlevel))
sub_mat_E5.5 <- GSE121690_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "GpC_level"

E5.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E5.5), GpC_level = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E5.5_GpC_df_result_common <- E5.5_GpC_df_result %>%
  filter(geneone %in% E5.5_common_data)
E5.5_GpC_df_result_common_unique <- E5.5_GpC_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)


#CpG GpC RNA
E5.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E5.5_common_data,
  RNA_Exp = E5.5_RNA_df_result_common_unique$Expression,
  CpG_level = E5.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E5.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E5.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "Enhancer E5.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer E5.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E5.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E5.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E5.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer E5.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()




#6.5----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[3])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[3])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[3])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[3])

E6.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E6.5_CpG_DEG_data <- E6.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E6.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E6.5_GpC_DEG_data <- E6.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E6.5_common_data <- unique(intersect(intersect(E6.5_CpG_DEG_data$geneone.x, E6.5_GpC_DEG_data$geneone.x),
                                     E6.5_RNA_df_result$gene_name))

E6.5_RNA_df_result_common <- E6.5_RNA_df_result %>%
  filter(gene_name %in% E6.5_common_data)
E6.5_RNA_df_result_common_unique <- E6.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E6.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E6.5"]
valid_cells <- intersect(CpG_cells_E6.5, colnames(GSE121690_CpG_enhancer_methlevel))
sub_mat_E6.5 <- GSE121690_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "CpG_level"

E6.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E6.5), CpG_level = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E6.5_CpG_df_result_common <- E6.5_CpG_df_result %>%
  filter(geneone %in% E6.5_common_data)
E6.5_CpG_df_result_common_unique <- E6.5_CpG_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)



#GpC
GpC_cells_E6.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E6.5"]
valid_cells <- intersect(GpC_cells_E6.5, colnames(GSE121690_GpC_enhancer_methlevel))
sub_mat_E6.5 <- GSE121690_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "GpC_level"

E6.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E6.5), GpC_level = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E6.5_GpC_df_result_common <- E6.5_GpC_df_result %>%
  filter(geneone %in% E6.5_common_data)
E6.5_GpC_df_result_common_unique <- E6.5_GpC_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)


#CpG GpC RNA
E6.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E6.5_common_data,
  RNA_Exp = E6.5_RNA_df_result_common_unique$Expression,
  CpG_level = E6.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E6.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E6.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "Enhancer E6.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer E6.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E6.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E6.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E6.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer E6.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()






#7.5----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[4])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[4])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[4])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[4])

E7.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E7.5_CpG_DEG_data <- E7.5_CpG_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)
E7.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chrdata",by.y = "chrdata",all.x = T)
E7.5_GpC_DEG_data <- E7.5_GpC_DEG_data %>%
  mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(numeric_logFC)) %>%
  filter(adj_p.value_fdr < 0.05) %>%
  slice_max(numeric_logFC, n = 25000)

E7.5_common_data <- unique(intersect(intersect(E7.5_CpG_DEG_data$geneone.x, E7.5_GpC_DEG_data$geneone.x),
                                     E7.5_RNA_df_result$gene_name))

E7.5_RNA_df_result_common <- E7.5_RNA_df_result %>%
  filter(gene_name %in% E7.5_common_data)
E7.5_RNA_df_result_common_unique <- E7.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E7.5 <- GSE121690_CpG_Epis_sample$Sample_title[GSE121690_CpG_Epis_sample$Developmental_stage == "E7.5"]
valid_cells <- intersect(CpG_cells_E7.5, colnames(GSE121690_CpG_enhancer_methlevel))
sub_mat_E7.5 <- GSE121690_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E7.5 <- as.data.frame(rowMeans(sub_mat_E7.5, na.rm = TRUE))
colnames(avg_E7.5) <- "CpG_level"

E7.5_CpG_df_result <- data.frame(chrdata = rownames(avg_E7.5), CpG_level = avg_E7.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E7.5_CpG_df_result_common <- E7.5_CpG_df_result %>%
  filter(geneone %in% E7.5_common_data)
E7.5_CpG_df_result_common_unique <- E7.5_CpG_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)



#GpC
GpC_cells_E7.5 <- GSE121690_GpC_Epis_sample$Sample_title[GSE121690_GpC_Epis_sample$Developmental_stage == "E7.5"]
valid_cells <- intersect(GpC_cells_E7.5, colnames(GSE121690_GpC_enhancer_methlevel))
sub_mat_E7.5 <- GSE121690_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
avg_E7.5 <- as.data.frame(rowMeans(sub_mat_E7.5, na.rm = TRUE))
colnames(avg_E7.5) <- "GpC_level"

E7.5_GpC_df_result <- data.frame(chrdata = rownames(avg_E7.5), GpC_level = avg_E7.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(enhancer_data_paste, by = c("chrdata" = "chrdata"))

E7.5_GpC_df_result_common <- E7.5_GpC_df_result %>%
  filter(geneone %in% E7.5_common_data)
E7.5_GpC_df_result_common_unique <- E7.5_GpC_df_result_common %>%
  distinct(geneone, .keep_all = TRUE)


#CpG GpC RNA
E7.5_CpG_GpC_RNA_data <- data.frame(
  GeneName = E7.5_common_data,
  RNA_Exp = E7.5_RNA_df_result_common_unique$Expression,
  CpG_level = E7.5_CpG_df_result_common_unique$CpG_level,
  GpC_level = E7.5_GpC_df_result_common_unique$GpC_level
)

plot_data <- E7.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
                      names_to = "Omics_Type",
                      values_to = "Value")
# 3. 设置因子水平，让X轴按顺序排列
plot_data$Omics_Type <- factor(plot_data$Omics_Type,
                               levels = c("RNA_Exp", "CpG_level", "GpC_level"))

# 定义颜色
my_colors <- c("#E7B800", "#2E9FDF", "#FC4E07")

ggplot(plot_data, aes(x = Omics_Type, y = Value, fill = Omics_Type)) +
  # 绘制小提琴图
  geom_violin(trim = FALSE, alpha = 0.5) +
  # 在内部叠加一个细窄的箱线图
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # 分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(title = "Enhancer E7.5 Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer E7.5 CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E7.5 GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E7.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer E7.5 CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E7.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer E7.5 Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()






#Joint Dimensionality Reduction----
#genetss2k
CpG_genetss2k_pca_data
CpG_genetss2k_umap_data
CpG_genetss2k_mds_data

GpC_genetss2k_pca_data
GpC_genetss2k_umap_data
GpC_genetss2k_mds_data

CpG_enhancer_pca_data
CpG_enhancer_umap_data
CpG_enhancer_mds_data

GpC_enhancer_pca_data
GpC_enhancer_umap_data
GpC_enhancer_mds_data

CpG_genebody_pca_data
CpG_genebody_umap_data
CpG_genebody_mds_data

GpC_genebody_pca_data
GpC_genebody_umap_data
GpC_genebody_mds_data


#CpG genetss2k
pca_scaled  <- scale(CpG_genetss2k_pca_data[,c(1:30)])
mds_scaled  <- scale(CpG_genetss2k_mds_data[,c(1:2)])
umap_scaled <- scale(CpG_genetss2k_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
CpG_genetss2k_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = CpG_genetss2k_final_umap[,1],
  Y = CpG_genetss2k_final_umap[,2],
  Stage = CpG_genetss2k_pca_data$Stage
  )

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "GeneTSS2k CpG Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))


#GpC genetss2k
pca_scaled  <- scale(GpC_genetss2k_pca_data[,c(1:30)])
mds_scaled  <- scale(GpC_genetss2k_mds_data[,c(1:2)])
umap_scaled <- scale(GpC_genetss2k_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
GpC_genetss2k_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = GpC_genetss2k_final_umap[,1],
  Y = GpC_genetss2k_final_umap[,2],
  Stage = GpC_genetss2k_pca_data$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "GeneTSS2k GpC Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))


#CpG enhancer
pca_scaled  <- scale(CpG_enhancer_pca_data[,c(1:30)])
mds_scaled  <- scale(CpG_enhancer_mds_data[,c(1:2)])
umap_scaled <- scale(CpG_enhancer_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
CpG_enhancer_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = CpG_enhancer_final_umap[,1],
  Y = CpG_enhancer_final_umap[,2],
  Stage = CpG_enhancer_pca_data$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "Enhancer CpG Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))


#GpC enhancer
pca_scaled  <- scale(GpC_enhancer_pca_data[,c(1:30)])
mds_scaled  <- scale(GpC_enhancer_mds_data[,c(1:2)])
umap_scaled <- scale(GpC_enhancer_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
GpC_enhancer_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = GpC_enhancer_final_umap[,1],
  Y = GpC_enhancer_final_umap[,2],
  Stage = GpC_enhancer_pca_data$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "Enhancer GpC Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))

#CpG genebody
pca_scaled  <- scale(CpG_genebody_pca_data[,c(1:30)])
mds_scaled  <- scale(CpG_genebody_mds_data[,c(1:2)])
umap_scaled <- scale(CpG_genebody_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
CpG_genebody_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = CpG_genebody_final_umap[,1],
  Y = CpG_genebody_final_umap[,2],
  Stage = CpG_genebody_pca_data$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "GeneBody CpG Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))


#GpC genebody
pca_scaled  <- scale(GpC_genebody_pca_data[,c(1:30)])
mds_scaled  <- scale(GpC_genebody_mds_data[,c(1:2)])
umap_scaled <- scale(GpC_genebody_umap_data[,c(1:2)])

PCA_join_data <- cbind(pca_scaled, mds_scaled, umap_scaled)
GpC_genebody_final_umap <- umap(PCA_join_data)

meta_plot_data <- data.frame(
  X = GpC_genebody_final_umap[,1],
  Y = GpC_genebody_final_umap[,2],
  Stage = GpC_genebody_pca_data$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "GeneBody GpC Integration of PCA + MDS + UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))






#Joint Dimensionality Reduction Multi Omics----
RNA_umap_data <- read.csv("plot/RNA/umap.csv")
RNA_umap_data_choose <- left_join(
  RNA_umap_data,
  GSE121650_RNA_sample %>% dplyr::select(Run, Sample_title), # <-- This is the key change!
  by = c("X" = "Run")
)
rownames(RNA_umap_data_choose) <- RNA_umap_data_choose$Sample_title

#以MDS为例
CpG_genetss2k_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_CpG_mds.csv",row.names = 1)
GpC_genetss2k_mds_data <- read.csv("plot/genetss2k/UMAP_DATA/genetss2k_GpC_mds.csv",row.names = 1)
CpG_enhancer_mds_data <- read.csv("plot/enhancer/UMAP_DATA/enhancer_CpG_mds.csv",row.names = 1)
GpC_enhancer_mds_data <- read.csv("plot/enhancer/UMAP_DATA/enhancer_GpC_mds.csv",row.names = 1)
CpG_genebody_mds_data <- read.csv("plot/genebody/UMAP_DATA/genebody_CpG_mds.csv",row.names = 1)
GpC_genebody_mds_data <- read.csv("plot/genebody/UMAP_DATA/genebody_GpC_mds.csv",row.names = 1)


sample_list <- list(
  CpG_genetss2k=CpG_genetss2k_mds_data,
  GpC_genetss2k=GpC_genetss2k_mds_data,
  CpG_enhancer=CpG_enhancer_mds_data,
  GpC_enhancer=GpC_enhancer_mds_data,
  CpG_genebody=CpG_genebody_mds_data,
  GpC_genebody=GpC_genebody_mds_data,
  RNA=RNA_umap_data_choose
  )

sample_name_list <- lapply(sample_list, rownames)
sample_common_names <- Reduce(intersect, sample_name_list)

sample_common_final_list <- lapply(sample_list, function(df) {
  return(df[sample_common_names, , drop = FALSE])
})

common_final_dim <- data.frame(
  CpG_genetss2k=sample_common_final_list[[1]][,c(1,2)],
  GpC_genetss2k=sample_common_final_list[[2]][,c(1,2)],
  CpG_enhancer=sample_common_final_list[[3]][,c(1,2)],
  GpC_enhancer=sample_common_final_list[[4]][,c(1,2)],
  CpG_genebody=sample_common_final_list[[5]][,c(1,2)],
  GpC_genebody=sample_common_final_list[[6]][,c(1,2)]
  # RNA=sample_common_final_list[[7]][,c("UMAP1","UMAP2")]
)
common_final_dim_umap <- umap(common_final_dim)

meta_plot_data <- data.frame(
  X = common_final_dim_umap[,1],
  Y = common_final_dim_umap[,2],
  Stage = sample_common_final_list[[1]]$Stage
)

ggplot(meta_plot_data, aes(x=X, y=Y, color=Stage)) +
  geom_point(meta_plot_data,mapping=aes(X,Y,color=Stage),size = 3.5,alpha = 0.85)+
  theme_bw()+
  labs(x = "X",y = "Y",title = "Integration of MDS")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))




