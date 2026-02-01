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
library(Seurat)
library(ggpubr)
library(patchwork)
library(ggalluvial)
library(ggridges)
library(GGally)
library(SNFtool)
library(readr)


# save.image("Multi_Omics_Data.RData")
load("data/human/Multi_Omics_Data.RData")
load("data/human/Level_Data.RData")
load("data/human/List_Data.RData")


GSE255339_CpG_Epis_sample <- read_xlsx("GSE255339_sample.xlsx",sheet = "CpG")
GSE255339_GpC_Epis_sample <- read_xlsx("GSE255339_sample.xlsx",sheet = "GpC")
GSE255339_RNA_sample <- read_xlsx("GSE255339_sample.xlsx",sheet = "RNA3")


#RNA----
GSE255339_RNA_data <- read.csv("data/GSE255339_RNA_counts_genename_sum.csv")

GSE255339_RNA_data_gene <- GSE255339_RNA_data[,c("Geneid","Length","gene_name")]
GSE255339_RNA_data_gene <-  unique(GSE255339_RNA_data_gene, by = c("Geneid", "Length","gene_name"))

GSE255339_RNA_data_counts <- GSE255339_RNA_data[,c(-3,-2)]
# any(duplicated(GSE255339_RNA_data_counts$Geneid))
# any(is.na(GSE255339_RNA_data$gene_name))

dedup_dt <- GSE255339_RNA_data_counts[!duplicated(GSE255339_RNA_data_counts$Geneid), ]
GSE255339_RNA_data_counts_filter <- dedup_dt
rownames(GSE255339_RNA_data_counts_filter) <- GSE255339_RNA_data_counts_filter$Geneid
GSE255339_RNA_data_counts_filter <- GSE255339_RNA_data_counts_filter[,-1]
all(colnames(GSE255339_RNA_data_counts_filter)==GSE255339_RNA_sample$Sample_Geo)
colnames(GSE255339_RNA_data_counts_filter) <- GSE255339_RNA_sample$sample

#counts
GSE255339_pbmc <- CreatRNAObject(GSE255339_RNA_data_counts_filter,
                                 min_cells = 50, min_features = 0)
rownames(GSE255339_pbmc) <- gsub("-", "_", rownames(GSE255339_pbmc))

GSE255339_RNA_sample_seurat_group <- GSE255339_RNA_sample %>%
  filter(sample %in% colnames(GSE255339_pbmc))
GSE255339_pbmc <- AddMetaData(object = GSE255339_pbmc,
                              metadata = GSE255339_RNA_sample$group1,
                              col.name = "group1")
GSE255339_pbmc_counts <- GetAssayData(GSE255339_pbmc,slot = "data")
GSE255339_pbmc_counts_frame <- as.data.frame(GSE255339_pbmc_counts)

RNA_unique_groups <- unique(GSE255339_RNA_sample$group1)

# target_order <- c("E4.5", "E5.5", "E6.5", "E7.5")

RNA_AZA_DAC_marker_data <- read.csv("plot/RNA/GSE255339_pbmc_markers_AZA_DAC.csv")
RNA_AZA_Unt_marker_data <- read.csv("plot/RNA/GSE255339_pbmc_markers_AZA.csv")
RNA_DAC_Unt_marker_data <- read.csv("plot/RNA/GSE255339_pbmc_markers_DAC.csv")



#AZA_DAC
RNA_marker_data_E4.5 <- RNA_AZA_DAC_marker_data

RNA_marker_data_E4.5 <- RNA_marker_data_E4.5 %>%
  filter(p_val < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E4.5_counts <- GSE255339_pbmc_counts_frame[RNA_marker_data_E4.5$gene,]

RNA_cells_E4.5 <- GSE255339_RNA_sample$sample[GSE255339_RNA_sample$group1 == "AZA" | GSE255339_RNA_sample$group1 == "DAC"]
valid_cells <- intersect(RNA_cells_E4.5, colnames(RNA_marker_data_E4.5_counts))
sub_mat_E4.5 <- RNA_marker_data_E4.5_counts[, valid_cells, drop = FALSE]
avg_E4.5 <- as.data.frame(rowMeans(sub_mat_E4.5, na.rm = TRUE))
colnames(avg_E4.5) <- "Expression"

E4.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E4.5), Expression = avg_E4.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE255339_RNA_data_gene, by = c("GeneID" = "Geneid"))



#AZA_Unt
RNA_marker_data_E5.5 <- RNA_AZA_Unt_marker_data

RNA_marker_data_E5.5 <- RNA_marker_data_E5.5 %>%
  filter(p_val < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E5.5_counts <- GSE255339_pbmc_counts_frame[RNA_marker_data_E5.5$gene,]

RNA_cells_E5.5 <- GSE255339_RNA_sample$sample[GSE255339_RNA_sample$group1 == "AZA" | GSE255339_RNA_sample$group1 == "Unt"]
valid_cells <- intersect(RNA_cells_E5.5, colnames(RNA_marker_data_E5.5_counts))
sub_mat_E5.5 <- RNA_marker_data_E5.5_counts[, valid_cells, drop = FALSE]
avg_E5.5 <- as.data.frame(rowMeans(sub_mat_E5.5, na.rm = TRUE))
colnames(avg_E5.5) <- "Expression"

E5.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E5.5), Expression = avg_E5.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE255339_RNA_data_gene, by = c("GeneID" = "Geneid"))



#DAC_Unt
RNA_marker_data_E6.5 <- RNA_DAC_Unt_marker_data

RNA_marker_data_E6.5 <- RNA_marker_data_E6.5 %>%
  filter(p_val < 0.05 & abs(avg_log2FC) > 1)
RNA_marker_data_E6.5_counts <- GSE255339_pbmc_counts_frame[RNA_marker_data_E6.5$gene,]

RNA_cells_E6.5 <- GSE255339_RNA_sample$sample[GSE255339_RNA_sample$group1 == "DAC" | GSE255339_RNA_sample$group1 == "Unt"]
valid_cells <- intersect(RNA_cells_E6.5, colnames(RNA_marker_data_E6.5_counts))
sub_mat_E6.5 <- RNA_marker_data_E6.5_counts[, valid_cells, drop = FALSE]
avg_E6.5 <- as.data.frame(rowMeans(sub_mat_E6.5, na.rm = TRUE))
colnames(avg_E6.5) <- "Expression"

E6.5_RNA_df_result <- data.frame(GeneID = rownames(avg_E6.5), Expression = avg_E6.5) %>%
  # 使用 left_join 合并对照表
  # by = c("你的ID列" = "对照表的ID列")
  left_join(GSE255339_RNA_data_gene, by = c("GeneID" = "Geneid"))



#genetss2k----
GSE255339_CpG_genetss2k_meth_level_data <- GSE255339_CpG_genetss2k_methlevel
GSE255339_GpC_genetss2k_meth_level_data <- GSE255339_GpC_genetss2k_methlevel

GSE255339_CpG_genetss2k_methlevel <- Read_file_colname(GSE255339_CpG_genetss2k_meth_level_data,"level")
GSE255339_GpC_genetss2k_methlevel <- Read_file_colname(GSE255339_GpC_genetss2k_meth_level_data,"level")

CpG_sample_title <- GSE255339_CpG_Epis_sample %>%
  filter(level %in% colnames(GSE255339_CpG_genetss2k_methlevel))
GpC_sample_title <- GSE255339_GpC_Epis_sample %>%
  filter(level %in% colnames(GSE255339_GpC_genetss2k_methlevel))

colnames(GSE255339_CpG_genetss2k_methlevel) <- CpG_sample_title$sample
colnames(GSE255339_GpC_genetss2k_methlevel) <- GpC_sample_title$sample

CpG_genetss2k_meth_DEG_files <- list.files("site_DEG/genetss2k/",full.names = T,pattern = "\\CpG.csv*")
GpC_genetss2k_meth_DEG_files <- list.files("site_DEG/genetss2k/",full.names = T,pattern = "\\GpC.csv*")

CpG_genetss2k_methlevel_DEG_files <- list.files("level_DEG/genetss2k/",full.names = T,pattern = "\\CpG.csv*")
GpC_genetss2k_methlevel_DEG_files <- list.files("level_DEG/genetss2k/",full.names = T,pattern = "\\GpC.csv*")



#AZA_DAC----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[1])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[1])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[1])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[1])

E4.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E4.5_CpG_DEG_data <- E4.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E4.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E4.5_GpC_DEG_data <- E4.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E4.5_common_data <- unique(intersect(intersect(E4.5_CpG_DEG_data$gene_name.x, E4.5_GpC_DEG_data$gene_name.x),
                                E4.5_RNA_df_result$gene_name))

E4.5_RNA_df_result_common <- E4.5_RNA_df_result %>%
  filter(gene_name %in% E4.5_common_data)
E4.5_RNA_df_result_common_unique <- E4.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E4.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "DAC"]
valid_cells <- intersect(CpG_cells_E4.5, colnames(GSE255339_CpG_genetss2k_methlevel))
sub_mat_E4.5 <- GSE255339_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E4.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "DAC"]
valid_cells <- intersect(GpC_cells_E4.5, colnames(GSE255339_GpC_genetss2k_methlevel))
sub_mat_E4.5 <- GSE255339_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E4.5_CpG_GpC_RNA_data,"a.csv")


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
  labs(title = "GeneTSS2k AZA vs DAC Violin Plot", y = "Value")


#相关性散点图
# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs DAC CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs DAC GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs DAC CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3


#桑基图 / 冲击图
df_cat <- E4.5_CpG_GpC_RNA_data %>%
  mutate(
    RNA_Group = ifelse(RNA_Exp > median(RNA_Exp), "High_Exp", "Low_Exp"),
    CpG_Group = ifelse(CpG_level > median(CpG_level), "High_CpG", "Low_CpG"),
    GpC_Group = ifelse(GpC_level > median(GpC_level), "High_GpC", "Low_GpC")
  ) %>%
  group_by(CpG_Group, GpC_Group, RNA_Group) %>%
  summarise(Freq = n()) # 统计每种组合有多少个基因

# 2. 绘图
ggplot(df_cat,
       aes(y = Freq, axis1 = CpG_Group, axis2 = GpC_Group, axis3 = RNA_Group)) +
  geom_alluvium(aes(fill = RNA_Group), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CpG", "GpC", "RNA"), expand = c(.05, .05)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(title = "Gene Regulation Flow", y = "Number of Genes")


#山峦图
plot_data <- E4.5_CpG_GpC_RNA_data %>%
  tidyr::pivot_longer(cols = c("RNA_Exp", "CpG_level", "GpC_level"),
               names_to = "Omics_Type", values_to = "Value")

ggplot(plot_data, aes(x = Value, y = Omics_Type, fill = Omics_Type)) +
  # scale控制山峰重叠程度
  geom_density_ridges(scale = 1.5, alpha = 0.6, quantile_lines = TRUE) +
  scale_fill_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  theme_ridges() +
  # 因为RNA值很大，建议横坐标分面或者对RNA取log，这里简单演示分面
  facet_wrap(~Omics_Type, scales = "free", nrow = 3) +
  labs(title = "Density Distribution (Ridge Plot)", x = "Value", y = "")


#相关性矩阵图
ggpairs(E4.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k AZA vs DAC Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()


#映射色彩的散点图
ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  # 用点的颜色代表 RNA_Exp
  geom_point(aes(color = RNA_Exp), alpha = 0.8, size = 1.5) +
  # 使用 viridis 色阶，色彩区分度高（黄色高表达，紫色低表达）
  scale_color_viridis_c(option = "C", name = "Log(RNA)") +
  # 添加平滑趋势线（可选），看看整体趋势
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  theme_bw() +
  labs(title = "Methylation Landscape & Gene Expression",
       x = "CpG Methylation Level (%)",
       y = "GpC Level (%)")


#分层箱线图
df_grouped <- E4.5_CpG_GpC_RNA_data %>%
  mutate(Exp_Group = cut(RNA_Exp,
                         breaks = quantile(RNA_Exp, probs = c(0, 0.33, 0.66, 1)),
                         labels = c("Low", "Medium", "High"),
                         include.lowest = TRUE)) %>%
  # 转长数据以便同时画 CpG 和 GpC
  tidyr::pivot_longer(cols = c("CpG_level", "GpC_level"),
               names_to = "Meth_Type", values_to = "Level")

# 2. 绘图
ggplot(df_grouped, aes(x = Exp_Group, y = Level, fill = Exp_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) + # 加上噪点看分布
  facet_wrap(~Meth_Type, scales = "free_y") + # 分面展示 CpG 和 GpC
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw() +
  # 添加统计检验 (Kruskal-Wallis 或 Wilcoxon)
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Low") +
  labs(title = "Methylation Levels across Expression Groups",
       x = "Gene Expression Group", y = "Methylation Level (%)")


#AZA_Unt----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[2])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[2])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[2])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[2])

E5.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E5.5_CpG_DEG_data <- E5.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E5.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E5.5_GpC_DEG_data <- E5.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E5.5_common_data <- unique(intersect(intersect(E5.5_CpG_DEG_data$gene_name.x, E5.5_GpC_DEG_data$gene_name.x),
                                     E5.5_RNA_df_result$gene_name))

E5.5_RNA_df_result_common <- E5.5_RNA_df_result %>%
  filter(gene_name %in% E5.5_common_data)

E5.5_RNA_df_result_common_unique <- E5.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)



#CpG
CpG_cells_E5.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(CpG_cells_E5.5, colnames(GSE255339_CpG_genetss2k_methlevel))
sub_mat_E5.5 <- GSE255339_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E5.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(GpC_cells_E5.5, colnames(GSE255339_GpC_genetss2k_methlevel))
sub_mat_E5.5 <- GSE255339_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E5.5_CpG_GpC_RNA_data,"a.csv")

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
  labs(title = "GeneTSS2k AZA vs Unt Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs Unt CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs Unt GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k AZA vs Unt CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3


#相关性矩阵图
ggpairs(E5.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k AZA vs Unt Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()



#DAC_Unt----
CpG_meth_DEG_data <- read.csv(CpG_genetss2k_meth_DEG_files[3])
GpC_meth_DEG_data <- read.csv(GpC_genetss2k_meth_DEG_files[3])

CpG_methlevel_DEG_data <- read.csv(CpG_genetss2k_methlevel_DEG_files[3])
GpC_methlevel_DEG_data <- read.csv(GpC_genetss2k_methlevel_DEG_files[3])

E6.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E6.5_CpG_DEG_data <- E6.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E6.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E6.5_GpC_DEG_data <- E6.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E6.5_common_data <- unique(intersect(intersect(E6.5_CpG_DEG_data$gene_name.x, E6.5_GpC_DEG_data$gene_name.x),
                                     E6.5_RNA_df_result$gene_name))

E6.5_RNA_df_result_common <- E6.5_RNA_df_result %>%
  filter(gene_name %in% E6.5_common_data)
E6.5_RNA_df_result_common_unique <- E6.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E6.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "DAC" |GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(CpG_cells_E6.5, colnames(GSE255339_CpG_genetss2k_methlevel))
sub_mat_E6.5 <- GSE255339_CpG_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E6.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "DAC" |GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(GpC_cells_E6.5, colnames(GSE255339_GpC_genetss2k_methlevel))
sub_mat_E6.5 <- GSE255339_GpC_genetss2k_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E6.5_CpG_GpC_RNA_data,"a.csv")

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
  labs(title = "GeneTSS2k DAC vs Unt Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "GeneTSS2k DAC vs Unt CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k DAC vs Unt GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "GeneTSS2k DAC vs Unt CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E6.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "GeneTSS2k DAC vs Unt Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()







#enhancer----
GSE255339_CpG_enhancer_meth_level_data <- GSE255339_CpG_enhancer_methlevel
GSE255339_GpC_enhancer_meth_level_data <- GSE255339_GpC_enhancer_methlevel

GSE255339_CpG_enhancer_methlevel <- Read_file_colname(GSE255339_CpG_enhancer_meth_level_data,"level")
GSE255339_GpC_enhancer_methlevel <- Read_file_colname(GSE255339_GpC_enhancer_meth_level_data,"level")

CpG_sample_title <- GSE255339_CpG_Epis_sample %>%
  filter(level %in% colnames(GSE255339_CpG_enhancer_methlevel))
GpC_sample_title <- GSE255339_GpC_Epis_sample %>%
  filter(level %in% colnames(GSE255339_GpC_enhancer_methlevel))

colnames(GSE255339_CpG_enhancer_methlevel) <- CpG_sample_title$sample
colnames(GSE255339_GpC_enhancer_methlevel) <- GpC_sample_title$sample

CpG_enhancer_meth_DEG_files <- list.files("site_DEG/enhancer/",full.names = T,pattern = "\\CpG.csv*")
GpC_enhancer_meth_DEG_files <- list.files("site_DEG/enhancer/",full.names = T,pattern = "\\GpC.csv*")

CpG_enhancer_methlevel_DEG_files <- list.files("level_DEG/enhancer/",full.names = T,pattern = "\\CpG.csv*")
GpC_enhancer_methlevel_DEG_files <- list.files("level_DEG/enhancer/",full.names = T,pattern = "\\GpC.csv*")


#AZA_DAC----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[1])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[1])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[1])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[1])

E4.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E4.5_CpG_DEG_data <- E4.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E4.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E4.5_GpC_DEG_data <- E4.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E4.5_common_data <- unique(intersect(intersect(E4.5_CpG_DEG_data$geneone.x, E4.5_GpC_DEG_data$geneone.x),
                                     E4.5_RNA_df_result$gene_name))

E4.5_RNA_df_result_common <- E4.5_RNA_df_result %>%
  filter(gene_name %in% E4.5_common_data)
E4.5_RNA_df_result_common_unique <- E4.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E4.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "DAC"]
valid_cells <- intersect(CpG_cells_E4.5, colnames(GSE255339_CpG_enhancer_methlevel))
sub_mat_E4.5 <- GSE255339_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E4.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "DAC"]
valid_cells <- intersect(GpC_cells_E4.5, colnames(GSE255339_GpC_enhancer_methlevel))
sub_mat_E4.5 <- GSE255339_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E4.5_CpG_GpC_RNA_data,"a.csv")


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
  labs(title = "Enhancer AZA vs DAC Violin Plot", y = "Value")


#相关性散点图
# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer AZA vs DAC CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer AZA vs DAC GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E4.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer AZA vs DAC CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E4.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer AZA vs DAC Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()



#AZA_Unt----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[2])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[2])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[2])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[2])

E5.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E5.5_CpG_DEG_data <- E5.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E5.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E5.5_GpC_DEG_data <- E5.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E5.5_common_data <- unique(intersect(intersect(E5.5_CpG_DEG_data$geneone.x, E5.5_GpC_DEG_data$geneone.x),
                                     E5.5_RNA_df_result$gene_name))

E5.5_RNA_df_result_common <- E5.5_RNA_df_result %>%
  filter(gene_name %in% E5.5_common_data)
E5.5_RNA_df_result_common_unique <- E5.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E5.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(CpG_cells_E5.5, colnames(GSE255339_CpG_enhancer_methlevel))
sub_mat_E5.5 <- GSE255339_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E5.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "AZA" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(GpC_cells_E5.5, colnames(GSE255339_GpC_enhancer_methlevel))
sub_mat_E5.5 <- GSE255339_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E5.5_CpG_GpC_RNA_data,"a.csv")

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
  labs(title = "Enhancer AZA vs Unt Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer AZA vs Unt CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer AZA vs Unt GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E5.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer AZA vs Unt CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E5.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer AZA vs Unt Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()




#DAC_Unt----
CpG_meth_DEG_data <- read.csv(CpG_enhancer_meth_DEG_files[3])
GpC_meth_DEG_data <- read.csv(GpC_enhancer_meth_DEG_files[3])

CpG_methlevel_DEG_data <- read.csv(CpG_enhancer_methlevel_DEG_files[3])
GpC_methlevel_DEG_data <- read.csv(GpC_enhancer_methlevel_DEG_files[3])

E6.5_CpG_DEG_data <- merge(CpG_methlevel_DEG_data,CpG_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E6.5_CpG_DEG_data <- E6.5_CpG_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)
E6.5_GpC_DEG_data <- merge(GpC_methlevel_DEG_data,GpC_meth_DEG_data,by.x = "chr",by.y = "chr",all.x = T)
E6.5_GpC_DEG_data <- E6.5_GpC_DEG_data %>%
  # mutate(numeric_logFC = parse_number(logFC)) %>%
  filter(!is.na(logFC)) %>%
  filter(p.value < 0.05) %>%
  slice_max(logFC, n = 25000)

E6.5_common_data <- unique(intersect(intersect(E6.5_CpG_DEG_data$geneone.x, E6.5_GpC_DEG_data$geneone.x),
                                     E6.5_RNA_df_result$gene_name))

E6.5_RNA_df_result_common <- E6.5_RNA_df_result %>%
  filter(gene_name %in% E6.5_common_data)
E6.5_RNA_df_result_common_unique <- E6.5_RNA_df_result_common %>%
  distinct(gene_name, .keep_all = TRUE)


#CpG
CpG_cells_E6.5 <- GSE255339_CpG_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "DAC" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(CpG_cells_E6.5, colnames(GSE255339_CpG_enhancer_methlevel))
sub_mat_E6.5 <- GSE255339_CpG_enhancer_methlevel[, valid_cells, drop = FALSE]
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
GpC_cells_E6.5 <- GSE255339_GpC_Epis_sample$sample[GSE255339_CpG_Epis_sample$group1 == "DAC" | GSE255339_CpG_Epis_sample$group1 == "Unt"]
valid_cells <- intersect(GpC_cells_E6.5, colnames(GSE255339_GpC_enhancer_methlevel))
sub_mat_E6.5 <- GSE255339_GpC_enhancer_methlevel[, valid_cells, drop = FALSE]
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
# write.csv(E6.5_CpG_GpC_RNA_data,"a.csv")

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
  labs(title = "Enhancer DAC vs Unt Violin Plot", y = "Value")


# 1. CpG vs RNA (查看启动子甲基化是否抑制表达)
p1 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#2E9FDF") +
  geom_smooth(method = "lm", color = "black") + # 添加回归线
  stat_cor(method = "spearman") + # 添加相关系数
  theme_bw() +
  labs(title = "Enhancer DAC vs Unt CpG vs RNA", x = "CpG Level (%)", y = "log(RNA Expression)")

# 2. GpC vs RNA
p2 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = GpC_level, y = RNA_Exp)) +
  geom_point(alpha = 0.6, color = "#FC4E07") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer DAC vs Unt GpC vs RNA", x = "GpC Level (%)", y = "log(RNA Expression)")

# 3. CpG vs GpC (查看两种甲基化的一致性)
p3 <- ggplot(E6.5_CpG_GpC_RNA_data, aes(x = CpG_level, y = GpC_level)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(method = "spearman") +
  theme_bw() +
  labs(title = "Enhancer DAC vs Unt CpG vs GpC", x = "CpG Level (%)", y = "GpC Level (%)")

# 拼图
p1 + p2 + p3

#相关性矩阵图
ggpairs(E6.5_CpG_GpC_RNA_data,
        columns = c("RNA_Exp", "CpG_level", "GpC_level"),
        title = "Enhancer DAC vs Unt Multi-Omics Correlation Matrix",
        # 设置下三角（散点图）的样式：透明度0.5，点大小0.5
        lower = list(continuous = wrap("points", alpha = 0.5, size = 0.5, color = "#2E9FDF")),
        # 设置对角线（密度图）的样式
        diag = list(continuous = wrap("densityDiag", fill = "gray", alpha = 0.5)),
        # 设置上三角（相关系数）的样式
        upper = list(continuous = wrap("cor", size = 5))
) + theme_bw()

