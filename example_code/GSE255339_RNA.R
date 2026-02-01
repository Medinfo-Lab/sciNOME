library(readxl)
library(dplyr)
library(Seurat)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
library(ggVolcano)
library(tidydr)
library(Rgraphviz)
library(GOplot)
library(uwot)
library(ComplexHeatmap)
library(circlize)
library(Rtsne)


load("data/human/RNA_Data.RData")


#object----
GSE255339_pbmc <- CreatRNAObject(counts = GSE255339_RNA_data_counts_filter,
                                     min.cells = 50, min.features = 0, project = "SmartSeq2_project")
# GSE255339_pbmc[["percent.mt"]] <- PercentageFeatureSet(GSE255339_pbmc, pattern = "^mt-")
# VlnPlot(GSE255339_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

GSE255339_RNA_sample_seurat_group <- GSE255339_RNA_sample_data %>%
  filter(sample %in% colnames(GSE255339_pbmc))
GSE255339_pbmc <- AddMetaData(object = GSE255339_pbmc,
                              metadata = GSE255339_RNA_sample_data$group1,
                              col.name = "group1")

rownames(GSE255339_pbmc) <- gsub("-", "_", rownames(GSE255339_pbmc))


#pca umap----
normalized_matrix <- GetAssayData(GSE255339_pbmc, slot = "data")
normalized_matrix_frame <- as.data.frame(normalized_matrix)
normalized_matrix_frame_t <- t(normalized_matrix_frame)

#pca
pca_res <- prcomp(normalized_matrix_frame_t, center = T, scale. = T, rank. = 30)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
pca_embedding <- as.data.frame(pca_res$x)

pca_embedding$Group <- GSE255339_RNA_sample_seurat_group$group1

ggplot(pca_embedding, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  theme_classic()+
  labs(x = pc1_lab,y = pc2_lab,title = "RNA Developmental Stage PCA")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("Unt","AZA", "DAC"),
                     values = c("#cccccc", "#E64B35", "#3C5488"))


#tSNE
tsne_out <- Rtsne(normalized_matrix_frame_t,
                  dims = 2,
                  perplexity = 30,
                  verbose = TRUE,
                  max_iter = 1000,
                  check_duplicates = F,
                  pca = TRUE)

tsne_plot_data <- data.frame(
  tSNE_1 = tsne_out$Y[, 1],
  tSNE_2 = tsne_out$Y[, 2],
  Group = GSE255339_RNA_sample_seurat_group$group1  # 将分组信息合并进来用于上色
)

ggplot(tsne_plot_data, aes(x = tSNE_1, y = tSNE_2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "RNA Developmental Stage tSNE")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("Unt","AZA", "DAC"),
                     values = c("#868686", "#E64B35", "#3C5488"))


#umap
umap_result <- umap(normalized_matrix_frame_t,n_neighbors = 100, min_dist = 0.5)
# umap_result <- umap(pca_embedding,n_neighbors = 100, min_dist = 0.5)

umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Group <- GSE255339_RNA_sample_seurat_group$group1

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "RNA Developmental Stage UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("Unt","AZA", "DAC"),
                     values = c("#868686", "#E64B35", "#3C5488"))


#DEG----
Idents(GSE255339_pbmc) <- GSE255339_pbmc$group1

GSE255339_pbmc_markers_AZA <- FindMarkers(GSE255339_pbmc,
                                          ident.1 = "AZA",      # 感兴趣的簇
                                          ident.2 = "Unt",      # 对照簇
)
GSE255339_pbmc_markers_AZA$gene <- rownames(GSE255339_pbmc_markers_AZA)

GSE255339_pbmc_markers_DAC <- FindMarkers(GSE255339_pbmc,
                                          ident.1 = "DAC",      # 感兴趣的簇
                                          ident.2 = "Unt",      # 对照簇
)
GSE255339_pbmc_markers_DAC$gene <- rownames(GSE255339_pbmc_markers_DAC)

GSE255339_pbmc_markers_AZA_DAC <- FindMarkers(GSE255339_pbmc,
                                              ident.1 = "DAC",      # 感兴趣的簇
                                              ident.2 = "AZA",      # 对照簇
)
GSE255339_pbmc_markers_AZA_DAC$gene <- rownames(GSE255339_pbmc_markers_AZA_DAC)


GSE255339_pbmc_markers_result  <- merge(
  x = GSE255339_pbmc_markers_AZA_DAC,
  y = GSE255339_RNA_data_gene,
  by.x = "gene",
  by.y = "Geneid",
  all.x = TRUE
)

GSE255339_pbmc_markers_P <- GSE255339_pbmc_markers_result %>%
  filter(p_val < 0.05)

GSE255339_pbmc_markers_P_logfc <- GSE255339_pbmc_markers_P
  # group_by(cluster) %>%
  # slice_max(n = 500, order_by = avg_log2FC)


GSE255339_pbmc_markers_volcano <- GSE255339_pbmc_markers_P_logfc %>%
  filter(p_val < 0.01)

volcano_data <- GSE255339_pbmc_markers_volcano %>%
  mutate(
    gene = GSE255339_pbmc_markers_volcano$gene_name,  # 确保有基因名列
    log_padj = -log10(p_val),  # 转换校正p值
    direction = case_when(         # 标记上下调
      avg_log2FC > 1 ~ "Up",
      avg_log2FC < -1 ~ "Down",
      TRUE ~ "Normal"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_padj))

top_genes_up <- volcano_data %>%
  filter(direction == "Up") %>%
  arrange(p_val) %>%
  head(5)

top_genes_down <- volcano_data %>%
  filter(direction == "Down") %>%
  arrange(p_val) %>%
  head(5)


ggplot(volcano_data, aes(x = avg_log2FC, y = log_padj, color = direction)) +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=c("#00AFBB", "#999999", "#FC4E07"))+
  # scale_color_manual(values=c("#007BC3FF", "#d2dae2", "#C70E7BFF"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10(P-value)")+
  theme_bw()+
  ggtitle("AZA vs DAC DEG Volcano")+
  # 图例
  theme(legend.position="right",
        # legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 12),
        # axis.text.y = element_text(size = 12)
  )+
  xlim(-8, 8)+
# ylim(0, 40)+
  geom_label_repel(data = top_genes_up, aes(label = gene_name),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)+             # 设置标签重叠的最大次数
  geom_label_repel(data = top_genes_down, aes(label = gene_name),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)               # 设置标签重叠的最大次数




#KEGG GO----
enrichment_genes_UP <- volcano_data %>%
  filter(direction == "Up")

enrichment_genes_DOWN <- volcano_data %>%
  filter(direction == "Down")


gene_ids <- bitr(enrichment_genes_DOWN$gene_name,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb="org.Hs.eg.db")


#GO
ego <- enrichGO(gene = gene_ids$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.5,
                readable = TRUE)         # 转换ID为基因名


ego_data <- ego@result


dotplot(ego,
        showCategory=10,         # 显示top15条目
        color = "pvalue",
        title="AZA vs DAC GO Enrichment") +
  # scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))+
  scale_color_gradientn(
    colours = c("#E64B35", "#4DBBD5", "#00A087"))




#KEGG
kk <- enrichKEGG(gene = gene_ids$ENTREZID,
                 organism = "hsa",       # 人类'hsa'，小鼠'mmu'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.5)

kegg_data <- kk@result

dotplot(kk,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="KEGG Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))


