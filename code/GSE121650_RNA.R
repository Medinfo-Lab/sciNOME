library(readxl)
library(dplyr)
library(Seurat)
library(ggplot2)
library(org.Mm.eg.db)
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



load("GSE232650_RNA_Data.RData")


#object----
GSE121650_pbmc <- CreateSeuratObject(counts = GSE121650_RNA_data_counts_filter,
                                     min.cells = 50, min.features = 0, project = "SmartSeq2_project")
# GSE121650_pbmc[["percent.mt"]] <- PercentageFeatureSet(GSE121650_pbmc, pattern = "^mt-")
# VlnPlot(GSE121650_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

GSE121650_RNA_sample_seurat_group <- GSE121650_RNA_sample_data %>%
  filter(Run %in% colnames(GSE121650_pbmc))
GSE121650_pbmc <- AddMetaData(object = GSE121650_pbmc,
                              metadata = GSE121650_RNA_sample_data$Developmental_stage,
                              col.name = "group1")

GSE121650_pbmc <- NormalizeData(GSE121650_pbmc, normalization.method = "LogNormalize")
GSE121650_pbmc <- FindVariableFeatures(GSE121650_pbmc, selection.method = "vst", nfeatures = 2000)
GSE121650_pbmc_variable_data <- VariableFeatures(GSE121650_pbmc)
GSE121650_pbmc <- ScaleData(GSE121650_pbmc,features = rownames(GSE121650_pbmc))
GSE121650_pbmc <- RunPCA(GSE121650_pbmc, features = rownames(GSE121650_pbmc))
GSE121650_pbmc <- FindNeighbors(GSE121650_pbmc)
GSE121650_pbmc <- FindClusters(GSE121650_pbmc,resolution = 0.5)
GSE121650_pbmc <- RunUMAP(GSE121650_pbmc, dims = 1:20)
# GSE121650_pbmc <- RunUMAP(GSE121650_pbmc,dims = 1:50,
#                           min.dist = 1,n.neighbors = 100)


#pca umap----
normalized_matrix <- GetAssayData(GSE121650_pbmc, slot = "data")
normalized_matrix_frame <- as.data.frame(normalized_matrix)
normalized_matrix_frame_t <- t(normalized_matrix_frame)

#pca
pca_res <- prcomp(normalized_matrix_frame_t, center = T, scale. = T, rank. = 30)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
pca_embedding <- as.data.frame(pca_res$x)

pca_embedding$Group <- GSE121650_RNA_sample_seurat_group$Developmental_stage

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
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))



#umap
normalized_scale_matrix <- GetAssayData(GSE121650_pbmc, slot = "data")
# umap_result <- umap(t(normalized_scale_matrix),n_neighbors = 100, min_dist = 1, metric = "cosine")
umap_result <- umap(pca_embedding,n_neighbors = 100, min_dist = 1, metric = "cosine")

umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Group <- GSE121650_RNA_sample_seurat_group$Developmental_stage

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "RNA Developmental Stage UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=18,face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
        # axis.text.x = element_text(size = 13),
        # axis.text.y = element_text(size = 13)
  )+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5","E7.5"),
                     values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))
# write.csv(umap_df,"plot/RNA/umap.csv")


#DEG----
Idents(GSE121650_pbmc) <- GSE121650_pbmc$group1
# Idents(GSE121650_pbmc) <- GSE121650_pbmc$seurat_clusters
GSE121650_pbmc_markers <- FindAllMarkers(GSE121650_pbmc)

GSE121650_pbmc_markers_result  <- merge(
  x = GSE121650_pbmc_markers,
  y = GSE121650_RNA_data_gene,
  by.x = "gene",
  by.y = "Geneid",
  all.x = TRUE
)
# write.csv(GSE121650_pbmc_markers_result,"data/RNA/GSE121650_markers_result.csv")

GSE121650_pbmc_markers_P <- GSE121650_pbmc_markers_result %>%
  filter(p_val_adj < 0.01)

GSE121650_pbmc_markers_P_logfc <- GSE121650_pbmc_markers_P %>%
  group_by(cluster) %>%
  slice_max(n = 200, order_by = avg_log2FC)

DoHeatmap(GSE121650_pbmc,features = GSE121650_pbmc_markers_P_logfc$gene, angle = 20,
          group.colors = c("E4.5"="#E69F00","E5.5"="#56B4E9","E6.5"="#009E73","E7.5"="#CC79A7"),
          label = T,slot = "scale.data")+
  scale_fill_gradient2(
    low = "#0571B0",
    mid = "#F7F7F7",
    high = "#CA0020",
    midpoint = 0)+
  theme(
    axis.text.y = element_blank(),  # 隐藏Y轴文本（基因名）
    axis.ticks.y = element_blank()   # 可选：隐藏Y轴刻度线
  )

mat <- GetAssayData(GSE121650_pbmc, layer = "scale.data")
mat <- mat[GSE121650_pbmc_markers_P_logfc$gene, ]
cluster_info <- Idents(GSE121650_pbmc)
col_fun = scales::hue_pal()(length(levels(cluster_info)))
names(col_fun) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(
  Cluster = cluster_info,
  # Group = group_info, # 如果想加第二条注释
  col = list(Cluster = col_fun), # 指定颜色
  show_annotation_name = TRUE
)
Heatmap(
  mat,
  name = "Expression", # 图例的名字

  # --- 颜色设置 ---
  # 经典的紫-黑-黄配色 (Seurat风格)
  col = colorRamp2(c(-2, 0, 2), c("#0571B0", "#F7F7F7", "#CA0020")),
  # 或者红-白-蓝: colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),

  # --- 行列设置 ---
  cluster_rows = FALSE,  # 是否聚类行（基因），通常 Marker 图不需要聚类，按 Top10 顺序排列
  cluster_columns = FALSE, # 是否聚类列（细胞），通常我们按 Cluster 分割，内部不聚类
  show_column_names = FALSE, # 细胞太多，不显示细胞名
  show_row_names = F,     # 显示基因名

  # --- 分割与注释 ---
  top_annotation = top_anno, # 加上刚才做的顶部注释条
  column_split = cluster_info, # 核心参数：按 Cluster 切分热图
  # row_split = top10$cluster,   # (可选) 按 Cluster 切分基因行，看起来更清晰

  # --- 美化细节 ---
  # row_names_gp = gpar(fontsize = 8), # 基因名字体大小
  column_title = NULL, # 去掉列的大标题
  use_raster = TRUE # 如果细胞很多，开启光栅化渲染，速度更快
)




volcano_data <- GSE121650_pbmc_markers_P %>%
  mutate(
    gene = GSE121650_pbmc_markers_P$gene_name,  # 确保有基因名列
    log_padj = -log10(p_val_adj),  # 转换校正p值
    direction = case_when(         # 标记上下调
      avg_log2FC > 1 & p_val_adj < 0.01 ~ "Up",
      avg_log2FC < -1 & p_val_adj < 0.01 ~ "Down",
      TRUE ~ "Normal"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_padj))

top_genes_up <- volcano_data %>%
  filter(direction == "Up") %>%
  arrange(p_val_adj) %>%
  head(5)

top_genes_down <- volcano_data %>%
  filter(direction == "Down") %>%
  arrange(p_val_adj) %>%
  head(5)


ggplot(volcano_data, aes(x = avg_log2FC, y = log_padj, color = direction)) +
  geom_point(alpha=0.9, size=2) +
  scale_color_manual(values=c("#00AFBB", "#999999", "#FC4E07"))+
  # scale_color_manual(values=c("#007BC3FF", "#d2dae2", "#C70E7BFF"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10(P-value)")+
  theme_bw()+
  # ggtitle("RNA DEG Volcano")+
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
  xlim(-10, 10)+
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


gene_ids <- bitr(enrichment_genes_UP$gene_name,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb="org.Mm.eg.db")
entrez_ids <- gene_ids$ENTREZID


#GO
ego <- enrichGO(gene = gene_ids$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.5,
                readable = TRUE)         # 转换ID为基因名

#提取并简化结果（避免重复）
# go_simplified <- simplify(
#   ego,
#   cutoff = 0.1,           # 相似性阈值
#   by = "p.adjust",        # 按调整后p值筛选
#   select_fun = min
# )

ego_data <- ego@result

ego_description <- c("mature B cell differentiation involved in immune response",
                     "regulation of B cell differentiation",
                     "mature B cell differentiation",
                     "B cell differentiation",
                     "marginal zone B cell differentiation",
                     "positive regulation of B cell differentiation",
                     "type B pancreatic cell development",
                     "stem cell development",
                     "stem cell differentiation",
                     "regulation of stem cell differentiation"
                     )

dotplot(ego,
        showCategory=ego_description,         # 显示top15条目
        color = "pvalue",
        title="GO Enrichment") +
  # scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))+
  scale_color_gradientn(
    colours = c("#E64B35", "#4DBBD5", "#00A087"))




#KEGG
kk <- enrichKEGG(gene = gene_ids$ENTREZID,
                 organism = "mmu",       # 人类'hsa'，小鼠'mmu'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

kegg_data <- kk@result

kegg_description <- c("mTOR signaling pathway - Mus musculus (house mouse)",
                     "Signaling pathways regulating pluripotency of stem cells - Mus musculus (house mouse)",
                     "B cell receptor signaling pathway - Mus musculus (house mouse)"
)

dotplot(kk,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="KEGG Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))


