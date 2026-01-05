library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)


GRCm38_gtf <- rtracklayer::import("/work/genomes/Gencode/GRCm38/gencode.vM25.annotation.gtf")
GRCm38_gtf_data <- as.data.frame(GRCm38_gtf)

GRCm38_gtf_data_gene <- GRCm38_gtf_data %>%
  filter(type=="gene")
# GRCm38_gtf_data_gene$seqnames
a <- unique(GRCm38_gtf_data_gene,by=c("seqnames","start","end"))

GRCm38_gtf_data_gene_DC <- GRCm38_gtf_data_gene %>%
  filter(strand=="+")
GRCm38_gtf_data_gene_RC <- GRCm38_gtf_data_gene %>%
  filter(strand=="-")

GRCm38_gtf_data_gene_DC_data <- data.frame(chr=GRCm38_gtf_data_gene_DC$seqnames,
                                           TSS=GRCm38_gtf_data_gene_DC$start,
                                           gene_id=GRCm38_gtf_data_gene_DC$gene_id,
                                           gene_name=GRCm38_gtf_data_gene_DC$gene_name)
for (i in 1:nrow(GRCm38_gtf_data_gene_DC_data)) {
  GRCm38_gtf_data_gene_DC_data$start[i] <- GRCm38_gtf_data_gene_DC_data$TSS[i]-2000
  GRCm38_gtf_data_gene_DC_data$end[i] <- GRCm38_gtf_data_gene_DC_data$TSS[i]+2000
}

GRCm38_gtf_data_gene_RC_data <- data.frame(chr=GRCm38_gtf_data_gene_RC$seqnames,
                                           TSS=GRCm38_gtf_data_gene_RC$end,
                                           gene_id=GRCm38_gtf_data_gene_RC$gene_id,
                                           gene_name=GRCm38_gtf_data_gene_RC$gene_name)
for (i in 1:nrow(GRCm38_gtf_data_gene_RC_data)) {
  GRCm38_gtf_data_gene_RC_data$start[i] <- GRCm38_gtf_data_gene_RC_data$TSS[i]-2000
  GRCm38_gtf_data_gene_RC_data$end[i] <- GRCm38_gtf_data_gene_RC_data$TSS[i]+2000
}

GRCm38_gtf_data_gene <- rbind(GRCm38_gtf_data_gene_DC_data,GRCm38_gtf_data_gene_RC_data)
GRCm38_gtf_data_gene_choose <- GRCm38_gtf_data_gene[,-2]

GRCm38_gtf_data_gene_choose_unique <- GRCm38_gtf_data_gene_choose %>%
  distinct(chr, start, end, .keep_all = TRUE)
result <- GRCm38_gtf_data_gene_choose_unique %>%
  group_by(chr) %>%            # 按group_col分组
  arrange(start, .by_group = TRUE)
# write.csv(result,"data/list/GRCm38_GeneTSS2k.csv",row.names = F)


genebody <- read.csv("data/list/GRCm38_Genebody.csv")
genetss2k <- read.csv("data/list/GRCm38_GeneTSS2k.csv")

genebody_result <- genebody %>%
  distinct(chr, start, end, .keep_all = TRUE)
genebody_result <- genebody_result %>%
  group_by(chr) %>%            # 按group_col分组
  arrange(start, .by_group = TRUE)
# write.csv(genebody_result,"data/list/mm10_genebody_order.csv",row.names = F)

genetss2k_result <- genetss2k %>%
  group_by(chr) %>%            # 按group_col分组
  arrange(start, .by_group = TRUE)
# write.csv(genetss2k_result,"data/list/mm10_genetss2k_order.csv",row.names = F)

screen_mm10_enhancer <- read.csv("data/list/SCREEN_mm10-ELS.csv")
screen_mm10_enhancer_result <- screen_mm10_enhancer %>%
  group_by(chr) %>%            # 按group_col分组
  arrange(start, .by_group = TRUE)
# write.csv(screen_mm10_enhancer_result,"data/list/screen_mm10_enhancer_order.csv",row.names = F)


screen_mm10_enhancer_result <- read.csv("data/list/screen_mm10_enhancer_order.csv")
enhancer_ranges <- GRanges(
  seqnames = screen_mm10_enhancer_result$chr,
  ranges = IRanges(start = screen_mm10_enhancer_result$start, end = screen_mm10_enhancer_result$end)
)
# export.bed(enhancer_ranges, "data/list/enhancer_order.bed")


screen_mm10_enhancer_result <- read.csv("data/list/screen_mm10_enhancer_order.csv")
encode_order_genes <- screen_mm10_enhancer_result %>%
  separate(gene,
           into = c("Gene1", "Position1", "Gene2", "Position2"),
           sep = "[(), ]+",
           remove = TRUE,
           fill = "right")
encode_order_genes$Position1 <- as.numeric(encode_order_genes$Position1)
encode_order_genes$Position2 <- as.numeric(encode_order_genes$Position2)

encode_order_genes$geneone <- ifelse(is.na(encode_order_genes$Position2),encode_order_genes$Gene1,
                                     ifelse(abs(encode_order_genes$Position1) < abs(encode_order_genes$Position2),
                                            encode_order_genes$Gene1,
                                            encode_order_genes$Gene2))
# write.csv(encode_order_genes,"data/list/enhancer_order.csv")
