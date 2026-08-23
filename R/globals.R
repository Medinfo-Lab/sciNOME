# Suppress R CMD check notes for non-standard evaluation (NSE) variables
utils::globalVariables(c(
  ".", "..final_cols", "..keep_cols", "..level_cols", "Associated_Regions",
  "Auto_Cluster", "Cluster", "CpG_level", "Dim1", "Dim2", "GeneID", "GeneID.x",
  "GeneName", "Global_Mean", "GpC_level", "Group", "PlotGroup", "Plot_Group",
  "Pseudotime", "RNA_Exp", "RegionID", "Significance", "Value", "avg_log2FC",
  "chr", "chrdata", "chrdata_nochr", "end", "final_gene_id", "final_gene_name",
  "gene", "meth", "p_val_adj", "region_id", "region_idx", "start", "unmeth"
))
