# sciNOME: Integrative Analysis of Nucleosome Occupancy, Genome Methylation and Expression across Single-Cells

# Installation instructions

***sciNOME**: A Region-Centric Framework for Integrative Analysis of Single-Cell Epigenomic and Transcriptomic Data.* 

***sciNOME*** is an R package for jointly analyzing transcriptomic, DNA methylation, and chromatin accessibility data. The package is designed to process sequencing data from scNOMe-seq and scRNA-seq experiments. Additionally, ***sciNOME*** converts sequencing data into region-based formats for efficient storage and subsequent analysis. The R package incorporates features such as differential analysis, dimensionality reduction analysis, and differential analysis of horizontal sites.

```R
# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sciNOME")

# Development version from GitHub
BiocManager::install("Medinfo-Lab/sciNOME")
```

If you want to install offline, please follow the steps below:

- Step 1: Download the package file

Download the latest ***sciNOME*** release from our [GitHub Releases](https://github.com/Medinfo-Lab/sciNOME/releases) page to your local machine.

- Step 2: Pre-install Dependencies

```R
# Install prerequisite packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

deps <- c("data.table", "dplyr", "ggplot2", "ggpubr", "ggridges", "patchwork", 
          "Matrix", "irlba", "Rtsne", "uwot", "parallel", "pbapply", "igraph", 
          "FNN", "NMF", "princurve", "tidyr", "ggrepel", "GenomicRanges", 
          "IRanges", "S4Vectors", "impute")
```

- Step 3: Install the Local Package

```R
install.packages("/path/to/your/download/sciNOME_latest.tar.gz", repos = NULL, type = "source")
```

# Usage

**The Workflow:**

![](https://imgur.com/qXWGDfb.png)

## How to use sciNOME

```R

library(sciNOME)

# Step 1: Enter the path to the apparent group coverage file and the path to the region file, and read the transcriptome expression matrix
cov_directoryCpG <- "data/cov/CpG/"
cov_directoryGpC <- "data/cov/GpC/"
bed_path <- "data/region/bed.csv"
region_data <- read.csv("data/region/bed.csv")
exper <- read.csv("data/RNA_counts.csv")
sample_data_all <- read.csv("data/sample_data.csv")

# Step 2: Generate aggregated apparent genomic data, and create transcriptomic objects
result_matrixCpG <- aggregate_epi_regions(
  cov_dir = cov_directoryCpG,
  bed_file = bed_path,
  n_cores = 8
)
result_matrixGpC <- aggregate_epi_regions(
  cov_dir = cov_directoryGpC,
  bed_file = bed_path,
  n_cores = 8
)

RNA_obj <- Build_RNAObject(
  expr_mat = exper,        			# Expression Matrix
  meta_data = sample_data_all,      # Grouped Data Frames
  meta_id_col = "Sample_Geo",       # The columns in the grouped data are used to match the column names in the matrix
  min_cells = 5,                    # (Optional) Remove genes that are completely non-expressed
  min_features = 200,               # (Optional) Remove completely empty samples
  project_name = "RNA_object"
)

# Step 3: Filter the epigenomic aggregation data and the transcriptomic objects
qc_result_matrixCpG <- QC_epiData(
  data = result_matrixCpG,
  top_n_rows = 5000,
  max_col_na_ratio = 0.8
)
qc_result_matrixGpC <- QC_epiData(
  data = result_matrixGpC,
  top_n_rows = 5000,
  max_col_na_ratio = 0.8
)

PlotQC_RNA(hg19_RNA_obj,"group")
RNA_obj <- ProcessQC_RNA(
  obj = RNA_obj,
  mt_pattern = "^MT-",          # Human mitochondrial prefix
  min_nCount = 100,             # At least 100 UMI
  max_nCount = 8000000,         # Up to 8,000,000 UMI
  min_nFeature = 5000,          # Express at least 5,000 genes
  max_nFeature = 12500,         # Up to 12,500 genes
  max_mt = 10,                  # The proportion of mitochondria must not exceed 10%
  norm_method = "LogNormalize", # Standardized Methods
  do_scale = TRUE               # Whether to scale
)

# RNA Dimensionality reduction analysis option (option)
hg19_RNA_obj <- RunDimReduction_RNA(
  obj = hg19_RNA_obj,
  method = "PCA",
  layer_name = "data",
  pca_rank = 50,
  umap_neighbors = 100,
  umap_mindist = 2
)
PlotDimRed_RNA(hg19_RNA_obj,"PCA","group",show_cluster = F)

# Step 3: Extract the level data for each sample
qc_result_matrixCpG_level <- Extract_epiData(qc_result_matrixCpG,".level")
qc_result_matrixGpC_level <- Extract_epiData(qc_result_matrixGpC,".level")

# Regional dimensionality reduction analysis (option)
MDS_df_CpG <- Reduce_epiDims(
  qc_result_matrixCpG_level, # choose GpC data (option)
  sample_data_all,
  "level","group",
  dr_method = "MDS",
  impute_method = "knn"
)
PlotDimRed_Epi(MDS_df_CpG,"group","MDS")

# Step 4: RNA, DNA Methylation, Chromatin Accessibility difference analysis
RNA_diff <- RunDEA_RNA(
  obj = mm10_RNA_obj,
  group_col = "group",
  ident_1 = "DAC"       
  # ident_2 = "Unt"
)

# Differential gene plot (option)
RNA_diff_choose <- RNA_diff %>%
  filter(p_val < 0.05 & abs(avg_log2FC) > 1)
PlotVolcano_RNA(RNA_diff_choose,
                fc_cut = 2,
                p_cut = 0.005)

CpG_diff <- Run_Diffanalysis(
    raw_mat = qc_result_matrixCpG,
    meta = sample_data_epiCpG,
    group_col = "group",
    target_group = "AZA",
    # control_group = "Unt",
    col_level = "level",
    col_meth = "meth",
    col_nonmeth = "nonmeth")
# DMRs plot (option)
CpG_diff_choose <- CpG_diff %>%
  filter(P.Value < 0.05)
PlotVolcano_Epi(CpG_diff_choose,th_effect = 0.05,title = "DAC vs other CpG Volcano") # Differential region analysis (DMRs)

GpC_diff <- Run_Diffanalysis(
    raw_mat = qc_result_matrixGpC,
    meta = sample_data_all,
    group_col = "group",
    target_group = "DAC",
    # control_group = "Unt",
    col_level = "level",
    col_meth = "meth",
    col_nonmeth = "nonmeth")
# DARs plot (option)
GpC_diff_choose <- GpC_diff %>%
  filter(P.Value < 0.05)
PlotVolcano_Epi(GpC_diff_choose,th_effect = 0.05,title = "DAC vs other GpC Volcano") # Differential region analysis (DARs)

# Step 5: Multi-omics data integration
result_df <- Integrate_MultiOmics(
  mode = "tri",                        # Integration mode
  target_group = "DAC",                # Target group
  meta_df = sample_data_all,           # Grouped Data Frames
  group_col = "group",                 # Columns that define groups
  region_df = region_data,             # Regional Information Table
  rna_obj = RNA_obj,              	   # Transcriptomic objects
  rna_id_col = "Sample_Geo_RNA",       # The column name containing the RNA groups
  # rna_diff_df = RNA_diff,			   # RNA differential data (option)
  cpg_mat = qc_result_matrixCpG_level, # The level data for each CpG sample
  cpg_id_col = "CpG_level",            # The column name containing the CpG groups CpG
  # cpg_diff_df = CpG_diff,			   # DNA methylation differential data  (option)
  gpc_mat = qc_result_matrixGpC_level, # The level data for each GpC sample
  gpc_id_col = "GpC_level"             # The column name containing the CpG groups GpC
  # gpc_diff_df = GpC_diff,            # Chromatin accessibility differential data  (option)
)

# Step 6: Creating multi-omics statistical plots
PlotOmicsMatrix(result_df)
PlotOmicsScatter(result_df)
```
