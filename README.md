# sciNOME: Integrative Analysis of Nucleosome Occupancy, Genome Methylation and Expression across Single-Cells

# Installation instructions

***sciNOME**: A Region-Centric Framework for Integrative Analysis of Single-Cell Epigenomic and Transcriptomic Data.* 

***sciNOME*** is an R package for jointly analyzing transcriptomic, DNA methylation, and chromatin accessibility data. The package is designed to process sequencing data from scNOMe-seq and scRNA-seq experiments. Additionally, ***sciNOME*** converts sequencing data into region-based formats for efficient storage and subsequent analysis. The R package incorporates features such as differential analysis, dimensionality reduction analysis, and differential analysis of horizontal sites.

```R
# Development version from GitHub
BiocManager::install("Medinfo-Lab/sciNOME")

# From Bioconductor (once accepted)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sciNOME")
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

## Quick Start

## Loading the package

```
library(sciNOME)
library(dplyr)
```

**Here, we will use the built-in data as an example.**

## Step 1 — load data

```R
# 1.1 Dynamically obtain the absolute path of test data included in the package
rna_path <- system.file("extdata", "RNA_counts.txt", package = "sciNOME")
meta_path <- system.file("extdata", "metadata.csv", package = "sciNOME")
bed_path <- system.file("extdata", "target_promoters.bed", package = "sciNOME")
cpg_dir <- system.file("extdata/cpg_cov", package = "sciNOME")
gpc_dir <- system.file("extdata/gpc_cov", package = "sciNOME")

# 1.2 Read basic expression profiles, metadata, and target interval files
rna_counts <- read.table(rna_path, sep="\t", header=TRUE)
metadata <- read.csv(meta_path)
bed_data <- read.table(bed_path,header = T)
```

## Step 2 — Transcriptome Analysis

```R
# 2.1 Build lightweight RNA analysis objects
RNA_obj <- Build_RNAObject(
  expr_mat = rna_counts,            # Expression Matrix
  meta_data = metadata,             # Grouped Data Frames
  meta_id_col = "RNA_ID",           # The columns in the grouped data are used to match the column names in the matrix
  min_cells = 0,                    # (Optional) Remove genes that are completely non-expressed
  min_features = 0,                 # (Optional) Remove completely empty samples
  project_name = "RNA_object"
)
PlotQC_RNA(RNA_obj,"CellType")

# 2.2 Data Filtering and Standardization
RNA_obj <- ProcessQC_RNA(
  obj = RNA_obj,
  mt_pattern = "^MT-",          # Human mitochondrial prefix
  min_nFeature = 0,             # Express at least 0 genes
  max_nFeature = 6000,          # Up to 6000 genes
  min_nCount = 0,               # At least 0 UMI
  max_nCount = 100000,          # Up to 100000 UMI
  max_mt = 10,                  # The proportion of mitochondria must not exceed 10%
  norm_method = "LogNormalize", # Standardized Methods
  do_scale = TRUE               # Whether to scale
)

# 2.3 Dimensionality Reduction Analysis (PCA)
RNA_obj <- RunDimReduction_RNA(
  obj = RNA_obj,
  method = "PCA",
  n_hvg = 2000,
  layer_name = "data",
  pca_rank = 10,
  umap_neighbors = 100,
  umap_mindist = 2
)
PlotDimRed_RNA(RNA_obj,"PCA","CellType",show_cluster = F)

# 2.4 Differential Expression Analysis (DEA)
RNA_diff <- RunDEA_RNA(
  obj = RNA_obj,
  group_col = "CellType",
  ident_1 = "CellType.A"
  # ident_2 = "CellType.B"
)
PlotVolcano_RNA(RNA_diff,
                fc_cut = 0.5,
                p_cut = 0.05)
```

## Step 3 — Epigenomic Analysis

```R
# 3.1 Aggregate apparent site data based on the given interval
result_matrixCpG <- Aggregate_epiRegions(
  cov_dir = gpc_dir,
  bed_file = bed_path,
  n_cores = 8
)
result_matrixGpC <- Aggregate_epiRegions(
  cov_dir = gpc_dir,
  bed_file = bed_path,
  n_cores = 8
)

# 3.2 Missing Value Quality Control and Filtering
qc_result_matrixCpG <- QC_epiData(
  data = result_matrixCpG,
  top_n_rows = 2000,
  max_col_na_ratio = 0.8
)
qc_result_matrixGpC <- QC_epiData(
  data = result_matrixGpC,
  top_n_rows = 2000,
  max_col_na_ratio = 0.8
)

# 3.3 Extract methylation/open chromatin level matrix for downstream dimensionality reduction
qc_result_matrixCpG_level <- Extract_epiData(qc_result_matrixCpG,".level")
qc_result_matrixGpC_level <- Extract_epiData(qc_result_matrixGpC,".level")

# 3.4 Multidimensional Scaling (MDS) for Phenotypic Data
MDS_df_CpG <- Reduce_epiDims(
  qc_result_matrixCpG_level, # choose GpC data (option)
  metadata,
  "CpG_level","CellType",
  dr_method = "MDS",
  impute_method = "knn"
)
PlotDimRed_Epi(MDS_df_CpG,"CellType","MDS")

# 3.5 Identifying Differentially Methylated Regions (DMRs)
CpG_diff <- Run_Diffanalysis(
  raw_mat = qc_result_matrixCpG,
  meta = metadata,
  group_col = "CellType",
  target_group = "CellType.A",
  # control_group = "Unt",
  col_level = "CpG_level",
  col_meth = "CpG_meth",
  col_nonmeth = "CpG_nonmeth")
# DMRs plot (option)
PlotVolcano_Epi(CpG_diff,th_effect = 0.05,title = "CellType.A vs other CpG Volcano") # Differential region analysis (DMRs)

# 3.6 Identifying Differentially Accessible Chromatin Regions (DARs)
GpC_diff <- Run_Diffanalysis(
  raw_mat = qc_result_matrixGpC,
  meta = metadata,
  group_col = "CellType",
  target_group = "CellType.A",
  # control_group = "Unt",
  col_level = "CpG_level",
  col_meth = "CpG_meth",
  col_nonmeth = "CpG_nonmeth")
# DARs plot (option)
PlotVolcano_Epi(GpC_diff,th_effect = 0.05,title = "CellType.A vs other GpC Volcano") # Differential region analysis (DARs)
```

## Step 4 — Multi-omics integrated analysis

```R
# 4.1 Multi-omics data integration
integrated_obj <- Integrate_MultiOmics(
  mode = "tri",                        	# Integration mode
  target_group = "ALL",              	# Target group
  meta_df = metadata,                  	# Grouped Data Frames
  group_col = "CellType",              	# Columns that define groups
  region_df = bed_data,                	# Regional Information Table
  rna_obj = RNA_obj,              	   	# Transcriptomic objects
  rna_id_col = "RNA_ID",               	# The column name containing the RNA groups
  # rna_diff_df = RNA_diff,			   	# RNA differential data (option)
  cpg_mat = qc_result_matrixCpG_level,  # The level data for each CpG sample
  cpg_id_col = "CpG_level",            	# The column name containing the CpG groups CpG
  # cpg_diff_df = CpG_diff,			   	# DNA methylation differential data  (option)
  gpc_mat = qc_result_matrixGpC_level,  # The level data for each GpC sample
  gpc_id_col = "GpC_level"             	# The column name containing the CpG groups GpC
  # gpc_diff_df = GpC_diff,            	# Chromatin accessibility differential data  (option)
)

# 4.2 Draw multi-omics level heatmaps / matrix plots
PlotOmicsMatrix(integrated_obj)

# 4.3 Draw multi-omics scatter correlation plots (e.g., Methylation vs Expression vs Accessibility)
PlotOmicsScatter(integrated_obj)
```

# License

MIT License — see [LICENSE](LICENSE) for details.
