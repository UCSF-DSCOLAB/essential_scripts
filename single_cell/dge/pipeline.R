setwd("/krummellab/data1/almonteloya/essential_scripts/")
source('single_cell/dge/fxns_load__SOURCE_ME.R', chdir = TRUE)
library(Seurat)
packageVersion("Seurat")
expr <- as.matrix(GetAssayData(pbmc_small, layer = 'counts')) 
head(expr)
expr <- as.matrix(pbmc_small@assays$RNA@counts)
head(expr)

expr <- as.matrix(GetAssayData(pbmc_small, slot = 'counts')) 
head(expr)

run_dge(
  counts = expr, metadata = pbmc_small@meta.data,
  method = 'edger', dge_by = 'groups', cell_by = 'RNA_snn_res.0.8',
  case_group = 'g1', reference_group = 'g2',
  cell_targets = valid_clusters,
  min_frac = 0.3
)

