# Code for running WNN 
# Usage:
#   RScript run_wnn.R ${in_dir_harmony} ${in_dir_rpca} ${out_dir_wnn}
#
# TODO:
# [ ] Missing options:
#     - npcs
#     - prune.SNN
# [ ] add additional print messages

library(Seurat)
library(tidyverse)
library(cowplot)
library(dittoSeq)

args = commandArgs(trailingOnly=T)
IN_DIR_HARMONY= args[1]
IN_DIR_RPCA=args[2]
OUT_DIR_WNN = args[3]

# you may want to adjust the number of PCs
npcs.rna = 30
npcs.adt = 18

# prune.SNN is smaller than the default (1/15) to avoid tiny clusters
# you may need to decrease this if you dataset is large and you get tiny clusters
# you can also try setting it back to the default
prune.SNN.param = 1/20 


# load the data
load(file=sprintf('%s/merged_processed.RData', IN_DIR_HARMONY)) 
sobj = merged_data
rm(merged_data); gc()

integ_data = readRDS(file = sprintf("%s/integ_ADT.rds", IN_DIR_RPCA))
olap_ids = intersect(rownames(sobj@meta.data), rownames(integ_data@meta.data))
sobj = subset(sobj, cells=olap_ids)
integ_data = subset(integ_data, cells=olap_ids)

# add the ADT
DefaultAssay(sobj) = "RNA"
adt_assay = integ_data@assays[["ADT"]] #sobj[['ADT']] # save this for later, remove for integration
# NOTE THAT sobj[["ADT"]] is only correct if the `processed` and not `filtered` .RDS was used for harmony
sobj[['ADT']] = NULL
sobj[["integrated.ADT"]] = integ_data@assays[["integrated.ADT"]]
DefaultAssay(sobj) = "integrated.ADT"

sobj[["a.pca"]] = integ_data@reductions$pca
sobj[["a.umap"]] = integ_data@reductions$umap

# run WNN on harmony + a.harmony reduction
sobj = FindMultiModalNeighbors(
  sobj, reduction.list=list("harmony", "a.pca"),
  dims.list=list(1:npcs.rna, 1:npcs.adt), 
  prune.SNN = prune.SNN.param, 
  modality.weight.name="RNA.weight"
)
sobj = RunUMAP(sobj, nn.name = "weighted.nn", 
               reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

for (res in c(0.8, 1.0, 1.2, 1.4, 1.6)){
  
  sobj = FindClusters(sobj, graph.name = "wsnn", 
                      algorithm = 4 , method="igraph", 
                      resolution = res, 
                      verbose = FALSE)
  sobj@meta.data[[paste0('wsnn_res.', res)]] <- sobj@meta.data$seurat_clusters
  
  png(filename=sprintf('%s/wsnn_%s.png', OUT_DIR_WNN, res), 
      width = 5, height = 5, units = "in", 
      res = 300)
  print(DimPlot(sobj, group.by=paste0('wsnn_res.', res), label=T, 
                reduction="wnn.umap") + NoLegend())
  dev.off()
}

sobj[['ADT']] = adt_assay # add back in
save(sobj, file=sprintf("%s/wnn_obj.RData", OUT_DIR_WNN))



