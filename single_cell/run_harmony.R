# code for running harmony on all files in a list
# adapted from code from Arjun Rao
#
# usage:
#  RScript run_harmony.R ${in_dir} ${out_dir} ${path_to_sample_list}
#
# The ${in_dir} should be the `processed` directory within a project
# The sample file should be a tsv formatted as one sample name per line:
#  sample1
#  sample2
# Note that it assumes standard immunox filesystem structure for inputs
#
# TODO: 
# [ ] Missing options:
#    - vars.to.regress
#    - filter.doublets
# [ ] add additional print messages
# [ ] set up split plots to scale to the number of libraries

library(Seurat)
library(harmony)
library(tidyverse)
library(cowplot)

args = commandArgs(trailingOnly=T)
IN_DIR= args[1]
OUT_DIR=args[2]
SAMPLE_FILE=args[3]
list_samples = read_tsv(SAMPLE_FILE, col_names=F)[,1]


if (!file.exists(sprintf('%s/merged_temp.RData', OUT_DIR))){
  dir.create(OUT_DIR)
  
  sobjs.list <- lapply(list_samples, FUN = function(SAMPLE_NAME) {

    sobj = readRDS(sprintf("%s/%s/automated_processing/%s_filtered.rds", IN_DIR,
                                   SAMPLE_NAME, SAMPLE_NAME))
    
    DefaultAssay(sobj) = "RNA" 
    sobj = NormalizeData(sobj)
    
    # remove doublets
    sobj = subset(sobj, DROPLET.TYPE.FINAL=="SNG")  

    sobj@meta.data$LIBRARY = sobj@meta.data$orig.ident
    
    return(sobj)
  })
  
  names(sobjs.list) = list_samples
  
  first_sobj <- sobjs.list[[names(sobjs.list)[1]]]
  sobjs.list[[names(sobjs.list)[1]]] <- NULL
  merged_data <- merge(x = first_sobj, y = unname(sobjs.list)) 
  
  merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
  merged_data <- ScaleData(merged_data, vars.to.regress=c('percent.mt', 'percent.ribo', 
  	      	 					'nCount_RNA', 'nFeature_RNA',
                                                          'S.Score', 'G2M.Score'), verbose = FALSE)
  merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)
  merged_data <- RunUMAP(merged_data, dims = 1:30)
  
  save(merged_data, file=sprintf('%s/merged_temp.RData', OUT_DIR))
} else {
  load(sprintf('%s/merged_temp.RData', OUT_DIR))
}

raw_metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  raw_metadata <- merge(raw_metadata,
                        cell_embeddings,
                        by=0)
  rownames(raw_metadata) <- raw_metadata$Row.names
  raw_metadata$Row.names <- NULL
}

write.table(raw_metadata,
            file=sprintf('%s/merged_metadata_raw.tsv', OUT_DIR),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)


png(filename=sprintf('%s/merged_raw_library_umap.png', OUT_DIR), 
    width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=sprintf('%s/merged_raw_split_library_umap.png', OUT_DIR), 
    width = 15, height = 30, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4, raster=F) + 
        theme(legend.position="none", axis.title=element_blank(), 
              axis.text=element_blank()))
dev.off()

pdf(sprintf('%s/merged_harmony_convergence.pdf', OUT_DIR), height=5, width=5)
merged_data <- RunHarmony(merged_data,
                          "LIBRARY",
                          assay.use='RNA',
                          plot_convergence = TRUE,
                          max.iter.harmony=20,
                          max.iter.cluster=30)
dev.off()

merged_data <- RunUMAP(merged_data,
                       dims = 1:30,  # Num PCs to use
                       reduction='harmony',
                       n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                       b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                       verbose = FALSE)

# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             dims = 1:30,  # Num PCs to use
                             reduction='harmony',
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
)

save(merged_data, file=sprintf('%s/merged_processed.RData', OUT_DIR))

png(filename=sprintf('%s/merged_harmony_library_umap.png', OUT_DIR), 
    width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=sprintf('%s/merged_harmony_split_library_umap.png', OUT_DIR), 
    width = 15, height = 30, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4, raster=F) + 
        theme(legend.position="none", axis.title=element_blank(), 
              axis.text=element_blank()))
dev.off()

# do it by highlighting cells - easier to viz
libraries = unique(merged_data@meta.data$LIBRARY)
list_cl = sapply(libraries, function(my_library) rownames(merged_data@meta.data %>% 
                                                            filter(LIBRARY==my_library)))
png(filename=sprintf('%s/merged_harmony_split_library_umap_hl.png', OUT_DIR), 
    width = 15, height = 30, units = "in", res = 300)
plots = lapply(list_cl, function(.x) DimPlot(merged_data, cells.highlight=.x, combine=T)+NoLegend())
plot_grid(plotlist=plots)
dev.off()

for (res in c( 0.8, 1.0, 1.2, 1.4)){
  if (paste0('louvain_res', res) %in% colnames(merged_data@meta.data)){
    next
  }
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                              algorithm = 1,
                              resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters
  
  png(filename=sprintf('%s/merged_louvain_res_%s.png', OUT_DIR, res), width = 5, height = 5, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T) + NoLegend())
  dev.off()
}

png(filename=sprintf('%s/merged_harmony_split_library_umap_clus.png', OUT_DIR), width = 15, height = 12, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()


metadata <- merged_data@meta.data
save(metadata, file=sprintf('%s/merged_clus_meta.RData', OUT_DIR))

for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  metadata <- merge(metadata,
                    cell_embeddings,
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}


write.table(metadata,
            file=sprintf('%s/merged_metadata.tsv', OUT_DIR),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

