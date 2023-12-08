# code for running RPCA on the ADT of all files in a list
# adapted from code from Ravi Patel
#
# usage:
#  RScript run_rpca.R ${in_dir} ${out_dir} ${path_to_sample_list}
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
#    - npcs
#    - remove isotype controls (and/or provide a list)
# [ ] add additional print messages
# [ ] set up split plots to scale to the number of libraries

library(Seurat)
library(tidyverse)

args = commandArgs(trailingOnly=T)

IN_DIR = args[1]
OUT_DIR=args[2]
SAMPLE_FILE=args[3]
list_samples = read_tsv(SAMPLE_FILE, col_names=F)[,1]


remove.isotype.ctls = TRUE
npcs=18


sobjs.list <- lapply(X = list_samples, FUN=function(SAMPLE_NAME){
  print(SAMPLE_NAME)
  sobj = readRDS(sprintf("%s/%s/automated_processing/%s_processed.rds", IN_DIR,
                             SAMPLE_NAME, SAMPLE_NAME))
  sobj = subset(sobj, DROPLET.TYPE.FINAL=="SNG")
  DefaultAssay(sobj) = "ADT"
  return(sobj)
})
names(sobjs.list) = list_samples

#Using all non-isotype markers for the following analyses.
 if (remove.isotype.ctls){
   features <- grep("isotype", rownames(sobjs.list[[1]][["ADT"]]), ignore.case = T, value = T, invert = T)
   # Remove the isotype controls
   features <- features[ features %in% grep("isotype", 
                                           rownames(sobjs.list[[1]][["ADT"]]), 
                                           ignore.case = T, value = T, 
                                           invert = T) ]
 }
 print(features)                                        

# pre-process before running integration
sobjs.list <- lapply(list_samples, FUN = function(SAMPLE_NAME) {
  x = sobjs.list[[SAMPLE_NAME]]
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = FALSE)

  x <- RunUMAP(x, dims=1:npcs)
  DimPlot(x)+NoLegend()
  ggsave(sprintf("%s/%s_pre_rpca.png", OUT_DIR, SAMPLE_NAME), height=5, width=5)
  
  return(x)
  
})
names(sobjs.list) = list_samples

# run integration
immune.anchors <- FindIntegrationAnchors(object.list = sobjs.list, 
                                             anchor.features = features, 
                                             reduction = "rpca" )
integ_data <- IntegrateData(anchorset = immune.anchors, 
  new.assay.name = "integrated.ADT")
  
DefaultAssay(integ_data) <- "integrated.ADT"

# make sure to regress out nCount/nFeature ADT prior to WNN
# S.Score, G2M.Score are optional
integ_data <- ScaleData(integ_data, vars.to.regress=c('nCount_ADT', 'nFeature_ADT', 
                                                      'S.Score', 'G2M.Score'), 
                                                      verbose = FALSE)

integ_data <- RunPCA(integ_data, npcs = npcs, verbose = FALSE)
integ_data <- RunUMAP(integ_data, reduction = "pca", dims = 1:npcs)
DimPlot(integ_data, split.by="orig.ident")+NoLegend()
ggsave(sprintf("%s/post_rpca.pdf", OUT_DIR), height=5, width=5*length(sobjs.list))

saveRDS(integ_data, file = sprintf("%s/integ_ADT.rds", OUT_DIR))
