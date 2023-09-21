# code for running find markers parallelized such that each job is a cluster
# usage:
#   sbatch --array=0-${num_clus} run_find_markers.sh ${input_dir} ${fname_sobj} ${output_dir} ${idents} ${test_use}
#
# to do: update so defaults to current idents and wilcoxon

library(Seurat)
library(tidyverse)
library(miceadds)
args = commandArgs(trailingOnly=T)
indir= args[[1]]
fname = args[[2]]
outdir = args[[3]]
idents = args[[4]]
test_use = args[[5]]
cluster_no = as.character(args[[6]])


sobj = load.Rdata2(path=indir, filename=fname) # --> sobj
DefaultAssay(sobj) = "RNA"
if (!dir.exists(outdir)){
   dir.create(outdir)
}
if (!idents %in% colnames(sobj@meta.data )){
  print(sprintf("Error. Ident %s specified is not in metadata", idents))

} else {
  Idents(sobj) = idents

  nclus = length(unique(sobj@active.ident))
  print(sprintf("Object has %s clusters", nclus))
  if (as.numeric(cluster_no) > (nclus-1)){
     print(sprintf("Max number of clusters is %s, %s exceeds that", idx, cluster_no))
  } else {
    if (!cluster_no %in% unique(sobj@active.ident)){
        cluster_idx = as.numeric(cluster_no)
    	cluster_no = unique(sobj@active.ident)[cluster_idx] 
    	print(sprintf("Warning, cluster number not in idents, using the number %s as an index, becomes cluster %s",cluster_idx, cluster_no))
    }
    markers = FindMarkers(sobj, ident.1=as.character(cluster_no), test.use=test_use, only.pos=T, min.pct=0.3)

  if (! dir.exists(sprintf("%s/markers/", outdir))){
    dir.create(sprintf("%s/markers/", outdir))
  }
  saveRDS(markers, file=sprintf("%s/markers/RNA_markers_%s_%s.RDS", outdir,  idents, cluster_no))
}
}




