require(Seurat)
require(tidyverse)
require(pals)
require(ggpubr)
# require(LaplacesDemon) # if suggest_multimodal = TRUE

get_plot_grid_layout <- function(no_of_plots) {
  sq = sqrt(no_of_plots)
  fl = floor(sq)
  if(sq > fl+0.5) {
    nrow = fl + 1
    ncol = fl + 1
  } else if(sq == fl) {
    nrow = fl
    ncol = fl
  } else {
    nrow = fl + 1
    ncol = fl
  }
  if(nrow * ncol < no_of_plots)
    print("The grid size doesn't fit all plots.")
  return(list("nrow"=nrow, "ncol"=ncol))
}



#' umap_by_cluster
#' @export
#' @author Ravi K. Patel
#' @description UMAP plots by cluster/group.
#' @param sobj	Seurat Object.
#' @param reduction Character Reduction name to be used from the Seurat object.
#' @param cluster_clm Character Name of a column in meta.data to use for grouping single-cells.
#' @param subset_n Integer Down-sample to \code{subset_n}. Default -1 (no down-sampling).
#' @param suggest_multimodal Logical Whether to indicate which clusters may have dispersed clusters. This is indicated using "*" in the top-right corner of the UMAP. The unimodal/multimodal calls are made using \code{LaplacesDemon::is.unimodal}.
#' @return A ggarrange object.
#' @examples
#' p = umap_by_cluster(seurat_obj, reduction="umap", cluster_clm="wsnn_res.2", subset_n = -1, suggest_multimodal = FALSE)
#' # It is recommended to save the plots as a PNG if working with large Seurat objects and not down-sampling the data using "subset_n"
#' png("umap_by_cluster.png", width=20, height=20, units="in", res=600)
#' print(p)
#' dev.off()

umap_by_cluster <- function(sobj, reduction, cluster_clm, subset_n = -1, suggest_multimodal = FALSE) {
  set.seed(1234)
  if(! reduction %in% Reductions(sobj)) {
  	stop(reduction, " is not in list of reductions. Available reductions are ", Reductions(sobj), "\n")
  }
  if(! cluster_clm %in% colnames(sobj@meta.data)) {
  	stop(cluster_clm, " is not in list of meta.data columns. Available columns in meta.data are ", colnames(sobj@meta.data), "\n")
  }
  reduction_xy = colnames(sobj@reductions[[reduction]]@cell.embeddings)
  cell_metadata = FetchData(sobj, vars=c(reduction_xy[1], reduction_xy[2], cluster_clm)) %>% setNames(c("UMAP1","UMAP2","cluster"))
  if(subset_n == -1) {
  	cell_metadata_sub = cell_metadata
  } else {
  	cell_metadata_sub = cell_metadata %>% sample_n(min(nrow(.), subset_n))
  }
  x_range = range(cell_metadata_sub$UMAP1) %>% scales::expand_range(add = 0.2)
  y_range = range(cell_metadata_sub$UMAP2) %>% scales::expand_range(add = 0.2)
  cm_tmp = cell_metadata_sub %>% group_by(cluster) %>% sample_n(50, replace=T)
  plot_list = list()
  colors = rep(pals::kelly(22)[-1], 20)
  groups = sort(unique(cell_metadata_sub$cluster))
  for(i in 1:length(groups) ) {
  	clust = groups[i]
    cm_tmp = cell_metadata_sub %>% filter(cluster == clust)
    is_unimodal = TRUE
    if(suggest_multimodal)
      is_unimodal = LaplacesDemon::is.unimodal(cm_tmp$UMAP1) & LaplacesDemon::is.unimodal(cm_tmp$UMAP2)
    p = ggplot(cm_tmp, aes(UMAP1, UMAP2)) +
      geom_rect(data=cm_tmp, aes(xmin=UMAP1-0.05,xmax=UMAP1+0.05,ymin=UMAP2-0.05,ymax=UMAP2+0.05), color="gray95", fill="gray95") +
      geom_point(size=0.1, alpha=0.3, color=colors[as.numeric(i)+1]) + theme_classic() + xlim(x_range) + ylim(y_range) + ggtitle(clust) +
      labs(x = reduction_xy[1], y = reduction_xy[2])
    if(! is_unimodal & suggest_multimodal)
      p = p + annotate("text", x=x_range[2]-1, y=y_range[2]-1, label="*")
    plot_list[[clust]] = ggExtra::ggMarginal(p, type = "histogram", color = colors[as.numeric(i)+1], fill = colors[as.numeric(i)+1], binwidth = 0.1)
  }
  grid_size = get_plot_grid_layout(length(groups))
  p = ggpubr::ggarrange(plotlist=plot_list, ncol=grid_size$ncol, nrow=grid_size$nrow)
  return(p)
}

