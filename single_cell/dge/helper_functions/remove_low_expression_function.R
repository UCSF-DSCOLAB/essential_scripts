### Remove lowly expressed genes from the count matrix.
#
#' Remove lowly expressed genes from the count matrix. The function first determines the size of the smallest group (denoted min_gsize), where groups are defined by \code{dge_by}. A gene is retaind only if it is expressed (CPM > 1) in at least a \code{min_frac} proportion of min_gsize samples overall (i.e., not within each group separately).
#' @param counts Feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param dge_groups A vector of group labels (from \code{dge_by}) to consider. Only samples belonging to these groups are used to determine the minimum number of samples required to have CPM > 1. If NULL, all groups in \code{dge_by} are considered.
#' @param min_frac Genes with counts per million (CPM) more than 1 in at least \code{min_frac} proportion of samples of the smaller group are retained, regardless of group membership. Default 0.6 (i.e., 60%).
#' @return \code{counts} after removing lowly expressed genes.
.remove_low_expression <- function(counts, metadata, dge_by, dge_groups=NULL, min_frac=0.6) {
  # Determine the number of sample in the smaller group (min_gsize) and keep features that are expressed with CPM > 1 in at least min_frac of min_gsize samples, regardless of the sample group identity.
  warning("Assuming that the order of samples (columns) in 'counts' matches the order of samples (rows) in 'metadata'.")
  if(is.null(dge_groups)) {
    warning("dge_groups is not provided. Using all groups in dge_by to determine the minimum number of samples required to have CPM > 1.")
    sample_idx = rep(TRUE, nrow(metadata))
  } else {
    sample_idx = metadata[, dge_by] %in% dge_groups
  }
  min_gsize = min(table(metadata[ sample_idx, dge_by]))
  feats = names(which(rowSums(edgeR::cpm(counts[, sample_idx]) > 1) > (min_gsize * min_frac)))
  cat(paste0("Number of final features: ", length(feats)), "\n")
  counts = counts[feats,]
  return(counts)
}


