### Remove lowly expressed genes from the count matrix.
#
#' Remove lowly expressed genes from the count matrix.
#' @param counts: feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}
#' @param metadata sample metadata, with rownames containing column names of \code{counts}
#' @param dge_by name of \code{metadata} column to use DGE comparion. Must have exactly two levels
#' @param min_frac Genes with counts per million (CPM) more than 1 in at least \code{min_frac} fraction of samples of the smaller group are retained. Default 0.6 (i.e., 60%)
#' @return \code{counts} after removing lowly expressed genes
remove_low_expression <- function(counts, metadata, dge_by, min_frac=0.6) {
  # Determine the number of sample in the smaller group (min_gsize) and keep features that are expressed with CPM > 1 in at least min_frac of min_gsize samples, regardless of the sample group identity.
  min_gsize = min(table(metadata[ , dge_by]))
  feats = names(which(rowSums(cpm(counts) > 1) > (min_gsize * min_frac)))
  cat(paste0("Number of final features: ", length(feats)), "\n")
  counts = counts[feats,]
  return(counts)
}

