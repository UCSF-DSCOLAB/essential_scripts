### Remove lowly expressed genes from the count matrix.
#
#' Remove lowly expressed genes from a count matrix
#'
#' Filters genes based on their expression in counts per million (CPM).
#' The minimum number of samples required for a gene to pass the filter is
#' calculated by multiplying \code{min_frac} by the number of samples in the
#' smallest group defined by \code{dge_by} and rounding to the nearest whole
#' number. A gene is retained if its CPM exceeds \code{min_exp} in at least
#' this many samples, regardless of group membership.
#' @param counts Feature-by-sample matrix of raw counts. Only samples intended
#'   for DGE analysis should be included.
#' @param metadata Sample metadata. Row names must correspond to the column
#'   names of \code{counts}.
#' @param dge_by Name of the column in \code{metadata} defining the groups used
#'   for DGE comparison.
#' @param dge_groups Optional vector of group labels from \code{dge_by} to
#'   consider when determining the smallest group size. Only samples belonging
#'   to these groups are used to determine \code{min_gsize}. If \code{NULL},
#'   all groups in \code{dge_by} are considered.
#' @param min_frac Fraction of the smallest group used to determine the minimum
#'   number of samples in which a gene must exceed \code{min_exp}. The resulting
#'   sample count is rounded to the nearest whole number. Default is \code{0.6}.
#' @param min_exp Minimum expression cutoff, in counts per million (CPM), used
#'   to define whether a gene is expressed in a sample. Genes must have
#'   \code{CPM > min_exp} in the required number of samples. Default is
#'   \code{1}.
#'
#' @return A count matrix containing only genes that pass the low-expression
#'   filter.
.remove_low_expression <- function(
    counts,
    metadata,
    dge_by,
    dge_groups = NULL,
    min_frac = 0.6,
    min_exp = 1
) {
  # warning("Assuming that the order of samples (columns) in 'counts' matches the order of samples (rows) in 'metadata'.")
  if(is.null(dge_groups)) {
    warning(paste0("dge_groups is not provided. Using all groups in dge_by to determine the minimum number of samples required to have CPM > ", min_exp, "."))
    sample_idx = rep(TRUE, nrow(metadata))
  } else {
    sample_idx = metadata[, dge_by] %in% dge_groups
  }
  min_gsize = min(table(metadata[ sample_idx, dge_by]))
  feats = names(which(rowSums(edgeR::cpm(counts[, sample_idx]) > min_exp) > round(min_gsize * min_frac, 0)))
  removed_feats = setdiff(rownames(counts), feats)
  cat(paste0("Number of final features: ", length(feats), " (removed ", length(removed_feats), " of ", nrow(counts) , " input features).\n"))
  counts = counts[feats,]
  return(counts)
}


