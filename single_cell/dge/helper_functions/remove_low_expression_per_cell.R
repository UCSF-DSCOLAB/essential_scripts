### Remove lowly expressed genes from the count matrix.
#
#' Remove lowly expressed genes from a count matrix
#'
#' Filters genes based on minimum numbers of cells, per dge_group, with non-zero
#' capture.
#' In order to be retained, genes must be captured (non-zero expression) in at
#' least \code{min_exp_frac} fraction of cells of any group defined by
#' \code{dge_by}.  If \code{dge_groups} is provided, only the expression
#' fraction among these groups will be assessed.
#' @param counts Feature-by-sample matrix of raw counts. Only samples intended
#'   for DGE analysis should be included.
#' @param metadata Sample metadata. Row names must correspond to the column
#'   names of \code{counts}.
#' @param dge_by Name of the column in \code{metadata} defining the groups used
#'   for DGE comparison.
#' @param dge_groups Optional vector of group labels from \code{dge_by} to
#'   consider. Only fraction of expression among these groups will be considered.
#' @param min_cells_frac_exp Decimal setting the minimum fraction filter.  To be
#'   genes must be captured in at least this fraction of cells of at east
#'   one group of cells defined by \code{dge_by}.
#'
#' @return A count matrix containing only genes that pass the low-expression
#'   filter.
.remove_low_expression_per_cells <- function(
    counts,
    metadata,
    dge_by,
    dge_groups = NULL,
    min_cells_frac_exp = 0.2
) {
    # (Assuming counts and metdata pre-trimmed to samples with sufficient numbers of cells to be trusted)
    if (!is.null(dge_groups)) {
      in_groups <- metadata[,dge_by] %in% dge_groups,
      metadata <- metadata[in_groups,]
      counts <- counts[,in_groups]
    }
    dge_vals <- metadata[,dge_by]
    expressed <- counts > 0
    frac_exp_per_group <- do.call(cbind, lapply(
        unique(dge_vals),
            function(grp) {
              rowSums(expressed[,dge_vals==grp] > 0) / length(dge_vals==grp)
            }
    ))
    keep <- apply(frac_exp_per_group, 1, function(fracs) { any(fracs >= min_exp_frac) })
    cat(paste0("Number of final features: ", sum(keep)), "\n")
    counts[keep,]
}
