### Remove samples with too few cells.
# QUESTIONS:
#  - should there be a toggle for removing entire samples if any dge_group cell set of the sample is too small?
#  - How to log...
#      - Messaging is ephemeral.
#      - Should it be included in the output itself?
#' Remove samples with too few cells.
#'
#' @details Filters cells coming from samples, or parts of samples, with fewer cells than
#' the 'min_cells' cutoff.
#' If 'dge_by' groups entire samples (e.g. dge across disease versus healthy subjects), then samples with fewer than 'min_cells' are filtered out. 
#' If 'dge_by' cuts within samples (e.g. dge across cell phenotypes), cells within each
#' dge_group are counted then filtered out for any sample & dge_group set that
#' does not meet the 'min_cells' minimum.
#' This function is run per each individual cell_target, so sample filter can
#' differ from cell type to cell type.
#'
#' @param metadata A data.frame. Sample metadata where rows are cells and columns
#'   corresponding to \code{dge_by} and \code{sample_by} must exist.
#' @param dge_by String. Name of the column in \code{metadata} defining the groups used
#'   for DGE comparison.
#' @param sample_by String. Name of the column in \code{metadata} defining sample
#'   groupings of cells.
#' @param min_cells Integer providing the minimum number of cells cutoff where
#'   samples, or sample subsets, which contain fewer cells will be ignored.
#' @param case_group,reference_group Strings representing the values of
#'   \code{metadata[,dge_by]} reflecting the groupings to be assessed. These are
#'   only used, here, for a final check of if all relevant cells for either
#'   group have been dropped.
#' @param ct String. For logging purposes, the cell_target being assessed.
#'
#' @return rownames of metadata that passing min_cells filtering.
.remove_low_cell_samples <- function(
    metadata,
    dge_by,
    sample_by,
    min_cells,
    case_group,
    reference_group,
    ct
) {
    # (Assuming metadata already cut to dge_groups only!)

    # For eliminations check at end
    orig_dge_groups <- unique(metadata[,dge_by])

    samp_ignore <- c()
    dge_ignore <- c()
    cell_nums_ignore <- c()
    # Check if dge cuts within samples
    counts_samp_dge <- table(metadata[,sample_by], metadata[,dge_by])
    if (if(any(rowSums(counts_samp_dge!=0)>1))) {
        # Count within dge_by groups per each sample
        for (samp in unique(metadata[,sample_by])) {
            # Assessment, within each sample
            too_small <- counts_samp_dge[samp,][counts_samp_dge[samp,]<min_cells]
            if (length(too_small)>0) {
                samp_ignore <- c(samp_ignore, rep(samp, length(too_small)))
                dge_ignore <- c(dge_ignore, names(too_small))
                cell_nums_ignore <- c(cell_nums_ignore, unname(too_small))
                # Removal, within each sample
                metadata <- metadata[!(metadata[,sample_by]==samp & metadata[,dge_by] %in% names(too_small)),]
            }
        }
        # Logging
        if (length(samp_ignore) > 0) {
            cat(
                "For cell type, ", ct,
                " Ignoring cells from ", length(samp_ignore), " 'sample_by' + 'group_by' cell groupings with fewer than 'min_cells' cells:\n\t",
                paste0(
                    paste0(samp_ignore," + ", dge_ignore, ": ", cell_nums_ignore, " cells"),
                    collapse = "\n\t"
                ),
                "\n"
            )
        }
    } else {
        # Count all cells per each sample
        counts_samp <- table(metadata[,sample_by])
        # Assessment
        samp_ignore <- names(counts_samp)[counts_samp < min_cells]
        cell_nums_ignore <- unname(counts_samp)[counts_samp < min_cells]
        # Removal
        metadata <- metadata[!(metadata[,sample_by] %in% samp_ignore),]
        # Logging
        if (length(samp_ignore) > 0) {
            cat(
                "For cell type, ", ct,
                " Ignoring cells from ", length(samp_ignore), " 'sample_by' groupings with fewer than 'min_cells' cells:\n\t",
                paste0(
                    paste0(samp_ignore, ": ", cell_nums_ignore, " cells"),
                    collapse = "\n\t"
                ),
                "\n"
            )
        }
    }

    # Ensure have not removed an entire important dge_group.  Return NA if so.
    lost_groups <- orig_dge_groups[!orig_dge_groups %in% unique(metadata[,dge_by])]
    if (length(lost_groups)>0) {
        non_key_lost <- lost_groups[!lost_groups %in% c(case_group, reference_group)]
        key_lost <- lost_groups[lost_groups %in% c(case_group, reference_group)]
        if (length(non_key_lost) > 0) {
            message(
                "No samples of group '",
                paste(non_key_lost, collapse = "' or '"),
                "' meet the 'min_cells', cutoff for celltype '", ct, "'.  But these groups are not the target case or reference.")
        }
        if (length(key_lost) > 0) {
            warning(
                "No samples of group '",
                paste(key_lost, collapse = "' or '"),
                "' meet the 'min_cells', cutoff for celltype '", ct, "'.  Returning NA, and moving on.")
            return(NA)
        }
    }
    rownames(metadata)
}
