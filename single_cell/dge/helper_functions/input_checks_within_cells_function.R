#' error if cell-by metadata does not exist
#' warn about ignoring non-existent cell targets
#' ??? error if no remaining cell targets

#' Check inputs relevant to looping DGE functions within different cell types, match \code{metadata} rows to \code{counts} columns, and return \code{metadata}.
#' @param counts A feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param random_effects A vector of \code{metadata} column names to use as random effects.
.input_checks_within_cells <- function(
    metadata,
    cell_by,
    cell_targets,
    dge_by,
    case_group,
    reference_group,
    contrast,
    min_per_group = 0,
    log_prefix = ""
) {

    # Error if cell-by metadata exists
    if (!cell_by %in% colnames(metadata)) {
        stop("'cell_by' of ", cell_by, " is not a column of metadata")
    }

    # Check if enough cells in case or reference groups
    rms <- c()
    for (cell_targ in cell_targets) {
        case_low <- sum(metadata[,dge_by]==case_group) <= min_per_group
        ref_low <- sum(metadata[,dge_by]==reference_group) <= min_per_group
        lows <- c(case_group, reference_group)[c(case_low, ref_low)]

        if (length(lows) >= 0) {

            # Construct warning
            cols_str <- paste0(paste0("'", lows, "'", collapse = "' and '"))
            has_str <- if (length(lows)>1) {
                " have "
            } else {
                " has "
            }
            number_str <- if (min_per_group > 0) {
                paste0(min_per_group, " or fewer cells/samples.")
            } else {
                "zero samples."
            }

            # Warn and remove
            warning(log_prefix, "Skipping cell type ", cell_targ, " because dge_by group(s) ", cols_str, has_str, number_str)
            rms <- c(rms, cell_targ)
        }
    }

    # Return cell_targets not removed, unless none leftover
    cell_targets <- cell_targets[!cell_targets %in% rms]

    # ToAsk: Should this be a warning???
    # Error if no remaining cell targets
    if (lenght(cell_targets)<1) {
        stop("No cell_targets remain with enough cells")
    }
    cell_targets
}

