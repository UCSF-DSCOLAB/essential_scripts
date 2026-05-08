#' Check inputs relevant to looping DGE functions within different cell types, and adjust 'cell_targets'.
#' @param metadata data.frame containing cell or sample metadata, with rownames containing column names of \code{counts}.
#' @param cell_by String naming a column of \code{metadata} which holds clusters, annotations, or other cell groupings to explore within
#' @param cell_targets NULL (to default to all options) or a string vector naming particular levels of the \code{cell_by}-data to run DGE within.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#' @param min_per_group To be changed? Sets a direct numeric minimum of number of cells/samples per group for keeping a given \code{cell_target}
#' @details Performs all checks necessary to looping accross cell types.
#'
#' Ensures that \code{cell_by} names an existing column in \code{metadata}, erroring if not.
#'
#' Fills 'cell_targets' with all options if it was left as NULL.
#'
#' Checks for any 'cell_targets' with fewer than 'min_per_group' number of samples within \code{case_group} or \code{reference_group}.
#' Any targets with low cell numbers are trimmed from the cell_targets output, with a warning noting the group with low numbers.
#'
#' Finally, ensures at least one 'cell_target' remains (erroring if not) before outputting the final set.
#' @examples
#' cell_targets <- .input_checks_within_cells(
#'     metadata, cell_by, cell_targets,
#'     dge_by, case_group, reference_group, contrasts,
#'     min_per_group = 0
#' )
#' @author Dan Bunis
#' @importFrom dittoViz colLevels
.input_checks_within_cells <- function(
    metadata,
    cell_by,
    cell_targets,
    dge_by,
    case_group,
    reference_group,
    contrasts,
    min_per_group = 0
) {

    # Error if cell-by metadata does not exist
    if (!cell_by %in% colnames(metadata)) {
        stop("'cell_by' of ", cell_by, " is not a column of metadata")
    }
    # check min-per-group is numeric
    if (!as.numeric(min_per_group) == min_per_group) {
        stop("'min_per_group' of ", min_per_group, " is not numeric")
    }
    
    # check cell targets are in the cell_by column
    if (!is.null(cell_targets) & !all(cell_targets %in% metadata[,cell_by])) {
        stop("The following 'cell_targets' are not levels of the 'cell_by' column: ", paste0(cell_targets[!cell_targets %in% metadata[,cell_by]], collapse = ", "))
    }

    # Establish cell targets if not given
    if (is.null(cell_targets)) {
        # dittoViz::colLevels will automatically drop any empty factor-levels, avoiding a potential warning later on
        cell_targets <- dittoViz::colLevels(cell_by, metadata)
    }

    # needed for min_per_group check below
    if (!is.null(case_group) | !is.null(reference_group)) {
        .input_check_dge_case_ref(metadata, dge_by, case_group, reference_group) 
    } else if (is.null(contrasts)) {
        stop("Error. Either case_group and reference_group or contrasts must be provided for DGE analysis.")  
    }

    # TODO add version for contrasts

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
            warning("Skipping cell type ", cell_targ, " because dge_by group(s) ", cols_str, has_str, number_str)
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

