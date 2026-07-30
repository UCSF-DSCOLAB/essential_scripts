### Check inputs to DGE functions
#
#' Check and adjust inputs for DGE functions.
#' @param counts A feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param random_effects A vector of \code{metadata} column names to use as random effects.
#' @details Under construction, but will: match \code{metadata} rows to \code{counts} columns; trim genes by expression minimums; more; and finally return adjusted \code{counts}, \code{metadata}, and \code{contrasts}.
.input_validation <- function(
    counts,
    metadata,
    dge_by,
    method,
    cell_by,
    case_group,
    reference_group,
    cell_targets,
    contrasts,
    dge_groups,
    fixed_effects,
    random_effects,
    min_frac,
    return_model
) {
    if (! (is.matrix(counts) || is.data.frame(counts) || is(counts, "sparseMatrix")) ) stop("Validation Error: Input data is not a matrix or a data.frame.")
    if (ncol(counts)==0 | nrow(counts)==0) stop("Validation Error: Count matrix does not have at least one column and row.")
    if(! dge_by %in% colnames(metadata) ) stop("Validation Error: dge_by column does not exist in metadata.")
    if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) stop("Validation Error: One or more fixed_effects do not exist in metadata.")
    if(! is.null(random_effects) & ! all(random_effects %in% colnames(metadata)) ) stop("Validation Error: One or more random_effects do not exist in metadata.")

    # Subset and order metadata according to the order of samples in counts
    if(! all( colnames(counts) %in% rownames(metadata)) ) stop("Validation Error: All column names of counts must exists as rownames of metadata.")
    metadata <- metadata[ colnames(counts), ]

    if (identical(cell_targets, NULL)) {
        cell_targets <- unique(metadata[,cell_by])
    }
    # STUB, to be replaces with 'min_frac' usage!
    if (any(table(metadata[,cell_by])[cell_targets] <= 5)) {
        cells <- table(metadata[,cell_by])
        cell_targets <- names(cells)[cells>5]
    }

    outs <- list(
        counts = counts,
        metadata = metadata,
        contrasts = contrasts,
        cell_targets = cell_targets
    )
    case_ref = FALSE
    cntrst = FALSE
    if( !is.null(case_group) & !is.null(reference_group) ) {
        case_ref = TRUE
        if( ! case_group %in% metadata[,dge_by] ) stop("Validation Error: Case group does not exist in dge_by column of metadata.")
        if( ! reference_group %in% metadata[,dge_by] ) stop("Validation Error: Reference group does not exist in dge_by column of metadata.")
    }
    if(! is.null(contrasts) )
        cntrst = TRUE
    if(case_ref & cntrst) {
        warning("Both contrasts and case-reference groups are given. The contrasts will be ignored for the analysis.")
        outs[['contrasts']] <- NULL
    }

    outs
}

