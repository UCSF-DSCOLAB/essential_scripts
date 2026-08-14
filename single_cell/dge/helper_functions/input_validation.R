### Check inputs to DGE functions
#
#' Check and adjust inputs for DGE functions.
#' @param counts A feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param method Name of DGE method to use. One of "deseq2", "dreamlet", "edger", or "mast".
#' @param cell_by Name of \code{metadata} column to use for cell grouping.
#' @param sample_by Name of \code{metadata} column to use for sample grouping. 
#' @param metadata_cell_count Metadata column for cell counts in pseudobulk. 
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param cell_targets NULL or vector of cell types to use for DGE analysis. If NULL, all cell types in \code{cell_by} will be used.
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#' @param dge_groups NULL or the vector of group names from dge_by column that should be used for filtering out lowly expressed genes
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param random_effects A vector of \code{metadata} column names to use as random effects.
#' @param min_exp Minimum expression cutoff, in counts per million (CPM), used to define whether a gene is expressed in a sample.
#' @param min_frac Minimum proportion of samples required to have at least MIN.COUNT counts per gene.
#' @param return_model Boolean indicating whether the function should return the model or the DGE result data frame.
#' @details Under construction, but will: match \code{metadata} rows to \code{counts} columns; trim genes by expression minimums; more; and finally return adjusted \code{counts}, \code{metadata}, and \code{contrasts}.
.input_validation <- function(
    counts,
    metadata,
    dge_by,
    method,
    cell_by,
    sample_by,
    metadata_cell_count,
    case_group,
    reference_group,
    cell_targets,
    contrasts,
    dge_groups,
    fixed_effects,
    random_effects,
    min_exp,
    min_frac,
    return_model
) {
    if (! (is.matrix(counts) || is.data.frame(counts) || is(counts, "sparseMatrix")) ) stop("Validation Error: Input data is not a matrix or a data.frame.")
    if (ncol(counts)==0 | nrow(counts)==0) stop("Validation Error: Count matrix does not have at least one column and row.")
    if (!is.data.frame(metadata)) stop("Error. metadata must be a data.frame.")
    if(! all( colnames(counts) == rownames(metadata)) ) stop("Validation Error: Column names of counts and rownames of metadata do not match.")
    # QUESTION: Do we want to swap to allowing extra columns to exist in the metadata?
    # if(! all( colnames(counts) %in% rownames(metadata)) ) stop("Validation Error: All column names of counts must exists as rownames of metadata.")
    # metadata <- metadata[ colnames(counts), ]
    metadata <- .remove_NA_cols(metadata)
    if (ncol(metadata) < 1) stop("Validation Error: Metadata must contain at least one usable column after removing all-NA columns.")
    if(! dge_by %in% colnames(metadata) ) stop("Validation Error: dge_by column does not exist in metadata.")
    if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) stop("Validation Error: One or more fixed_effects do not exist in metadata.")
    if(! is.null(random_effects) & ! all(random_effects %in% colnames(metadata)) ) stop("Validation Error: One or more random_effects do not exist in metadata.")

    # check dge, fixed_effects and random_effects dont have spaces or special characters
    vars <- unique(c(dge_by, fixed_effects, random_effects))
    if(is.null(vars)) stop("Validation Error: No variables provided for DGE analysis.")
    if(any(grepl(" ", vars)) | any(grepl("[^a-zA-Z0-9_]", vars))) stop("Validation Error: One or more of dge_by, fixed_effects, random_effects elements contain spaces or special characters.")

    # If cell_by is NULL, assume the data is bulk RNA-seq. In this case, create a new cell_by column with a same value for all rows in metadata
    if(is.null(cell_by)) {
	ts_log( paste0("Warning! `cell_by` is NULL, therefore, assuming the input data represents bulk RNA-seq. If the input data represents pseudobulk (with multiple cell types) or single-cell data, provide 'cell_by' and run_dge() again.") )
        cell_by = "new_cell_by"
        metadata[,cell_by] = "bulk_data"
    }

    # check cell_by, sample_by, and metadata_cell_count are in metadata
    if(! is.null(cell_by) & ! cell_by %in% colnames(metadata)) stop("Validation Error: cell_by column does not exist in metadata.")
    if(! is.null(sample_by) & ! sample_by %in% colnames(metadata)) stop("Validation Error: sample_by column does not exist in metadata.")
    if (method=="dreamlet" & ! metadata_cell_count %in% colnames(metadata)) stop("Validation Error: metadata_cell_count column does not exist in metadata and is required for dreamlet.")
    if (method=="dreamlet" & is.null(cell_by)) stop("Validation Error: cell_by must be provided for dreamlet.")
    if (method=="dreamlet" & is.null(sample_by)) stop("Validation Error: sample_by must be provided for dreamlet.")

    # check min_frac is numeric and between 0 and 1
    if (!is.numeric(min_frac) | min_frac < 0 | min_frac > 1) stop("Error. min_frac of ", min_frac, " is not a numeric value between 0 and 1.")
    # check min_exp is numeric and a positive value
    if (!is.numeric(min_exp) | min_exp < 0) stop("Error. min_exp of ", min_exp, " is not a numeric value greater than 0.")
    # check return_model is boolean
    if (!is.logical(return_model)) stop("Validation Error: return_model of ", return_model, " is not a boolean value.")

    # Check dge_groups are in dge_by column
    if (!is.null(dge_groups) & !all(dge_groups %in% metadata[,dge_by])) {
        stop("Validation Error: The following 'dge_groups' are not levels of the 'dge_by' column: ", paste0(dge_groups[!dge_groups %in% metadata[,dge_by]], collapse = ", "))
    }

    if (!all(counts == as.integer(counts))) warning("Warning. The count matrix contains non-integer entries")
    if (!is.null(dge_by) & is.numeric(metadata[[dge_by]])){
        warning("Warning. You have provided a numeric variable as your grouping of interest. This will be converted to a factor when calculating group size.")
    }
    count_level <- "sample"
    if( any(duplicated( paste(metadata[,sample_by], metadata[,cell_by]) )) ) count_level <- "cell"
    warning(
        paste0(
            "Count data are assumed to be at the ", count_level, " level because ",
            if (count_level == "cell") {
                paste0("duplicate '", sample_by, "' values were found within ", "'", cell_by, "' groups.")
            } else {
                paste0("no duplicate '", sample_by, "' values were found within ", "'", cell_by, "' groups.")
            }
        )
    )


    if (!is.null(case_group) & !is.null(reference_group)) {
        .input_check_dge_case_ref(metadata, dge_by, case_group, reference_group)
        if (!is.null(contrasts)) {
            warning("Both contrasts and case-reference groups are given. 'contrasts' will be ignored.")
            contrasts <- NULL
        }
    } else if (is.null(contrasts)) {
        stop("Validation Error: Either case_group and reference_group or contrasts must be provided for DGE analysis.")
    } else {
        if (method %in% c('mast', 'deseq2')) {
            stop("Validation Error: 'contrasts' usage is not implemented for the requested 'method'.  Use 'case_group' and 'reference_group' instead.")
        }
        ### ToDo Utilize .input_checks_contrasts()
    }

    if (count_level == "sample" && method %in% c("mast", "memento")) {
        stop(
            method, " is intended for single-cell data, but the count data appear to be ",
            "at the sample level. Check that `", sample_by,
            "` correctly identifies samples and that the count matrix contains ",
            "single-cell observations."
        )
    }
    if (count_level == "cell" && ! (method %in% c("mast", "memento")) ) {
        stop(
            method, " is intended for sample-level data, but the count data appear to be ",
            "single-cell data. Check that `", sample_by,
            "` correctly identifies samples and that the count matrix contains ",

            "sample-level observations."
        )
    }


    # Make sure cell_targets include only those cell types that are in cell_by column
    cell_targets <- intersect( unique(cell_targets), unique(metadata[ , cell_by]) )
    if(length(cell_targets) == 0) {
        warning(paste0("None of the cell_targets values are in '", cell_by, "' column of metadata. Using all cell types in '", cell_by, "' column."))
        cell_targets <- NULL
    }


    if (identical(cell_targets, NULL)) {
        cell_targets <- unique(metadata[,cell_by])
    }

    ####### STUB, to be replaced with 'min_frac' usage!  And perhaps other per-cell checks.
    if (any(table(metadata[,cell_by])[cell_targets] <= 5)) {
        cells <- table(metadata[,cell_by])
        cell_targets <- names(cells)[cells>5]
    }


    # Set updated values in outs for sending returning
    outs <- list(
        counts = counts,
        metadata = metadata,
        contrasts = contrasts,
        cell_targets = cell_targets,
        cell_by = cell_by
    )


    outs
}

#' Check that case and reference groups for DGE analysis are valid.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as  numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @details Checks that dge_by column exists in metadata and that case and reference groups exist in the dge_by column. Errors if any of these conditions are not met.
.input_check_dge_case_ref <- function(metadata, dge_by, case_group, reference_group) {
    if(! dge_by %in% colnames(metadata) ) stop("Error. dge_by column does not exist in metadata.")
    if( ! case_group %in% metadata[,dge_by] ) stop("Error. Case group does not exist in dge_by column of metadata.")
    if( ! reference_group %in% metadata[,dge_by] ) stop("Error. Reference group does not exist in dge_by column of metadata.")
    if (case_group==reference_group) stop("Error. Case and reference groups cannot be the same.")
}

.remove_NA_cols <- function(metadata) {
    metadata[,!apply(metadata, 2, function(col) {all(is.na(col))})]
}

