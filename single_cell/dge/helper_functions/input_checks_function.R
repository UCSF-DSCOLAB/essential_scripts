### Check inputs to DGE functions
#
#' Check inputs to DGE functions, match \code{metadata} rows to \code{counts} columns, and return \code{metadata}.
#' @param counts A feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param random_effects A vector of \code{metadata} column names to use as random effects.
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#' @param min_frac Minimum proportion of samples required to have at least MIN.COUNT counts per gene.
#' @param return_model Boolean indicating whether the function should return the model or the DGE result data frame.
.input_checks <- function(counts, metadata, dge_by, case_group=NULL, reference_group=NULL, contrasts=NULL, fixed_effects=NULL, random_effects=NULL, min_frac=0.6, return_model=FALSE) {
  if (!is.matrix(counts) & !is.data.frame(counts) ) stop("Error. Input data is not a matrix or a data.frame.")
  if (ncol(counts)==0 | nrow(counts)==0) stop("Error. Count matrix does not have at least one column and row.")
  if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) stop("Error. One or more fixed_effects do not exist in metadata.")
  if(! is.null(random_effects) & ! all(random_effects %in% colnames(metadata)) ) stop("Error. One or more random_effects do not exist in metadata.")
  if(! all( colnames(counts) == rownames(metadata)) ) stop("Error. Column names of counts and rownames of metadata do not match.")

  # check min_frac is numeric and between 0 and 1
  if (!is.numeric(min_frac) | min_frac < 0 | min_frac > 1) stop("Error. min_frac of ", min_frac, " is not a numeric value between 0 and 1.")
  # check return_model is boolean
  if (!is.logical(return_model)) stop("Error. return_model of ", return_model, " is not a boolean value.")

  # Subset and order metadata according to the order of samples in counts
  #metadata = metadata[ match(colnames(counts), rownames(metadata)), ]
  if (!is.null(case_group) | !is.null(reference_group)) {
    .input_check_dge_case_ref(metadata, dge_by, case_group, reference_group) 
  } else if (is.null(contrasts)) {
    stop("Error. Either case_group and reference_group or contrasts must be provided for DGE analysis.")  
  }

  if(! is.null(contrasts) & (!is.null(case_group) & !is.null(reference_group) ))  warning("Both contrasts and case-reference groups are given. The contrasts will be ignored for the analysis.")

  if (!all(counts == as.integer(counts))) warning("Warning. The count matrix contains non-integer entries")
  if (!is.null(dge_by) & is.numeric(metadata[[dge_by]])){
    warning("Warning. You have provided a numeric variable as your grouping of interest. This will be converted to a factor when calculating group size.")
  }
  
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
}

