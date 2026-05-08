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
.input_checks <- function(counts, metadata, dge_by, case_group=NULL, reference_group=NULL, contrasts=NULL, fixed_effects=NULL, random_effects=NULL) {
  if (!is.matrix(counts) & !is.data.frame(counts) ) stop("Error. Input data is not a matrix or a data.frame.")
  if (ncol(counts)==0 | nrow(counts)==0) stop("Error. Count matrix does not have at least one column and row.")
  if(! dge_by %in% colnames(metadata) ) stop("Error. dge_by column does not exist in metadata.")
  if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) stop("Error. One or more fixed_effects do not exist in metadata.")
  if(! is.null(random_effects) & ! all(random_effects %in% colnames(metadata)) ) stop("Error. One or more random_effects do not exist in metadata.")
  if(! all( colnames(counts) == rownames(metadata)) ) stop("Error. Column names of counts and rownames of metadata do not match.")

  # Subset and order metadata according to the order of samples in counts
  #metadata = metadata[ match(colnames(counts), rownames(metadata)), ]

  case_ref = FALSE
  cntrst = FALSE
  if( !is.null(case_group) & !is.null(reference_group) ) {
    case_ref = TRUE
    if( ! case_group %in% metadata[,dge_by] ) stop("Error. Case group does not exist in dge_by column of metadata.")
    if( ! reference_group %in% metadata[,dge_by] ) stop("Error. Reference group does not exist in dge_by column of metadata.")
  }
  if(! is.null(contrasts) )
    cntrst = TRUE
  if(case_ref & cntrst) warning("Both contrasts and case-reference groups are given. The contrasts will be ignored for the analysis.")
}

