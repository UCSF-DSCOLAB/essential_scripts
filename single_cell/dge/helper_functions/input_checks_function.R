### Check inputs to DGE functions
#
#' Check inputs to DGE functions, match \code{metadata} rows to \code{counts} columns, and return \code{metadata}.
#' @param counts: feature (row) x sample (column) raw count matrix.  Only samples intended for DGE analysis should be included in \code{counts}
#' @param metadata: sample metadata, with rownames containing column names of \code{counts}
#' @param dge_by: name of \code{metadata} column to use DGE comparion. Must have exactly two levels
#' @param case_group: level to be used as numerator in log2FC calculation
#' @param reference_group: level to be used as denominator in log2FC calculation
#' @param fixed_effects: a vector of \code{metadata} column names to use as fixed effects
#' @param random_effects: a vector of \code{metadata} column names to use as random effects
#' @return \code{metadata} after matching the rows to the columns (samples) in \code{counts}
input_checks <- function(counts, metadata, dge_by, case_group, reference_group, fixed_effects=NULL, random_effects=NULL) {
  if (!is.matrix(counts) & !is.data.frame(counts) ) stop("Error. Input data is not a matrix or a data.frame.")
  if (ncol(counts)==0 | nrow(counts)==0) stop("Error. Count matrix does not have at least one column and row.")
  if(! dge_by %in% colnames(metadata) ) stop("Error. dge_by colum does not exist in metadata.")
  if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) stop("Error. One or more fixed_effects do not exist in metadata.")
  if(! is.null(random_effects) & ! all(random_effects %in% colnames(metadata)) ) stop("Error. One or more random_effects do not exist in metadata.")
  if(! all( colnames(counts) %in% rownames(metadata)) ) stop("Error. One or more samples in counts are not present in metadata rownames.")

  # Subset and order metadata according to the order of samples in counts
  metadata = metadata[ match(colnames(counts), rownames(metadata)) , ]

  if( length(unique(as.vector(metadata[,dge_by]))) != 2) stop("Error. dge_by column of metadata does not have exactly two levels (unique values).")
  if( ! case_group %in% metadata[,dge_by] ) stop("Case group does not exist in dge_by column of metadata.")
  if( ! reference_group %in% metadata[,dge_by] ) stop("Reference group does not exist in dge_by column of metadata.")
  return(metadata)
}

