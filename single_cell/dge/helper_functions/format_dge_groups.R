#' Format DGE groups
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use DGE comparison. Must have exactly two levels.
#' @param case_group Group to be used as  numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @details Converts the dge_by column to a factor if it is not already and orders case/ref group. 
.format_dge_groups = function(metadata, dge_by, case_group, reference_group) {

  # TODO - need to confirm that all methods use the same ordering of case and ref groups
  if (!is.null(dge_by) ) {
    if (!is.null(case_group) & !is.null(reference_group)) {
      curr_levels = unique(metadata[[dge_by]])
      metadata[[dge_by]] = factor(metadata[[dge_by]], levels = c(case_group, reference_group, setdiff(curr_levels, c(case_group, reference_group))))
    } else {
      metadata[[dge_by]] = as.factor(metadata[[dge_by]])
    }
  }
  return(metadata)
}