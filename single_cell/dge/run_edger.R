library(edgeR)
library(tidyverse)
### DGE analysis using EdgeR using pseudobulk data
#
#' @param counts A feature (row) x sample (column) raw count matrix.  Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use for DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param return_model A boolean indicating whether the function should return the model (TRUE) or the data.frame of differential expression analysis results.#'
#'
#' @return DGE result data frame or model
run_edger <- function(counts, metadata, dge_by,
                      case_group, reference_group,
                      fixed_effects = NULL, return_model = FALSE) 
{
  y <- DGEList(counts, samples = metadata)
  y$samples[[dge_by]] <- factor(y$samples[[dge_by]], levels = c(reference_group, case_group))
  y <- calcNormFactors(y)
  
  formula_str = paste("~", dge_by)
  if(! is.null(fixed_effects)) {
    for(fe in fixed_effects) {
      formula_str = paste(formula_str, "+", fe)
    }
  }
  
  design <- model.matrix(as.formula(formula_str), data = y$samples)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  
  if(return_model) {
    return(fit)
  }
  
  res <- glmQLFTest(fit, coef = paste0(dge_by, case_group))
  res = topTags(res, n = Inf)$table
  res = arrange(res, FDR)
  
  res$AveExp=2^(res$logCPM)
  res = dplyr::rename(res, "log2FC"="logFC", "pval"="PValue", "padj"="FDR")
  return(res)
}
