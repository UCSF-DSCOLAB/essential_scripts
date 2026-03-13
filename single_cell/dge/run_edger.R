#-- DGE analysis using EdgeR

#' @param counts: feature (row) x sample (column) raw count matrix.  Only samples intended for DGE analysis should be included in \code{counts}
#' @param metadata: sample metadata, with rownames containing column names of \code{counts}
#' @param dge_by: name of \code{metadata} column to use DGE comparion. Must have exactly two levels
#' @param case_group: level to be used as numerator in log2FC calculation
#' @param reference_group: level to be used as denominator in log2FC calculation
#' @param fixed_effects: a vector of \code{metadata} column names to use as fixed effects


library(edgeR)
library(tidyverse)


run_edger <- function(counts, metadata, dge_by,
                      case_group, reference_group,
                      fixed_effects = NULL) 
{
  metadata = input_checks(counts, metadata, dge_by, case_group, reference_group, fixed_effects)
  
  ## Create object
  y <- DGEList(counts, samples = metadata)
  ## Order case vs reference
  
  y$samples[[dge_by]]=as.factor(y$samples[[dge_by]])
  levels(y$samples[[dge_by]]) = c(reference_group, case_group)    
  
  # Filter low expression
  keep_genes <- filterByExpr(y, group = metadata[[dge_by]], min.prop=.6)
  y <- y[keep_genes, ]
  y <- calcNormFactors(y)
  
  ## Make formula string
  formula_str = paste("~", dge_by)
  if(! is.null(fixed_effects)) {
    for(fe in fixed_effects) {
      formula_str = paste(formula_str, "+", fe)
    }
  }
  
  ## Make contrast and desing matrix
  design <- model.matrix(as.formula(formula_str), data = y$samples)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, coef = paste0(dge_by, case_group))
  return(topTags(res, n = Inf))
}
