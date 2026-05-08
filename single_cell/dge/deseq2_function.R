library(edgeR)
library(DESeq2)
library(tidyverse)
### DGE analysis using DESeq2
#
#' Run DGE between 2 groups using DESeq2.
#' @param counts A feature (row) x sample (column) raw count matrix.  Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param dge_by Name of \code{metadata} column to use for DGE comparion. Must have exactly two levels.
#' @param case_group Group to be used as numerator in log2FC calculation.
#' @param reference_group Group to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' @param return_model A boolean indicating whether the function should return the model (TRUE) or the data.frame of differential expression analysis results.
run_deseq2 <- function(counts, metadata, dge_by, case_group, reference_group, fixed_effects=NULL, return_model = FALSE) {

  # Commong input checks
  contrasts = .input_checks(counts, metadata, dge_by, case_group, reference_group, fixed_effects)

  # Method specific input checks
  # None for DESeq2

  counts = .remove_low_expression(counts=counts, metadata=metadata, dge_by=dge_by, dge_groups=c(case_group, reference_group), min_frac=0.6)
  
  # Make formula string
  formula_str = paste("~", dge_by)
  if(! is.null(fixed_effects)) {
    for(fe in fixed_effects) {
      formula_str = paste(formula_str, "+", fe)
    }
  }
  
  # Perform DGE analysis
  dds = DESeqDataSetFromMatrix(counts, colData = metadata, design = as.formula( formula_str ))
  dds = DESeq(dds)
  
  # Extract and format DGE results
  res = results(dds, contrast = c(dge_by, case_group, reference_group)) %>% as.data.frame()
  res = arrange(res, padj)
  
  res = dplyr::rename(res, "log2FC"="log2FoldChange", "aveExpr"="baseMean", "pval"="pvalue", "padj"="padj")

  if(return_model) {
    return(dds)
  } else {
    return(res)
  }
}


