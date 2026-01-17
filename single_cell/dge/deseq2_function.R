library(edgeR)
library(DESeq2)
library(tidyverse)
##### DGE analysis using DESeq2
#' @param mat: feature (row) x sample (column) raw count matrix
#' @param meta: sample metadata, with rownames containing column names of mat
#' @param comparison: name of meta column to use DGE comparion. Must have exactly two levels
#' @param numerator: level to be used as numerator in log2FC calculation
#' @param denominator: level to be used as denominator in log2FC calculation
#' @param fixed_effects: a vector of meta column names to use as fixed effects
run_deseq2 <- function(mat, meta, comparison, numerator, denominator, fixed_effects=NULL) {
  if(! all( colnames(mat) %in% rownames(meta)) ) {
    stop("One or more samples in mat are not present in meta.")
  }
  # Subset and order metadata according to the order of samples in mat
  meta = meta[ match(colnames(mat), rownames(meta)) , ]
  
  if( length(unique(as.vector(meta[,comparison]))) != 2) {
    stop("Comparison column of meta does not have exactly two unique values.")
  }
  if( ! numerator %in% meta[,comparison] ) {
    stop("Numerator value does not exist in comparison column of meta.")
  }
  if( ! denominator %in% meta[,comparison] ) {
    stop("Denominator value does not exist in comparison column of meta.")
  }
  if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(meta)) ) {
    stop("One or more fixed_effects do not exist in meta.")
  }

  # Keep features that are expressed with CPM > 1 in more than 50% of samples in at least one comparison group.
  # feats = names(which(apply(cpm(mat) > 1, 1, function(x) { 
  #                                            any(tapply(x, meta[,comparison], mean) > 0.5) 
  #                                          } )
  #                    )
  #              )

  # Determine the number of sample in the smaller group (min_gsize) and keep features that are expressed with CPM > 1 in at least 60% of min_gsize samples, regardless of the sample group identity.
  min_gsize = min(table(meta[,comparison]))
  feats = names(which(rowSums(cpm(mat) > 1) > (min_gsize * 0.6)))
  print(paste0("Number of final features: ", length(feats)))
  mat = mat[feats,]
  
  # Make formula string
  formula_str = paste("~", comparison)
  if(! is.null(fixed_effects)) {
    for(fe in fixed_effects) {
      formula_str = paste(formula_str, "+", fe)
    }
  }
  
  # Perform DGE analysis
  dds = DESeqDataSetFromMatrix(mat, colData = meta, design = as.formula( formula_str ))
  dds = DESeq(dds)
  
  # Extract and format DGE results
  res = results(dds, contrast = c(comparison, numerator, denominator)) %>% as.data.frame()
  res = arrange(res, padj)
  
  res = dplyr::rename(res, "logFC"="log2FoldChange", "AveExpr"="baseMean", "P.Value"="pvalue", "adj.P.Val"="padj")
  return(res)
}


