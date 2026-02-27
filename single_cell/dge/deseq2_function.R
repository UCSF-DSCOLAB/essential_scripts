library(edgeR)
library(DESeq2)
library(tidyverse)
##### DGE analysis using DESeq2
#' @param counts: feature (row) x sample (column) raw count matrix
#' @param metadata: sample metadata, with rownames containing column names of \code{counts}
#' @param dge_by: name of \code{metadata} column to use DGE comparion. Must have exactly two levels
#' @param case_group: level to be used as numerator in log2FC calculation
#' @param reference_group: level to be used as denominator in log2FC calculation
#' @param fixed_effects: a vector of \code{metadata} column names to use as fixed effects
run_deseq2 <- function(counts, metadata, dge_by, case_group, reference_group, fixed_effects=NULL) {
  if(! all( colnames(counts) %in% rownames(metadata)) ) {
    stop("One or more samples in counts are not present in metadata.")
  }
  # Subset and order metadata according to the order of samples in counts
  metadata = metadata[ match(colnames(counts), rownames(metadata)) , ]
  
  if( length(unique(as.vector(metadata[,dge_by]))) != 2) {
    stop("Comparison column of metadata does not have exactly two unique values.")
  }
  if( ! case_group %in% metadata[,dge_by] ) {
    stop("Case group does not exist in dge_by column of metadata.")
  }
  if( ! reference_group %in% metadata[,dge_by] ) {
    stop("Reference group does not exist in dge_by column of metadata.")
  }
  if(! is.null(fixed_effects) & ! all(fixed_effects %in% colnames(metadata)) ) {
    stop("One or more fixed_effects do not exist in metadata.")
  }

  # Keep features that are expressed with CPM > 1 in more than 50% of samples in at least one comparison group.
  # feats = names(which(apply(cpm(counts) > 1, 1, function(x) { 
  #                                            any(tapply(x, metadata[,dge_by], mean) > 0.5) 
  #                                          } )
  #                    )
  #              )

  # Determine the number of sample in the smaller group (min_gsize) and keep features that are expressed with CPM > 1 in at least 60% of min_gsize samples, regardless of the sample group identity.
  min_gsize = min(table(metadata[,dge_by]))
  feats = names(which(rowSums(cpm(counts) > 1) > (min_gsize * 0.6)))
  print(paste0("Number of final features: ", length(feats)))
  counts = counts[feats,]
  
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
  
  res = dplyr::rename(res, "logFC"="log2FoldChange", "AveExpr"="baseMean", "P.Value"="pvalue", "adj.P.Val"="padj")
  return(res)
}


