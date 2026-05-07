library(edgeR)
library(tidyverse)

#' Run DGE between 2 groups using edgeR, iterating over cell types.
#' @param counts A feature (row) x sample (column) raw count matrix. Only samples intended for DGE analysis should be included in \code{counts}.
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}. Must contain a 'celltype' column.
#' @param dge_by Name of \code{metadata} column to use for DGE comparison. Must have exactly two levels.
#' @param case_group Level to be used as numerator in log2FC calculation.
#' @param reference_group Level to be used as denominator in log2FC calculation.
#' @param fixed_effects A vector of \code{metadata} column names to use as fixed effects.
#' 
run_edger <- function(counts, metadata, dge_by,
                      case_group, reference_group,
                      fixed_effects = NULL)
{
  metadata <- input_checks(counts, metadata, dge_by, case_group, reference_group, fixed_effects)
  
  # Ensure 'celltype' column exists
  if (!"celltype" %in% colnames(metadata)) {
    stop("metadata must contain a 'celltype' column.")
  }
  
  celltypes <- unique(metadata$celltype)
  
  results_list <- lapply(celltypes, function(ct) {
    # Subset metadata and counts to this cell type
    meta_sub <- metadata[metadata$celltype == ct, , drop = FALSE]
    counts_sub <- counts[, rownames(meta_sub), drop = FALSE]
    
    ## Create object
    y <- DGEList(counts_sub, samples = meta_sub)
    
    ## Order case vs reference
    y$samples[[dge_by]] <- as.factor(y$samples[[dge_by]])
    levels(y$samples[[dge_by]]) <- c(reference_group, case_group)
    
    # Filter low expression
    keep_genes <- filterByExpr(y, group = meta_sub[[dge_by]], min.prop = 0.6)
    y <- y[keep_genes, ]
    y <- calcNormFactors(y)
    
    ## Make formula string
    formula_str <- paste("~", dge_by)
    if (!is.null(fixed_effects)) {
      for (fe in fixed_effects) {
        formula_str <- paste(formula_str, "+", fe)
      }
    }
    
    ## Make contrast and design matrix
    design <- model.matrix(as.formula(formula_str), data = y$samples)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    res <- glmQLFTest(fit, coef = paste0(dge_by, case_group))
    
    # Extract and format DGE results
    res <- topTags(res, n = Inf)$table
    res <- arrange(res, FDR)
    
    ## Approximate AveMean
    res$AveExp <- 2^(res$logCPM)
    res <- dplyr::rename(res, "log2FC" = "logFC", "pval" = "PValue", "padj" = "FDR")
    
    res$celltype <- ct
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    return(res)
  })
  
  # Combine all cell type results into a single data frame
  combined_results <- bind_rows(results_list)
  
  return(combined_results)
}