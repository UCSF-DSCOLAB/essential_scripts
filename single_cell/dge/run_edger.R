#-- DGE analysis using EdgeR


## To do:
##--- Add make contrast
## --- pseudobulk here vs not 
##---- Single_cell experiment, cold accept a seurat objects
##--- All categorical values, could accept integers

#' @param object: a singlecellexperiment
#' @param metadata: sample metadata, with rownames containing column names of \code{counts_matrix}
#' @param pseudobulk_by:  vector of \code{metadata}  metadata name holding clusters or annotations to pseudibulk with
#' @param cells_group_by: a vector of \code{metadata}  metadata name holding clusters or annotations to explore within
#' @param skip_cell_type: a vector of \code{metadata}  cells.group.by to skip
#' @param dge_by: name of \code{metadata} column to use DGE comparion. Must have exactly two levels. 
#' @param fixed_effects: a vector of \code{metadata} column names to use as fixed effects


library(edgeR)
library(scater)
library(scuttle)

run_edger <- function(object, metadata, pseudobulk_by,
                      cells_group_by, dge_by, case_group, reference_group,
                      fixed_effects = NULL, skip_cells = NULL) {
  
  celltypes <- unique(metadata[[cells_group_by]])
  if (!is.null(skip_cells)) {
    celltypes <- setdiff(celltypes, skip_cells)
  }
  
  counts_matrix <- scuttle::aggregateAcrossCells(object, 
                                                 id = colData(object)[, c(cells_group_by, pseudobulk_by)])
  
  deg_list <- list()
  
  for (label in celltypes) {
    message("Processing: ", label)
    
    counts_subset <- counts_matrix[, metadata[[cells_group_by]] == label]
    y <- DGEList(counts(counts_subset), samples = colData(counts_subset))
    
    ncells <- colData(counts_subset)$ncells
    keep_cells <- ncells >= 10
    y <- y[, keep_cells]
    
    keep_genes <- filterByExpr(y, group = colData(y)[[dge_by]])
    y <- y[keep_genes, ]
    
    y <- calcNormFactors(y)
    
    formula_str <- paste("~", paste(c(paste0("factor(", dge_by, ")"), 
                                      paste0("factor(", fixed_effects, ")")), collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = y$samples)
    
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    
    coef_index <- which(colnames(design) == paste0("factor(", dge_by, ")1"))
    res <- glmQLFTest(fit, coef = coef_index)
    
    deg_list[[label]] <- topTags(res, n = Inf)
  }
  
  return(deg_list)
}
