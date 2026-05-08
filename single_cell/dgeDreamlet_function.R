source('helpers_load.R')

### This script provides a function to run dreamlet based on pseudobulked data.
#'
#' @param counts matrix of pseudobulked count data where columns are samples and rows are genes
#' @param metadata metadata where rows are the samples from the count matrix and columns are variables
#' @param dge_by NULL or a String naming the metadata column that represents the main grouping of interest. 
#'    This is used to identify the size of the smallest group for calculating which genes to retain.
#'    It is also included in the model as a fixed effect if it was not already listed.
#' @param case_group NULL or name of the case group to be used as the numerator in logFC
#' @param reference_group NULL or name of the reference group to be used as the denom in logFC
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#'    Example: list(c(contr_ssc = "diseaseSSC - diseaseHC"), c(contr_cl="diseaseCL - diseaseHC"))
#' @param dge_groups NULL or the vector of group names from dge_by column that should be used for filtering out
#     lowly expressed genes
#' @param fixed_effects String vectors naming metadata columns to treat as a fixed effects. this cannot be NULL.
#' @param random_effects NULL or String vectors naming metadata columns to treat as a random effects. 
#'    We recommend including batch as a random effect.
#' @param min_frac minimum proportion of samples required to have at least MIN.COUNT counts per gene (default=0.6)
#'    If dge.group.by is not NULL, then the minimum proportion of the smallest dge.group.by group will be used. 
#' @param MIN.COUNT minimum number of counts required per sample when calculating whether a gene should be retained (default = 5)
#' @param return_model boolean indicating whether the function should return the model or the DGE result
#'     data frame DEFAULT TO FALSE?
#'
#' @return DGE result data frame or model
run_dreamlet = function(counts, metadata, dge_by=NULL, case_group=NULL, reference_group=NULL,
	     contrasts=NULL, dge_groups=NULL,
	    fixed_effects=NULL, random_effects=NULL, 
       min_frac=0.6, MIN.COUNT=5, return_model=F){

  .input_checks(counts, metadata, dge_by, case_group, reference_group, contrasts, fixed_effects, random_effects)
 
  # TODO .input_checks_function.R
  # - add these
  # - add an input check for contrasts
  # - add an input check for return_model, min_frac
  # - MIN.COUNT?
  if (!all(counts == as.integer(counts))) warning("Warning. The count matrix contains non-integer entries")
  if (!is.null(dge_by) & is.numeric(metadata[[dge_by]])){
    warning("Warning. You have provided a numeric variable as your grouping of interest. This will be converted to a factor when calculating group size.")
  }
  
  counts = .remove_low_expression_genes(counts, metadata, dge_by, dge_groups, min_frac)

  # create an object from the pb counts and metadata
  pb = SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(counts)),
                                                    colData=metadata)

  # construct the fixed effects portion of the model
  model="~"
  if (!is.null(contrasts)){
    model= paste0(model, "0")
  } 
  for (fe in unique(c(fixed_effects, dge_by))) {
    model <- paste0(model, " + ", fe)
  }

  # dreamlet-required: drop samples that do not have sufficient cells and variables that are collinear
  res.proc = dreamlet::processAssays(pb, as.formula(model), min.count=MIN.COUNT, min.prop=min_frac)

  # add random effects
  for (rand in random_effects) {
      model <- paste0(model, " + (1 | ", rand, ")")
  }
  
  # TODO add option for running dreamlet w/o contrasts

  # fit model
  print(sprintf("Fitting model %s", model))
  res.dl = dreamlet::dreamlet(res.proc, as.formula(model), contrasts=contrasts)

  if (return_model) {
    return(res.dl)
  }

  # pull out and format DGEs
  dge.df = do.call(rbind,
               lapply(setdiff(dreamlet::coefNames(res.dl), "(Intercept)"), 
               function(x){
                 coef.df = tibble::as_tibble(limma::topTable(res.dl, x, number=nrow(counts)))
                 coef.df$contrast = x
                 return(coef.df)
               }))
  

  ## reformat the DGE output
  # log2FC, avgExpr, pval, padj, celltype, fracCase, fracRef
  return(dge.df)
}

run_dreamlet_within_cells = function(counts, metadata, cell_by, cell_targets, min_per_group=10, ...){
  # this function will run dreamlet within each cell type and add a column for cell type to the output

  cell_targets = .input_checks_within_cells(
        metadata,
        cell_by, cell_targets,
        dge_by, case_group, reference_group, contrasts, 
        min_per_group
    ) 
  # TODO:
  # .input_checks_within_cells needs to:
  #  - check that cell_targets are part of the cell_by column
  #  - also work if it's contrasts
  #  - check min_per_grp is numeric


  dge_list = lapply(cell_targets, function(ct){
    print(paste0("Running dreamlet for cell type ", ct))
    ct_idx = which(metadata[[cell_by]] == ct)
    ct_counts = counts[, ct_idx]
    ct_metadata = metadata[ct_idx, ]
    ct_dge = run_dreamlet(ct_counts, ct_metadata, ...)
    return(ct_dge)
  })
  names(dge_list) = cell_targets
  if (return_model) {
    return(dge_list)
  }

  dge.df = do.call(rbind, lapply(names(dge_list), function(x) dge_list[[x]] %>% mutate(cell_type = x)))
  return(dge.df)

}


