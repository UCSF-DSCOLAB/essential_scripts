### This script provides a function to run dreamlet based on pseudobulked data.
#'
#' @param counts matrix of pseudobulked count data where columns are samples and rows are genes
#' @param metadata metadata where rows are the samples from the count matrix and columns are variables
#' @param cell_by String naming the metadata column that represents the cell type or cluster.
#' @param sample_by String naming the metadata column that represents the sample.
#' @param cell_target String naming the cell type or cluster to be used for DGE analysis.
#' @param metadata_cell_count String naming the metadata column that contains the number of cells in each pseudobulk.
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
#' @param return_model boolean indicating whether the function should return the model or the DGE result
#'     data frame DEFAULT TO FALSE?
#'
#' @return DGE result data frame or model
run_dreamlet = function(counts, metadata, cell_by, sample_by, 
       cell_target, metadata_cell_count ,
       dge_by=NULL,  case_group=NULL, reference_group=NULL,
	     contrasts=NULL, dge_groups=NULL,
	     fixed_effects=NULL, random_effects=NULL,return_model=F){


  # create an object from the pb counts and metadata
  assay_list = list(as.matrix(counts))
  assay_list = setNames(assay_list, cell_target)

  pb = SingleCellExperiment::SingleCellExperiment(assays=assay_list,colData=metadata)

  int_colData(pb)$n_cells = lapply(colData(pb)[,metadata_cell_count], function(x) {
    y=list(x); y=setNames(y, cell_target); return(y)})
  names(int_colData(pb)$n_cells) = rownames(colData(pb))
  metadata(pb)$agg_pars <- list(assay = as.character(cell_target), by = c(cell_by, sample_by), fun = "sum", scale = FALSE)
  metadata(pb)$aggr_means = 
        tidyr::as_tibble(colData(pb)[,c(cell_by, sample_by)]) 
  metadata(pb)$aggr_means[,sample_by]=rownames(colData(pb))

  # construct the fixed effects portion of the model
  model="~"
  if (!is.null(contrasts)){
    model= paste0(model, "0")
  } 
  for (fe in unique(c(fixed_effects, dge_by))) {
    model <- paste0(model, " + ", fe)
  }

  # process the pseudobulked data for dreamlet
  # do not perform filtering here as this has already been done
  res.proc = dreamlet::processAssays(pb, as.formula(model), min.prop=0, min.cells=1, min.count=1, min.samples=1)

  # add random effects
  for (rand in random_effects) {
      model <- paste0(model, " + (1 | ", rand, ")")
  }
  
  print(sprintf("Fitting model %s", model))
  if (!is.null(case_group) & !is.null(reference_group)) {
    res.dl = dreamlet::dreamlet(res.proc, as.formula(model))
    if (return_model) return(res.dl)
    
    var_str = paste0(dge_by, case_group)
    dge.df = dreamlet::topTable(res.dl, var_str, number=nrow(counts)) |>
      tidyr::as_tibble() |>
      dplyr::select(ID, logFC, AveExpr, P.Value, adj.P.Val) 
      
  } else {
    res.dl = dreamlet::dreamlet(res.proc, as.formula(model), contrasts=contrasts)
    if (return_model) return(res.dl)
    
    dge.df = do.call(rbind, lapply(contrasts, function(x){
                 coef.df = tibble::as_tibble(dreamlet::topTable(res.dl, x, number=nrow(counts)))
                 coef.df$contrast = x
                 return(coef.df)})) |> tidyr::as_tibble() |>
      dplyr::select(ID, logFC, AveExpr, P.Value, adj.P.Val, contrast) 
  }

  dge.df = dplyr::rename(dge.df, gene=ID, log2FC = logFC, avgExpr = AveExpr, pval = P.Value, padj = adj.P.Val)
  return(dge.df)
}
