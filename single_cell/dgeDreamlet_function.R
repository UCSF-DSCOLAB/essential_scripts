### This script provides a function to run dreamlet based on pseudobulked data.
#'
#' @param pb.count.matrix matrix of pseudobulked count data where columns are samples and rows are genes
#' @param metadata metadata where rows are the samples from the count matrix
#' @param fixed.effects String vectors naming metadata columns to treat as a fixed effects. this cannot be NULL.
#' @param dge.group.by NULL or a String naming the metadata column that represents the main grouping of interest. 
#'    This is used to identify the size of the smallest group for calculating which genes to retain.
#'    It is also included in the model as a fixed effect if it was not already listed.
#' @param random.effects NULL or String vectors naming metadata columns to treat as a random effects. 
#'    We recommend including batch as a random effect.
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#'    Example: list(c(contr_ssc = "diseaseSSC - diseaseHC"), c(contr_cl="diseaseCL - diseaseHC"))
#' @param MIN.PROP minimum proportion of samples required to have at least MIN.COUNT counts per gene (default=0.6)
#'    If dge.group.by is not NULL, then the minimum proportion of the smallest dge.group.by group will be used. 
#' @param MIN.COUNT minimum number of counts required per sample when calculating whether a gene should be retained (default = 5)
#'
#' @return a tibble with the DGE results for all genes, coefficients, and contrasts
run_dreamlet = function(pb.count.matrix, metadata, 
	    fixed.effects, dge.group.by=NULL, random.effects=NULL, 
      contrasts=NULL, MIN.PROP=0.6, MIN.COUNT=5){

  # input checks
  if (!is.matrix(pb.count.matrix) ) stop("Error. Input data is not a matrix.")
  if (ncol(matrix)==0 | nrow(matrix)==0) stop("Error. Matrix does not have at least one column and row.")
  if (!all(rownames(metadata)==colnames(pb.count.matrix))) stop("Error. The column names of the count matrix do not match the metadata")
  if (!all(c(dge.group.by, fixed.effects, random.effects) %in% colnames(metadata))) stop("Error. The listed variables are not available in the metadata.")
  if (is.null(fixed.effects) & is.null(dge.group.by)) stop("No fixed effects or grouping variable provided, so model cannot be fit")

  # warnings
  if (!all(pb.count.matrix == as.integer(pb.count.matrix))) warning("Warning. The count matrix contains non-integer entries")
  if (!is.null(dge.group.by) & is.numeric(metadata$dge.group.by)){
    warning("Warning. You have provided a numeric variable as your grouping of interest. This will be converted to a factor when calculating group size.")
  }

  # create a pseudobulked object
  pb = SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(pb.count.matrix)),
                                                    colData=metadata)

  # construct the fixed effects portion of the model
  model="~"
  if (!is.null(contrasts)){
    model= paste0(model, "0")
  } 
  for (fe in unique(c(fixed.effects, dge.group.by))) {
    model <- paste0(model, " + ", fe)
  }

  # calculate the minimum proportion to keep
  if (!is.null(dge.group.by)){
    print("Calculating the minimum proportion of samples required to keep a gene based on the size of the smallest group")
    MIN.PROP = MIN.PROP*(min(table(as.factor(metadata$dge.group.by)))/nrow(metadata))
  }

  # drop samples that do not have sufficient cells and variables that are collinear
  res.proc = dreamlet::processAssays(pb, as.formula(model), min.count=MIN.COUNT, min.prop=MIN.PROP)

  # add random effects
  for (rand in random.effects) {
      model <- paste0(model, " + (1 | ", rand, ")")
  }
  
  # fit model
  print(sprintf("Fitting model %s", model))
  res.dl = dreamlet::dreamlet(res.proc, as.formula(model), contrasts=contrasts)

  # pull out and format DGEs
  dge.df = do.call(rbind,
               lapply(setdiff(dreamlet::coefNames(res.dl), "(Intercept)"), 
               function(x){
                 coef.df = tibble::as_tibble(limma::topTable(res.dl, x, number=nrow(pb.count.matrix)))
                 coef.df$contrast = x
                 return(coef.df)
               }))

  return(dge.df)
}


