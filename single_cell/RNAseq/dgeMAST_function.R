##### This script initializes a function called 'dgeMAST' for running differential gene expression
# This is an adaptation by Dan Bunis of code originally formulated by Ravi Patel
#
# The function loop across multiple clusters running DGE between 2 groups, gathering eveything into a single data.frame.
# Think like FindAllMarkers, but for a biological question per each cluster.
#
# Some key notes:
# - This is a method meant for droplet-based single-cell transcriptome data.
# - Raw counts (integers) are expected to exist and to be:
#   - in the 'counts' slot of a Seurat object's 'RNA' assay.
#   - in an assay named 'counts' if using a SingleCellExperiment object.
#
# To use, `source("path/to/this/file")`
# Then, something like:
# dge <- dgeMAST(
#   object = sobj, # Seurat or SCE object
#   cells.group.by = 'full_annots', # metadata name holding clusters or annotations to explore within
#   dge.group.by = 'sex_at_birth', # metadata name holding the samples'/cell's groupings to compare across
#   dge.group.1 = 'M', # value of 'dge.group.by' metadata value of cells in the primary group.
#   dge.group.2 = 'F', # value of 'dge.group.by' metadata value of cells in the secondary / control group.
#   cells.group.targets = c('T4 naive', 'T8 em') # Without this, it will run within ALL cell types / clusters (a.k.a. all 'full_annots' per this example)
# )
# See roxygen notes below for additional details!

suppressPackageStartupMessages({
  library(assertthat)
  library(dplyr)
  library(Seurat)
  library(SingleCellExperiment)
  library(MAST)
  library(data.table)
  library(dittoSeq)
})

#' Run DGE between 2 groups while looping through cell types
#' @param object A Seurat or SingleCellExperiment (SCE) object
#' @param cells.group.by String naming a metadata column of \code{object} which holds clusters, annotations, or other cell groupings to explore within
#' @param cells.group.targets NULL or a string vector naming particular levels of the \code{cell.group.by}-data to run DGE within.  When left as \code{NULL}, all cell groupings will be targeted
#' @param dge.group.by String naming a metadata column of \code{object} which holds the sample or cell groupings that you wish to compare between
#' @param dge.group.1,dge.group.2 Strings naming the groups of \code{dge.group.by}-data to compare between. Directionality: Positive log2FC (foldChange) will mean upregulation in \code{dge.group.1} with respect to \code{dge.group.2}.
#' @param cells.use NULL or a Logical vector (of length equal to the number of cells in \code{object}) indicating cells to keep (TRUE) or remove (FALSE) before calculating DGE. Can be useful for ignoring specific samples.
#' @param log.prefix String to append at the beginning of log messages. Useful when wrapping this function inside an external loop.
#' @param mast.freq.expressed.min A fraction between 0 and 1, default = 0.2, which sets the minimal percent of cells that must express a gene for it to be considered. MAST performs less well for genes with low percent expression. This cutoff is run per each considered \code{cell.group.targets} cell grouping, and only expression among cells of the targeted \code{dge.group.1} and \code{dge.group.2} groups are considered.
#' @param add.random.effect,add.fixed.effect # NULL or String vectors naming metadata columns of \code{object} to treat as a random or fixed effects, respectively, in the mixed effect model used for DGE calculation. Additional details:
#' \itemize{
#' \item NOTE: \code{add.random.effect} defaults to "orig.ident" with the assumption that this column generally 1) exists and 2) holds batch information that is normally desired to be modeled as a random effect. Set to NULL, or something different, to turn this off.
#' \item \code{add.random.effect} must target discrete metadata and will be turned into a factor before use.
#' \item 'ngeneson' (the number of expressed genes per cell, which will be calculated internally) and \code{dge.group.by} are automatically added as fixed effects and do not need to be added here.
#' }
#' 
#' @details This function loops across multiple clusters running DGE between 2 groups, gathering eveything into a single data.frame.
#' Think like FindAllMarkers, but for a biological question per each cluster.
#' 
#' Raw counts (integers) are expected to exist and to be in the 'counts' slot of the 'RNA' assay if \code{object} is a Seurat object, or in an assay named 'counts' if \code{object} is a SingleCellExperiment object.
#' 
#' In preparation:
#' 
#' Seurat data is converted to a SingleCellObject (SCE).
#' The SCE is trimmed to only the cells indicated by \code{cells.use}, and only those where \code{dge.group.by} is either \code{dge.group.1} or \code{dge.group.2}.
#' And the set of cell targets given by \code{cells.group.targets} is trimmed of any targets that don't have any cells of code{dge.group.1} or \code{dge.group.2}
#' 
#' Then, the function loops through remaining \code{cells.group.targets}, performing these steps on just that subset of cells:
#' 
#' Log and read-depth normalization of the counts data,
#' trimming genes expressed in fewer than \code{mast.freq.expressed.min} proportion of cells,
#' calculation of the number of genes expressed per cell ('ngeneson'),
#' the DGE calculation by MAST with its zlm function with \code{method='glmer', ebayes = F, strictConvergence = FALSE, fitArgsD = list(nAGQ = 0)} and a \code{formula} based on the contents of \code{add.random.effect} and \code{add.fixed.effect}, followed by running the liklihood ratio test on the contrast of interest.
#' By default, the formula ammounts to \code{~ <dge.group.by> + ngeneson + (1 | orig.ident)}.
#' FDR-correction is applied to p-values at this stage -- per each cell type.
#' 
#' It then renames columns, converts the column holding fold change information to be base2, and adds an "up_in" column to help explain which code{dge.group.1} or \code{dge.group.2} was found to have higher expression for each gene.
#' @return a data.frame
#' @author Daniel Bunis and Ravi Patel
dgeMAST <- function(
  object,
  cells.group.by, # metadata name holding clusters or annotations to explore within
  cells.group.targets = NULL, # Optional but RECOMMENDED. String ("T4 naive" / "0") or strings ('c("T4 naive", "T8 naive")' / 'c(1,3)'). If not changed from the default here, will loop through all annotations/clusters given by 'cells.group.by'.
  dge.group.by, # metadata column holding samples'/cell's groupings
  dge.group.1, # The "primary" group = when foldChange is positive, it means it's up in this population. 
  dge.group.2, # The "other" group
  cells.use = NULL, # NULL or Logical vector (length equal to the number of cells in object) indicating cells to keep (TRUE) or remove (FALSE) before calculating DGE.
  log.prefix = '', # String to append at the beginning of log messages.  Useful when wrapping this function inside an external loop.
  mast.freq.expressed.min = 0.2, # Fraction. A gene filter before dge calculation: Minimal percent of cells (among targets groups, within each targetted cell set) that must express a gene for it to be considered.
  add.random.effect = "orig.ident", # NULL or a String vector naming a metadata to treat as a random effects in the mixed effect model used for DGE calculation. These metadata must contain discrete data and will be turned into a factor before use. Providing batch information is highly recommended, and is what the default of "orig.ident" aims to capture. 
  add.fixed.effect = NULL # NULL or a String vector naming a metadata to treat as fixed effects in the mixed effect model used for DGE calculation. Note: You do not need to provide 'dge.group.by' here, and the number of expressed genes will be calculated and stored/used as 'ngeneson'; Both will be added as fixed effects and do not need to be named here.
) {

  log_message <- function(...) {
    cat("[", format(Sys.time()), "]", ..., "\n")
  }

  if (!isMeta(cells.group.by, object)) stop("'cells.group.by', ", cells.group.by,", is not a metadata of 'object'")
  if (!isMeta(dge.group.by, object)) stop("'dge.group.by', ", dge.group.by,", is not a metadata of 'object'")
  if (!dge.group.1 %in% object[[dge.group.by, drop = TRUE]]) stop("No cells have 'dge.group.1', ", dge.group.1, ", as their 'dge.group.by' value.")
  if (!dge.group.2 %in% object[[dge.group.by, drop = TRUE]]) stop("No cells have 'dge.group.2', ", dge.group.2, ", as their 'dge.group.by' value.")
  if (!is.null(cells.group.targets) && !all(cells.group.targets %in% object[[cells.group.by, drop = TRUE]])) stop("At least one target of 'cells.group.targets' does not exist within 'cells.group.by' values.")

  # Tranform minimal bits into an SCE
  if (!is(object, "SingleCellExperiment")) {
    sce <- as.SingleCellExperiment(DietSeurat(object, assays="RNA"), assay="RNA")
  } else {
    sce <- object
  }

  # Trim to only cells passing cells.use
  if (!is.null(cells.use)) {
    sce <- sce[, dittoSeq:::.which_cells(
      cells.use,
      sce)]
  }

  # Trim to only target dge.groups, and use factor levels to have it target the right groups.
  sce <- sce[, dittoSeq:::.which_cells(
    meta(dge.group.by, sce) %in% c(dge.group.1, dge.group.2),
    sce)]
  sce[[dge.group.by]] <- factor(
    meta(dge.group.by, sce),
    levels = c(dge.group.2, dge.group.1)
  )

  # Establish cell targets if not given
  if (is.null(cells.group.targets)) {
    cell.group.targets <- metaLevels(cells.group.by, sce)
  }

  # Check for any lack of cell groups, and trim cell targets if so
  rms <- c()
  grps <- dittoSeq::meta(dge.group.by, sce)
  cts <- as.character(dittoSeq::meta(cells.group.by, sce))
  for (cell_targ in cells.group.targets) {
    if (any(table(grps[cts==cell_targ])==0)) {
      log_message(log.prefix, "Skipping cell type ", cell_targ, " because a dge.group has no cells")
      rms <- c(rms, cell_targ)
    }
  }
  cells.group.targets <- cells.group.targets[!cells.group.targets %in% rms]
  rm(cell_targ)

  # Loop through cell types, building dge output for each
  dge_all_list <- lapply(
    cells.group.targets,
    function(cell_targ) {

      log_message(log.prefix, "Running MAST dge for ", cell_targ, " cells")
      # Trim to just these cells
      sce_ct <- sce[, dittoSeq:::.which_cells(
        meta(cells.group.by, sce) == cell_targ,
        sce)]

      counts <- assay(sce_ct, "counts")
      libsizes <- colSums(counts)
      size.factors <- libsizes/mean(libsizes)
      logcounts(sce_ct) <- log2(t(t(counts)/size.factors) + 1)

      ### Ensure logcounts comes first!  ...by simply removing all others.
      for (assay in assayNames(sce_ct)) {
        if (assay!="logcounts") {
          assay(sce_ct, assay) <- NULL
        }
      }
      stopifnot(assayNames(sce_ct)[1]=="logcounts")

      # Transform into MAST's required SingleCellAssay structure.
      sca_ct <- SceToSingleCellAssay(sce_ct)

      # Trim to only genes expressed in enough cells
      expressed_genes <- freq(sca_ct) > mast.freq.expressed.min
      sca_ct <- sca_ct[expressed_genes,]

      # Add meta data for num genes captured per cell
      cdr2_ct <- colSums(assay(sca_ct)>0)
      colData(sca_ct)$ngeneson <- scale(cdr2_ct)

      # Build model
      model <- paste0("~ ngeneson + ", dge.group.by)
      # Ensure random effects, like batch, are factors too
      for (rand in add.random.effect) {
        colData(sca_ct)[,rand] <- as.factor(colData(sca_ct)[,rand])
        model <- paste0(model, " + (1 | ", rand, ")")
      }
      for (fix in add.fixed.effect) {
        model <- paste0(model, " + fixed")
      }

      # Run MAST
      contrast_val <- paste0(dge.group.by, dge.group.1)
      zlmCond <- zlm(
        as.formula(model),
        sca_ct, exprs_value = 'logcounts',
        method='glmer', ebayes = F,
        strictConvergence = FALSE,
        fitArgsD = list(nAGQ = 0))
      summaryCond <- summary(zlmCond, doLRT=contrast_val)$datatable
      fcHurdle <- merge(
        summaryCond[contrast==contrast_val & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryCond[contrast==contrast_val & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
      fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

      # Rebase coef toward log2() instead of ln()
      fcHurdle$coef <- fcHurdle$coef / log(2, base = exp(1))

      # Rename columns
      names(fcHurdle)[1:3] <- c("gene", "p", "log2FC")

      fcHurdle$cell_group <- cell_targ
      fcHurdle$up_in <- ifelse(
          fcHurdle$log2FC>0, dge.group.1,
          ifelse(fcHurdle$log2FC<0, dge.group.2, NA)
      )

      # Return
      fcHurdle
    }
  )

  # Collapse tables from all cell types into one, and return!
  do.call(rbind, dge_all_list)
}

