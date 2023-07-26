##### This script initializes a function called 'dgeMAST' for running differential gene expression
# This is an adaptation by Dan Bunis of code originally formulated by Ravi Patel
#
# The function runs DGE between 2 groups, and can loop across multiple clusters, running DGE within each one.
#
# Some key notes:
# - This is a method meant for droplet-based single-cell transcriptome data.
# - Raw counts (integers) are expected to exist and to be:
#   - in the 'counts' slot of a Seurat object's 'RNA' assay.
#   - in an assay named 'counts' if using a SingleCellExperiment object.
#
# To use, `source("path/to/this/file")`
# Then, something like:
# dge <- dgeMAST_allImmStates(
#   object = cd45_data, # Seurat or SCE object
#   cells.group.by = 'full_annots', # metadata name holding clusters or annotations to explore within
#   dge.group.by = 'sex_at_birth', # metadata name holding the samples'/cell's groupings to compare across
#   dge.group.1 = 'M', # value of 'dge.group.by' metadata for cells in the primary group = when foldChange is positive, it means it's up in this population.
#   dge.group.2 = 'F', # value of 'dge.group.by' metadata for cells in the secondary / control group 
#   # Optional, but important to note!
#   cells.group.targets = c('T4 naive', 'T8 em') # Without this, it will run within ALL cell types / clusters (a.k.a. all 'full_annots' per this example)
# )
# See comments in the function definition for details on remaining inputs -- batch.by, cells.use, log.prefix, mast.freq.expressed.min
#
# DGE methodology: To be filled before merge.
#
# Output: a data.frame containing the following columns... TBD by merge!
#
# Authors: Dan Bunis, Ravi Patel

suppressPackageStartupMessages({
  library(assertthat)
  library(dplyr)
  library(Seurat)
  library(SingleCellExperiment)
  library(MAST)
  library(data.table)
  library(dittoSeq)
})

log_message <- function(...) {
  cat("[", format(Sys.time()), "]", ..., "\n")
}

##### This function loops through cell types, performing a single DGE comparison within each one. #####
# You *might* never need to call this fucntion directly.  See it as just an "inner" function that's called by the main function that I want you to call. 
dgeMAST <- function(
  object, # Seurat or SCE object
  cells.group.by, # metadata name holding clusters or annotations to explore within
  cells.group.targets = NULL, # Optional but RECOMMENDED. String ("T4 naive" / "0") or strings ('c("T4 naive", "T8 naive")' / 'c(1,3)'). If not changed from the default here, will loop through all annotations/clusters given by 'cells.group.by'.
  dge.group.by = "ImmState", # metadata column holding samples'/cell's groupings
  dge.group.1, # The "primary" group = when foldChange is positive, it means it's up in this population. 
  dge.group.2, # The "other" group
  batch.by = "orig.ident", # String name of a metadata to treat as an independent variable in the DGE calculation.
  cells.use = NULL, # NULL or Logical vector (length equal to the number of cells in object) indicating cells to keep (TRUE) or remove (FALSE) before calculating DGE.
  # 'Inner' function optimization
  log.prefix = '',
  # MAST method tweaks
  mast.freq.expressed.min = 0.2 # Fraction. A gene filter before dge calculation: Minimal percent of cells (among targets groups, within each targetted cell set) that must express a gene for it to be considered.
) {

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
      
      # Add meta data for num genes captured per cell, and processing batch
      cdr2_ct <- colSums(assay(sca_ct)>0)
      colData(sca_ct)$ngeneson <- scale(cdr2_ct)
      colData(sca_ct)$batch <- factor(meta(batch.by, sce_ct))
      
      # Run MAST
      contrast_val <- paste0(dge.group.by, dge.group.1)
      zlmCond <- zlm(
        as.formula(paste0("~ ngeneson + ", dge.group.by, " + (1 | batch)")),
        sca_ct, exprs_value = 'logcounts',
        method='glmer', ebayes = F,
        strictConvergence = FALSE,
        fitArgsD = list(nAGQ = 0))
      summaryCond <- summary(zlmCond, doLRT=contrast_val)$datatable
      fcHurdle <- merge(
        summaryCond[contrast==contrast_val & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryCond[contrast==contrast_val & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
      fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

      fcHurdle$cell_group <- cell_targ
      fcHurdle$up_in_pos_FC <- dge.group.1

      # Return (Functions in R return the output of their last line, but the `return()` function can also be used!)
      fcHurdle
    }
  )

  # Collapse tables from all cell types into one, and return!
  do.call(rbind, dge_all_list)
}

