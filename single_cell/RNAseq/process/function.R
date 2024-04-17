##### This script initializes a function called 'process_normalized_rnaseq_data' which aims to achieve full "pre-processing" / "clustering" / "sub-clustering" of a Seurat object containing RNAseq data
# To use:
# - 'source("path/to/this/file")'
# - Subset / merge / prepare the data you plan to process.
# - Adjust inputs as needed.
# - Invoke the 'process_normalized_rnaseq_data()' function!

### Primary steps of the function:
# 1. Takes in a Seurat object which should contain only the set of cells you wish to cluster or "sub"cluster, then runs pre-processing steps from right after the data normalization stage.
#    **Object / Processing Plan assumptions**:
#      - Gene expression data in "RNA" assay. (Additional assays can exist, but will be lost unless given to 'assays_keep', or re-added after the function.)
#      - Desired batch correction methodology.
#      - Aside from "S.Score" and "G2M.Score" which will be (re-)created internally via CellCycleScoring on the RNA assay, any metadata intended to be regressed out in scaling steps and given to 'rna_vars_to_regress' should already exist in the object.
# 2. Runs DietSeurat to clean old dimensionality reductions and possible SCT assays
# 3. Runs CellCycleScoring, highly variable gene selection, scaling, pca
# 4. Uses RunHarmony to batch correct.
# 7. Runs umap and clustering based on the RNA assay only
# 8. Returns this newly pre-processed object.

### R Package Requirements:
# - Seurat v4 or later, originally tested in 4.1.1
# - harmony

### Inputs:
# Documented within the function definition below!

### Outputs:
# Two output methods can be used:
#  - The re-processed Seurat object can be saved to a .Rds file by giving a filepath to the "rds_path" input.
#  - The re-processed Seurat object can also be output directly as the result of the function call when "output" is set to TRUE, which is the default.
# Components of the output Seurat object:
#  - All elements of the input object except for metadata & RNA and ADT assay counts and normalized counts matrices are cleared at the start. Everything else is re-created in the function:
#  - Assays:
#    - RNA
#    - any additional assays given to 'assays_keep', default = ADT, *with only 'counts' and 'data' slots retained.*
#  - Dimensionality reduction names:
#    - pca: PCA, pre-batch correction
#    - harmony: harmonized PCA
#    - umap (or alternative name given to 'umap_name'): umap built from harmonized PCA
#    - umap_no_harmony(or '{umap_name}_no_harmony'): umap built pre-batch correction PCA
#  - Clustering Prefixes:
#    - RNA_snn_: RNA-based clustering, built from nearest neighbors determined on harmonized PCA
#    - Note that clustering is stored as metadata, and the DietSeurat cleanup step does not affect metadata. Thus, *previous clustering metadata will remain unless overwritten via the prefix above.*

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

process_normalized_rnaseq_data <- function(
  ### Primary Inputs
  object, # The Seurat object to target, which should already be subset to your cells of interest
  log_prefix = "", # String which will be added to the beginning of any timestamped log lines. E.g.: "T cells: ". Useful when multiple subsets will be processed with the same log file.
  rds_path = NULL, # String or NULL. File path naming where to output the processed Seurat as an Rds file.  If left NULL, the object will not be written to a file.
  output = TRUE, # Boolean. Whether to return the Seurat object, or not, at the end.
  ### Batch Correction Control
  harmony_group_by = "orig.ident", # String (or String vector) naming your batch metadata (and any other metadata) intended to be given to the harmony::RunHarmony()'s 'group.by.vars' input.
  ### Secondary Inputs
  vars_to_regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), # Value is passed to Seurat::ScaleData()'s 'vars.to.regress' input in RNA assay re-scaling.
  pca_dims = 1:30, # Number vector. The RNA PCs to utilize in WNN (and optional RNAonly) clustering and UMAP calculations.
  clustering_resolutions = seq(0.1, 2.0, 0.1), # Number vector used for the resolution parameter of Seurat::FindClusters() calls.
  assays_keep = c("RNA", "ADT"), # String vector giving names of all assays whose raw and normalized counts you wish to retain.  Can be a larger set than exists, but if you have data in another assay that you will want to access later, it must be named here or re-added later.
  umap_name = "umap", # String giving the name to use for the umap dimensionalty reduction.
  warn_if_no_rds = TRUE # Logical. Controls whether or not a message-level warning will be given from the start if the re-processed object will only be returned as output, and not first saved to a file.
  ) {

  print_message <- function(...) {
    cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
  }

  warn_start <- ifelse(log_prefix=="", "", paste0("For prefix '", log_prefix, ": "))
  if (is.null(rds_path)) {
    if (!output) {
      stop(warn_start, "No 'rds_path' given && 'output' set to FALSE. Nothing would be returned!")
    }
    if (warn_if_no_rds) {
      message(warn_start, "No 'rds_path' given. Object will be returned as output but not saved to disk within this function. (Set 'warn_if_no_rds = FALSE' to turn off this message.)")
    }
  }
  
  print_message(log_prefix, "Starting processing")
  DefaultAssay(object) <- "RNA"
  assays_keep <- assays_keep[assays_keep %in% Assays(object)]
  object <- DietSeurat(object, counts = TRUE, data = TRUE, scale.data = FALSE, assays = assays_keep, dimreducs = NULL, graphs = NULL, misc = FALSE)

  print_message(log_prefix, "Running Cell Cycle Scoring")
  object <- CellCycleScoring(object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, nbin = 12)
  print_message(log_prefix, "Selecting HVGs, regressing metrics and scaling, and PCA")
  object <- FindVariableFeatures(object, verbose = FALSE)
  object <- ScaleData(object, vars.to.regress = vars_to_regress, verbose = FALSE)
  object <- RunPCA(object, reduction.name = "pca", verbose = FALSE)
  print_message(log_prefix, "Batch Correcting PCA with Harmony")
  object <- harmony::RunHarmony(object, group.by.vars = harmony_group_by, reduction = "pca", reduction.save = "harmony", verbose = FALSE)
  
  print_message(log_prefix, "Nearest Neighbors")
  object <- FindNeighbors(object, reduction = "harmony", dims = pca_dims, verbose = FALSE)
  print_message(log_prefix, "UMAP")
  object <- RunUMAP(object, reduction = "harmony", dims = pca_dims, reduction.name = umap_name, reduction.key = "UMAP_", verbose = FALSE)
  object <- RunUMAP(object, reduction = "pca", dims = pca_dims, reduction.name = paste0(umap_name, "_no_harmony"), reduction.key = "UMAPnoHarmony_", verbose = FALSE)
  print_message(log_prefix, "Clustering")
  object <- FindClusters(object, algorithm = 2, resolution = clustering_resolutions, graph.name = "RNA_snn", verbose = FALSE)

  print_message(log_prefix, "Object Processing Complete!")
  
  if (!is.null(rds_path)) {
    print_message(log_prefix, "Saving to '", rds_path, "' ...")
    saveRDS(object, file = rds_path)
    print_message(log_prefix, "Saved")
  }

  if (output) {
    print_message(log_prefix, "Processed Object:")
    object
  }
}
