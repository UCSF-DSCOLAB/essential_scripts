##### This script initializes a function called 'process_normalized_citeseq_data' which aims to achieve full "sub-clustering" of a Seurat object containing CITEseq data
# To use:
# - 'source("path/to/this/file")'
# - Optional but often necessary, details below: 'options(future.globals.maxSize = 8000 * 1024^2)'
# - Subset / prepare the data you plan to process.
# - Adjust inputs as needed.
# - Invoke the 'process_normalized_citeseq_data()' function!

### Primary steps of the function:
# 1. Takes in a Seurat object which should be the target subset of your original / full object, then re-runs pre-processing steps from right after the data normalization stage.
#    **Object / Processing Plan assumptions**:
#      - Gene expression data in an "RNA" assay
#      - Protein expression data in an "ADT" assay
#      - batch correction methodologies are harmony for RNA and either harmony or Seurat's rpca-based integration approach for ADT.
#      - Aside from "S.Score" and "G2M.Score" which will be (re-)created internally via CellCycleScoring on the RNA assay, any metadata intended to be regressed out in scaling steps and given to 'rna_vars_to_regress' and/or 'adt_vars_to_regress' should already exist in the object.
#      - When using Seurat's rpca-based integration approach for ADT batch correction, a metadata holding library identities exists and is a factor with levels ordered by batch. (See 'integration_references' input details for why.)
# 2. Runs CellCycleScoring, highly variable gene selection, scaling, pca, and harmony batch correction for the RNA assay
# 3a. If using Seurat's rpca-based integration for ADT batch correction: Runs scaling on PCA in a per-library fashion and uses Seurat's 'IntegrateData' methodology for batch corretion of ADT data, pulling a PCA run on the integrated data as the batch corrected dimensionality reduction assay here.
# 3b. If using harmony for ADT batch correction: Runs scaling, pca, and harmony batch correction on all non-isotype markers of the ADT assay
# 4. Runs WNN and umap and clustering based on the WNN.
# 5. Optionally also runs umap and clustering based on the RNA assay only
# 6. Returns this re-processed object which has been, in part, re-clustered within its chosen subset of cells.

### R Package Requirements:
# - Seurat v4 or later, originally tested in 4.1.1
# - harmony
# - dittoSeq (due to a single getGenes call)

### Inputs:
# Documented within the function definition below!

### Outputs:
# Two output methods can be used:
#  - The re-processed Seurat object can be saved to a .Rds file by giving a filepath to the "rds_path" input.
#  - The re-processed Seurat object can also be output directly as the result of the function call when "output" is set to TRUE, which is the default.

### Usage notes:
# You may need to run 'options(future.globals.maxSize = 8000 * 1024^2)' prior to calling the function in order to avoid memory usage complaint/error within the ADT reprocessing loop. Per https://github.com/satijalab/seurat/issues/1845

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dittoSeq)
})

process_normalized_citeseq_data <- function(
  ### Primary Inputs
  object, # The Seurat object to target, which should already be subset to your cells of interest
  log_prefix = "", # String which will be added to the beginning of any timestamped log lines. E.g.: "T cells: ". Useful when multiple subsets will be processed with the same log file.
  rds_path = NULL, # String or NULL. File path naming where to output the processed Seurat as an Rds file.  If left NULL, the object will not be written to a file.
  output = TRUE, # Boolean. Whether to return the Seurat object, or not, at the end.
  ### Secondary Inputs
  rna_vars_to_regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), # Value is passed to Seurat::ScaleData()'s 'vars.to.regress' input in RNA assay re-scaling.
  adt_vars_to_regress = c("nCount_ADT", "S.Score", "G2M.Score"), # Value is passed to Seurat::ScaleData()'s 'vars.to.regress' input in ADT assay re-scaling after integration.
  harmony_group_by = "orig.ident", # String (or String vector) naming your batch metadata (and any other metadata) intended to be given to the harmony::RunHarmony()'s 'group.by.vars' input. 
  adt_batch_correction_method = c("seurat-integration", "harmony"), # One of the two values to the left.  Sets the methodology used for ADT batch correction.
  integration_split_by = "orig.ident", # String naming a metadata which holds original sequencing library ids. Only used when 'adt_batch_correction_method = "seurat-integration"'. Used for splitting the ADT data into 1 object per library for the integration batch correction approach.
  # NOTE: This integration_split_by metadata should be a factor so that you can know that reference indices indicated with 'integration_references' are the first of each batch!
  integration_references = NULL, # NULL to ignore, but highly recommended to give this the indexes of your first library per sequencing / preparation batch among the integration_split_by-metadata's levels
  rna_pca_dims = 1:30, # Number vector. The RNA PCs to utilize in WNN (and optional RNAonly) clustering and UMAP calculations.
  adt_pca_dims = 1:20, # Number vector. The ADT PCs to utilize in WNN, and subsequent clustering and UMAP calculations.
  clustering_resolutions = seq(0.1, 2.0, 0.1), # Number vector used for the resolution parameter of Seurat::FindClusters() calls.
  rna_only_umap_and_clustering = TRUE, # Boolean, whether or not to calculate a umap and clusterings based only on the RNA assay.
  warn_if_no_rds = TRUE # Logical. Controls whether or not a message-level warning will be given from the start if the re-processed object will only be returned as output, and not first saved to a file.
  ) {

  print_message <- function(...) {
    cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
  }

  adt_batch_correction_method <- match.arg(adt_batch_correction_method)

  if (is.null(rds_path)) {
    warn_start <- ifelse(log_prefix=="", "", paste0("For prefix '", log_prefix, ": "))
    if (!output) {
      stop(warn_start, "No 'rds_path' given && 'output' set to FALSE. Nothing would be returned!")
    }
    if (warn_when_no_rds) {
      message(warn_start, "No 'rds_path' given. Object will be returned as output but not saved to disk within this function. (Set 'warn_if_no_rds = FALSE' to turn off this message.)")
    }
  }
  
  print_message(log_prefix, "Starting processing")

  print_message(log_prefix, "RNA: Running Cell Cycle Scoring")
  object <- CellCycleScoring(object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, nbin = 12)
  print_message(log_prefix, "RNA: Selecting HVGs, regressing metrics and scaling, and PCA")
  object <- FindVariableFeatures(object, verbose = FALSE)
  object <- ScaleData(object, vars.to.regress = rna_vars_to_regress, verbose = FALSE)
  object <- RunPCA(object, reduction.name = "pca", verbose = FALSE)
  print_message(log_prefix, "RNA: Batch Correcting PCA with Harmony")
  object <- RunHarmony(object, group.by.vars = harmony_group_by, reduction.save = "harmony", verbose = FALSE)

  if (adt_batch_correction_method == "harmony") {
    print_message(log_prefix, "ADT: Selecting non-isotypes as \"HVGs\", regressing metrics and scaling, and PCA")
    VariableFeatures(object, assay = "ADT") <- grep("isotype", getGenes(object, assay = "ADT"), ignore.case = TRUE, value = TRUE, invert = TRUE)
    object <- ScaleData(object, vars.to.regress = adt_vars_to_regress, verbose = FALSE, assay = "ADT")
    object <- RunPCA(object, reduction.name = "adt.pca", verbose = FALSE, assay = "ADT")
    print_message(log_prefix, "ADT: Batch Correcting PCA with Harmony")
    object <- RunHarmony(object, group.by.vars = harmony_group_by, reduction.save = "adt.harmony", verbose = FALSE)
    adt_reduction <- "adt.harmony"
  } else {
    print_message(log_prefix, "ADT: Selecting non-isotypes, then subsetting to per-library objects")
    features <- grep("isotype", getGenes(object, assay = "ADT"), ignore.case = TRUE, value = TRUE, invert = TRUE)
    DefaultAssay(object) <- "ADT"
    diet_object <- DietSeurat(object, counts = TRUE, data = TRUE, scale.data = FALSE, assays = "ADT", dimreducs = NULL, graphs = NULL, misc = FALSE)
    DefaultAssay(object) <- "RNA"
    libs <- levels(as.factor(diet_object[[integration_split_by, drop = TRUE]]))
    sobjs.list <- lapply(libs, function(lib) {
      diet_object[, diet_object[[integration_split_by, drop = TRUE]]==lib]
    })
    print_message(log_prefix, "ADT: Scaling and PCA per library\n\tWorking on:")
    sobjs.list <- lapply(seq_along(libs), function(x) {
      cat("\t", libs[x], "\n", sep = "")
      x <- sobjs.list[[x]]
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    names(sobjs.list) <- libs
    print_message(log_prefix, "ADT: Finding integration anchors")
    immune.anchors <- FindIntegrationAnchors(
      object.list = sobjs.list,
      anchor.features = features,
      scale = FALSE,
      reference = integration_references,
      reduction = "rpca",
      verbose = FALSE)
    print_message(log_prefix, "ADT: Integrating")
    adt_int <- IntegrateData(anchorset = immune.anchors, new.assay.name = "integrated.ADT", verbose = FALSE)
    DefaultAssay(adt_int) <- "integrated.ADT"
    print_message(log_prefix, "ADT: Scaling and PCA of Integrated ADT data")
    adt_int <- ScaleData(adt_int, features = features, vars.to.regress = adt_vars_to_regress, verbose = FALSE)
    adt_int <- RunPCA(adt_int, reduction.name = "int.adt.pca", reduction.key = "integratedADTpca_", verbose = FALSE)
    print_message(log_prefix, "ADT: Pulling integrated ADT assay and PCA to in-progress object")
    stopifnot(all(colnames(adt_int) %in% colnames(object)))
    stopifnot(all(colnames(object) %in% colnames(adt_int)))
    adt_int <- adt_int[, colnames(object)]
    object[["integrated.ADT"]] <- adt_int[["integrated.ADT"]]
    object@reductions$int.adt.pca <- adt_int@reductions$int.adt.pca
    adt_reduction <- "int.adt.pca"
  }

    
  
  print_message(log_prefix, "WNN: Weighted Nearest Neighbors")
  object <- FindMultiModalNeighbors(
    object, reduction.list = list("harmony", adt_reduction), 
    dims.list = list(rna_pca_dims, adt_pca_dims), modality.weight.name = c("RNA.weight", "ADT.weight"),
    knn.graph.name = "wknn",
    snn.graph.name = "wsnn",
    weighted.nn.name = "wnn",
    prune.SNN = 1/20,
    verbose = FALSE)
  print_message(log_prefix, "WNN: UMAP")
  object <- RunUMAP(object, nn.name = "wnn", reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)
  print_message(log_prefix, "WNN: Clustering")
  object <- FindClusters(object, graph.name = "wsnn", algorithm = 2, resolution = clustering_resolutions, verbose = FALSE)
  
  if (rna_only_umap_and_clustering) {
    print_message(log_prefix, "RNAonly: Nearest Neighbors")
    object <- FindNeighbors(object, reduction = "harmony", dims = rna_pca_dims, verbose = FALSE)
    print_message(log_prefix, "RNAonly: UMAP")
    object <- RunUMAP(object, reduction = "harmony", dims = rna_pca_dims, reduction.name = "umap_rna", reduction.key = "UMAPrna_", verbose = FALSE)
    print_message(log_prefix, "RNAonly: Clustering")
    object <- FindClusters(object, algorithm = 2, resolution = clustering_resolutions, graph.name = "RNA_snn", verbose = FALSE)
  }

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
