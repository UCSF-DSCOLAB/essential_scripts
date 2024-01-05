### This script provides a function for importing annotations from a csv file
# Athough a relatively simple task, the goal is standardization and interfacing with annotation methodolgy of the Data Library
#
# (Documenting in roxygen syntax to be "library"-ready)
#' Import annotations from a csv file into a Seurat or SCE object
#' @param object A Seurat or SingleCellExperiment object, or a String representing the file path of an RDS file containing such an object
#' @param annots_file String representing the file path of the annotations csv you wish to import from.
#' See the details section for file structure requirements.
#' @param annots_regex A String representing a regex that imported annotation columns will be expected to match.
#' Alternatively, set to NULL to import from all columns of the annotation file.
#' @param sep A String denoting the field separator in the \code{annots_file}
#' @param stop_if_overwrite Logical which sets whether the function should warn (FALSE, default) versus error (TRUE) if annotation import would overwrite metadata that already exist in the \code{object}
#' @param verbose Logical which controls whether log messages should be output
#' @return A Seurat or SCE object, with clustering-mapped annotations added to the per-cell metadata
#' @details
#' The primary role of this function is to map annotations from the \code{annots_file} to clustering in the single-cell \code{object}.
#' 
#' It expects \code{annots_file} to target a csv file containing named columns where:\itemize{
#' \item the first column is named after the cell-metadata of \code{object} holding the target clustering.
#' This metadata must exist in the \code{object}, and the values of the column should match to all of the unique values (cluster names) contained in that metadata.
#' \item all other columns either represent annotations of said clusters, or notes (a.k.a. free text).
#' \item when columns containing notes exist, columns containing annotations must be distinguishable via a regex-definable naming scheme.
#' The default expectation is that annotation column names start with "annot", while notes column names will not.
#' }
#' 
#' How it works:\enumerate{
#' \item If \code{object} is given an RDS file path rather than the Seurat or SCE object itself, the object is first loaded in
#' \item The annotation csv is then read in with \code{read.csv(annots_file, header = TRUE, sep = sep)}
#' \item The cell-metadata matching the first column name (clustering) of the annotations csv is extracted from the \code{object}
#' \item An index map is built between values of this clustering metadata and values in this first column
#' \item Columns matching the \code{annots_regex} regex are pulled in to the object making use of 1) the column name as the metadata name, 2) the pre-built index map.
#' }
#' 
#' @author Daniel Bunis
#' @importFrom dittoSeq getMetas
#' @importFrom dittoSeq meta
#' @examples
#' # We'll use the Seurat example dataset for this example
#' sobj <- SeuratObject::pbmc_small
#' 
#' # Our annots_file will be a csv version of:
#' example_annotations <- data.frame(
#'   RNA_snn_res.1 = 0:2,
#'   annots_broad = c("T", "Myeloid", "B"),
#'   notes = c("has CD3E expression", "has MS4A1 expression", "has CD14 expression")
#' )
#' print(example_annotations)
#' 
#' if (!file.exists("annotation_import_example.csv")) {
#'   write.csv(
#'     example_annotations,
#'     file = "annotation_import_example.csv",
#'     row.names = FALSE
#'   )
#' }
#' 
#' # By default, only columns whose name starts with "annot" will be pulled in
#' sobj <- importAnnotations(
#'   object = sobj,
#'   annots_file = "annotation_import_example.csv")
#'   
#' # As we can see, the annotations are now pulled in:
#' table(sobj$RNA_snn_res.1, sobj$annots_broad)
#'
#' # Adjust 'annots_regex' if needed, OR set it to NULL to pull in from all columns
#' #   of the 'annots_file'.
#' sobj <- importAnnotations(
#'   object = sobj,
#'   annots_file = "annotation_import_example.csv",
#'   annots_regex = NULL)
#'
#' # Notice that in this second run, we are warned when the 'annots_broad' metadata
#' #   which was already added in the first run, is being brought in again.
#' # If you want to have the function error out instead of overwriting such
#' #   metadata, set 'stop_if_overwrite' to TRUE.
#' \dontrun{
#' sobj <- importAnnotations(
#'   object = sobj,
#'   annots_file = "annotation_import_example.csv",
#'   stop_if_overwrite = TRUE)
#' }
#' 
importAnnotations <- function(
    object,
    annots_file,
    annots_regex = "^annot",
    sep = ",",
    stop_if_overwrite = FALSE,
    verbose = TRUE) {

    # Until we decide a better way to handle this function I use all the time,
    # simply adding this here, again
    timestamped_message <- function(...) {
        cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
    }
    log <- function(..., do = verbose, timestamp = TRUE) {
      if (do) {
        if (timestamp) {
          timestamped_message(...)
        } else {
          cat(..., "\n", sep = "")
        }
      }
    }
    
    # Optionally load object from file
    if (is.character(object)) {
        log("Loading 'object' from file")
        object <- readRDS(object)
        log("'object' loaded")
    }
    
    # Load annots_file as data.frame
    log("Reading in ", annots_file)
    annots_df <- read.csv(annots_file, header = TRUE, sep = sep)
    
    # Extract clustering metadata from object
    object_metas <- dittoSeq::getMetas(object)
    cluster_meta <- colnames(annots_df)[1]
    if (!cluster_meta %in% object_metas) {
        stop("'object' does not contain a cell-metadata matching the name of the first column of 'annots_file'.")
    }
    clustering <- meta(cluster_meta, object)
    
    # Map cells to rows of annots_df matching their cluster
    log("Mapping 'object' cells' '", cluster_meta, "' identities to cluster names in 'annots_file'")
    ind_map <- match(clustering, annots_df[,1])
    
    # Trim annots_df to only columns we want to pull in
    regex_log_bit <- ""
    if (!identical(annots_regex, NULL)) {
        cols_keep <- c(TRUE, str_detect(colnames(annots_df)[-1], annots_regex))
        annots_df <- annots_df[,cols_keep]
        regex_log_bit <- " matching with the 'annots_regex'"
    }
    if (ncol(annots_df)<2) {
        stop("The 'annots_file' does not contain any annotation columns", regex_log_bit, ".")
    }
    
    # Create annotation metadata
    log("Mapping and adding columns from 'annots_file'", regex_log_bit," to 'object' as cell metadata.")
    for (col in colnames(annots_df)[2:ncol(annots_df)]) {
        log("\tAdding ", col, timestamp=FALSE)
        if (col %in% object_metas) {
            if (stop_if_overwrite) {
                stop("'stop_if_overwrite' set to TRUE and the '", col, "' metadata of 'object' would be overwritten.")
            }
            warning("A metadata named '", col, "' already existed in 'object' and will be overwritten.")
        }
        # At least for the moment, this syntax works for adding cell-metadata for both Seurats and SCEs
        object[[col]] <- annots_df[ind_map,col]
    }
    
    # Output
    log("Annotation adds completed")
    object
}
