### This script provides a function for pseudobulking single-cell RNAseq data, with a focus towards DGE functions also being planned
#
#' Standardized Pseudobulk Implementations
#' @param object A Seurat object
#' @param sample.by String or string vector naming metadata columns of \code{object} to use for assigning sample identities of cells.
#' A biospecimen, timepoint, or subject name is the common target here, and a batch identifier may be desired as well.
#' @param cell.by String naming the metadata column of \code{object} to use for assigning annotation or cluster identities of cells.
#' @param min.cells Number, 10 by default, which sets the minimum number of cells that a pseudobulk should represent in order to be retained.
#' @param cell.targets Optionally, a string vector naming particular \code{cell.by} values to target. Only these cell identities will be retained and pseudobulked.
#' @param features Optionally, a string vector naming a subset of genes to target. Only counts of these genes will be retained in pseudobulking.
#' @param assay String naming the assay of \code{object} to target. Default is "RNA".
#' @param meta.targets Optionally, a string vector naming particular metadata columns of \code{object} to target for retention. By default, as much a possible will be retained.
#' @param meta.num.summary.method Function like \code{\link[base]{mean}} or \code{\link[base]{median}} (with an \code{na.rm = TRUE} option) describing how to summarize numeric metadata from all cells of each pseudobulk
#' @param method "Seurat", "scuttle", or "dreamlet". Determines which pseudobulking tool to rely on. \itemize{
#' \item "Seurat" = \code{\link[Seurat]{AggregateExpression}}
#' \item Others planned but not yet implemented
#' }
#' @param output.style "Seurat" or "raw". Determines how data should be returned. Either as a Seurat object, or as a raw list of the 'counts' (matrix) and 'metadata' (data.frame)
#' @param output.metadata.cell.count String denoting the metadata column name added to hold the number of original cells going into each pseudobulk
#' @param verbose Logical which controls whether timestamped log messages should be shown
#' @return A Seurat object or, if \code{output.style = 'raw'} a named list with elements 'counts' and 'metadata'.
#' @details
#' The primary role of this function is standardizing the pseudobulking process.
#'
#' \emph{It is under construction and subject to change!}
#'
#' Currently, a Seurat object is expected and pseudobulking is performed per the groupings given by \code{cell.by} and any number of \code{sample.by} metadata.
#' Pseudobulk counts reflect the sum of counts across all cells that they represent.
#'
#' For the \code{method = "Seurat"} (the only method currently implemented) this is done using \code{Seurat::\link[Seurat]{AggregateExpression}} function with parameters:
#' \code{return.seurat = TRUE, group.by = c(sample.by, cell.by), features = features, assays = assay}.
#'
#' The number of cells that each pseudbulk represents are added to a metadata named 'cells_in_pseudobulk' by default.
#' Use \code{output.metadata.cell.count} to use a different column name.
#' \emph{Pseudobulks representing fewer than \code{min.cells} cells are removed.}
#'
#' \code{meta.targets}, or all metadata or \code{object} are then extracted for each pseudobulk.
#' For discrete metadata, the column will be ignored if data are not consistent within ALL pseudobulks.
#' For numeric metadata, values from each cell in the pseudobulks are summarized by the \code{meta.num.summary.method} function (\code{\link[base]{median}} by default).
#'
#' Finally, the output structure is determined by the \code{output.style} input.
#'
#' Prior to pseudobulking, particular cell identities or genes can be targeted using the \code{cell.targets} and \code{features} inputs, respectively.
#'
#' @author Daniel Bunis
#' @examples
#' library(Seurat)
#' dsco_pseudobulk(
#'     object = pbmc_small,
#'     sample.by = "groups",
#'     cell.by = "RNA_snn_res.0.8"
#' )
#'
dsco_pseudobulk <- function(
    object, sample.by, cell.by,
    min.cells = 10,
    cell.targets = NULL, features = NULL,
    meta.targets = NULL,
    assay = "RNA",
    method = c('Seurat', 'scuttle', 'dreamlet'),
    meta.num.summary.method = median,
    output.style = c('Seurat', 'raw'),
    output.metadata.cell.count = 'cells_in_pseudobulk',
    verbose = TRUE
) {
    method <- match.arg(method)
    output.style <- match.arg(output.style)

    orig_metas <- colnames(object@meta.data)

    # Ensure ts_log exists, or source it.
    if (verbose) {
        if (!exists('ts_log')) {
            stop("'ts_log' function required when 'verbose = TRUE'. Set 'verbose = FALSE' or run 'source('<path_to_this_essential_scripts_repo>/R_utils/ts_log.R', chdir = TRUE, local = TRUE)' and then try again.")
        }
        msg_if <- ts_log
    } else {
        msg_if <- function(...) {}
    }

    # Ensure all of cell.by and sample.by exist as metadata
    if (!all(c(cell.by, sample.by) %in% orig_metas)) {
        stop("A 'cell.by' or 'sample.by' element does not exist as a column the 'object' metadata.")
    }
    # Ensure assay exists
    if (!assay %in% Seurat::Assays(object)) {
        stop("'assay' is not a valid assay of the 'object'.")
    }

    # Trim to 'cell.targets'
    if (!is.null(cell.targets)) {
        object <- object[,object@meta.data[,cell.by] %in% cell.targets]
        retained_cell_num <- ncol(object)
        retained_cell_idents <- unique(object@meta.data[,cell.by])
        msg_if(
            "Trimming to 'cell.targets' retained ", retained_cell_num,
            " cells of identities:\n", paste0(retained_cell_idents, collapse=", ")
        )
    }

    ### Pseudobulk
    # (Trim to features internally)
    if (method == 'Seurat') {
        print("Initiating pseudobulking with Seurat's AggregateExpression...")
        group.metas <- c(sample.by, cell.by)
        if (packageVersion('Seurat')>'5.0') {
            psobject <- Seurat::AggregateExpression(
                object = object, return.seurat = TRUE, group.by = group.metas,
                features = features,
                assays = assay)
        } else {
            psobject <- Seurat::AggregateExpression(
                object = object, return.seurat = TRUE, group.by = group.metas,
                features = features,
                assays = assay,
                slot = "counts")
            # Prior versions do not retain group.by metadata
            expected_name <- apply(object@meta.data[,group.metas], 1, FUN = function(x) {paste(x, collapse = "_")})
            ind_match <- match(colnames(psobject), expected_name)
            for (col in group.metas) {
                psobject@meta.data[col] <- object@meta.data[ind_match,col]
            }
        }
        
        find_cells_in_pseudo <- function(object, psobject, col_targs, i) {
            in_pseudo <- rep(TRUE, ncol(object))
            for (col in col_targs) {
                in_pseudo <- in_pseudo & object@meta.data[,col]==as.vector(psobject@meta.data[i,col])
            }
            in_pseudo
        }

        ### Add cell counts and Trim too small pseudobulks
        print("Adding cell counts as metadata")
        msg_if("Adding cell counts as metadata '", output.metadata.cell.count, "'.")
        psobject@meta.data[,output.metadata.cell.count] <- vapply(
            seq_len(ncol(psobject)),
            function(i) {
                sum(find_cells_in_pseudo(object, psobject, group.metas, i))
            },
            numeric(1)
        )
        num_small <- sum(psobject@meta.data[,output.metadata.cell.count] < min.cells)
        if (num_small == ncol(psobject)) {
            warning(paste0("Skipping triming pseudobulks smaller than 'min_cells' as NONE were built from more than ", min_cells, " cells."))
        } else if (num_small > 0) {
            msg_if("\tTrimming ", num_small, " pseudobulks built from fewer than ", min_cells, " cells.")
            psobject <- psobject[,psobject@meta.data[,output.metadata.cell.count] >= min.cells]
        }

        ### Ensure meta.targets, or as many metadata as possible, are retained.
        meta_ignored <- c()
        meta_kept <- colnames(psobject@meta.data)
        if (is.null(meta.targets)) {
            meta.targets <- orig_metas
        }
        meta.targets <- meta.targets[!meta.targets %in% meta_kept]
        msg_if("Grabbing additional metadata '", paste0(meta.targets, collapse = "', '"), "'.")
        meta_collapse <- function(df_col, name) {
            if (is.numeric(df_col)) {
                meta.num.summary.method(df_col, na.rm = TRUE)
            } else {
                if (name %in% meta_ignored) {
                    NA
                } else {
                    if (length(unique(df_col))>1) {
                        msg_if("\tignoring discrete metadata column '", name, "' because it is not consistent within all pseudobulks")
                        assign('meta_ignored', c(meta_ignored, name), inherits = TRUE)
                        NA
                    } else {
                        as.vector(df_col[1])
                    }
                }
            }
        }
        for (i in seq_len(ncol(psobject))) {
            in_pseudo <- find_cells_in_pseudo(object, psobject, group.metas, i)
            for (meta in meta.targets) {
                if (i == 1) {
                    psobject@meta.data[,meta] <- NA
                }
                psobject@meta.data[i,meta] <- meta_collapse(object@meta.data[in_pseudo,meta], meta)
            }
        }
        psobject@meta.data <- psobject@meta.data[,!colnames(psobject@meta.data)%in% meta_ignored]
    } else {
        stop('Alternative Pseudobulking methods not yet implemented.')
    }

    # Output in requested style
    msg_if("Pseudobulk function COMPLETE.")
    if (output.style == 'Seurat') {
        psobject
    } else {
        exp <- if (packageVersion("Seurat")>='5.0') {
            Seurat::GetAssayData(psobject, layer = 'counts')
        } else {
            Seurat::GetAssayData(psobject, slot = 'counts')
        }
        list(counts = exp, metadata = psobject@meta.data)
    }
}
