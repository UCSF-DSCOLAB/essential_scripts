#' Run differential gene expression (DGE) analysis across cell types
#'
#' Performs differential gene expression analysis on single-cell or
#' pseudobulk count data, iterating independently over each cell-type or
#' cluster present in the input. Supports multiple DGE methods and flexible
#' model specification via fixed and random effects, pairwise group
#' comparisons, or custom contrasts.
#'
#' @param counts A raw gene x sample (or gene x cell) count matrix. Can
#'   contain data for a single cell-type, or be the row-bound
#'   (\code{rbind}) combination of counts across multiple cell-types, in
#'   which case \code{cell_by} must be supplied to indicate how to split
#'   them.
#' @param metadata A data.frame of sample-level (or cell-level) metadata
#'   corresponding to the columns of \code{counts}. Should match the
#'   structure of \code{counts} (single cell-type, or combined across
#'   multiple cell-types).
#' @param dge_by Character string giving the name of the metadata column
#'   containing the sample group labels used to define the DGE comparison
#'   (e.g., case vs. reference).
#' @param method Character string specifying the DGE method to use. One of
#'   \code{"deseq2"}, \code{"dreamlet"}, \code{"edger"}, or \code{"mast"}.
#' @param cell_by Character string giving the name of the metadata column
#'   containing cell-type or cluster annotations, used to iterate the DGE
#'   analysis over each cell-type separately. If \code{NULL}, \code{counts}
#'   and \code{metadata} are assumed to represent a single cell-type or
#'   dataset.
#' @param case_group Character string giving the name of the group (a level
#'   of \code{dge_by}) to use as the case/numerator in log2 fold-change
#'   calculations.
#' @param reference_group Character string giving the name of the group (a
#'   level of \code{dge_by}) to use as the reference/denominator in log2
#'   fold-change calculations.
#' @param contrasts Character vector of explicit model contrasts to test,
#'   e.g. \code{c("varCase - varRef", "varCase2 - varRef")}. Used as an
#'   alternative to specifying \code{case_group}/\code{reference_group} when
#'   more complex or multiple comparisons are needed.
#' @param dge_groups Character vector of group names (levels of
#'   \code{dge_by}) to include when filtering lowly expressed genes prior to
#'   DGE testing.
#' @param fixed_effects Character vector of metadata column name(s) to
#'   include as fixed effects in the DGE model.
#' @param random_effects Character vector of metadata column name(s) to
#'   include as random effects in the DGE model (e.g., for mixed models
#'   such as \code{dreamlet}).
#' @param min_frac Numeric threshold (between 0 and 1) specifying the
#'   minimum fraction of samples that must have CPM > 1 for a gene to be
#'   retained for DGE testing. Default is \code{0.6}.
#' @param return_model Logical indicating whether to return the fitted
#'   model object instead of the DGE results data.frame. Default is
#'   \code{FALSE}.
#'
#' @return If \code{return_model = FALSE}, a data.frame of DGE results (one
#'   row per gene, per cell-type if \code{cell_by} is specified) including
#'   log2 fold change, p-values, and adjusted p-values. If
#'   \code{return_model = TRUE}, the fitted model object(s) instead.
#'
#' @export
#' @examples
#' source('single_cell/dge/fxns_load__SOURCE_ME.R', chdir = TRUE)
#' library(Seurat)
#' expr <- GetAssayData(pbmc_small, layer = 'counts')
#' run_dge(
#'     counts = expr, metadata = pbmc_small@meta.data,
#'     method = 'mast',
#'     dge_by = 'groups',
#'     cell_by = 'RNA_snn_res.0.8',
#'     case_group = 'g1',
#'     reference_group = 'g2'
#' )
#'
run_dge <- function(counts,
                    metadata,
                    dge_by,
                    method = c('deseq2', 'dreamlet', 'edger', 'mast'),
                    cell_by = NULL,
                    case_group = NULL,
                    reference_group = NULL,
                    cell_targets = NULL,
                    contrasts = NULL,
                    dge_groups = c(case_group, reference_group),
                    fixed_effects = NULL,
                    random_effects = NULL,
                    min_frac = 0.6,
                    return_model = FALSE) {

    # Ensures the value of method is among the options given as default, or sets to the first, deseq2, if none was set.
    method <- match.arg(method)

    # Obtain method's function
    dge_fxn <- switch(
        method,
        'deseq2' = run_deseq2,
        'dreamlet' = run_dreamlet,
        'edger' = run_edger,
        'mast' = run_mast
    )

    # Collect input arguments in a list for passing arguments to downstream functions easily
    input_args = list(
        counts = counts,
        metadata = metadata,
        dge_by = dge_by,
        method = method,
        cell_by = cell_by,
        case_group = case_group,
        reference_group = reference_group,
        cell_targets = cell_targets,
        contrasts = contrasts,
        dge_groups = dge_groups,
        fixed_effects = fixed_effects,
        random_effects = random_effects,
        min_frac = min_frac,
        return_model = return_model
    )

    # Perform input checks by invoking ".input_checks()"
    validated_data = do.call( .input_validation, input_args)
    input_args[['counts']] = validated_data[['counts']]
    input_args[['metadata']] = validated_data[['metadata']]
    input_args[['contrasts']] = validated_data[['contrasts']]

    # Iterate over cell-types or clusters
    cell_types = validated_data[['cell_targets']]
    dge_results = list()

    for(ct in cell_types) {
        # Subset to data specific to a cell-type
        input_args[['counts']] = input_args$counts[, input_args$metadata[, input_args$cell_by] == ct ]
        input_args[['metadata']] = input_args$metadata[ input_args$metadata[, input_args$cell_by] == ct, ]

        # Remove lowly expressed genes
        input_args[['counts']] =
            do.call(
                .remove_low_expression,
                input_args[ names(input_args) %in% names(formals(.remove_low_expression)) ]
            )

        # Perform DGE analysis
        dge_results[[ct]] <- do.call(
            dge_fxn,
            input_args[ names(input_args) %in% names(formals(dge_fxn)) ]
        )
    }

    # Return the DGE results
    dge_results
}
