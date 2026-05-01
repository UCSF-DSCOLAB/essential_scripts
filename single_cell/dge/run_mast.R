## TODO: Remove/reduce these requirements!
suppressPackageStartupMessages({
  library(dplyr)
  library(MAST)
  library(data.table)
  library(dittoSeq)
})
source('helpers_load.R')

#' Run DGE between 2 groups using the MAST method
#' @param counts A feature (row) x sample (column) raw count matrix.  Values should be integers.  Only samples intended for DGE analysis should be included.
#' @param metadata data.frame containing cell metadata, with rownames containing column names of \code{counts}.
#' @param dge_by String naming a column of \code{metadata} which holds the sample or cell groupings that you wish to compare between
#' @param case_group,reference_group Strings naming the groups of \code{dge_by}-data to compare between. Directionality: Positive log2FC (foldChange) will mean upregulation in \code{case_group} with respect to \code{reference_group}.
#' @param log.prefix String to append at the beginning of log messages. Useful when wrapping this function inside an external loop (such us the '_within_cells' version of this function).
#' @param mast.freq.expressed.min A fraction between 0 and 1, default = 0.2, which sets the minimal percent of cells that must express a gene for it to be considered.
#' MAST performs less well for genes with low percent expression.
#' The calculation run only among cells/samples of the targeted \code{case_group} and \code{reference_group} groups.
#' @param random_effects,fixed_effects # NULL or String vectors naming metadata columns of \code{object} to treat as a random or fixed effects, respectively, in the mixed effect model used for DGE calculation. Additional details:
#' \itemize{
#' \item NOTE: \code{random_effects} defaults to "orig.ident" with the assumption that this column generally 1) exists and 2) holds batch information that is normally desired to be modeled as a random effect. Set to NULL, or something different, to turn this off.
#' \item \code{random_effects} must target discrete metadata and will be turned into a factor before use.
#' \item 'ngeneson' (the number of expressed genes per cell, which will be calculated internally) and \code{dge_by} are automatically added as fixed effects and do not need to be added here.
#' }
#' @param contrast Not yet implemented
#' @param min_frac Not yet implemented
#' @param return_model Not yet impemented
#' @return a data.table
#' @details Need to re-write
#' @section DGE methodology details:
#' DGE is run by first modeling the data using MAST's zlm function, followed by determining p-values using its likelihood ratio test methodology on the contrast of interest.
#'
#' Specifically, \code{zlm} is run with settings: \itemize{
#' \item Always: \code{ebayes = FALSE, strictConvergence = FALSE}
#' \item Additionally:\itemize{
#'   \item when there is a random effect given (default) to \code{random_effects}: \code{method='glmer', fitArgsD = list(nAGQ = 0)}
#'   \item when there is no random effect given to \code{random_effects}: \code{method='glm'}
#' }
#' }
#' @author Daniel Bunis and Ravi Patel
run_mast <- function(
    counts,
    metadata,
    dge_by,
    case_group,
    reference_group,
    contrast = NULL,
    min_frac, # IMPLEMENT
    log.prefix = '',
    mast.freq.expressed.min = 0.2,
    random_effects = "orig.ident",
    fixed_effects = NULL,
    return_model = FALSE
) {

    ### ToDo: Checks will soon come in from elsewhere.  Use those!
    ### TODO: Implement contrasts

    # Trim to case/ref groups only
    in_grps <- metadata[,dge_by] %in% c(case_group, reference_group)
    counts <- counts[,in_grps]
    metadata <- metadata[in_grps,]
    metadata[,dge_by] <- factor(
        metadata[,dge_by],
        levels = c(reference_group, case_group)
    )

    # logcounts normalization
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts <- log2(t(t(counts)/size.factors) + 1)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            logcounts = logcounts
        ),
        colData = metadata
    )
    # if (!dittoSeq::isMeta(dge_by, object)) stop("'dge_by', ", dge_by,", is not a metadata of 'object'")
    # if (!case_group %in% object[[dge_by, drop = TRUE]]) stop("No cells have 'case_group', ", case_group, ", as their 'dge_by' value.")
    # if (!reference_group %in% object[[dge_by, drop = TRUE]]) stop("No cells have 'reference_group', ", reference_group, ", as their 'dge_by' value.")

    # Transform into MAST's required SingleCellAssay structure.
    sca <- SceToSingleCellAssay(sce)

    # Trim to only genes expressed in enough cells
    expressed_genes <- freq(sca) > mast.freq.expressed.min
    sca <- sca[expressed_genes,]

    # Add meta data for num genes captured per cell
    cdr2_ct <- colSums(assay(sca)>0)
    colData(sca)$ngeneson <- scale(cdr2_ct)

    # Build model
    model <- paste0("~ ngeneson + ", dge_by)
    # Ensure random effects, like batch, are factors too
    for (rand in random_effects) {
        colData(sca)[,rand] <- as.factor(colData(sca)[,rand])
        model <- paste0(model, " + (1 | ", rand, ")")
    }
    for (fix in fixed_effects) {
        model <- paste0(model, " + " fix)
    }

    # Run MAST
    zlm_args <- list(
        formula = as.formula(model),
        sca = sca_ct,
        exprs_values = 'logcounts',
        ebayes = FALSE,
        strictConvergence = FALSE,
        method = 'glm'
    )
    if (!identical(random_effects, NULL)) {
        zlm_args$method = 'glmer'
        zlm_args$fitArgsD = list(nAGQ = 0)
    }

    ts_log(log.prefix, "\twith model formula: ", model)
    contrast_val <- paste0(dge_by, case_group)
    zlmCond <- do.call(zlm, zlm_args)
    summaryCond <- summary(zlmCond, doLRT=contrast_val)$datatable
    fcHurdle <- merge(
        summaryCond[contrast==contrast_val & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryCond[contrast==contrast_val & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

    # Rebase coef toward log2() instead of ln()
    fcHurdle$coef <- fcHurdle$coef / log(2, base = exp(1))

    # Rename columns
    names(fcHurdle)[1:3] <- c("gene", "pval", "log2FC")

    # Return
    fcHurdle
}

#' Run DGE between 2 groups using the MAST method, looping through cell types
#' @param counts Integer counts matrix
#' @param metadata 
#' @param cell_by String naming a column of \code{metadata} which holds clusters, annotations, or other cell groupings to explore within
#' @param cell_targets NULL (to default to all options) or a string vector naming particular levels of the \code{cell_by}-data to run DGE within.
#' @param dge_by String naming a metadata column of \code{object} which holds the sample or cell groupings that you wish to compare between
#' @param case_group,reference_group Strings naming the groups of \code{dge_by}-data to compare between. Directionality: Positive log2FC (foldChange) will mean upregulation in \code{case_group} with respect to \code{reference_group}.
#' @param log.prefix String to append at the beginning of log messages. Useful when wrapping this function inside an external loop.
#' @param mast.freq.expressed.min A fraction between 0 and 1, default = 0.2, which sets the minimal percent of cells that must express a gene for it to be considered. MAST performs less well for genes with low percent expression. This cutoff is run per each considered \code{cell.group.targets} cell grouping, and only expression among cells of the targeted \code{case_group} and \code{reference_group} groups are considered.
#' @param random_effects,fixed_effects # NULL or String vectors naming metadata columns of \code{object} to treat as a random or fixed effects, respectively, in the mixed effect model used for DGE calculation. Additional details:
#' \itemize{
#' \item NOTE: \code{random_effects} defaults to "orig.ident" with the assumption that this column generally 1) exists and 2) holds batch information that is normally desired to be modeled as a random effect. Set to NULL, or something different, to turn this off.
#' \item \code{random_effects} must target discrete metadata and will be turned into a factor before use.
#' \item 'ngeneson' (the number of expressed genes per cell, which will be calculated internally) and \code{dge_by} are automatically added as fixed effects and do not need to be added here.
#' }
#' @param min_per_cell Number giving the minimum number of samples or cells to require exist in all groups for each celltype.  Otherwise, the celltype will be skipped.
#' @param contrast Not yet implemented
#' @param min_frac Not yet implemented
#' @param return_model Not yet impemented
#' @return a data.table
#' @details Need to re-write
#' @section DGE methodology details:
#' DGE is run by first modeling the data using MAST's zlm function, followed by determining p-values using its likelihood ratio test methodology on the contrast of interest.
#'
#' Specifically, \code{zlm} is run with settings: \itemize{
#' \item Always: \code{ebayes = FALSE, strictConvergence = FALSE}
#' \item Additionally:\itemize{
#'   \item when there is a random effect given (default) to \code{random_effects}: \code{method='glmer', fitArgsD = list(nAGQ = 0)}
#'   \item when there is no random effect given to \code{random_effects}: \code{method='glm'}
#' }
#' }
#' @author Daniel Bunis and Ravi Patel
run_mast_within_cells <- function(
    counts,
    metadata,
    cell_by,
    cell_targets = NULL,
    dge_by,
    case_group,
    reference_group,
    contrast = NULL,
    min_frac, # IMPLEMENT
    min_per_group = 10, # IMPLEMENT elsewhere too?
    log.prefix = '',
    mast.freq.expressed.min = 0.2,
    random_effects = "orig.ident",
    fixed_effects = NULL,
    return_model = FALSE
) {

    cell_targets <- .input_checks_within_cells(
        metadata,
        cell_by, cell_targets,
        dge_by, case_group, reference_group, contrast,
        min_per_group
    )

    # Loop through cell types, building dge output for each
    dge_all_list <- lapply(
        cell_targets,
        function(cell_targ) {

            ts_log(log.prefix, "Running MAST dge for ", cell_targ, " cells")
            out <- run_mast(
                counts = counts[,cts==cell_targ],
                metadata = metadata[cts==cell_targ,],
                dge_by = dge_by,
                case_group = case_group,
                reference_group = reference_group,
                contrast = contrast,
                min_frac = min_frac,
                log.prefix = log.prefix,
                mast.freq.expressed.min = mast.freq.expressed.min,
                random_effects = random_effects,
                fixed_effects = fixed_effects,
                return_model = return_model
            )
            if (!return_model) {
                out$celltype <- cell_targ
            }
            out
        }
    )

    if (return_model) {
        names(dge_all_list) <- cell_targets
        dge_all_list
    } else {
        # Collapse tables from all cell types into one, and return!
        do.call(rbind, dge_all_list)
    }
}

