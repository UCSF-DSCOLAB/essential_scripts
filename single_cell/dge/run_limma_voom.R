#' Run differential gene expression analysis with limma-voom
#'
#' Fits a limma-voom model to pseudobulk RNA-seq counts. Input validation,
#' sample filtering, cell-type iteration, and low-expression filtering are
#' expected to be performed by \code{run_dge()} before this function is called.
#'
#' @param counts A gene-by-sample matrix or data.frame of raw pseudobulk counts.
#' @param metadata A data.frame with one row per sample, in the same order as
#'   the columns of \code{counts}.
#' @param dge_by Character scalar naming the metadata column that defines the
#'   groups being compared. It is always included as a fixed effect.
#' @param case_group Character scalar naming the numerator group for a pairwise
#'   comparison.
#' @param reference_group Character scalar naming the denominator group for a
#'   pairwise comparison.
#' @param contrasts A named list or character vector of limma contrast
#'   expressions. This is used only when \code{case_group} and
#'   \code{reference_group} are not both supplied. Design coefficient names are
#'   available in the returned model's \code{design} component.
#' @param fixed_effects Optional character vector naming additional metadata
#'   columns to include as fixed effects.
#' @param random_effects Optional character scalar naming a repeated-measures
#'   blocking variable. limma supports one blocking variable; supplying more
#'   than one is an error.
#' @param return_model Logical. If \code{TRUE}, return the fitted, moderated
#'   \code{MArrayLM} model instead of a results data.frame.
#'
#' @return If \code{return_model = FALSE}, a data.frame with columns
#'   \code{gene}, \code{log2FC}, \code{avgExpr}, \code{pval}, and \code{padj}.
#'   Results from user-supplied contrasts also include a \code{contrast} column.
#'   If \code{return_model = TRUE}, a moderated \code{MArrayLM} object is
#'   returned.
#'
#' @details Counts are TMM-normalized before the voom transformation. When
#'   \code{random_effects} is supplied, the within-block correlation is
#'   estimated with \code{duplicateCorrelation}; voom and the linear model are
#'   then refit using that consensus correlation. P-values are adjusted with
#'   the Benjamini-Hochberg method separately for each contrast.
#'
#' @export
run_limma_voom <- function(
    counts,
    metadata,
    dge_by,
    case_group = NULL,
    reference_group = NULL,
    contrasts = NULL,
    fixed_effects = NULL,
    random_effects = NULL,
    return_model = FALSE
) {
    # Check and validate dependencies
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("The 'limma' package is required to run run_limma_voom().", call. = FALSE)
    }
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("The 'edgeR' package is required to run run_limma_voom().", call. = FALSE)
    }

    if (length(random_effects) > 1L) {
        stop(
            "limma-voom accepts only one random effect as a blocking variable.",
            call. = FALSE
        )
    }
    # Setup the design matrix (fixed effects only; random effects added later)
    metadata <- as.data.frame(metadata)
    metadata[[dge_by]] <- factor(metadata[[dge_by]])

    model_variables <- unique(c(dge_by, fixed_effects))
    design_formula <- stats::reformulate(model_variables, intercept = FALSE)
    design <- stats::model.matrix(design_formula, data = metadata)
    colnames(design) <- make.names(colnames(design), unique = TRUE)
    # Check whether every coefficient in the design matrix can be niquely estimated
    non_estimable <- limma::nonEstimable(design)
    if (!is.null(non_estimable)) {
        stop(
            "The design matrix is not full rank. Non-estimable coefficient(s): ",
            paste(non_estimable, collapse = ", "),
            call. = FALSE
        )
    }
    # Automatically construct comparison framework
    pairwise <- !is.null(case_group) && !is.null(reference_group)
    if (pairwise) {
        group_columns <- which(attr(design, "assign") == 1L)
        group_coefficients <- stats::setNames(
            colnames(design)[group_columns],
            levels(metadata[[dge_by]])
        )

        case_coefficient <- unname(group_coefficients[as.character(case_group)])
        reference_coefficient <- unname(
            group_coefficients[as.character(reference_group)]
        )
        if (is.na(case_coefficient) || is.na(reference_coefficient)) {
            stop(
                "Case and reference groups must be levels of metadata[['",
                dge_by,
                "']].",
                call. = FALSE
            )
        }

        contrast_matrix <- matrix(
            0,
            nrow = ncol(design),
            ncol = 1L,
            dimnames = list(
                colnames(design),
                paste0(case_group, "_vs_", reference_group)
            )
        )
        contrast_matrix[case_coefficient, 1L] <- 1
        contrast_matrix[reference_coefficient, 1L] <- -1
        include_contrast_column <- FALSE
    } else {
        contrast_expressions <- .limma_contrast_expressions(contrasts)
        contrast_matrix <- limma::makeContrasts(
            contrasts = unname(contrast_expressions),
            levels = design
        )
        colnames(contrast_matrix) <- names(contrast_expressions)
        include_contrast_column <- TRUE
    }
    # store count data in a DGEList object
    dge <- edgeR::DGEList(counts = counts)
    # Compute normalization factors
    dge <- edgeR::calcNormFactors(dge)

    block <- NULL
    consensus_correlation <- NULL
    # Model random effects using (consensus) duplicated correlation
    if (length(random_effects) == 1L) {
        block <- factor(metadata[[random_effects]])
        if (anyNA(block)) {
            stop("The limma blocking variable cannot contain missing values.", call. = FALSE)
        }
        if (!anyDuplicated(block)) {
            stop(
                "The limma blocking variable has no repeated samples, so a ",
                "within-block correlation cannot be estimated.",
                call. = FALSE
            )
        }

        voom_fit <- limma::voom(dge, design, plot = FALSE)
        correlation_fit <- limma::duplicateCorrelation(
            voom_fit,
            design,
            block = block
        )
        consensus_correlation <- correlation_fit$consensus.correlation
        if (!is.finite(consensus_correlation)) {
            stop(
                "limma could not estimate a finite within-block correlation.",
                call. = FALSE
            )
        }
        # Apply voom again, now specifying the correlation structure
        voom_fit <- limma::voom(
            dge,
            design,
            plot = FALSE,
            block = block,
            correlation = consensus_correlation
        )
    } else {
        voom_fit <- limma::voom(dge, design, plot = FALSE)
    }
    # Fit the Linear Model
    fit <- limma::lmFit(
        voom_fit,
        design,
        block = block,
        correlation = consensus_correlation
    )
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    # Apply empirical Bayes moderation
    fit <- limma::eBayes(fit)

    if (return_model) {
        return(fit)
    }

    result_list <- lapply(seq_len(ncol(contrast_matrix)), function(coefficient) {
        result <- limma::topTable(
            fit,
            coef = coefficient,
            number = Inf,
            adjust.method = "BH",
            sort.by = "P"
        )
        output <- data.frame(
            gene = rownames(result),
            log2FC = result$logFC,
            avgExpr = result$AveExpr,
            pval = result$P.Value,
            padj = result$adj.P.Val,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        if (include_contrast_column) {
            output$contrast <- colnames(contrast_matrix)[coefficient]
        }
        output
    })

    do.call(rbind, result_list)
}

.limma_contrast_expressions <- function(contrasts) {
    if (is.null(contrasts) || length(contrasts) == 0L) {
        stop(
            "Provide both case_group and reference_group, or provide contrasts.",
            call. = FALSE
        )
    }

    if (is.list(contrasts)) {
        valid <- vapply(
            contrasts,
            function(contrast) is.character(contrast) && length(contrast) == 1L,
            logical(1)
        )
        if (!all(valid)) {
            stop("Each contrast must be a single character string.", call. = FALSE)
        }
        unlisted_contrasts <- unlist(contrasts, use.names = TRUE)
        contrast_expressions <- unname(unlisted_contrasts)
        contrast_names <- names(contrasts)
        if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
            contrast_names <- names(unlisted_contrasts)
        }
    } else if (is.character(contrasts)) {
        contrast_expressions <- contrasts
        contrast_names <- names(contrasts)
    } else {
        stop("Contrasts must be a list or character vector.", call. = FALSE)
    }

    if (is.null(contrast_names) || any(!nzchar(contrast_names))) {
        contrast_names <- make.unique(contrast_expressions)
    }
    stats::setNames(contrast_expressions, contrast_names)
}
