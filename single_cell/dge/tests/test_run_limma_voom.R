source(file.path("single_cell", "dge", "run_limma_voom.R"))
source(file.path("R_utils", "ts_log.R"))
for (helper in list.files(
    file.path("single_cell", "dge", "helper_functions"),
    pattern = "[.]R$",
    full.names = TRUE
)) {
    source(helper)
}
source(file.path("single_cell", "dge", "run_dge.R"))

set.seed(101)
sample_ids <- paste0("sample_", seq_len(12))
metadata <- data.frame(
    condition = rep(c("reference", "case"), each = 6),
    batch = rep(c("batch_1", "batch_2"), 6),
    patient = rep(paste0("patient_", seq_len(6)), 2),
    row.names = sample_ids
)

counts <- matrix(
    stats::rnbinom(40 * 12, mu = 100, size = 10),
    nrow = 40,
    dimnames = list(paste0("gene_", seq_len(40)), sample_ids)
)
counts[seq_len(5), metadata$condition == "case"] <-
    counts[seq_len(5), metadata$condition == "case"] + 150L

pairwise_results <- run_limma_voom(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    case_group = "case",
    reference_group = "reference",
    fixed_effects = "batch"
)
stopifnot(
    identical(
        colnames(pairwise_results),
        c("gene", "log2FC", "avgExpr", "pval", "padj")
    ),
    nrow(pairwise_results) == nrow(counts),
    pairwise_results$log2FC[pairwise_results$gene == "gene_1"] > 0
)

blocked_results <- run_limma_voom(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    case_group = "case",
    reference_group = "reference",
    fixed_effects = "batch",
    random_effects = "patient"
)
stopifnot(nrow(blocked_results) == nrow(counts))

design_model <- run_limma_voom(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    case_group = "case",
    reference_group = "reference",
    fixed_effects = "batch",
    return_model = TRUE
)
stopifnot(inherits(design_model, "MArrayLM"))

custom_results <- run_limma_voom(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    contrasts = c(case_vs_reference = "conditioncase-conditionreference"),
    fixed_effects = "batch"
)
stopifnot(
    identical(
        colnames(custom_results),
        c("gene", "log2FC", "avgExpr", "pval", "padj", "contrast")
    ),
    all(custom_results$contrast == "case_vs_reference")
)

multiple_contrast_results <- run_limma_voom(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    contrasts = list(
        c(up = "conditioncase-conditionreference"),
        c(down = "conditionreference-conditioncase")
    ),
    fixed_effects = "batch"
)
stopifnot(
    nrow(multiple_contrast_results) == 2 * nrow(counts),
    identical(unique(multiple_contrast_results$contrast), c("up", "down"))
)

multiple_random_effects_error <- tryCatch(
    {
        run_limma_voom(
            counts = counts,
            metadata = metadata,
            dge_by = "condition",
            case_group = "case",
            reference_group = "reference",
            random_effects = c("patient", "batch")
        )
        FALSE
    },
    error = function(error) TRUE
)
stopifnot(multiple_random_effects_error)

# Verify that the unified wrapper dispatches the limma method and forwards
# fixed and random effects to run_limma_voom().
metadata$cell_type <- "test_cell"
pipeline_results <- suppressWarnings(run_dge(
    counts = counts,
    metadata = metadata,
    dge_by = "condition",
    method = "limma",
    cell_by = "cell_type",
    case_group = "case",
    reference_group = "reference",
    cell_targets = "test_cell",
    fixed_effects = "batch",
    random_effects = "patient",
    min_frac = 0.5,
    return_model = FALSE
))
stopifnot(
    identical(names(pipeline_results), "test_cell"),
    is.data.frame(pipeline_results[["test_cell"]]),
    identical(
        colnames(pipeline_results[["test_cell"]]),
        c("gene", "log2FC", "avgExpr", "pval", "padj")
    )
)
