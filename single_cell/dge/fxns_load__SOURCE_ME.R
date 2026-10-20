### This script loads in all functions required for running unified-dge functions
# 'source' this file first to use the 'run_dge' function!
#
# Note: the script would nolonger be required if we packaged up these functions!


.source_dir <- dirname(normalizePath(sys.frame(1)$ofile))

.to_source <- c(
    list.files(
        file.path(.source_dir, "helper_functions"),
        pattern = "\\.R$",
        full.names = TRUE
    ),
    list.files(
        .source_dir,
        pattern = "^run_.*\\.R$",
        full.names = TRUE
    ),
    file.path(.source_dir, "../../R_utils/ts_log.R"),
    file.path(.source_dir, "../pseudobulk_function.R")
)

.to_source <- normalizePath(.to_source, mustWork = TRUE)

for (.file in .to_source) {
    source(.file)
}

rm(.source_dir, .to_source, .file)

