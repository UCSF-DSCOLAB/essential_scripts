### This script loads in all functions required for running unified-dge functions
# 'source' the file first to use 'run_<method>' or 'run_<method>_within_cells' functions.
#
# Note: the script would nolonger be required if we packaged up these functions!

.to_source <- c(
    list.files('helper_functions', pattern = '.R', full.names = TRUE),
    '../../R_utils/ts_log.R'
)
for (.file in .to_source) {
    source(.file)
}
rm(.to_source, .file)