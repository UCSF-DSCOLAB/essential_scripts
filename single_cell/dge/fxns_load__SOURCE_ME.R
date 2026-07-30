### This script loads in all functions required for running unified-dge functions
# 'source' this file first to use the 'run_dge' function!
#
# Note: the script would nolonger be required if we packaged up these functions!

.to_source <- c(
    list.files('helper_functions', pattern = '.R', full.names = TRUE),
    list.files(pattern = 'run_.*\\.R', full.names = TRUE),
    '../../R_utils/ts_log.R'
)
for (.file in .to_source) {
    source(.file)
}
rm(.to_source, .file)
