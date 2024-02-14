##### This script provides example usage of the 'process_normalized_rnaseq_data' function from the adjacent script.
### This script is intended to be copied elsewhere and edited for direct use
### It generates lots of timestamped log messages to always know what's going on.
### In this example, I'm subsetting on broad annotations Tcell calls.

# Replace this path with the path to your own data
full_object_rds_file <- "/path/to/full_data.Rds"
source("/path/to/single_cell/RNAseq/cluster/function.R")
# I recommend doing a DRYRUN for double-checking proper subseting and reference setup.
# Set this to FALSE to run processing!
DRYRUN <- TRUE

suppressPackageStartupMessages({
  # These 3 are loaded in sourcing the function's script.
  # library(Seurat)
  # library(harmony)
  library(stringr)
})

print_message <- function(...) {
  cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
}

### Load Original Object / Full Dataset
print_message("Reading previously processed rds object")
all_data <- readRDS(full_object_rds_file)
print_message("Finished reading in full data, summary:")
all_data

## Subset to T cells
regex <- c("^T4|^T8|^Tgd|^MAIT")
keep <- grepl(regex, all_data$annots_broad)
print_message(sum(keep), " Total cells being kept")

print_message("KEEPING annotations:\n\t", paste0(unique(all_data$annots_broad[keep]), collapse=", "))
table(all_data$annots_broad[keep])
print_message("REMOVING annotations:\n\t", paste0(unique(all_data$annots_broad[!keep]), collapse=", "))
table(all_data$annots_broad[!keep])

if (DRYRUN) {
  print_message("DRYRUN=TRUE. Check the above logs for correctness. If all looks good re-run with DRYRUN set to FALSE.")
} else {
  print_message("DRYRUN=FALSE. Moving on with to processing the subset data:")
  process_normalized_rnaseq_data(
    all_data[,keep],
    rds_name = 't_data.Rds'
    # NOTE: There are lots more control-points than just these! See the function definition for details.
  )
}

print_message("ALL DONE")
