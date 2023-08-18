##### This script provides example usage of the 'process_normalized_citeseq_data' function from the adjacent script.
### This script is intended to be copied elsewhere and edited for direct use
### It generates lots of timestamped log messages to always know what's going on.
### Here I'm subsetting on broad annotations that aren't yet in the object, so will be loading those in, and then subsetting to all the Tcell calls

# Replace this path with the path to your own data
full_object_rds_file <- "/path/to/full_data.Rds"
source("/path/to/subcluster_function.R")
# I recommend doing a DRYRUN for double-checking proper subseting and reference setup.
# Set this to FALSE to run processing!
DRYRUN <- TRUE

suppressPackageStartupMessages({
  # These 3 are loaded in sourcing the function's script.
  # library(Seurat)
  # library(harmony)
  # library(dittoSeq)
  library(stringr)
})

options(future.globals.maxSize = 8000 * 1024^2)
print_message <- function(...) {
  cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
}

### Load Original Object / Full Dataset
print_message("Reading previously processed rds object")
all_data <- readRDS(full_object_rds_file)
print_message("Finished reading in full data, summary:")
all_data

### Pull in annotations for subsetting on them.
print_message("Pulling in broad annotations")
annots <- readxl::read_xlsx("annotation_materials__wnn_res.2/annotation.xlsx", sheet = 2)
all_data$annots_broad <- factor(
  all_data$wsnn.harmony.int_res.2,
  levels = annots$cluster,
  labels = annots$broad_annotation
)
print_message("All Broad Annotation Counts")
print(table(all_data$annots_broad))

## Subset to T cells
regex <- c("^T4|^T8|^Tgd|^MAIT")
keep <- grepl(regex, all_data$annots_broad)
print_message(sum(keep), " Total cells being kept")

print_message("KEEPING annotations:\n\t", paste0(unique(all_data$annots_broad[keep]), collapse=", "))
table(all_data$annots_broad[keep])
print_message("REMOVING annotations:\n\t", paste0(unique(all_data$annots_broad[!keep]), collapse=", "))
table(all_data$annots_broad[!keep])

## Make orig.ident a factor to check order & references
all_data$orig.ident <- as.factor(all_data$orig.ident)
print_message("Libraries, ordered:\n\t", paste0(levels(all_data$orig.ident), collapse = "\n\t"))
ref_libs <- grep("SCG1|SCG5", levels(all_data$orig.ident))
print_message("References for ADT integration:\n\t", paste0(paste0(ref_libs, ": ", levels(all_data$orig.ident)[ref_libs]), collapse = "\n\t"))

if (DRYRUN) {
  print_message("DRYRUN=TRUE. Check the above logs for correctness. If all loks good re-run with DRYRUN set to FALSE.")
} else {
  print_message("DRYRUN=FALSE. Moving on with to processing the subset data:")
  process_normalized_citeseq_data(
    all_data[,keep],
    rds_name = 't_data.Rds',
    integration_split_by = "orig.ident",
    integration_references = ref_libs
    # NOTE: There are lots more control-points than just these! See the function definition for details.
  )
}

print_message("ALL DONE")
