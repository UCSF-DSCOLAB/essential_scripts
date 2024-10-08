# essential_scripts
Essential Scripts utilized by members and collaborators of DSCoLab

# Using this repo
It is expected that there will be continued development to many of the scripts contained in this repositiory. Thus, to maintain certainty of what versions were used for your particular projects, it is recommended to create a personal clone of this repo. Then, you can personally control timing of any pulls.

# Contributing
This repo is currently under construction and in a rapid population phase.

During this early phase, if you are a DSCoLab member and there are files that
you would like to add:

- If functionality of your files does not yet exist in the repo:
  - Use the standard branch & PR process to add them, but feel free to merge without waiting for review.
  - Add your files to the Files Index at the end of this README.md
- If functionality of your files replicates / adds to something already in the repo:
  - Start a PR where details can be scoped & reviewed with the current code's author.
  - Feel free to start out by just adding your own version of the code and using the PR's space to scope out the path to merging functionalities.
  - Be sure to update the Files Index below if needed.
  - Final form should coalesce around a single code path.
- Documentation Expectations:
  - For bash functions = create a documentation readout if the function is called without required parameters (or with --help)
  - For R / Python / other scripts and functions = A comment chunk at the top of the script.

# Expectations for all 'Library' scripts

- (This list will likely grows in the future.)
- Running bash script with no arguments, -h, or --help should print usage and purpose of script.
- Programs should have no hard-coded paths and allow users to pass paths in.
- Scripts should not contain passwords or any encryption keys.

# Index
| File | Purpose | Maintainer | Branch |
| --- | --- | --- | --- |
| single_cell/add_module_score_from_excel_gene_sets.R | R code for reading gene sets from an excel file, running Seurat::AddModuleScore, and visualizing the results | Dan | main |
| single_cell/density_plotter.R | When `source()`'d, defines an R function that plots density of clusters across the umap space | Dan | main |
| single_cell/annotation_import.R | When `source()`'d, defines an R function for pulling annotations into Seurat or SCE objects from a csv. An [example 'annots_file'](single_cell/annotation_import_example.csv) and [txt version of the function documentation](single_cell/annotation_import.txt) is also included. | Dan | main | 
| single_cell/CITEseq/process/function.R | When `source()`'d, defines an R function for subclustering Seurat CITEseq data. A script.R is also included to provide example usage. | Dan | db/scrna-process |
| single_cell/RNAseq/process/function.R | When `source()`'d, defines an R function for subclustering Seurat RNAseq data. A script.R is also included to provide example usage. | Dan | db/scrna-process |
| single_cell/RNAseq/dgeMAST_function.R | When `source()`'d, defines an R function dgeMAST() for performing MAST DGE across categories, within cell types, with flexible modeling. | Dan | db/sc-dge-mast |
| count_cores | a command line executable that allows a user to 1) self-monitor their active cores on `krummellab` nodes (default) or 2) use optional flags to query all DSCoLab active jobs to test for core monopoly | Rebecca | main |
| seff | a command line util that will collect time, core, and memory usages stats for a given job; dependency for `core_count` | Rebecca | main |
