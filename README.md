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
| File | Purpose | Maintainer | Library |
| --- | --- | --- | --- |
| single_cell/add_module_score_from_excel_gene_sets.R | R code for reading gene sets from an excel file, running Seurat::AddModuleScore, and visualizing the results | Dan | No |
