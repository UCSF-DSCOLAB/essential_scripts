# essential_scripts
Essential Scripts utilized by members and collaborators of DSCoLab

# Using this repo
It is expected that there will be continued development to many of the scripts contained in this repositiory. Thus, to maintain certainty of what versions were used for your particular projects, it is recommended to create a personal clone of this repo. Then, you can personally control timing of any pulls.

When loading functions from an R script, in order to ensure internal cross-repo-file references can work, please use `source("path/to/function/file", chdir = TRUE)`.  (For python, no special equivalent is needed because standard `import` methodology follows relative pathing from the target script's location by defaut!) 

# Contributing
This repo is meant to be used for sharing, in a version controlled way, either heavily-tested and validated code, or code that is simply useful for a collaborator.

- If functionality of your files does not yet exist in the repo:
  - Use the standard branch & commit process to add them.
  - No matter what branch, add your files to the Index table at the end of this README.md.
  - For widely useful items, merging into the 'main' branch is desirable, but we do ask that all functions and script in the main brach meet requirements listed below. Open a PR into 'main' at any time, but only merge after these requirements are met.
- If functionality of your files replicates / adds to something already in the repo:
  - Start a PR where details can be scoped & reviewed with the current code's author.
  - Feel free to start out by just adding your own version of the code and using the PR's space to scope out the path to merging functionalities.
  - Be sure to update the Files Index below if needed.
  - Final form should coalesce around a single code path, and must still meet all 'main'-branch requirements if the code is on the 'main'-branch.

# Expectations for 'main' branch scripts

- Programs should never have hard-coded paths and instead allow users to pass paths in.
  - Only relative paths should be used, specifically for referencing other files within this repo 
- **Scripts should never contain passwords or any encryption keys.**
- **All scripts must be tested by a user other than their creator before merging into 'main'.**
- **All scripts should be documentated**:
  - The more documentation the better!
  - For scripts focused around achieving some function rather than defining a function:
    - Code comment at top which describes the goals
    - Occasional interspersed comments as needed
    - If possible, consider collecting all variable requiring editing from run to run into a chunk at the top of the script
  - For scripts defining functions:
    - Minimally:
      - Overall goal described
      - Each inputs should be described (both type + function)
    - Ideally also:
      - Typing / structure of outputs
      - General description of how the function uses inputs to arrive at outputs
      - Examples
    - Per language requirements:
      - Bash: calling the function without required parameters, or with -h or --help, should yield display of the documentation requirements
        - Example reference: [count_cores](c4_utils/count_cores)
      - R: documentation in roxygen syntax is heavily recommended, as well as publishing via document::document() to create a .txt file rendering due to lack of language support for retaining documentation alongside od source()-acquired functions.
        - Example reference: [essential_scripts/single_cell/density_plotter.R](essential_scripts/single_cell/density_plotter.R)
      - Python: documentation in sphinx syntax is heavily recommended.
          <details>

          <summary>
            Example reference (no python functions to point toward directly yet!)
          </summary>
            
            ```
            def calculate_monthly_payment(principal: float, annual_rate: float, years: int) -> float:
              """
              Calculates the monthly payments of a fixed-rate loan.
          
              Args:
                  principal (float): The initial amount of the loan (the principal).
                  annual_rate (float): The annual interest rate in percentage (e.g., 5.5 for 5.5%).
                  years (int): The total duration of the loan in years.
          
              Returns:
                  float: The amount to be paid monthly.
          
              Raises:
                  ValueError: If the annual_rate is negative or if the years are not positive.
          
              Examples:
                  >>> monthly_payment = calculate_monthly_payment(10000, 5.5, 10)
                  >>> print(f"Monthly Payment: {monthly_payment:.2f}")
                  Monthly Payment: 108.85
          
              """
            ```
          
          </details>


# Index

Note: All files deemed share-worthy can be referenced in this table, but only files in the 'main' branch are guaranteed to have been tested by someone other than their original creator. 

| File | Purpose | Maintainer | Branch |
| --- | --- | --- | --- |
| single_cell/add_module_score_from_excel_gene_sets.R | R code for reading gene sets from an excel file, running Seurat::AddModuleScore, and visualizing the results | Dan | main |
| single_cell/density_plotter.R | When `source()`'d, defines an R function that plots density of clusters across the umap space | Dan | main |
| single_cell/annotation_import.R | When `source()`'d, defines an R function for pulling annotations into Seurat or SCE objects from a csv. An [example 'annots_file'](single_cell/annotation_import_example.csv) and [txt version of the function documentation](single_cell/annotation_import.txt) is also included. | Dan | main | 
| single_cell/CITEseq/subcluster/function.R | When `source()`'d, defines an R function for subclustering Seurat CITEseq data. A script.R is also included to provide example usage. | Dan | main |
| c4_utils/count_cores | a command line executable that allows a user to 1) self-monitor their active cores on `krummellab` nodes (default) or 2) use optional flags to query all DSCoLab active jobs to test for core monopoly | Rebecca | main |
| c4_utils/seff | a command line util that will collect time, core, and memory usages stats for a given job; dependency for `core_count` | Rebecca | main |
