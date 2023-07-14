### Scoring & visualizing Seurat object cells' by expression of Gene Signatures
# Originall by Dan Bunis
# This script gives *example* code for:
# - Reading in sets of genes, a.k.a. gene signatures / modules, from excel files
# - Using Seurat::AddModuleScore() to calculate expression scores per module
# - Renaming the score metadata (because Seurat's number method is meh to me.)
# - A few example visualizations using dittoSeq
# It does not:
# - Give full code for all loading in all signatures, just an example for each excel-tab structure type.
# - Give code for loading in the target Seurat object.
# Also note:
# - treats the target Seurat object as being named 'sobj'

library(Seurat)
library(readxl)
library(dittoSeq)
timestamped_message <- function(...) {
  cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
}

##### Reading in genesets #####
# Structure of the target gene set-containing excel file is 1 tab per reference
# and those tabs are of 2 distinct structures, so example code below for both:
#
# The plan is to build a single named list object holding all genesets so that
# renaming step later will be easiest.
gene_sets <- list()

### Structure1: "multi cell type DGE table", Important features:
#  - Has a column holding cell-group annotations
#  - Has another column holding gene names
#  - All genes are positively associated with the labeled cell-groups (upregulated in, correlated with, etc.)
this_tab <- read_xlsx("Gene_lists_reformatted.xlsx", sheet = 1)
View(this_tab)
# Remove some rows that are all NA at the bottom...
this_tab <- this_tab[!is.na(this_tab$Associated.Cell.Population), ]

# Column 'Gene.ID' contains the genes
# Column 'Associated.Cell.Population' contains the associated annotations
this_set_names <- unique(this_tab$Associated.Cell.Population)
this_set <- lapply(
  this_set_names,
  function(celltype) {
    this_tab$Gene.ID[this_tab$Associated.Cell.Population==celltype]
})

# Let's add an extra tag to these gene_set names to help us know what they are...
names(this_set) <- paste0("vilani2017_DCs__", this_set_names)

# Add to the collection
gene_sets <- c(gene_sets, this_set)

### Structure2: "collection of sets", Important features:
#  - Collection of columns each holding a single geneset
#  - Top column is the set names
#  - All rows below contain the genes in those sets
#  - Shorter sets get NA values from empty excel cells
this_tab <- read_xlsx("Gene_lists_reformatted.xlsx", sheet = 2)
View(this_tab)
# Remove the 'Notes' column which is not a geneset 
this_tab <- this_tab[, grep("^Notes", colnames(this_tab), invert = TRUE, value = TRUE)]

this_set <- lapply(
  names(this_tab),
  function(celltype) {
    set <- this_tab[, celltype]
    set[!is.na(set)]
  })
names(this_set) <- paste0("vilani2017_cytochemo__", names(this_tab))

# Add to the collection
gene_sets <- c(gene_sets, this_set)

### Trim all sets to genes inside our object
# IMPORTANT TO NOTE WARNINGS HERE.
# If too many genes are trimmed form a particular set, you may not want to use it.
gene_set_names <- names(gene_sets)
gene_sets <- lapply(
  gene_set_names,
  function(set) {
    genes <- gene_sets[[set]]
    missing <- genes[!isGene(genes, object = sobj, assay = "RNA")]
    if (length(missing)>0) {
      warning(length(missing), " genes out of ", length(genes), " from ", set, " are missing from the Seurat object: ", paste0(missing, collapse=", "))
    }
    isGene(genes, object = sobj, assay = "RNA", return.values = TRUE)
  }
)
names(gene_sets) <- gene_set_names

##### Run Scoring #####
# Seurat's AddModuleScore function does this for us, so it's pretty simple.
# The one thing I don't like about the system is just its naming scheme, so we'll fix that in the next section
# Per the docs, ?AddModuleScore, it calculates average expression for genes of each geneset,
# then subtracts the aggregated expression of a set of (internally determined) control features.

# Let's first decide on a prefix for the temporary metadata names.
# It can be anything as long as we won't end up overwriting old metadata.
# This name choice is probably low risk, but you can check getMetas(sobj) to be sure!
module_prefix <- "score_for_module_"

# This will take some time. Took 3 minutes for 8 sets in my own testing with 146k cells.
timestamped_message("Starting AddModuleScore() on ", length(gene_sets), " modules.")
sobj <- AddModuleScore(
  sobj,
  features = gene_sets,
  name = module_prefix # We'll rename this later.
)
timestamped_message("AddModuleScore()'ing complete.")

##### Rename #####

# Determine current names
temp_module_metas <- grep(module_prefix, getMetas(sobj), value = TRUE)
names(temp_module_metas) <- paste0("mod", str_extract(temp_module_metas, "(\\d)+$"))

# Build final names
# This will get a bit long, but to be have a bit more description inside, let's add "module_score__" to the start of our gene_set names.
gene_set_names_final <-  paste0("modulescore__", names(gene_sets))

# **Perhaps necessary**: Removing special characters.
# dittoHeatmap's underlying heatmap plotters can error from special characters
# Feel free to skip if not using those and the plotters work without this fix.
# We'll replace any special characters with dots here.
gene_set_names_final <- str_replace_all(gene_set_names_final, "[^[:alnum:]_]", ".")

# Rename the metadata column matching each temporary name
for (i in seq_along(gene_sets)) {
  names(sobj@meta.data)[which(names(sobj@meta.data)==temp_module_metas[paste0("mod",i)])] <- gene_set_names_final[i]
}
rm(temp_module_metas)

# Now, each module should be named in the way we want!
# Let's look. They should be at the end of this output.
getMetas(sobj)

##### Visualizations #####
# Note: swap out the "score__vilani2017_DCs__CD141/CLEC9A" used here for any of your own gene_set names!
# 
# Note2: For the plotters which show multiple gene set scores, you can also make a subset of your gene_set_names to use instead of showing all of them.
#   2 versions of optional code for this Note2.
#   Use/modify however you need, knowing that you can the variable whatever you like.
#   Afterwards, just use your own variable instead of the 'gene_set_names_final' in the visualization funciton calls!
# gene_set_subset <- gene_set_names_final[1:3]
# gene_set_subset <- c(
#   "modulescore__vilani2017_DCs__CD141/CLEC9A",
#   "modulescore__vilani2017_DCs__CD1C_A",
#   "modulescore__vilani2017_DCs__CD1C_B"
# )
#
# Note3: I use clustering and a clinical metadata as my groupings metadata examples here,
#   but annotations are another option you might want to target!

### On the umap:
# For the sobj I know will be used, the final umap is called "umap_harmony", but in many cases
# you will be able to skip the 'reduction.use = ' line and the function will use your 'umap'-named dimensionality reduction by default.

# Single
# Also a 'order' focus which can be applied to other dittoDimPlot code as well!
dittoDimPlot(
  sobj, "modulescore__vilani2017_DCs__CD141.CLEC9A",
  # This next bit will put the cells with highest scores in front.
  # It's good to look without this as well to be sure doing so does not bias the plot too much!
  order = "increasing"
)

# Or multiple!
dittoDimPlot(
  sobj, gene_set_names_final, # <- Just give the vector in place of the single one here!
)

### By violin plot, grouped by some important grouping
# Useful groupings might be clinical or clusters
# You can group by two things by making use of 'color.by' as well

# Single, one grouping variable
dittoPlot(
  sobj, "modulescore__vilani2017_DCs__CD141.CLEC9A",
  # Use this to set your groupings.  Using clustering here.
  group.by = "harmonized_res.2"
)

# Single, tweaked to ensure violin visibility with MANY groups
dittoPlot(
  sobj, "modulescore__vilani2017_DCs__CD141.CLEC9A",
  group.by = "harmonized_res.2",
  # The tweak:
  vlnplot.scaling = "width"
)

# Single, two grouping variables, subset of data
# Also a 'cells.use' focus which can be applied in any of the other examples as well!
dittoPlot(
  sobj, "modulescore__vilani2017_DCs__CD141.CLEC9A",
  group.by = "harmonized_res.2",
  color.by = "Iatrogenic.spontaneous",
  # Optional subsetting to certain clusters
  cells.use = sobj$harmonized_res.2 %in% 1:5
)

# Multiple
# (With same options extras as 2nd example above.  Simplest version can be just the first line of inputs!)
dittoPlot(
  sobj, gene_set_names_final, # <- Just give the vector in place of the single one here!
  group.by = "harmonized_res.2",
  color.by = "Iatrogenic.spontaneous",
  # Optional subsetting to certain clusters
  cells.use = sobj$harmonized_res.2 %in% 1:5
)

# DotPlot
# Note: Makes sense only when you want to look at multiple modules together
# Also, this is really a quick stand-in for making a heatmap as there is no per-group summary heatmap in dittoSeq at the moment...

# Scaled
# Scaling is nice because it makes differences between groups easier to see.
dittoDotPlot(
  sobj, gene_set_names_final,
  group.by = "harmonized_res.2",
  scale = TRUE
)

# Unscaled
# I recommend looking at the unscaled version too, because sometimes scaling blow tiny differences WAY out of proportion!
dittoDotPlot(
  sobj, gene_set_names_final,
  group.by = "harmonized_res.2",
  scale = FALSE
)

# Multiple groupings
# There is no color.by for this function.  Instead, use 'split.by' for secondary cell groupings & 'group.by' for what you want to directly compare between.
dittoDotPlot(
  sobj, gene_set_names_final,
  split.by = "harmonized_res.2",
  group.by = "Iatrogenic.spontaneous",
  # (I still recommend running with both TRUE / FALSE below.)
  scale = FALSE
)
