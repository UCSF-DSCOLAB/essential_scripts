### This script provides a visualization function that is helpful for assessing cluster dispersion
# (Documenting in roxygen syntax to be "library"-ready)
#' Visualize the density of each cluster over the umap space, wraps dittoDimHex
#' @param object A Seurat or SingleCellExperiment object
#' @param clustering String, the name of the metadata to facet by, especially one holding cells' cluster identities. Extra detail: This value will be passed to \code{\link[dittoSeq]{dittoDimHex}} in its \code{split.by} argument.
#' @param scale_legend Number, a scale factor for the size of the density legend. By default, this legend is expanded to 3x the normal legend size in order to be able to reveal the multiple color changes in its lower end. 
#' @param main,data.out,... Parameters passed on to \code{\link[dittoSeq]{dittoDimHex}}
#' @param split.by NULL or a Single String giving the name of a additional metadata to facet by. Not recommended normally, but not blocked as there may be valid use-cases.
#' @param color.var Ignored. Usage of this parameter of \code{\link[dittoSeq]{dittoDimHex}} changes density to be plotted via opacity instead of color, which would defeat the entire purpose of this wrapper function.
#' @param dont_warn Logical. Set to TRUE to block output of \code{color.var} and \code{split.by} input-related warnings.
#' @return A ggplot object where colored hexagonal bins are used to summarize cell density of \code{clustering}-identities across the UMAP (or dimensionality reduction) space.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing
#' @author Daniel Bunis
#' @examples
#' # We'll use the Seurat example dataset for this example
#' sobj <- SeuratObject::pbmc_small
#' 
#' # We'll also make use of the fact that dittoSeq allows use of the string
#' #  "ident" to refererence the 'Idents(seurat_object)' because Seurat
#' #  could change the metadata name sof 'pbmc_small' at any time.
#' 
#' density_plotter(
#'     sobj,
#'     "ident" # or the name of your desired clustering metadata
#' )
#' 
#' # You make the density legend bigger or smaller by adjusting the
#' #  'scale_legend' input.
#' density_plotter(sobj, "ident",
#'     scale_legend = 1
#' )
#' 
#' # You can also make use of additional dittoDimHex inputs.
#' density_plotter(sobj, "ident",
#'     bins = 5,
#'     reduction.use = "pca",
#'     main = "Cluster Densities on PCA"
#' )
#' 
#' # You can also turn off the warning from supplying an extra 'split.by'
#' #  argument by adding 'dont_warn = TRUE'.
#' density_plotter(sobj, "ident",
#'     split.by = "groups",
#'     dont_warn = TRUE)
#'
density_plotter <- function(
        object,
        clustering,
        scale_legend=3,
        main = clustering,
        data.out = FALSE,
        dont_warn = FALSE,
        split.by = NULL,
        color.var = NULL,
        ...) {
    
    warn_if <- function(..., do = !dont_warn) {
        if (do) warning(...)
    }
    
    if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        stop('The dittoSeq package is required to run this function.')
    }
    
    if (!identical(color.var, NULL)) {
        warn_if("Ignoring the requested 'color.var' as density would then be plotted by dittoDimHex using opacity instead of color, and thus defeat the entire purpose of this function.")
    }
    
    if (!identical(split.by, NULL)) {
        warn_if("This function uses splitting to show density of each individual cluster. Additional use of the 'split.by' input is allowed but limited and not recommended.")
        if (length(split.by)>1) {
            warn_if("In favor of faceting by clustering, split.by element '", split.by[2], "' will be ignored.")
            split.by <- split.by[1]
        }
    }
    split.by <- c(clustering, split.by)
    
    legend_breaks <- c(
        0, 0.001, 0.005,
        0.01, 0.05,
        0.1, 0.2, 1)
    
    lightblue <- dittoSeq::dittoColors()[2]
    yellow <- dittoSeq::dittoColors()[4]
    orangered <- dittoSeq::dittoColors()[6]
    legend_colors <- c(
        "grey95",
        dittoSeq::Lighten(lightblue, 0.6),
        dittoSeq::Lighten(lightblue, 0.3),
        dittoSeq::Lighten(yellow, 0.3),
        dittoSeq::Lighten(orangered, 0.3),
        orangered,
        dittoSeq::Darken(orangered, 0.3),
        dittoSeq::Darken(orangered, 0.75))
    mod_plot <- function(p, expand_legend=scale_legend, breaks=legend_breaks, colors=legend_colors) {
        p + 
            ggplot2::scale_fill_gradientn(
                name = "Cells",
                values = breaks, colors = colors,
                breaks = function(limits) {c(1, 50, 200, 500, scales::breaks_extended(5)(limits))}
            ) +
            ggplot2::theme(legend.key.height = ggplot2::unit(1.2*expand_legend, 'lines'))
    }
    
    ditto_out <- dittoSeq::dittoDimHex(
        object, split.by = split.by, main = main, data.out = data.out, ...)
    
    if (data.out) {
        ditto_out$plot <- mod_plot(ditto_out$p)
        ditto_out
    } else {
        mod_plot(ditto_out)
    }
}
