density_plotter <- function(
    object,
    clustering,
    split.by = NULL,
    main = clustering,
    scale_legend=3,
    color.var = NULL,
    ...) {

    if (!requireNamespace("dittoSeq")) {
        stop('The dittoSeq package is required to run this function.')
    }

    if (!identical(color.var, NULL)) {
        warning("Ignoring the requested 'color.var' as density would then be plotted by dittoDimHex using opacity instead of color, and thus defeat the entire purpose of this function.")
    }

    if (!identical(split.by, NULL)) {
        warning("This function uses splitting to show density of each individual cluster. Additional use of the 'split.by' input is allowed but limited and not recommended.")
        if (length(split.by)>1) {
            warning("In favor of faceting by clustering, split.by element '", split.by[2], "' will be ignored.")
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
        Lighten(lightblue, 0.6),
        Lighten(lightblue, 0.3),
        Lighten(yellow, 0.3),
        Lighten(orangered, 0.3),
        orangered,
        Darken(orangered, 0.3),
        Darken(orangered, 0.75))

    dittoSeq::dittoDimHex(
        object, split.by = split.by, main = main, ...) +
        scale_fill_gradientn(
            name = "Cells",
            values = legend_breaks, colors = legend_colors,
            breaks = function(limits) {c(1, 50, 200, 500, scales::breaks_extended(5)(limits))}
        ) +
        theme(legend.key.height = unit(1.2*scale_legend, 'lines'))
}
