library(Seurat)
# Plot RNA and ADT on the same plot with a bar next to it indicating
# which modality it is from.
# Adapted from FuncPctDotPlot.R code
#
# TODO: 
#   - combine with that code to reduce redundancy

library(tidyverse)
library(dittoSeq)
library(cowplot)
library(gridExtra)

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

FuncDotPlotBar <- function (object, assay = NULL, features, cols = c("lightgrey",
                                                                     "blue"),
                            col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                            idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
                            scale = TRUE, scale.by = "radius",
                            scale.min = NA, scale.max = NA, saturate.func=mean, ...)
{

  
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- (if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  })
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  saturated.features = apply(data.features[, 1:(ncol(x = data.features) - 1), drop = FALSE],
                             MARGIN=2,
                             FUN=function(x){
                               ifelse(is.null(saturate.func), 0, saturate.func(x))
                             })
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = x))  #return(mean(x = expm1(x = x)))
    })
    pct.exp <- sapply(X = colnames(data.use), FUN=function(x) {
      return(Seurat:::PercentAbove(data.use[[x]], threshold = saturated.features[x]))
    })
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  # repair the names!
  data.plot2 = data.plot %>%
    mutate(features.plot=str_replace_all(features.plot, "rna_", ""),
           features.plot=str_replace_all(features.plot, "\\.1$", ""))
  features_repaired = str_replace_all(str_replace_all(features, "\\.1$", ""), "rna_", "")
  data.plot2$features.plot <- factor(x = data.plot2$features.plot, 
                                     levels = features_repaired)
  # create the bar
  longest.id = id.levels[which.max(sapply(id.levels, str_length))]
  tib = tibble("features" = features) %>%
    mutate(assay=ifelse(str_detect(features, "rna"), "RNA", "ADT")) %>%
    mutate(features=features_repaired) %>%
    mutate(id=longest.id)
  tib$features <- factor(tib$features, levels=tib$features)
  bar = ggplot(tib, aes(x=id, y=features, fill=assay))+
    geom_tile()+
    scale_fill_manual(values=c(dittoColors(1), dittoColors()[4]))+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, 
                                   vjust=0.1, color="white"),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_blank(),
          axis.line.x =element_line(color="white"))+
    theme(plot.margin = unit(c(7, 0, 7, 7), "pt"))
  
  bar_legend_only = get_legend(bar)
  
  plot <- ggplot(data = data.plot2, mapping = aes_string(x = "features.plot", 
                                                         y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale),  limits = c(scale.min, scale.max)) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + 
    theme_cowplot()+
    coord_flip()+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.25),
          axis.text.y=element_blank()) +
    theme(plot.margin = unit(c(7, 7, 7, 0), "pt"))
  
  plot <- plot + scale_color_distiller(palette = cols)
  plot_legend_only = get_legend(plot)
  
  
  blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
    cowplot::theme_nothing()
  legend_plot = arrangeGrob(bar_legend_only,
                             plot_legend_only, blankPlot, ncol=1, 
                             heights=c(0.7, 2, 1.5))
  
  p = arrangeGrob(bar+ theme(legend.position="none"), 
               plot+theme(legend.position="none"),
               legend_plot,
               ncol=3, nrow=1,
               widths=c(1, 4, 1)) 
  return(p )
}


### EXAMPLE USAGE ###

# p = FuncDotPlotBar(sobj, # seurat object
#                    assay="ADT", # set to integrated.ADT if you are using that
#                    features=genes_of_interest,  # list of genes
#                    cols="RdYlBu",
#                    cluster.idents=F)
# 
# # note you must save like this
# ggsave("adt_rna_plot.pdf", p, height=12, width=10)
