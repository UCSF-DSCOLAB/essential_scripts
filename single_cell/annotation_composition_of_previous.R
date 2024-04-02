#' Add previous annotation composition statistics to a cluster annotation dataframe
#' @param annotationDF A data.frame-like object or String, giving the file path of a csv or tsv encoding such an object, in which each line represents a given cluster.
#' @param object A Seurat or SingleCellExperiment object, or data.frame representing per-cell metadata
#' @param transfer_target String naming a metadata column of \code{object} to holding the previous annotations you wish to quantify within each cluster
#' @param clustering NA or String naming a metadata column of \code{object} which holds the cluster-identities of each cell.
#' If left as NA, the name of the first column of \code{annotationDF} will be used and must exist in \code{object}.

add_previous_annotation_composition <- function(
    # Primary inputs
    annotationDF, object, transfer_target,
    # Optional inputs
    clustering = NA, decimals = 4, next_tops = 4,
    # For reading in from file
    read_sep = NA, read_check.names = FALSE, ...,
    # Logging
    verbose = TRUE, timestamps = FALSE) {

    # Until we decide a better way to handle this function I use all the time,
    # simply adding this here, again
    timestamped_message <- function(...) {
        cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
    }
    log <- function(..., do = verbose, timestamp = timestamps) {
        if (do) {
            if (timestamp) {
                timestamped_message(...)
            } else {
                cat(..., "\n", sep = "")
            }
        }
    }

    if (!requireNamespace("dittoViz", quietly = TRUE)) {
        stop("dittoViz installation required for this function. `remotes::install_github('dtm2451/dittoViz')`")
    }

    # Get annotationDF
    if (is(annotationDF, "character")) {
        if (identical(read_sep, NA)) {
            read_sep <- ifelse(endsWith(annotationDF, ".csv"), ",", "\t")
        }
        log("Reading in annotationDF from ", annotationDF)
        annotationDF <- read.csv(
            annotationDF, sep = read_sep, check.names = read_check.names, ...
        )
    }

    # Get metadata, if not given directly
    if (is(object,"SingleCellExperiment") || is(object,"Seurat")) {
        if (!requireNamespace("dittoSeq")) {
            stop("dittoSeq installation required for starting from a SingleCellExperiment or Seurat object")
        }
        meta <- dittoSeq::getMetas(object, FALSE)
    }

    # Interpret clustering
    if (identical(clustering, NA)) {
        clusts <- annotationDF[,1]
        clustering <- names(annotationDF)[1]
        if (! clustering %in% names(meta)) {
            stop("The first column name of annotationDF does not exist in the given cell metadata. Do you need to give 'clustering'? Or is it accidentally missing from 'object'?")
        }
    } else {
        clusts <- dittoViz::colLevels(clustering, meta)
    }

    # Determine Compositions
    log(
        "Calculating compositions of ", transfer_target, " per ", clustering,
        " clusters, using dittoViz::barPlot(..., data.only=TRUE)")
    data <- dittoViz::barPlot(meta, transfer_target, clustering, data.only = TRUE)

    log("Summarizing top compositions")
    transfer_counts <- table(meta[,transfer_target])

    transfer_summary <- do.call(rbind, lapply(
        clusts, function(this_clust) {
            this_data <- data[data$grouping==this_clust,]

            max_ind <- which.max(this_data$count)
            max <- as.character(this_data$label[max_ind])
            max.percent <- round(this_data$percent[max_ind], decimals)
            percent.of.ref <- round(this_data$count[max_ind]/transfer_counts[max], decimals)

            next_tops_str <- ""
            num_top <- 2
            this_data <- this_data[this_data$label != max,]
            while (next_tops > 0) {
                new_max_ind <- which.max(this_data$count)
                new <- as.character(this_data$label[new_max_ind])
                new.percent <- round(this_data$percent[new_max_ind], decimals)
                new.percent.of.ref <- round(this_data$count[new_max_ind]/transfer_counts[new], decimals)
                next_tops_str <- paste0(
                    next_tops_str, ifelse(num_top>2, " | ", ""),
                    num_top, ": ",
                    new, "; ",
                    new.percent, " of cluster; ",
                    new.percent.of.ref, " of ", new
                )
                num_top <- num_top + 1
                next_tops <- next_tops - 1
                this_data <- this_data[this_data$label != new,]
            }

            data.frame(
                cluster = this_clust,
                transfer_max = max,
                percent_max_of_this_cluster = max.percent,
                percent_max_here_of_total_max = percent.of.ref,
                transfer_next_tops = next_tops_str
            )
        }
    ))

    log("Merging with annotationDF, then outputting")
    merge(annotationDF, transfer_summary, by.x = 1, by.y = "cluster")
}
