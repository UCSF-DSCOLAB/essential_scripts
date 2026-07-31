library(tidyverse)
library(plotly)
library(htmlwidgets)

#' Create a summary heatmap of enrichment analysis results
#'
#' Generates an interactive heatmap from filtered enrichment results and 
#' exports the visualization to a standalone HTML file. Heatmap is colored by 
#' \code{signed -log10(padj)}, with red or blue values indicating positive or 
#' negative enrichment, respectively. The heatmap is generated using ployly; 
#' hovering a cursor over a heatmap cell shows enrichment statistics and 
#' leading-edge features. The plot is printed in an HTML file to retain its
#' interactive functionality.
#' 
#' @author Ravi K. Patel
#' @param df Data frame of filtered (e.g., padj < 0.1) enrichment analysis 
#'   results with the following columns:
#'   \itemize{
#'     \item \code{pathway}: Name of the pathway.
#'     \item \code{padj}: Adjusted p-value.
#'     \item \code{NES}: Normalized enrichment score or direction of enrichment (e.g., -1 vs 1).
#'     \item \code{leadingEdge}: String of leading-edge features, separated by \code{le_separator}.
#'     \item \code{dataset}: Name of the source dataset (name of the celltype/comparison used for 
#'     calculating the input log2FC for the enrichment analysis). 
#'     The \code{dataset} will become columns of the heatmap.
#'   }
#' @param order_by Cluster the heatmap by adjusted p-values (\code{padj}) or
#'   by the overlap between leading-edge features (combined across rows for 
#'   clustering columns and vice versa) between pathways (\code{leadingEdge}).
#'   Clustering by \code{padj} does put similar pathways together due to a high 
#'   chance for similar p-values between pathways with highly overlapping features. 
#'   But clustering by leadingEdge tries to do that explicitly, though this 
#'   functionality is not fully tested yet.
#' @param le_separator Separator used to separate leading-edge features in \code{leadingEdge} column of \code{df}. Default ", " (a command and a space).
#' @param html_output_path File path for the output HTML. (e.g., \code{/path/to/my.html})
#' @param plot_title Plot title displayed above the plot.
#' @param return_ggplot A boolean indicating if the plot should be returned or be printed to an HTML file.
#'
#' @return The function has no R return value; the heatmap is exported 
#'   directly to the specified HTML file.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # An example data frame to demonstrate the structure of the input data frame.
#' my_gsea_results = data.frame(
#'     pathway = c('LABCD8_IFN','LABCD8_IFN_Module','LABCD8_ISG','KEGG_PROTEASOME','REACTOME_ABC_TRANSPORTER_DISORDERS','REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION'),
#'     padj = c('0.00221604027912997','0.0717798958531707','0.000731840566458411','0.0467085946901833','0.0571030113172619','0.0746338661273512'),
#'     NES = c('1.91638152494958','1.74199063631165','2.03991528803013','1.91334841630672','1.62183620773234','1.53659018078308'),
#'     leadingEdge = c('GBP5, TRIM21, OAS2, STAT1, EPSTI1, GBP4', 'GBP2, TAP1, STAT1, EPSTI1, IFI16, APOBEC3D, GBP4',  'PSMC4, PSMA5, PSMD4, PSME2, PSMB2, PSMB10, PSMB9, PSMB4, SEM1',  'RNF5, PSMC4, PSMA5, PSMD4, PSME2, PSMB2', 'PSMB10, PSMD10, PSMB9, PSMB4', 'TAP1, IKBKG, LY96, PSMC4, PSMA5, TAP2, PSMD4, PSME2, PSMB2, PSMB10, PSMD10K'),
#'     dataset = c('CD4_CTL','B_mem_CD11c','B_mem_not_switched','B_mem_switched','B_naive','CD4_CTL')
#' )
#' 
#' # Make a heatmap
#' make_enrichment_heatmap(
#'   df = my_gsea_results, 
#'   order_by = "padj",
#'   html_output_path = "my_gsea_heatmap.html",
#'   plot_title = "Pathway Enrichment",
#'   return_ggplot = FALSE
#' )
#' } 
#' 
#' plot = make_enrichment_heatmap(
#'   df = my_gsea_results, 
#'   order_by = "leadingEdge",
#'   le_separator = ", ",
#'   plot_title = "Pathway Enrichment",
#'   return_ggplot = TRUE
#' )
#' } 
make_enrichment_heatmap <- function(df, order_by = c("leadingEdge", "padj"), 
                                    le_separator = ", ",
                                    html_output_path = "../analysis/10x/dge_analysis/tmp/interactive_heatmap.html", 
                                    plot_title="Enrichment analysis results",
                                    return_ggplot = FALSE) {
  
  # Calculate signed padj
  df$padj_signed = sign(df$NES) * -log10(df$padj)
  
  if(order_by == "leadingEdge") {
    orders = order_by_le_overlap(df, le_separator)
  } else if(order_by == "padj") {
    orders = order_by_padj(df)
  }
  
  row_order = orders$row_order
  col_order = orders$col_order
  
  df = df %>% 
    select(c("pathway","padj","NES","leadingEdge","dataset","padj_signed")) %>% 
    mutate(padj_lbl = format.pval(padj, 1),
      leadingEdge = paste0(dataset, "<br>", pathway, "<br>padj: ", padj_lbl, "|NES: ", round(NES,2), "<br>", leadingEdge))
  df$leadingEdge = sapply(df$leadingEdge, function(x) {
    genes <- strsplit(x, le_separator)[[1]]
    groups <- split(genes, ceiling(seq_along(genes) / 8))
    lines <- sapply(groups, paste, collapse = le_separator)
    result <- paste(lines, collapse = paste0(le_separator, "<br>"))
    result
  })
  
  
  df <- df %>%
    mutate(
      pathway = factor(pathway, levels = rev(row_order)),   # reverse for top-to-bottom in ggplot
      dataset   = factor(dataset, levels = col_order)
    )
  
  
  gg <- ggplot(df, aes(x = dataset, y = pathway, fill = padj_signed, text = leadingEdge)) +
    geom_tile() +
    labs(x = "Cell types", y = "Pathways", title = plot_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
    scale_fill_distiller(palette = "RdBu")
  
  
  # assume your interactive object is:
  p <- ggplotly(gg, tooltip = "text") %>% layout(
    height = length(unique(df$pathway))*12 + (max(nchar(as.character(df$dataset))) * 10),
    width = length(unique(df$dataset))*14 + (max(nchar(as.character(df$pathway))) * 10)
    )
  
  if(return_ggplot)
    return(gg)
  
  # save as self-contained HTML
  saveWidget(
    p,
    file = html_output_path,
    selfcontained = TRUE
  )
  
  
  html_lines <- readLines(html_output_path, warn = FALSE)
  
  # defining <div> with scolling ability
  div_style <- c(
    '<div id="htmlwidget_container" style="height: 90vh; overflow-y: auto; overflow-x: auto;">'
  )
  
  # Finding the <div> to be replaced
  div_idx <- grep('<div id="htmlwidget_container">', html_lines, ignore.case = TRUE)[1]
  
  if (! is.na(div_idx)) {
    
    # Insert <div> with style to enable scrolling
    new_html <- c(
      html_lines[1:(div_idx - 1)],
      div_style,
      html_lines[div_idx:length(html_lines)]
    )
    
    # Write back to file
    writeLines(new_html, html_output_path)
  }
  cat("Interactive heatmap is exported in ", html_output_path)
}


#' Calculate order of rows and columns of a heatmap based on the overlap between leading-edge features.
#' For calculating the order of rows, the leading-edge features across columns are combined, pair-wise overlap matrix is calculated, and 
order_by_le_overlap <- function(df, le_separator) {
  le = sapply(df$leadingEdge, function(x) unlist(str_split(x, le_separator)) )
  
  # Order of rows / pathways
  le_by = list()
  for(i in 1:nrow(df)) {
    le_by[[ df[i,"pathway"] ]] = unique(c(le_by[[ df[i,"pathway"] ]], le[[i]]))
  }
  feats = unique(df[,"pathway"])
  sim_mat = matrix(0, nrow=length(feats), ncol=length(feats))
  for( i in 1:(length(feats)) ) {
    for(j in (i):length(feats) ) {
      sim_mat[i,j] = length( intersect( le_by[[feats[i]]], le_by[[feats[j]]] ) )
    }
  }
  sim_mat[lower.tri(sim_mat)] = t(sim_mat)[lower.tri(sim_mat)]
  
  # sim_mat has similarity values, i.e., the number of leading edge features shared by two pathways. We need to convert it to a distance matrix for hclust.
  # I will normalize by the Maximum Possible Overlap, i.e., divide each overlap score by the smaller size of the two vectors being compared.
  # 1. Extract the diagonal (self-overlap reveals the size of each vector)
  vector_sizes <- diag(sim_mat)
  # 2. Create a matrix of the minimum possible size for each pair
  min_size_matrix <- outer(vector_sizes, vector_sizes, pmin)
  # 3. Normalize similarity to a 0-1 scale
  norm_sim <- sim_mat / min_size_matrix
  # 4. Convert to a distance object
  dist_mat <- as.dist(1 - norm_sim)
  
  row_hc = hclust(dist_mat, method = "ward.D2")
  row_order = feats[row_hc$order]
  
  # Order of cols / datasets
  le_by = list()
  for(i in 1:nrow(df)) {
    le_by[[ df[i,"dataset"] ]] = unique(c(le_by[[ df[i,"dataset"] ]], le[[i]]))
  }
  feats = unique(df[,"dataset"])
  sim_mat = matrix(0, nrow=length(feats), ncol=length(feats))
  for( i in 1:(length(feats)) ) {
    for(j in (i):length(feats) ) {
      sim_mat[i,j] = length( intersect( le_by[[feats[i]]], le_by[[feats[j]]] ) )
    }
  }
  sim_mat[lower.tri(sim_mat)] = t(sim_mat)[lower.tri(sim_mat)]
  
  # sim_mat has similarity values, i.e., the number of leading edge features shared by two pathways. We need to convert it to a distance matrix for hclust.
  # I will normalize by the Maximum Possible Overlap, i.e., divide each overlap score by the smaller size of the two vectors being compared.
  # 1. Extract the diagonal (self-overlap reveals the size of each vector)
  vector_sizes <- diag(sim_mat)
  # 2. Create a matrix of the minimum possible size for each pair
  min_size_matrix <- outer(vector_sizes, vector_sizes, pmin)
  # 3. Normalize similarity to a 0-1 scale
  norm_sim <- sim_mat / min_size_matrix
  # 4. Convert to a distance object
  dist_mat <- as.dist(1 - norm_sim)
  
  col_hc = hclust(dist_mat, method = "ward.D2")
  col_order = feats[col_hc$order]
  
  return(list(row_order=row_order, col_order=col_order))
}

order_by_padj <- function(df) {
  mat = pivot_wider(df[,c("pathway","dataset","padj_signed")], names_from = "dataset", values_from = "padj_signed") %>% column_to_rownames("pathway")
  mat[is.na(mat)] = 0
  row_dist <- dist(mat, method = "euclidean")
  col_dist <- dist(t(mat), method = "euclidean")
  
  row_hc <- hclust(row_dist, method = "ward.D2")
  col_hc <- hclust(col_dist, method = "ward.D2")
  
  # get orderings
  row_order <- rownames(mat)[row_hc$order]   # from top-to-bottom in dendrogram leaves order
  col_order <- colnames(mat)[col_hc$order]
  
  return(list(row_order=row_order, col_order=col_order))
}
