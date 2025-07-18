#' @name gl.grm.network
#' @title Represents a similarity matrix as a network
#' @description
#' This script takes any similarity matrix and represents
#' the relationship among the specimens as a network diagram. 
#'
#' @param G A similarity matrix [required].
#' @param x A genlight object from which the matrix was generated [required].
#' @param standardise Whether to standardise matrix using Goudet et al method, 
#' see details [default FALSE].
#' @param method One of 'fr', 'kk', 'gh' or 'mds' [default 'fr'].
#' @param node.size Size of the symbols for the network nodes [default 8].
#' @param node.label TRUE to display node labels [default TRUE].
#' @param node.label.size Size of the node labels [default 3].
#' @param node.label.color Color of the text of the node labels
#' [default 'black'].
#' @param link.color  Colors for links [default gl.select.colors].
#' @param link.size Size of the links [default 2].
#' @param relatedness_factor Factor of relatedness [default 0.125].
#' @param title Title for the plot [default 'Network of similarity matrix'].
#' @param legend.title Title for the legend [default "Populations"].
#' @param palette_discrete A discrete set of colors
#'  with as many colors as there are populations in the dataset
#'  [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir]
#' @param plot.file Name for the RDS binary file to save (base name only,
#'  exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @details
#' The gl.grm.network function creates a network diagram that represents 
#' genetic relationships among individuals in a dataset. 
#' 
#'\strong{Layout options}
#'  
#' Four layout options are implemented in this function:
#' \itemize{
#' \item 'fr' Fruchterman-Reingold layout  \link[igraph]{layout_with_fr}
#' (package igraph)
#' \item 'kk' Kamada-Kawai layout \link[igraph]{layout_with_kk} (package igraph)
#' \item 'gh' Graphopt layout \link[igraph]{layout_with_graphopt}
#' (package igraph)
#' \item 'mds' Multidimensional scaling layout \link[igraph]{layout_with_mds}
#' (package igraph)
#' }
#' 
#'  \strong{Standardise matrix using Goudet et al method}
#'  
#' Choosing meaningful thresholds to represent relationships between 
#' individuals can be challenging because kinship and inbreeding coefficients 
#' are relative measures. To standardize a genomic relationship matrix (GRM),
#' such as the one produced by the function gl.grm, and facilitate interpretation, 
#' the function adjusts the matrix through the following steps:
#' 
#' 1. Centering Inbreeding Coefficients: Subtract 1 from the mean of the 
#' diagonal elements to calculate the average inbreeding coefficient. This 
#' centers the inbreeding coefficients around zero, providing a reference point 
#' relative to the population's average inbreeding level.
#' 
#' 2. Calculating Kinship Coefficients: Divide the off-diagonal elements by 2 
#' to obtain the kinship coefficients. This conversion reflects the probability 
#' of sharing alleles IBD between pairs of individuals.
#' 
#' 3. Centering Kinship Coefficients: Subtract the adjusted mean inbreeding 
#' coefficient (from step 1) from each kinship coefficient (from step 2). This 
#' centers the kinship coefficients relative to the population average, 
#' allowing for meaningful comparisons.
#' 
#' This adjustment method aligns with the approach used by Goudet et al. (2018), 
#' enabling the relationships to be interpreted in the context of the overall 
#' genetic relatedness within the population.
#'
#' Below is a table modified from Speed & Balding (2015) showing kinship values,
#'  and their confidence intervals (CI), for different relationships that could
#'  be used to guide the choosing of the relatedness threshold in the function.
#'
#' \tabular{lll}{
#'   \strong{Relationship} \tab \strong{Kinship} \tab \strong{95\% CI} \cr
#'   Identical twins / clones / same individual \tab 0.5   \tab –              \cr
#'   Sibling / Parent–Offspring                \tab 0.25  \tab (0.204, 0.296)\cr
#'   Half‑sibling                              \tab 0.125 \tab (0.092, 0.158)\cr
#'   First cousin                              \tab 0.062 \tab (0.038, 0.089)\cr
#'   Half‑cousin                               \tab 0.031 \tab (0.012, 0.055)\cr
#'   Second cousin                             \tab 0.016 \tab (0.004, 0.031)\cr
#'   Half‑second cousin                        \tab 0.008 \tab (0.001, 0.020)\cr
#'   Third cousin                              \tab 0.004 \tab (0.000, 0.012)\cr
#'   Unrelated                                 \tab 0     \tab –              \cr
#' }
#' 
#'
#' @return A network plot showing relatedness between individuals
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE) & requireNamespace("rrBLUP",
#'   quietly = TRUE
#' ) & requireNamespace("fields", quietly = TRUE)) {
#'   t1 <- possums.gl
#'   # filtering on call rate
#'   t1 <- gl.filter.callrate(t1)
#'   t1 <- gl.subsample.loc(t1, n = 100)
#'   # relatedness matrix
#'   res <- gl.grm(t1, plotheatmap = FALSE)
#'   # relatedness network
#'   res2 <- gl.grm.network(res, t1, relatedness_factor = 0.125)
#' }
#' @references
#' \itemize{
#' \item Endelman, J. B. , Jannink, J.-L. (2012). Shrinkage estimation of the
#' realized relationship matrix. G3: Genes, Genomics, Genetics 2, 1405.
#' \item Goudet, J., Kay, T., & Weir, B. S. (2018). How to estimate kinship.
#' Molecular Ecology, 27(20), 4121-4135.
#' \item Speed, D., & Balding, D. J. (2015). Relatedness in the post-genomic era:
#' is it still useful?. Nature Reviews Genetics, 16(1), 33-44.
#'  }
#' @seealso \code{\link{gl.grm}}
#' @family inbreeding functions
#' @export

gl.grm.network <- function(G,
                           x,
                           standardise = FALSE,
                           method_relatedness = "",
                           method = "fr",
                           node.size = 8,
                           node.label = TRUE,
                           node.label.size = 2,
                           node.label.color = "black",
                           link.color = NULL,
                           link.size = 2,
                           relatedness_factor = 0.125,
                           title = "Network of a similarity matrix",
                           legend.title = "Populations",
                           palette_discrete = gl.select.colors(x, 
                                                               library = "brewer", 
                                                               palette = "PuOr", 
                                                               ncolors = nPop(x), 
                                                               verbose = 0),
                           plot.dir = NULL,
                           plot.file = NULL,
                           verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING Set a population if none is specified
  # (such as if the genlight object has been generated manually)
  if (is.null(pop(x)) |
      is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2) {
      cat(
        important(
          "  Population assignments not detected, individuals assigned
                    to a single population labelled 'pop1'\n"
        )
      )
    }
    pop(x) <- array("pop1", dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  }
  
  # check if package is installed
  pkg <- "igraph"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  if (!(method == "fr" ||
        method == "kk" ||
        method == "gh" || 
        method == "mds")) {
    if(verbose>0) cat(warn(
      "Warning: Layout method must be one of fr, or kk, gh or mds, set to fr\n"
    ))
    method <- "fr"
  }
  
  # DO THE JOB
  G2 <- as.matrix(G)
  G2[upper.tri(G2, diag = TRUE)] <- NA
  links <- as.data.frame(as.table(G2))
  links <- links[which(!is.na(links$Freq)), ]
  
  colnames(links) <- c("from", "to", "weight")
  
  if(isTRUE(standardise)){
    # using the average inbreeding coefficient (1-f) of the diagonal elements as
    #the reference value
    MS <- mean(diag(G) - 1)
    # the result of the GRM is the summation of the IBD of each allele .
    links$kinship <- (links$weight / 2) - MS
    links_tmp <- links[,c(1,2,4)]
  }
  
  if(isFALSE(standardise)){
    links_tmp <- links
    colnames(links_tmp) <- c("from","to","kinship")
    links_tmp$kinship <- as.character(links_tmp$kinship)
    colnames(links) <- c("from","to","kinship")
    links_tmp <- links
  }
  
  links_tmp <- rbind(links_tmp,cbind(from=indNames(x),
                                     to=indNames(x),kinship=0))
  
  links_matrix <- as.matrix(reshape2::acast(links_tmp, 
                                            from~to, value.var="kinship"))
  links_matrix <- apply(links_matrix, 2, as.numeric)
  rownames(links_matrix) <- colnames(links_matrix)
  
  links_plot_tmp <- links_tmp[links_tmp$kinship>relatedness_factor,]
  links_plot_2 <- links_plot_tmp[,c("from","kinship")]
  colnames(links_plot_2) <- c("label.node","kinship")
  links_plot_3 <- links_plot_tmp[,c("to","kinship")]
  colnames(links_plot_3) <- c("label.node","kinship")
  links_plot <- rbind(links_plot_2,links_plot_3)
  
  nodes <- data.frame(cbind(x$ind.names, as.character(pop(x))))
  colnames(nodes) <- c("name", "pop")
  
  network <- igraph::graph_from_data_frame(d = links,
                                           vertices = nodes,
                                           directed = FALSE)
  
  q <- relatedness_factor
  network.FS <-
    igraph::delete_edges(network, igraph::E(network)[links$kinship < q])
  
  if (method == "fr") {
    layout.name <- "Fruchterman-Reingold layout"
    plotcord <-
      data.frame(igraph::layout_with_fr(network.FS))
  }
  
  if (method == "kk") {
    layout.name <- "Kamada-Kawai layout"
    plotcord <-
      data.frame(igraph::layout_with_kk(network.FS))
  }
  
  if (method == "gh") {
    layout.name <- "Graphopt layout"
    plotcord <-
      data.frame(igraph::layout_with_graphopt(network.FS))
  }
  
  if (method == "mds") {
    layout.name <- "Multidimensional scaling layout"
    plotcord <-
      data.frame(igraph::layout_with_mds(network.FS))
  }
  
  # get edges, which are pairs of node IDs
  edgelist <- igraph::as_edgelist(network.FS, names = F)
  # convert to a four column edge data frame with source and destination
  # coordinates
  edges <- data.frame(plotcord[edgelist[, 1], ], plotcord[edgelist[, 2], ])
  # using kinship for the size of the edges
  edges$size <- links[links$kinship >= q, "kinship"]
  X1 <- X2 <- Y1 <- Y2 <- label.node <- NA
  colnames(edges) <- c("X1", "Y1", "X2", "Y2", "size")
  
  # node labels
  plotcord$label.node <- igraph::V(network.FS)$name
  
  # adding populations
  pop_df <-
    as.data.frame(cbind(indNames(x), as.character(pop(x))))
  colnames(pop_df) <- c("label.node", "pop")
  plotcord <- merge(plotcord, pop_df, by = "label.node")
  plotcord$pop <- as.factor(plotcord$pop)
  
  plotcord <- merge(plotcord,links_plot,by="label.node",all.x = TRUE)
  plotcord[is.na(plotcord$kinship),"kinship"] <- 0
  plotcord$kinship <- as.numeric(plotcord$kinship)
  plotcord$kinship <- scales::rescale(plotcord$kinship, to = c(0.1, 1))
  
  # assigning colors to populations
  if(is.null(palette_discrete)){
    # function to replicate defaults colors of ggplot
    discrete_palette <- gl.select.colors(x, 
                                         library = "brewer", 
                                         palette = "PuOr", 
                                         ncolors = nPop(x), 
                                         verbose = 0)
    
    palette_discrete <- discrete_palette
  }
  
  if (is(palette_discrete, "function")) {
    colors_pops <- palette_discrete(length(levels(pop(x))))
  }
  
  if (!is(palette_discrete, "function")) {
    colors_pops <- palette_discrete
  }
  
  if(is.null(link.color)){
    # taken from wes_palette::Zissou1
    diverging_palette <-
      colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))
    # creating convergent to 0 palette
    link.color <- diverging_palette
  }
  
  names(colors_pops) <- as.character(levels(x$pop))
  pal <- link.color(length(levels(x$pop)))
  size <- NULL
  new_levels <- unique(plotcord$pop[order(as.character(plotcord$pop))])
  # Rebuild the factor
  plotcord$pop <- factor(plotcord$pop, levels = new_levels)
  
  breaks <- round(seq(min(edges$size),
                      max(edges$size),
                      0.1),1)
  
  p1 <-
    ggplot() + 
    geom_segment(data = edges,
                 aes( x = X1, y = Y1, xend = X2,yend = Y2,color = size),
                 linewidth = link.size) +
    scale_colour_gradientn(
      name = "Relatedness",
      colours = pal,
      breaks = breaks,
      labels = as.character(breaks)
    ) +
    geom_point(data = plotcord,aes(x = X1,y = X2, fill = pop), 
               pch = 21,
               alpha=plotcord$kinship,
               size = node.size) +
    scale_fill_manual(name = legend.title, values = colors_pops)+
    coord_fixed(ratio = 1) + 
    theme_void() +
    ggtitle(paste(title, "\n[", layout.name, "]")) + 
    theme(legend.position = "bottom",
          plot.title = element_text( hjust = 0.5, face = "bold",size = 14)) 
  
  if (node.label == T) {
    p1 <-
      p1 + geom_text(
        data = plotcord,
        aes(
          x = X1,
          y = X2,
          label = label.node
        ),
        size = node.label.size,
        show.legend = FALSE,
        color = node.label.color
      )
  }
  
  # PRINTING OUTPUTS
  print(p1)
  
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose
    )
  }
  
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  
  return(invisible(list(p1, links_matrix)))
}
