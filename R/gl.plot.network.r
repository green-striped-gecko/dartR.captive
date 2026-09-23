#' @name gl.plot.network
#' @title Represents a similarity or distance matrix as a network
#' @description
#' This script takes a similarity matrix (such as the output of gl.grm) or a
#' distance matrix (such as the output of dist() or gl.dist.ind()) and
#' represents the relationship among the specimens as a network diagram. In
#' order to use this script, a decision is required on a threshold for
#' relatedness to be represented as link in the network, and on the layout
#' used to create the diagram.
#'
#' The threshold for relatedness to be represented as a link in the network is
#' specified as the proportion alpha of all pairs to draw: for a similarity
#' matrix the alpha proportion of pairs with the highest values, for a
#' distance matrix the alpha proportion with the lowest values. Often you are
#' looking for relatedness outliers in comparison with the overall relatedness
#' among individuals, so a very conservative value is used (e.g. 0.004), but
#' ultimately, this decision is made as a matter of trial and error. One way
#' to approach this trial and error is to try to achieve a sparse set of links
#' between unrelated 'background' individuals so that the stronger links are
#' preferentially shown.
#'
#' Link widths are proportional to the strength of the relationship (highest
#' similarity or lowest distance widest), rescaled to between 0.5 and 4.
#'
#' There are several layouts from which to choose. The most popular are given as
#' options in this script.
#' \itemize{
#' \item fr -- Fruchterman, T.M.J. and Reingold, E.M. (1991). Graph Drawing by
#' Force-directed Placement. Software -- Practice and Experience 21:1129-1164.
#' \item kk -- Kamada, T. and Kawai, S.: An Algorithm for Drawing General
#' Undirected Graphs. Information Processing Letters 31:7-15, 1989.
#' \item drl -- Martin, S., Brown, W.M., Klavans, R., Boyack, K.W., DrL:
#' Distributed Recursive (Graph) Layout. SAND Reports 2936:1-10, 2008.
#' }
#'
#' Colors of node symbols are those of the rainbow.
#'
#' @param D A similarity matrix (e.g. from gl.grm) or a distance matrix (e.g.
#' from dist() or gl.dist.ind()); set type accordingly. Row and column names
#' must be the individual names of x when x is given [required].
#' @param x A genlight object from which the D matrix was generated; used to
#' colour nodes by population [default NULL].
#' @param method One of "fr", "kk" or "drl" [default "fr"].
#' @param node.size Size of the symbols for the network nodes [default 3].
#' @param node.label TRUE to display node labels [default FALSE].
#' @param node.label.size Size of the node labels [default 0.7].
#' @param node.label.color Color of the text of the node labels
#' [default 'black'].
#' @param alpha Proportion of all pairs drawn as links: those with the highest
#' similarity, or the lowest distance [default 0.005].
#' @param title Title for the plot
#' [default "Network based on genetic distance"].
#' @param type Whether D holds similarities ("similarity", high values =
#' closely related, e.g. gl.grm) or distances ("distance", low values =
#' closely related, e.g. dist) [default "similarity"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return Invisibly, the network that was drawn (an igraph object with the
#' drawn links; the edge attribute weight holds the value from D).
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' if ((requireNamespace("rrBLUP", quietly = TRUE)) & (requireNamespace("gplots", quietly = TRUE))) {
#'   if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#'   test <- gl.subsample.loc(platypus.gl, n = 100)
#'   test <- gl.keep.ind(test, ind.list = indNames(test)[1:10])
#'   D <- gl.grm(test, legendx = 0.04)
#'   gl.plot.network(D, test)
#' }
#' @family captive management
#' @importFrom grDevices rgb
#' @importFrom graphics legend
#' @export

gl.plot.network <- function(D,
                            x = NULL,
                            method = "fr",
                            node.size = 3,
                            node.label = FALSE,
                            node.label.size = 0.7,
                            node.label.color = "black",
                            alpha = 0.005,
                            title = "Network based on genetic distance",
                            type = "similarity",
                            verbose = NULL) {
  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "igraph"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  # x is optional; without it nodes are drawn in one colour
  if (!is.null(x)) {
    datatype <- utils.check.datatype(x, verbose = verbose)
  } else if (verbose >= 2) {
    cat(warn(
      "  Note: genlight object not specified, population assignments not",
      "available for plotting\n"
    ))
  }

  # FUNCTION SPECIFIC ERROR CHECKING

  method <- match.arg(method, c("fr", "kk", "drl"))
  type <- match.arg(type, c("similarity", "distance"))

  if (!is(D, "dist") & !is(D, "matrix")) {
    stop(error(
      "  D must be a similarity or distance matrix (matrix or dist object)\n"
    ))
  }

  m <- as.matrix(D)
  if (nrow(m) != ncol(m) || is.null(rownames(m))) {
    stop(error("  D must be a square matrix with row and column names\n"))
  }
  if (!is.null(x) && !setequal(rownames(m), indNames(x))) {
    stop(error(
      "  The row and column names of D must be the individual names of x",
      "(indNames(x))\n"
    ))
  }

  # DO THE JOB

  # one row per pair of individuals
  pairs.idx <- which(upper.tri(m), arr.ind = TRUE)
  links <- data.frame(
    from = rownames(m)[pairs.idx[, "row"]],
    to = colnames(m)[pairs.idx[, "col"]],
    weight = m[pairs.idx],
    stringsAsFactors = FALSE
  )

  if (!is.null(x)) {
    nodes <- data.frame(name = indNames(x), pop = as.character(pop(x)))
  } else {
    nodes <- data.frame(name = rownames(m))
  }

  network <-
    igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )

  if (!is.null(x)) {
    colors <- rainbow(nlevels(pop(x)))
    my_colors <- colors[pop(x)]
  } else {
    my_colors <- "red"
  }

  # keep the alpha proportion of the closest pairs: the highest similarities
  # or the lowest distances
  if (type == "similarity") {
    strength <- links$weight
  } else {
    strength <- -links$weight
  }
  q <- stats::quantile(strength, p = 1 - alpha, na.rm = TRUE)
  keep <- !is.na(strength) & strength >= q
  network.FS <- igraph::delete_edges(network, which(!keep))

  # widths from the drawn links' own strength, strongest widest
  kept.strength <- strength[keep]
  if (length(kept.strength) > 0 &&
      diff(range(kept.strength)) > 0) {
    edge.widths <- 0.5 + 3.5 * (kept.strength - min(kept.strength)) /
      diff(range(kept.strength))
  } else {
    edge.widths <- 2
  }

  if (method == "fr") {
    layout.name <- "Fruchterman-Reingold layout"
    l <- igraph::layout_with_fr(network.FS)
  }
  if (method == "kk") {
    layout.name <- "Kamada-Kawai layout"
    l <- igraph::layout_with_kk(network.FS)
  }
  if (method == "drl") {
    layout.name <- "DrL Graph layout"
    l <- igraph::layout_with_drl(network.FS)
  }
  title <- paste(title, "\n[", layout.name, "]")

  if (node.label) {
    node.label <- igraph::V(network)$name
  } else {
    node.label <- NA
    node.label.size <- NA
    node.label.color <- NA
  }

  plot(
    network.FS,
    edge.arrow.size = 0,
    edge.curved = 0,
    edge.width = edge.widths,
    vertex.size = node.size,
    vertex.color = my_colors,
    vertex.frame.color = "#555555",
    vertex.label = node.label,
    vertex.label.color = node.label.color,
    vertex.label.cex = node.label.size,
    layout = l,
    main = title
  )

  if (!is.null(x)) {
    legend(
      "bottomleft",
      legend = levels(pop(x)),
      col = colors,
      bty = "n",
      pch = 20,
      pt.cex = 3,
      cex = 1,
      text.col = colors,
      horiz = FALSE,
      inset = c(0.1, 0.1)
    )
  }

  if (verbose >= 3) {
    cat(report("  Links drawn:", igraph::ecount(network.FS), "of",
               nrow(links), "pairs\n"))
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(invisible(network.FS))
}
