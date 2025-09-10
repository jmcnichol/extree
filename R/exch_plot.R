#' Plot Exchangeability Test Results on a Phylogenetic Tree
#'
#' Visualizes the results of an exchangeability test by mapping node-level
#' entropy values onto a phylogenetic tree using `ggtree`.
#' Entropy is computed as \eqn{\log_2(1/p)}
#' for each node's branching probability,
#' with higher values indicating greater uncertainty.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @param result The output of \code{exch_test}.
#'
#' @return A `ggplot` object displaying the tree with nodes that
#'    contribute to the test statistic coloured by entropy values.
#'
#' @details
#' - Only nodes with `s == 1` are plotted (nodes whose child node is
#' chronologically next to branch).
#' - Colors are scaled using the "plasma" palette from `viridis`.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' library(ggtree)
#' tree <- rtree(10)
#' result <- exch_test(tree)
#' exch_plot(tree, result)
#' }
#'
#' @importFrom ggtree ggtree %<+% geom_nodepoint theme_tree2
#' @importFrom ggplot2 scale_colour_viridis_c expand_limits
#' scale_x_continuous aes guides
#' @importFrom scales pretty_breaks
#' @export
exch_plot <- function(tree, result) {
  res <- result$details
  colnames(res)[colnames(res) == "nodeids"] <- "node"
  res$entropy <- log2(1 / res$probs)
  res <- res[res$s == 1, ]
  nodedists <- node.depth.edgelength(tree) 
  distdf <- data.frame(
    node = seq_along(nodedists),
    disttoroot = nodedists
  )
  res <- merge(res, distdf, by = "nodeids", all.x = TRUE)
  ggtree(tree,
         color = "grey",
         layout = "rectangular") %<+% res +
    geom_nodepoint(aes(color = res$entropy,
                       subset = res$s,
                       size = 3)) +
    scale_colour_viridis_c(option = "plasma",
                           direction = 1,
                           name = "Entropy") +
theme_tree2() +
expand_limits(y = -5) +
scale_x_continuous(name = " ",
                   breaks = scales::pretty_breaks(10)) +
guides(size = "none")

}