#' Plot phylogenetic tree with nodes whose child is "next to branch"
#'
#' Plots a tree, using `ggtree`, where nodes that contribute to the
#' value of the exchnageability test statistc (i.e., their child
#' node is chronolocially next to branch) are labelled with a dot.
#' To add to the visualization the nodes are coloured in according to
#' their value of the Shannon entropy, \eqn{\log_2(1/p)}.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @param result The output of \code{exch_test}.
#' @param showInfo Whether to show the information content, -p log(p), for each node. Defaults to FALSE.
#' @param rootDate For a plot in calendar time, add the date of the root in Date format. Defaults to numeric format with value 0.  
#' @param breaks The breaks parameter for the x axis  
#' @param dotsize The size of the dots
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
exch_plot <- function(tree, result,showInfo=FALSE, rootDate=0,breaks=6,dotsize=3) {
  res <- result$details
  colnames(res)[colnames(res) == "nodeids"] <- "node"
  res$entropy <- log2(1 / res$probs)
  res <- res[res$s == 1, ]
  nodedists <- node.depth.edgelength(tree) 
  distdf <- data.frame(
    node = seq_along(nodedists),
    disttoroot = nodedists)
  res <- merge(res, distdf, by = "node", all.x = TRUE)
  if (showInfo) { 
   if (class(rootDate)=="Date") {
     ggtree(tree,
             color = "grey",
             layout = "rectangular") %<+% res +
          geom_nodepoint(aes(color = entropy,
                             subset = s ), size=dotsize) +
          scale_colour_viridis_c(option = "plasma",
                                 direction = 1,
                                 name = "Information") +
          theme_tree2() +
          expand_limits(y = -5) +
          scale_x_continuous(name = " ",
                             breaks = scales::pretty_breaks(10),
                             labels = function(x) format(rootDate + 365.25*x, "%Y-%m")) +
          guides(size = "none")
   } else { 
     ggtree(tree,
            color = "grey",
            layout = "rectangular") %<+% res +
       geom_nodepoint(aes(color = entropy,
                          subset = s ), size=dotsize) +
       scale_colour_viridis_c(option = "plasma",
                              direction = 1,
                              name = "Information") +
       theme_tree2() +
       expand_limits(y = -5) +
       scale_x_continuous(name = " ",
                          breaks = scales::pretty_breaks(10),
                          labels = function(x) (rootDate + x)) +
       guides(size = "none")
     
   }
      } 
  if (!showInfo) { 
    if (class(rootDate)=="Date") {
      ggtree(tree,
             color = "grey",
             layout = "rectangular") %<+% res +
          geom_nodepoint(aes(subset = s),color="black",size = dotsize) +
                  theme_tree2() +
          expand_limits(y = -5) +
          scale_x_continuous(name = " ",
                             breaks = scales::pretty_breaks(5),
                             labels = function(x) format(rootDate + 365.25*x, "%Y-%m")) +
          guides(size = "none")
      } else {
        ggtree(tree,
               color = "grey",
               layout = "rectangular") %<+% res +
          geom_nodepoint(aes(subset = s),color="black",size = dotsize) +
          theme_tree2() +
          expand_limits(y = -5) +
          scale_x_continuous(name = " ",
                             breaks = scales::pretty_breaks(5),
                             labels = function(x) (rootDate + x)) +
          guides(size = "none")
        
      }
  }
  }

