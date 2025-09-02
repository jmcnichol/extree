#' Identify sequential nodes in a phylogenetic tree
#'
#' Internal function to calculate which internal nodes in a phylogenetic tree 
#' have child lineages that branch sequentially (i.e., differ in rank by 1).
#'
#' @param tree An object of class \code{phylo}.
#'
#' @return An integer vector of length \code{2 * nTip - 1}, where the first 
#' \code{ntip} entries (corresponding to tips) are 0, and the remaining 
#' entries correspond to nodes. A value of 1 indicates a sequential node, 0 otherwise.
#'
#' @keywords internal
#'

snodes = function(tree) { 
  ntip = length(tree$tip.label)
  nd.edges =tree$edge[ which(tree$edge[,2] > ntip),]
  heights = node.depth.edgelength(tree)[(ntip+1):(2*ntip-1)]
  ranks = order(order(heights)) 
  rank.version.of.nd.edges = 0*nd.edges
  rank.version.of.nd.edges[,1]=ranks[nd.edges[,1]-ntip]
  rank.version.of.nd.edges[,2]=ranks[nd.edges[,2]-ntip]
  diffs = rank.version.of.nd.edges[,2]-rank.version.of.nd.edges[,1]
  node.ss = ifelse(diffs ==1, 1, 0)
  svals = c(0*(1:ntip), node.ss) 
  return(svals) 
} 

