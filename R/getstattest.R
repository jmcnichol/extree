#' Test statistic
#'
#' Internal function that evaluates whether each node in a phylogenetic tree
#' has a child lineage that branches "next". Computes node-level summaries
#' and an overall tree-level test statistic.
#'
#' @param tree An object of class \code{phylo}.
#'
#' @return A list with two elements:
#' \enumerate{
#'   \item A data frame with the following columns:
#'     \describe{
#'       \item{\code{node}}{The node IDs (internal nodes).}
#'       \item{\code{s}}{Indicator (1/0) for whether the node’s child lineage branched sequentially.}
#'       \item{\code{nlins}}{The number of lineages immediately after the node.}
#'       \item{\code{probs}}{The probability (\eqn{2/k}) that one of the node's child lineages would be next.}
#'       \item{\code{depth}}{The topological depth of the node.}
#'       \item{\code{time}}{The branching time of the node.}
#'     }
#'   \item A numeric value: the test statistic \eqn{W} for the tree. Under the null,
#'         \eqn{W \sim N(0,1)}.
#' }
#'
#' @keywords internal
#'
getstattest <- function(tree, usePoly = FALSE) {
  if (!inherits(tree, "phylo")) {
    tree <- as(tree, "phylo")
  }
    if (usePoly) {
        lttb <- lttpoly(tree)
    } else {
        lttb <- lttbelow(tree)
    }
    IND <- which(lttb$nodeids > length(tree$tip.label))
    intonly <- lttb$nodeids[IND]
    nlins <- lttb$lttbelow[IND]
    nexttobranch <- c(intonly[-1], 0)
    if (usePoly) {
    declist <- lapply(intonly, function(x) {
    tree$edge[which(tree$edge[, 1] == x), 2]
    })

        names(declist) <- as.character(intonly)
        numdesc <- vapply(intonly, function(x) {
            length(declist[[as.character(x)]])
        }, FUN.VALUE = 1)
        sresult <- vapply(1:length(intonly), function(x) {
            is.element(nexttobranch[x], declist[[as.character(intonly[x])]])
        }, FUN.VALUE = 0.5)
    } else {
        declist <- vapply(intonly, function(x) {
            tree$edge[which(tree$edge[, 1] == x), 2]
        }, FUN.VALUE = c(1, 1))
        colnames(declist) <- as.character(intonly)
        numdesc <- 2 + 0 * (1:length(intonly))
        sresult <- vapply(1:length(intonly), function(x) {
            is.element(nexttobranch[x], declist[, as.character(intonly[x])])
        }, FUN.VALUE = 0.5)
    }
    probs <- numdesc / nlins
    probs[length(probs)] <- 0
    W <- sum(sresult)
    L <- length(probs)
    return(list(details = data.frame(
        nodeids = intonly, s = sresult, nlins = nlins,
        probs = probs, topdepths = lttb$topdepths[IND],
        times = lttb$times[IND]
    ), W = W))
}
