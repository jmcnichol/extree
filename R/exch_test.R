#' Test for exchangeability
#'
#' Carries out a test for exchangeability on a timed
#' phylogenetic tree. The test is based on the
#' Poisson-binomial distribution.
#'
#' @param tree An object of class \code{phylo}.
#' @param poly If the tree has polytomies.
#'
#' @return A list with three elements:
#' \enumerate{
#'   \item \code{statistic} the test statistic for the tree.
#'   \item \code{pvalue} the corresponding p-value for the test.
#'   \item A data frame with the following columns:
#'     \describe{
#'       \item{\code{node}}{The node IDs (internal nodes).}
#'       \item{\code{s}}{Indicator (1/0) for whether the node’s child
#'            lineage branched sequentially.}
#'       \item{\code{nlins}}{The number of lineages
#'             immediately after the node.}
#'       \item{\code{probs}}{The probability \eqn{\frac{m}{k}} that one
#'             of the node's child lineages would be next.}
#'       \item{\code{depth}}{The topological depth of the node.}
#'       \item{\code{time}}{The branching time of the node.}
#'     }
#' }
#'
#' @examples
#' library(ape)
#'
#' # Simulate a random phylogenetic tree with 10 tips
#' set.seed(123)
#' tree <- rcoal(10)
#'
#' # Run the exchangeability test (upper tail by default)
#' result <- exch_test(tree)
#'
#' # Print results
#' print(result)
#'
#' # Access components
#' result$statistic
#' result$pvalue
#' head(result$details)
#'
#' # Tree with a polytomies
#' newick <- "((t1:2,t2:2,t3:2,t4:2):4,
#'            (((t5:1,t6:1):3,(t7:2,(t8:1,t9:1):1):2):2,
#'              ((t10:1,t11:1):3,((t12:1,t13:1):2,(t14:1,t15:1):2):2):2):2);"
#' tree <- read.tree(text=newick)
#' # Run the exchangeability test with poly=TRUE
#' result <- exch_test(tree, poly=TRUE)
#'
#' # Print results
#' print(result)
#'
#' # Access components
#' result$statistic
#' result$pvalue
#' head(result$details)
#'
#' @export

exch_test <- function(tree, poly = FALSE) {
  if (!inherits(tree, "phylo"))
    stop("Tree should be an object of class \"phylo\".")
  if (!inherits(tree, "phylo")) {
    tree <- as(tree, "phylo")
  }
  # add a check that it is timed
  if (poly) {
    lttb <- lttpoly(tree)
  } else {
    lttb <- lttbelow(tree)
  }
  ind <- which(lttb$nodeids > length(tree$tip.label))
  intonly <- lttb$nodeids[ind]
  nlins <- lttb$lttbelow[ind]
  nexttobranch <- c(intonly[-1], 0)
  if (poly) {
    declist <- lapply(intonly, function(x) {
      tree$edge[which(tree$edge[, 1] == x), 2]
    })
    names(declist) <- as.character(intonly)
    numdesc <- vapply(intonly, function(x) {
      length(declist[[as.character(x)]])
    }, 1)
    sresult <- vapply(seq_along(intonly), function(x) {
      is.element(nexttobranch[x], declist[[as.character(intonly[x])]])
    }, 0.5)
    edges <- data.frame(tree$edge)
    counts <- as.data.frame(table(edges$X1))
    colnames(counts) <- c("nodeids", "count")
    counts$countpoly <- ifelse(counts$count > 3, 1, 0)
    edges <- counts[, c("nodeids", "countpoly")]
    sresult <- ifelse(sresult | edges$countpoly, 1, 0)
  } else {
    declist <- vapply(intonly, function(x) {
      tree$edge[which(tree$edge[, 1] == x), 2]
    }, c(1, 1))
    colnames(declist) <- as.character(intonly)
    numdesc <- 2 + 0 * seq_along(intonly)
    sresult <- vapply(seq_along(intonly), function(x) {
      is.element(nexttobranch[x], declist[, as.character(intonly[x])])
    }, 0.5)
  }
  probs <- numdesc / nlins
  probs[length(probs)] <- 0
  W <- sum(sresult)
  # prep returns
  details <- data.frame(nodeids = intonly,
                        s = sresult,
                        nlins = nlins,
                        probs = probs,
                        topdepths = lttb$topdepths[ind],
                        times = lttb$times[ind])
  # carry out test
  lowertailp <- poibin::ppoibin(W, pp = probs)
  uppertailp <- 1 - poibin::ppoibin(W - 1, pp = probs)
  twotailp <- min(1, 2 * min(lowertailp, uppertailp))
  pval <- data.frame(upper = uppertailp,
                     twosided = twotailp,
                     lower = lowertailp)
  return(list(statistic = W, pvalue = pval, details = details))
}