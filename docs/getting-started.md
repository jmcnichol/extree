---
layout: default
title: Getting started
---

# Getting started

[Home]({{ '/' | relative_url }}) · [Interpreting results]({{ '/interpretation.html' | relative_url }}) · [Reference]({{ '/reference.html' | relative_url }})

## Install the R package

From R, install the package from its GitHub repository (replace `jmcnichol` with the repository owner):

```r
install.packages("remotes")
remotes::install_github("jmcnichol/extree")
```

This command assumes `extree` has a standard installable R-package structure. If it is not installable with `remotes`, document the repository's actual installation command here.

## Load a tree

`extree` works with phylogenetic trees represented as `ape` objects. For example:

```r
library(ape)
library(extree)

tree <- ape::read.nexus("path/to/your_dated_tree.nex")

# A NEXUS file may contain multiple trees.
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]

class(tree)
ape::Ntip(tree)
ape::is.rooted(tree)
ape::is.ultrametric(tree)
```

For a Newick file, use `ape::read.tree("path/to/tree.newick")` instead.

## Check branch lengths

```r
stopifnot(
  inherits(tree, "phylo"),
  !is.null(tree$edge.length),
  all(is.finite(tree$edge.length)),
  all(tree$edge.length >= 0)
)
```

For an analysis based on event timing, use a dated tree and review the time calibration. If the source is a posterior distribution of trees, decide in advance whether to analyze one representative tree or repeat the test over trees.

## Run the test

```r
result <- extree::exch_test(tree)

result$statistic
result$pvalue
head(result$details)
```

The p-value object includes `lower`, `upper`, and `twosided` values. See [Interpreting results]({{ '/interpretation.html' | relative_url }}) for their meanings.

## Visualize results

```r
extree::exch_plot2(tree, result)
```

The plotting function has additional display options in versions used in analyses, including `showInfo`, `rootDate`, `breaks`, `dotsize`, and `showTipLabels`. Check the version installed on your machine for exact argument requirements:

```r
args(extree::exch_plot2)
```

## Next steps

- Try a small, well-dated tree first.
- Check both tails of the test rather than relying on the two-sided p-value alone.
- Record the tree source, dating method, sampling scheme, and package version alongside the result.
