---
layout: default
title: FAQ
---

# Frequently asked questions

[Home]({{ '/' | relative_url }}) · [Getting started]({{ '/getting-started.html' | relative_url }})

## Can I use a tree with no branch lengths?

Not for an analysis that depends on the timing of branching events. Use a time-calibrated phylogeny with meaningful branch lengths; do not create artificial lengths just to make the function run.

## What if my file contains multiple trees?

`ape::read.nexus()` can return a `multiPhylo` object. Select a tree explicitly or plan an analysis across the posterior sample:

```r
trees <- ape::read.nexus("posterior.trees")
if (inherits(trees, "multiPhylo")) {
  tree <- trees[[1]]
} else {
  tree <- trees
}
```

Do not treat a single selected tree as automatically representing phylogenetic uncertainty.

## What does a small lower-tail p-value mean?

It indicates that the observed statistic falls unusually low under the exchangeability null. A small upper-tail p-value indicates the opposite direction.

## Does a nonsignificant test prove exchangeability?

No. Failure to reject the null is not proof that the model is correct. Power can depend on the number of tips and the form and magnitude of the departure.

## Is there a way to locate the contributing nodes?

Yes. Inspect `result$details` and use the package's plotting functions:

```r
result$details
extree::exch_plot2(tree, result)
```

## Why does GitHub Pages look different from my local preview?

Confirm that `docs/_config.yml` has the correct theme and `baseurl`. For a project repository, the base URL is typically `/<repository-name>`, for example `/extree`. Internal documentation links in these pages use Jekyll's `relative_url` filter so they work under that prefix.
