---
layout: default
title: extree
---

# extree

**Tools for testing exchangeability in time-calibrated phylogenetic trees.**

Phylogenetic models often assume that contemporaneous lineages are exchangeable: at a given time, the labels of lineages do not determine which lineage experiences the next branching event. The `extree` R package provides a way to investigate this assumption using information contained in the branching structure and timing of an observed phylogeny.

The principal function, `exch_test()`, summarizes the observed branching-event pattern and returns lower-tail, upper-tail, and two-sided p-values. The associated plotting functions help visualize where relevant branching events occur on the tree.

## Quick example

```r
library(ape)
library(extree)

tree <- read.nexus("my_dated_tree.nex")
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]

result <- exch_test(tree)
result$statistic
result$pvalue

# Inspect the branching events underlying the test
head(result$details)
```

Use a rooted, time-calibrated tree with meaningful branch lengths. Do not assign artificial branch lengths to a topology-only tree.

## Documentation

- [Getting started]({{ '/getting-started.html' | relative_url }}) — installation, loading a tree, and running your first test.
- [Interpreting results]({{ '/interpretation.html' | relative_url }}) — the statistic, p-values, and what departures from exchangeability mean.
- [Function reference]({{ '/reference.html' | relative_url }}) — the core functions and returned objects.
- [Frequently asked questions]({{ '/faq.html' | relative_url }}) — common input and plotting issues.

## Why test exchangeability?

Exchangeability underlies many commonly used stochastic models for phylogenetic trees, including several coalescent and branching-process models. A departure from exchangeability may point to lineage-specific differences in the processes generating the tree, but the test alone does not identify a unique biological mechanism. Interpret results in light of tree reconstruction, sampling, and other model assumptions.

## Citation

If you use `extree` in published work, cite the associated software or methods paper when a formal citation is available. **[Add the definitive package or paper citation here.]**
