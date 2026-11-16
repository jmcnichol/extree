---
layout: default
title: Function reference
---

# Function reference

[Home]({{ '/' | relative_url }}) · [Getting started]({{ '/getting-started.html' | relative_url }}) · [Interpreting results]({{ '/interpretation.html' | relative_url }})

This page documents functions and output fields used in current analyses. Consult the installed function signatures for the exact argument names in your package version.

## `exch_test()`

Runs the exchangeability test on a phylogenetic tree.

```r
result <- extree::exch_test(tree)
```

**Input:** An `ape` `phylo` tree with meaningful branch lengths and timing.

**Output:** A list with the following components:

| Component | Content |
| --- | --- |
| `statistic` | Observed test statistic W. |
| `pvalue` | Data frame with `lower`, `upper`, and `twosided` p-values. |
| `details` | Node-level data including `nodeids`, `s`, `nlins`, `probs`, `topdepths`, and `times`. |

```r
names(result)
str(result, max.level = 2)
```

## `exch_plot2()`

Visualizes the tree and information associated with a fitted exchangeability test.

```r
extree::exch_plot2(tree, result)
args(extree::exch_plot2)
```

Depending on package version, display controls may include `showInfo`, `rootDate`, `breaks`, `dotsize`, and `showTipLabels`.

## Inspect the installed package

```r
packageVersion("extree")
getNamespaceExports("extree")
help(package = "extree")
```

**To complete this reference:** Replace these brief summaries with function signatures, parameter definitions, return values, and worked examples drawn directly from the current R source or roxygen documentation.
