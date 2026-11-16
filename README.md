# extree

**Testing exchangeability in time-calibrated phylogenetic trees**

`extree` is an R package for assessing whether the branching patterns in a phylogenetic tree are consistent with **exchangeability of lineages**.

Exchangeability is an assumption underlying many models of phylogenetic tree generation: conditional on the relevant state of the process, contemporaneous lineages are treated symmetrically. The package uses the distribution of branching events across lineages to test for departures from that assumption.

The goal is to provide a diagnostic for the assumptions of phylogenetic models, rather than to assign a particular evolutionary mechanism to a significant result.

## Features

- Test exchangeability in rooted, time-calibrated phylogenetic trees.
- Summarize the branching-event pattern with a test statistic.
- Obtain lower-tail, upper-tail and two-sided p-values.
- Inspect node-level contributions to the analysis.
- Visualize branching events on a phylogeny.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("OWNER/extree")
```

Replace `OWNER` with the GitHub account or organization hosting this repository.

If you have already cloned the repository, you can install the package from the repository root:

```r
remotes::install_local(".")
```

## Quick start

The following example generates a tree using `ape` and runs the exchangeability test.

```r
library(ape)
library(extree)

set.seed(123)
tree <- ape::rcoal(50)

result <- exch_test(tree)

# Test statistic
result$statistic

# One-sided and two-sided p-values
result$pvalue

# Node-level details
head(result$details)
```

To analyze your own tree, load a rooted, time-calibrated phylogeny with meaningful branch lengths:

```r
tree <- ape::read.nexus("my_dated_tree.nex")
result <- exch_test(tree)
result$pvalue
```

If your file contains multiple trees, select the tree you intend to analyze before calling `exch_test()`.

## Interpreting the results

`exch_test()` returns a list with three main components:

| Component | Description |
| --- | --- |
| `statistic` | Observed test statistic describing the branching-event pattern. |
| `pvalue` | Data frame containing `lower`, `upper`, and `twosided` p-values. |
| `details` | Node-level information used to inspect the test calculation. |

A small **lower-tail** or **upper-tail** p-value indicates that the observed branching-event pattern lies in the corresponding tail of the reference distribution under exchangeability. The **two-sided** p-value assesses departures in either direction.

A significant result is evidence against the exchangeability assumption for the tree and test under consideration. It does **not**, by itself, identify the cause of the departure. Sampling, tree reconstruction, and evolutionary or epidemiological processes may all affect interpretation.

## Documentation

Additional documentation is available in [`docs/`](docs/index.md):

- [Getting started](docs/getting-started.md)
- [Interpreting results](docs/interpretation.md)
- [Function reference](docs/reference.md)
- [Frequently asked questions](docs/faq.md)

## Citation

A formal citation for the package or associated methods paper will be added when available. Until then, please reference the repository and the version or commit used in your analysis.

## License

Please consult the repository's license file, if present, before reusing or redistributing the code.
