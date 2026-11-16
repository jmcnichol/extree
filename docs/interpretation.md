---
layout: default
title: Interpreting results
---

# Interpreting results

[Home]({{ '/' | relative_url }}) · [Getting started]({{ '/getting-started.html' | relative_url }}) · [Reference]({{ '/reference.html' | relative_url }})

The exchangeability test is based on an observed branching-event statistic, usually denoted **W**. The result includes a test statistic, tail probabilities, and information about the contributing nodes.

## The test output

```r
result <- extree::exch_test(tree)

result$statistic          # observed statistic W
result$pvalue$lower       # lower-tail p-value
result$pvalue$upper       # upper-tail p-value
result$pvalue$twosided    # two-sided p-value
```

| Quantity | Interpretation |
| --- | --- |
| `statistic` | Observed test statistic W. |
| `pvalue$lower` | Evidence for fewer relevant branching events than expected under the exchangeability null. |
| `pvalue$upper` | Evidence for more relevant branching events than expected under the exchangeability null. |
| `pvalue$twosided` | Evidence for a departure in either direction. |

A small p-value indicates a discrepancy with the null model, **not** proof of a particular evolutionary mechanism. The two tails distinguish different directions of deviation.

## Node-level information

The `details` data frame can include:

| Column | Description |
| --- | --- |
| `nodeids` | Tree node identifiers. |
| `s` | Observed event indicator, typically 0 or 1. |
| `probs` | Null-model probability associated with that event. |
| `nlins` | Lineage-count information. |
| `topdepths` | Associated node-depth information. |
| `times` | Associated times. |

You can inspect the observed event count and its null expectation directly:

```r
observed <- sum(result$details$s == 1, na.rm = TRUE)
expected <- sum(result$details$probs, na.rm = TRUE)

c(observed = observed, expected = expected)
```

## How to report a result

> We applied the exchangeability test to a dated phylogeny and obtained an observed statistic of W = **[value]**. The lower-tail, upper-tail, and two-sided p-values were **[values]**, respectively. These results **[did / did not]** provide evidence against exchangeability, with the direction of deviation indicating **[fewer / more]** observed branching events than predicted under the null.

For biological interpretation, also report the number of tips, tree source, whether you analyzed a posterior sample or a single tree, and any relevant subsampling or sensitivity analyses.

## Important cautions

- A significant result does not by itself distinguish diversification-rate heterogeneity, selective sampling, tree-estimation error, or other causes.
- Results from a subclade answer a different question than results from the full phylogeny.
- Treat p-values from many exploratory clade tests as multiple comparisons, where appropriate.
