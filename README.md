# phonJSD

**phonJSD** is an R package for measuring phonological category separation using **Jensen–Shannon Divergence (JSD)**.  
It is designed for researchers working in sociophonetics, laboratory phonology, bilingualism, and speech perception who need a principled, distributional metric of category overlap in acoustic space.

Version **0.5.1** is a research release focused on stable core metrics, reproducible uncertainty estimates, and comparison with classical overlap measures.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20816586.svg)](https://doi.org/10.5281/zenodo.20816586)

---

## Motivation

Traditional measures of phonological contrast (e.g., Euclidean distance between means, Pillai scores) capture only limited aspects of category structure.  
**phonJSD** treats phonological categories as **probability distributions** in acoustic space and quantifies their separation using Jensen–Shannon Divergence—a symmetric, bounded, and interpretable information-theoretic metric.

This approach is especially useful when:
- Categories differ in **shape, variance, or multimodality**
- Token counts are unbalanced
- Category overlap is gradient rather than categorical
- You want a measure that generalizes naturally to **multidimensional features** (e.g., MFCCs)

---

## Core Features (v0.5.1)

- Jensen–Shannon Divergence for phonological category comparison  
- Kernel density–based estimation of acoustic distributions  
- Support for **1D and n-dimensional acoustic features**
- Global and group-level bootstrap summaries
- Comparison metrics including Pillai-Bartlett trace, Bhattacharyya distance/affinity, Mahalanobis distance, and percent overlap
- Reproducible pipelines compatible with tidyverse workflows
- Designed for integration with forced alignment and acoustic extraction tools

---

## Quick Start

Most users should start with `compare_overlap_metrics()`. It takes a token-level
vowel table, compares two category labels in a numeric acoustic space, and
returns the major overlap/separation metrics side by side.

```r
library(phonJSD)

set.seed(2026)
vowels <- data.frame(
  speaker = rep(c("s01", "s02"), each = 80),
  vowel = rep(rep(c("ih", "eh"), each = 40), 2),
  f1 = c(
    rnorm(40, 500, 55), rnorm(40, 560, 60),
    rnorm(40, 510, 60), rnorm(40, 575, 65)
  ),
  f2 = c(
    rnorm(40, 1980, 150), rnorm(40, 1880, 155),
    rnorm(40, 1960, 160), rnorm(40, 1840, 165)
  )
)

compare_overlap_metrics(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker",
  output = "wide"
)
```

For plotting or ranking, use the long output:

```r
metrics_long <- compare_overlap_metrics(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker",
  output = "long"
)

metrics_long[, c("group", "metric", "estimate", "orientation",
                 "separation_value", "separation_rank")]
```

If `ggplot2` is installed, the same workflow can be visualized directly:

```r
plot_category_space(
  data = vowels,
  features = c("f2", "f1"),
  category_col = "vowel",
  group_col = "speaker",
  reverse_x = TRUE,
  reverse_y = TRUE
)

plot_overlap_metrics(metrics_long)
```

If you only need the package's main information-theoretic estimate, use
`estimate_jsd()`:

```r
estimate_jsd(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker",
  do_boot = TRUE,
  n_boot = 100
)
```

Use lower-level helpers such as `jsd_kde_nd()`, `percent_overlap_kde()`,
`pillai_overlap()`, and `bhattacharyya_mvnorm()` when you are validating a
method, debugging one contrast, or need direct control over one metric.

---

## Conceptual Overview

Given two phonological categories (e.g., vowels /ɪ/ and /ɛ/), phonJSD:

1. Represents each category as a probability distribution over acoustic space  
2. Estimates densities using kernel density estimation (KDE)  
3. Computes Jensen–Shannon Divergence between the distributions  

KDE is the package's default density-estimation strategy for continuous
acoustic spaces. JSD itself, however, operates on probability distributions and
is not intrinsically tied to KDE. Manuscript sensitivity analyses show that the
same substantive contrast pattern is recovered when JSD is computed from KDE,
histogram-based, empirical-binned, and Gaussian-mixture distribution estimates.

JSD values:
- **0** → complete overlap (no separation)
- **Higher values** → greater distributional separation  
- **Bounded and symmetric**, making cross-study comparisons more interpretable

## Interpreting Metrics

The metrics are not all oriented in the same direction:

| Metric | Direction | Notes |
| --- | --- | --- |
| JSD | Higher = more separation | Bounded from 0 to 1 in bits |
| Jensen-Shannon distance | Higher = more separation | `sqrt(JSD)` |
| Pillai trace | Higher = more separation | Classical MANOVA-based comparison |
| Bhattacharyya distance | Higher = more separation | Assumes multivariate normality |
| Mahalanobis distance | Higher = more separation | Mean separation relative to pooled covariance |
| Percent overlap | Higher = more overlap | 1 means near-complete shared density |
| Bhattacharyya affinity | Higher = more overlap | `exp(-Bhattacharyya distance)` |

For side-by-side comparison, `compare_overlap_metrics(output = "long")` adds
`orientation`, `separation_value`, and `separation_rank` columns so overlap
metrics can be read on the same separation-oriented scale as distance metrics.

---

## Installation

This package is currently in early development and not yet on CRAN.

```r
# install.packages("remotes")
remotes::install_github("berrygrant/phonJSD")
```

## LabPhon 2026 Poster

The final LabPhon 2026 poster and reproducibility bundle are maintained on the
`labphon_2026` branch so the main package branch can stay lightweight.

- [Download the final poster PDF](https://raw.githubusercontent.com/berrygrant/phonJSD/labphon_2026/analysis/labphon_2026/poster/LabPhon2026_Poster_Berry_Final.pdf)
- [Browse the poster reproducibility bundle](https://github.com/berrygrant/phonJSD/tree/labphon_2026/analysis/labphon_2026)
