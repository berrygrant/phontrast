# phonJSD

**phonJSD** is an R package for measuring phonological category separation using **Jensen–Shannon Divergence (JSD)**.  
It is designed for researchers working in sociophonetics, laboratory phonology, bilingualism, and speech perception who need a principled, distributional metric of category overlap in acoustic space.

Version **1.0.0** was the first stable release of phonJSD, focused on core overlap/separation metrics, reproducible uncertainty estimates, visualization, and comparison with classical overlap measures. Version **1.2.0** replaces the KDE-based JSD and percent-overlap estimator with a consistent Monte-Carlo plug-in (now the default; results differ from 1.0.0, and `method = "legacy"` reproduces them) and adds opt-in high-dimensional KDE speed controls.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21465469.svg)](https://doi.org/10.5281/zenodo.21465469)

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

## Core Features

- Jensen–Shannon Divergence for phonological category comparison  
- Kernel density–based estimation of acoustic distributions  
- Support for **1D and n-dimensional acoustic features**
- Optional high-dimensional KDE speed controls, including diagonal Scott
  bandwidths, sampled evaluation points, and a fast diagonal-Gaussian engine
- Global and group-level bootstrap summaries
- Comparison metrics including Pillai-Bartlett trace, Bhattacharyya distance/affinity, Mahalanobis distance, and percent overlap
- ggplot2-backed visualizations for metric tables, category spaces, and PCA projections
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

For multidimensional features such as MFCCs or acoustic embeddings, compute
metrics on the full feature set and use PCA plots as a diagnostic projection:

```r
features <- paste0("mfcc", 1:13)

metrics_13d <- compare_overlap_metrics(
  data = vowel_tokens,
  features = features,
  category_col = "vowel",
  output = "long"
)

plot_category_pca(
  data = vowel_tokens,
  features = features,
  category_col = "vowel"
)
```

`plot_category_pca()` helps visualize high-dimensional structure, but it does
not change the estimand: the metric table above is still estimated in all 13
dimensions.

For larger high-dimensional datasets, the default KDE path is conservative:
it uses `ks::Hpi()` and evaluates on all pooled observations. You can opt into
a faster nonparametric KDE path by using a diagonal Scott bandwidth, a sampled
set of evaluation points, and the chunked diagonal evaluator:

```r
metrics_13d_fast <- compare_overlap_metrics(
  data = vowel_tokens,
  features = paste0("mfcc", 1:13),
  category_col = "vowel",
  bw = "scott.diag",
  eval_on = "pooled_sample",
  eval_n = 200,
  eval_seed = 2026,
  engine = "fast_diag",
  output = "long"
)
```

This remains KDE: one kernel is still centered on each observed token. The
speedup comes from a cheaper diagonal bandwidth rule, evaluating densities on a
bounded sampled support, and using a vectorized diagonal-Gaussian evaluator.
The package default is unchanged for backward compatibility.

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

By default (`method = "mc"`) the divergence is estimated with a Monte-Carlo
plug-in: each category's KDE is evaluated at that category's own tokens and the
log density ratio against the mixture is averaged, giving a consistent estimate
of the continuous JSD in any number of dimensions. The pre-1.2.0 self-normalized
sample-point estimate — a bounded relative separation index rather than the JSD
integral — remains available as `method = "legacy"` for reproducing 1.0.0
results.

JSD values:
- **0** → complete overlap (no separation)
- **Higher values** → greater distributional separation  
- **Bounded and symmetric**, making cross-study comparisons more interpretable

## Choosing and Interpreting Metrics

The metrics are not all oriented in the same direction:

| Metric | Direction | Scale | Best used when | Main caveat |
| --- | --- | --- | --- | --- |
| JSD | Higher = more separation | 0-1 bits | You want a bounded distributional separation metric | KDE estimates need adequate sample size per category |
| Jensen-Shannon distance | Higher = more separation | 0-1 | You want a distance transform of JSD | Same estimator caveats as JSD |
| Pillai trace | Higher = more separation | 0-1 | You want a classical MANOVA-style comparison | Primarily mean-based; less sensitive to distribution shape |
| Bhattacharyya distance | Higher = more separation | 0 to infinity | You want a parametric distribution-distance comparison | Assumes approximately multivariate normal categories |
| Mahalanobis distance | Higher = more separation | 0 to infinity | You want mean separation scaled by covariance | Sensitive to covariance estimation and small samples |
| Percent overlap | Higher = more overlap | 0-1 proportion | You want a directly interpretable shared-density estimate | Despite the name, output is a proportion, not 0-100 |
| Bhattacharyya affinity | Higher = more overlap | 0-1 | You want a parametric overlap analogue | Assumes approximately multivariate normal categories |

For side-by-side comparison, `compare_overlap_metrics(output = "long")` adds
`orientation`, `separation_value`, and `separation_rank` columns so overlap
metrics can be read on the same separation-oriented scale as distance metrics.
Confidence interval columns are named `ci_lower` and `ci_upper`; legacy
JSD-specific aliases `jsd_low` and `jsd_high` are retained for compatibility.

## PB52 Example Note

The Peterson and Barney 1952 data in `phonTools::pb52` are useful for global
vowel contrasts. After filtering the factor-valued `vowel` column, phonJSD
counts observed levels only, so manual `droplevels()` is not required.

```r
data(pb52, package = "phonTools")
pb_i <- subset(pb52, as.character(vowel) %in% c("I", "i"))

compare_overlap_metrics(
  data = pb_i,
  features = c("f1", "f2"),
  category_col = "vowel"
)
```

Do not use per-speaker F1/F2 `I/i` overlap metrics in PB52 without a coarser
grouping or different design: each speaker has only two repetitions per vowel,
which is too few for KDE-based two-dimensional overlap estimates.

---

## Installation

This package is not yet on CRAN. Install the latest tagged release (v1.2.0) with:

```r
# install.packages("remotes")
remotes::install_github("berrygrant/phonJSD@v1.2.0")
```

For the current development version, use:

```r
remotes::install_github("berrygrant/phonJSD")
```

## AI Use Disclosure

OpenAI ChatGPT/Codex and Anthropic Claude Code assisted with code review, documentation drafting,
test generation, and implementation support during development of this package.
All package design decisions, analyses, validation, and release decisions were
reviewed and directed by the maintainer.
