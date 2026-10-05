# phontrast

**phontrast** is an R package for measuring the **separation and overlap of phonological categories** in acoustic space, using multiple complementary contrast metrics through one interface.
It is designed for researchers in sociophonetics, laboratory phonology, bilingualism, and speech perception who need principled, distributional measures of category contrast.

The package's entry point, `phontrast()`, computes and compares a family of metrics for a two-category contrast — Jensen–Shannon divergence and distance, the Pillai–Bartlett trace, Bhattacharyya distance and affinity, Mahalanobis distance, and proportional overlap — globally or by group, with optional bootstrap confidence intervals. Because the metrics differ in what they capture (distribution shape vs. mean separation vs. overlap), reporting several together gives a fuller picture of a contrast than any one alone.

`rank_contrasts()` turns the recommendations of the simulation study behind the package — rank speakers by Jensen–Shannon distance, check the ranking against Pillai, and read ranks and flags only where the sample licenses them — into one call with a plot and an inspection view; see [Ranking speakers and checking Pillai agreement](#ranking-speakers-and-checking-pillai-agreement).

> **Formerly `phonJSD`.** phontrast is the continuation of the `phonJSD` package (through v1.2.0), broadened from a Jensen–Shannon-divergence focus to a general multi-metric contrast toolkit. The estimators are unchanged; see [`NEWS.md`](NEWS.md) for migration notes, and `ROADMAP.md` in the repository for what's planned next.

[![CRAN status](https://www.r-pkg.org/badges/version/phontrast)](https://CRAN.R-project.org/package=phontrast)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20816585.svg)](https://doi.org/10.5281/zenodo.20816585)

---

## Motivation

Traditional measures of phonological contrast (e.g., Euclidean distance between means, Pillai scores) each capture only limited aspects of category structure.
**phontrast** treats phonological categories as **probability distributions** in acoustic space and brings several contrast measures — information-theoretic, parametric, and mean-based — under one roof, so you can compare them on a common footing instead of committing to a single metric up front.

This approach is especially useful when:
- Categories differ in **shape, variance, or multimodality**
- Token counts are unbalanced
- Category overlap is gradient rather than categorical
- You want measures that generalize naturally to **multidimensional features** (e.g., MFCCs)

---

## Core Features

- **One call, many metrics:** `phontrast()` computes Jensen–Shannon divergence and distance, Pillai–Bartlett trace, Bhattacharyya distance/affinity, Mahalanobis distance, and percent overlap side by side; opt-in metrics add total variation, matched-kernel Bhattacharyya and Hellinger distances, and the Euclidean distance between standardized means
- Pick any subset with the `metrics` argument
- **The ranking protocol in one call:** `rank_contrasts()` computes √JSD, Pillai, and shared mass on the same tokens, ranks speakers on the percentile scale, flags Pillai/√JSD disagreements of 0.25 or more, applies the sample-size floors, and runs a bandwidth check; `plot()` draws the agreement plot and `inspect_contrast()` a flagged speaker across bandwidths
- Kernel density–based estimation of acoustic distributions for the distributional metrics, with an optional multivariate-normal backend (`density = "mvnorm"`)
- Support for **1D and n-dimensional acoustic features**
- Optional high-dimensional KDE speed controls, including diagonal Scott
  bandwidths, sampled evaluation points, and a fast diagonal-Gaussian engine
- Global and group-level bootstrap summaries
- A common `orientation` / `separation_value` scale so overlap and separation metrics can be read together
- **Distribution-aware, accountable visualization:** `plot_contrast()` draws the same density model the metrics use (KDE highest-density regions or fitted-Gaussian ellipses), shades the overlap itself, and annotates panels with the metric values; plus metric-table, category-space, and PCA plots, all on a colorblind-safe Okabe–Ito palette (`theme_phontrast()`)
- Reproducible pipelines compatible with tidyverse workflows

---

## Quick Start

Most users should start with `phontrast()`. It takes a token-level vowel table,
compares two category labels in a numeric acoustic space, and returns the
contrast metrics side by side.

```r
library(phontrast)

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

phontrast(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker"
)
```

Select specific metrics, or switch to tidy long output for plotting and ranking:

```r
metrics_long <- phontrast(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker",
  metrics = c("jsd", "pillai", "overlap"),
  output = "long"
)

metrics_long[, c("group", "metric", "estimate", "orientation",
                 "separation_value", "separation_rank")]
```

If `ggplot2` is installed, the same workflow can be visualized directly. The
flagship visualization is `plot_contrast()`, a distribution-aware,
accountable view of one contrast: it draws the *same density model the
distributional metrics are computed from* (KDE highest-density regions by
default, coverage ellipses of the fitted Gaussians under
`density = "mvnorm"`), shades the density overlap itself, and labels each
panel with the Jensen–Shannon divergence and proportional overlap computed
under the plotted model. The caption records the estimator settings, and the
full metric table is attached as `attr(p, "contrast_metrics")`.

```r
plot_contrast(
  data = vowels,
  features = c("f2", "f1"),
  category_col = "vowel",
  group_col = "speaker",
  reverse_x = TRUE,
  reverse_y = TRUE
)
```

`phontrast()` results plot themselves — `plot(metrics_long)` (or
`ggplot2::autoplot()`) draws the metric comparison via
`plot_overlap_metrics()`:

```r
plot(metrics_long)
```

All plots share a colorblind-safe Okabe–Ito palette and a publication theme;
`theme_phontrast()`, `scale_color_phontrast()`, `scale_fill_phontrast()`, and
`phontrast_palette()` are exported so your own figures can match. Simple
exploratory views remain available via `plot_category_space()`.

For multidimensional features such as MFCCs or acoustic embeddings, compute
metrics on the full feature set and use PCA plots as a diagnostic projection:

```r
features <- paste0("mfcc", 1:13)

metrics_13d <- phontrast(
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
metrics_13d_fast <- phontrast(
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

These KDE controls affect only the distributional metrics (JSD and percent
overlap): one kernel is still centered on each observed token, and the speedup
comes from a cheaper diagonal bandwidth rule, evaluating densities on a bounded
sampled support, and a vectorized diagonal-Gaussian evaluator. The package
default is unchanged for backward compatibility.

### Choosing the density backend

The distributional metrics (Jensen–Shannon divergence and percent overlap) are
computed from a density estimate for each category. By default that estimate is
a kernel density (`density = "kde"`). You can instead fit one multivariate
normal per category with `density = "mvnorm"`:

```r
phontrast(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  density = "mvnorm"
)
```

Under `density = "mvnorm"` the two categories are modelled as Gaussians and the
JSD and overlap between them are estimated by Monte-Carlo (Jensen–Shannon
divergence between two Gaussians has no closed form); `mc_n` sets the number of
Monte-Carlo samples and `eval_seed` makes the estimate reproducible. This
decouples the *density estimator* from the *metric*, so the estimator behind
JSD and overlap can be matched to the same multivariate-normal assumptions that
the Pillai, Bhattacharyya, and Mahalanobis columns already make — useful when
comparing metrics on a common footing, or for higher-dimensional feature spaces
where multivariate KDE is impractical. The Pillai, Bhattacharyya, and
Mahalanobis columns are parametric by construction and are unaffected by this
argument. `density` is available on `phontrast()`, `estimate_jsd()`,
`estimate_overlap()`, `jsd_summary()`, `global_boot_jsd()`, `jsd_kde_nd()`, and
`percent_overlap_kde()`.

### Ranking speakers and checking Pillai agreement

The simulation study behind phontrast ("Estimand or estimator? Comparing
vowel overlap measures against a known ground truth", Berry, under review)
closes with a measurement protocol: compute √JSD and Pillai on the same tokens
and report shared mass beside them; rank speakers by √JSD; flag speakers whose
Pillai percentile rank differs from their √JSD percentile rank by 0.25 of the
ordering or more, and inspect them. `rank_contrasts()` runs the three steps
and applies their conditions — the √JSD ceiling, the sample-size floors per
speaker and dimensionality, and a bandwidth check at half and twice the
diagonal Scott bandwidth — on the bundled `vowel_cohort` or your own data:

```r
ranking <- rank_contrasts(
  data = vowel_cohort,
  features = c("f1", "f2"),
  category_col = "vowel",
  group_col = "speaker"
)
ranking                      # the protocol table, with a header that states the conditions
plot(ranking)                # Pillai rank vs. sqrt(JSD) rank, with the +/- 0.25 band
inspect_contrast(ranking, "spk09", reverse_x = TRUE, reverse_y = TRUE)
```

`protocol_floors(d)` returns the licensing floors, `recommended_estimator(d)`
the kernel settings the study calibrated at each dimensionality, and
`percentile_rank()` the rank scale. The vignette *Ranking speakers by
Jensen–Shannon distance and checking Pillai agreement* walks through the whole
protocol.

If you only need the package's information-theoretic estimate on its own, use
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

### Proportion-standardized Pillai estimates

For an exactly two-category design, `pillai_overlap()` can optionally estimate
the squared Mahalanobis separation and the Pillai score the same estimated
system would have under a balanced category split:

```r
pillai_overlap(
  data = vowels,
  features = c("f1", "f2"),
  category_col = "vowel",
  proportion_standardized = TRUE
)
```

The default remains `FALSE`, so existing calls still return only `pillai` and
`p_value`. Raw Pillai describes the categories in the realized speech data;
the opt-in `pillai_eq` estimate characterizes the corresponding balanced-design
system. The extended result also reports `d2_plugin`, `d2_unbiased`, the class
counts and harmonic mean, and the finite-sample bias term.

When the unbiased separation estimate is negative, `pillai_eq` is `NA` and the
result sets `pillai_eq_fallback = TRUE`. The nonnegative fallback is reported
only as `d2_fallback`, alongside its known upward bias; it is never silently
mapped onto the corrected-score scale. This is a common outcome near merger,
not an exceptional error case. These estimates assume multivariate normality
within each category and a common within-category covariance, and no covariates
or multicategory extension is performed.

> `compare_overlap_metrics()` from phonJSD still works but is **deprecated**: it
> now calls `phontrast()` with `output = "wide"`. Switch calls to `phontrast()`.

---

## Conceptual Overview

Several of phontrast's metrics are distributional. For the Jensen–Shannon
family, given two phonological categories (e.g., vowels /ɪ/ and /ɛ/), phontrast:

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
of the continuous JSD in any number of dimensions. A sample-size-scaled partial
leave-one-out correction reduces resubstitution bias while keeping small real
divergences as small positive values rather than flooring them to exactly 0.
The pre-1.2.0 self-normalized sample-point estimate — a bounded relative
separation index rather than the JSD integral — remains available as
`method = "legacy"` for reproducing 1.0.0 results.

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
| Total variation (`tv`, opt-in) | Higher = more separation | 0-1 | You want the separation complement of percent overlap | Same kernel estimator as percent overlap |
| Bhattacharyya distance / affinity, kernel (`bhattacharyya_kde`, opt-in) | Distance: higher = more separation; affinity: higher = more overlap | 0 to infinity / 0-1 | You want the Bhattacharyya quantity under the same kernel densities as JSD, to separate measure from estimator | Kernel estimates move with the bandwidth; compare with the closed-form columns |
| Hellinger distance (`bhattacharyya_kde`, opt-in) | Higher = more separation | 0-1 | You want a bounded proper metric from the same kernel densities | As above |
| Euclidean distance of standardized means (`euclidean`, opt-in) | Higher = more separation | 0 to 2√d | You want the classic centroid distance on a comparable scale | Mean-based; blind to shape and variance |

For side-by-side comparison, `phontrast(output = "long")` adds `orientation`,
`separation_value`, and `separation_rank` columns so overlap metrics can be read
on the same separation-oriented scale as distance metrics. Confidence interval
columns are named `ci_lower` and `ci_upper`; legacy JSD-specific aliases
`jsd_low` and `jsd_high` are retained for compatibility.

## PB52 Example Note

The Peterson and Barney 1952 data in `phonTools::pb52` are useful for global
vowel contrasts. After filtering the factor-valued `vowel` column, phontrast
counts observed levels only, so manual `droplevels()` is not required.

```r
data(pb52, package = "phonTools")
pb_i <- subset(pb52, as.character(vowel) %in% c("I", "i"))

phontrast(
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

Install the current stable release from CRAN:

```r
install.packages("phontrast")
```

To install the latest tagged GitHub release instead, use:

```r
# install.packages("remotes")
remotes::install_github("berrygrant/phontrast@v2.4.1")
```

For the current development version, use:

```r
remotes::install_github("berrygrant/phontrast")
```

## Citation

To obtain the citation and BibTeX entry for the installed release, run:

```r
citation("phontrast")
toBibtex(citation("phontrast"))
```

## AI Use Disclosure

OpenAI ChatGPT/Codex and Anthropic Claude Code assisted with code review, documentation drafting,
test generation, and implementation support during development of this package.
All package design decisions, analyses, validation, and release decisions were
reviewed and directed by the maintainer.
