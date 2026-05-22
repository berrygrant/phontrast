# phonJSD

**phonJSD** is an R package for measuring phonological category separation using **Jensen–Shannon Divergence (JSD)**.  
It is designed for researchers working in sociophonetics, laboratory phonology, bilingualism, and speech perception who need a principled, distributional metric of category overlap in acoustic space.

Version **0.5.0** is a research release focused on stable core metrics, reproducible uncertainty estimates, and comparison with classical overlap measures.
[![DOI](https://zenodo.org/badge/1192114337.svg)](https://doi.org/10.5281/zenodo.20346257)

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

## Core Features (v0.5.0)

- Jensen–Shannon Divergence for phonological category comparison  
- Kernel density–based estimation of acoustic distributions  
- Support for **1D and n-dimensional acoustic features**
- Global and group-level bootstrap summaries
- Comparison metrics including Pillai-Bartlett trace, Bhattacharyya distance, and percent overlap
- Reproducible pipelines compatible with tidyverse workflows
- Designed for integration with forced alignment and acoustic extraction tools

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

---

## Installation

This package is currently in early development and not yet on CRAN.

```r
# install.packages("remotes")
remotes::install_github("berrygrant/phonJSD")
