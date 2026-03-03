# phonJSD

**phonJSD** is an R package for measuring phonological category separation using **Jensen–Shannon Divergence (JSD)**.  
It is designed for researchers working in sociophonetics, laboratory phonology, bilingualism, and speech perception who need a principled, distributional metric of category overlap in acoustic space.

Version **0.1.0** is an early research release focused on core functionality, transparency, and reproducibility.

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

## Core Features (v0.1.0)

- Jensen–Shannon Divergence for phonological category comparison  
- Kernel density–based estimation of acoustic distributions  
- Support for **1D and n-dimensional acoustic features**
- Reproducible pipelines compatible with tidyverse workflows
- Designed for integration with forced alignment and acoustic extraction tools

---

## Conceptual Overview

Given two phonological categories (e.g., vowels /ɪ/ and /ɛ/), phonJSD:

1. Represents each category as a probability distribution over acoustic space  
2. Estimates densities using kernel density estimation (KDE)  
3. Computes Jensen–Shannon Divergence between the distributions  

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
```

## Boundary Testing (Consonants + Tones)

To test JSD on commonly studied non-vowel boundaries (stop voicing, fricative distinctions, and tone contrasts), run:

```bash
BOUNDARY_DATA_PATH="/absolute/path/to/your_data.csv" \
BOUNDARY_GROUP_COL="speaker" \
BOUNDARY_SEGMENT_COL="segment" \
BOUNDARY_TONE_COL="tone" \
BOUNDARY_FEATURES="mfcc1,mfcc2,mfcc3,mfcc4,mfcc5,mfcc6,mfcc7,mfcc8,mfcc9,mfcc10,mfcc11,mfcc12,mfcc13" \
Rscript analysis/boundary_jsd_tests.R
```

This writes:
- `analysis/boundary_outputs/boundary_jsd_global.csv`
- `analysis/boundary_outputs/boundary_jsd_by_group.csv`
- `analysis/boundary_outputs/boundary_jsd_skipped.csv`

Optional tuning parameters:
- `BOUNDARY_N_BOOT_GLOBAL` (default `300`; set `0` for point-estimate only)
- `BOUNDARY_N_BOOT_GROUP` (default `150`; set `0` for point-estimate only)
- `BOUNDARY_MIN_TOKENS_GLOBAL` (default `80`)
- `BOUNDARY_MIN_TOKENS_GROUP` (default `30`)
- `BOUNDARY_MIN_PER_CATEGORY` (default `15`)
- `BOUNDARY_BW` (default `Hpi.diag`)
- `BOUNDARY_EVAL_ON` (default `pooled`)
- `BOUNDARY_OUT_DIR` (default `analysis/boundary_outputs`)
- `BOUNDARY_LOAD_LOCAL` (default `false`; set `true` to force `devtools::load_all()` from local source)

Notes:
- The script includes preset aliases for contrasts such as `b~p`, `d~t`, `g~k`, `s~sh`, `z~zh`, `f~th`, high~low tone, and rising~falling tone.
- If your labels use a different scheme, edit `analysis/boundary_jsd_tests.R` alias lists for each contrast.
