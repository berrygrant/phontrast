# Classical Boundary Validation Summary (Feb 24, 2026)

## Goal
Validate that `phonJSD` is useful for phonologists using classical consonant cues, specifically:
1. stop voicing with VOT,
2. fricative place with spectral measures,
3. comparison of JSD against classical metrics (Pillai, Bhattacharyya, overlap, distance, classifier accuracy),
4. relationship between JSD and Pillai on MFCC-based contrasts.

## Data and scripts

### Corpus sources
- Aligned TextGrids: `analysis/sbcae_mfa_segments_aligned/`
- Segment audio: `analysis/sbcae_mfa_segments/`

### New analysis scripts used
- `analysis/sbcae_extract_classical_consonant_metrics.py`
- `analysis/validate_classical_boundaries.R`
- `analysis/compare_classical_metrics.R`
- `analysis/boundary_jsd_tests.R`

### Derived datasets and outputs
- Classical cues table: `analysis/sbcae_classical_consonant_metrics.csv`
- Classical JSD summary: `analysis/classical_outputs/classical_contrast_jsd_summary.csv`
- Classical feature descriptives: `analysis/classical_outputs/classical_feature_descriptives.csv`
- Cross-metric comparison: `analysis/classical_outputs/classical_metric_comparison.csv`
- Cross-metric ranking: `analysis/classical_outputs/classical_metric_rankings.csv`
- Existing MFCC comparison table: `comparison_outputs/pairwise_metrics_mfcc13.csv`

## 1) Real-consonant boundary smoke test (JSD runner)

A quick end-to-end run on real consonant tokens (40 aligned files; 2,763 tokens) showed:
- Global JSD estimates succeeded for `b~p`, `d~t`, `g~k`, `s~sh`, `f~th`.
- `z~zh` skipped due sparse `zh`.
- tone contrasts skipped because no `tone` column in this corpus.

This confirmed the boundary pipeline works on actual aligned consonant data and correctly reports skipped conditions.

## 2) Classical cue extraction and validation

### Classical cue extraction
Extraction over 300 aligned files yielded 19,167 consonant tokens.
Counts: `b 1199, d 3330, f 995, g 792, k 1783, p 1099, s 3319, sh 454, t 3666, th 358, z 2154, zh 18`.

### Stop voicing (VOT): expected low JSD
From `classical_contrast_jsd_summary.csv`:
- `b~p`: JSD = **0.0156** (n=362 vs 350)
- `d~t`: JSD = **0.0241** (n=700 vs 700)
- `g~k`: JSD = **0.0286** (n=262 vs 351)

All three are below the low-separation threshold (0.2) and flagged `jsd_is_low = TRUE`.

VOT descriptives (`classical_feature_descriptives.csv`):
- `b~p`: mean difference = +6.62 ms, Cohen's d = 0.225
- `d~t`: mean difference = +11.89 ms, Cohen's d = 0.389
- `g~k`: mean difference = +13.31 ms, Cohen's d = 0.419

(Voiceless > voiced as expected in this extraction.)

### Fricative place: expected high JSD
From `classical_contrast_jsd_summary.csv`:
- `s~sh`: JSD = **0.9769** (n=700 vs 454), Pillai = 0.2267
- `f~th`: JSD = **0.9229** (n=700 vs 358), Pillai = 0.0403

Top `s~sh` cues by effect size:
- `band_ratio_hi_lo_db` (d = -0.959)
- `cog_hz` (d = -0.659)
- `spec_slope_db_per_khz` (d = -0.652)

`z~zh` was skipped due low `zh` tokens (18).

## 3) JSD vs classical metrics on the same contrasts

From `classical_metric_comparison.csv`:

| Contrast | JSD | Pillai | Bhatt | 1-Overlap | Euclid | Mahalanobis | LDA acc |
|---|---:|---:|---:|---:|---:|---:|---:|
| b~p (VOT) | 0.0156 | NA | 0.0071 | 0.1231 | 6.62 | 0.2237 | 0.549 |
| d~t (VOT) | 0.0275 | NA | 0.0217 | 0.1504 | 12.63 | 0.4059 | 0.571 |
| g~k (VOT) | 0.0286 | NA | 0.0228 | 0.1796 | 13.31 | 0.4109 | 0.639 |
| s~sh (8D spectral) | 0.9881 | 0.2370 | 0.3824 | 0.9965 | 1481.28 | 0.9961 | 0.728 |
| f~th (8D spectral) | 0.9386 | 0.0331 | 0.1014 | 0.9738 | 305.96 | 0.3844 | 0.625 |

Domain means:
- stop voicing: mean JSD = 0.0239
- fricative place: mean JSD = 0.9634

Ranking against expected separation labels (stop=low, fricative=high), Spearman:
- JSD: 0.866
- Bhattacharyya: 0.866
- 1-overlap: 0.866
- Euclidean: 0.866
- LDA acc: 0.577
- Mahalanobis: 0.289
- Pillai: NA (insufficient comparable points; unavailable for 1D stop contrasts)

## 4) MFCC-specific JSD vs Pillai comparison

Using existing `comparison_outputs/pairwise_metrics_mfcc13.csv` (351 vowel-pair contrasts):
- Pearson corr(JSD, Pillai) = **0.8739**
- Spearman corr(JSD, Pillai) = **0.9114**

So they track strongly overall, but are not interchangeable.

Using existing overlap-correlation summary (`comparison_outputs/metric_overlap_correlations.csv`):
- corr(JSD, overlap): Pearson -0.9512, Spearman -0.9624
- corr(Pillai, overlap): Pearson -0.8109, Spearman -0.8980

In this MFCC13 set, JSD tracks overlap more tightly than Pillai.

## Interpretation

1. **Stop voicing (VOT):** JSD is consistently low for voiced/unvoiced stop pairs, aligning with the intended "small but meaningful" separation on this single classical cue.
2. **Fricative place:** JSD is very high for `s~sh` and `f~th` in classical spectral space, matching expected place-based spectral separation.
3. **Metric comparison:** JSD performs at least as well as standard alternatives in separating low-separation (stop voicing) from high-separation (fricative place) contrasts.
4. **MFCC context:** Pillai and JSD are strongly correlated but JSD appears more sensitive to distributional overlap in these data.

## Caveats

- Tone contrasts were not testable in this SBCSAE pipeline because no tone-label column exists in the aligned segment table.
- `z~zh` remains underpowered due sparse `zh` tokens.
- Pillai is naturally unavailable for 1D MANOVA-style use in the stop-VOT setup, so it appears as `NA` in that subset.

## Bottom line

The tests support that `phonJSD` gives value to users working with classical phonological measurements:
- it behaves correctly on VOT-based stop voicing,
- strongly captures fricative place distinctions with standard spectral cues,
- and compares favorably with common metrics while providing a distribution-focused separation measure.

## 2026-03-03 Full-Corpus Re-Alignment Update

### Updated corpus and alignment basis
- Full-file alignments are now available for all 60 recordings in `analysis/SBCSAE_MFA_Aligned_Final/`.
- Classical cue extraction was rerun from the repaired local corpus mirror in `analysis/sbcae_mfa_corpus/` after relinking the stale SBCSAE WAV symlinks.
- The new full-corpus cue table is `analysis/sbcae_classical_consonant_metrics_full_alignments.csv`.
- The new full-corpus outputs are in `analysis/classical_outputs_full_alignments/`.

### Canonical deliverables
- Alignments: `analysis/SBCSAE_MFA_Aligned_Final/`
- Full-corpus cue table: `analysis/sbcae_classical_consonant_metrics_full_alignments.csv`
- Final summary outputs: `analysis/classical_outputs_full_alignments/`

### Full-corpus extraction coverage
Phone-level alignments were obtained for all 60 SBCSAE recordings; classical consonant metrics were extracted from **157,346** target consonant tokens, and KDE-based comparisons used a cap of **300** tokens per category for tractability.
Counts: `b 10510, d 25305, f 8518, g 8171, k 14674, p 9650, s 26759, sh 3545, t 30441, th 2719, z 16934, zh 120`.

### Full-corpus JSD validation
The final reproducibility rerun retained the expected pattern:

- Stop voicing via VOT remained low:
  - `b~p`: **0.0154**
  - `d~t`: **0.0446**
  - `g~k`: **0.0645**

- Fricative place remained very high in the 8D spectral cue space:
  - `s~sh`: **0.9927**
  - `f~th`: **0.9540**
  - `z~zh`: **0.9875**

This is the same qualitative result as the earlier segment-based run, but now on the full reconstructed alignment set and reproduced from the repaired local corpus path.

### Full-corpus comparison against classical metrics
From `analysis/classical_outputs_full_alignments/classical_metric_comparison.csv`:

| Contrast | JSD | Pillai | Bhatt | 1-Overlap | Euclid | Mahalanobis | LDA acc |
|---|---:|---:|---:|---:|---:|---:|---:|
| b~p (VOT) | 0.0154 | NA | 0.0113 | 0.1153 | 9.30 | 0.2926 | 0.583 |
| d~t (VOT) | 0.0349 | NA | 0.0314 | 0.1755 | 15.48 | 0.4824 | 0.594 |
| g~k (VOT) | 0.0382 | NA | 0.0302 | 0.1839 | 14.18 | 0.4482 | 0.539 |
| s~sh (8D spectral) | 0.9907 | 0.1833 | 0.8858 | 0.9973 | 1918.20 | 0.8555 | 0.678 |
| f~th (8D spectral) | 0.9256 | 0.0143 | 0.1358 | 0.9647 | 68.65 | 0.2390 | 0.483 |
| z~zh (8D spectral) | 0.9957 | 0.2768 | 0.5076 | 0.9988 | 2206.33 | 1.1631 | 0.722 |

Ranking against expected low vs high separation labels, Spearman:
- JSD: **0.8783**
- Bhattacharyya: **0.8783**
- 1-overlap: **0.8783**
- Euclidean: **0.8783**
- Mahalanobis: **0.2928**
- LDA accuracy: **0.2928**
- Pillai: **NA**

### Interpretation of the full-corpus rerun
1. The full-file alignment rerun reproduces the original conclusion: JSD is low for stop voicing under VOT and very high for fricative place under standard spectral cues.
2. The result is now grounded in the complete aligned corpus rather than the earlier partial segment subset.
3. On these full alignments, JSD remains in the top-performing group relative to standard comparison metrics and preserves the desired contrast ordering.
