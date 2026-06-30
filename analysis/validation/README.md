# Common contrast validation scaffold

This directory contains the shared P0 validation layer for new phonJSD contrast
families. It assumes corpus-specific scripts have already produced a token-level
feature table. The scaffold then enumerates binary contrasts, computes the
standard phonJSD metric bundle, runs feature-set ablations, and adds an
independent held-out classifier reference.

## Token table contract

Required columns:

- `token_id`
- `domain`: for example `tone`, `stop_voicing`, or `nasal_vowels`
- `language`
- `source_corpus`
- `speaker`
- `category`: the category being contrasted
- `control_group`: optional matched comparison bin such as stop place, rime,
  vowel quality, or word position
- one or more numeric acoustic feature columns

Recommended optional columns:

- `file`, `start`, `end`
- `word`, `phone`, `syllable`, `rime`
- `preceding_phone`, `following_phone`
- `measurement_method`
- `quality_flag`

The runner accepts CSV, TSV, or RDS input. Rows with missing or non-finite values
in the selected feature set are removed per contrast before estimation.

## Feature-set file

Feature sets are supplied as CSV/TSV with at least these columns:

- `feature_set`
- `features`

The `features` cell may separate columns with spaces, commas, or semicolons.
See [feature_sets_template.csv](feature_sets_template.csv).

If no feature-set file is supplied, pass `--features` with a comma/semicolon/space
separated list. If neither is supplied, the runner uses all numeric columns after
excluding metadata columns such as `start` and `end`.

## Runner

Example:

```sh
Rscript analysis/validation/run_validation_metrics.R \
  --input path/to/token_features.csv \
  --out-dir outputs/validation_stop_voicing \
  --domain stop_voicing \
  --feature-sets analysis/validation/feature_sets_template.csv \
  --category-col category \
  --control-col control_group \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

Outputs:

- `validation_metrics.csv`: one successful estimate per contrast and feature set.
- `validation_skipped_contrasts.csv`: contrasts or feature sets that were not
  estimable, with a reason.
- `validation_run_summary.csv`: compact run counts and configuration.

## Scope

This scaffold does not extract acoustic features or perform forced alignment.
Those steps remain contrast-specific:

- stop voicing: event and spectral cue extraction around stop release;
- lexical tone: vowel/syllable nucleus alignment and F0 trajectory extraction;
- nasal/oral vowels: vowel segmentation plus formant, MFCC, and nasality-proxy
  extraction.

The point of P0 is to keep everything after token-feature extraction identical
across these domains.

## Corpus adapters in progress

- `prepare_sbcsae_stop_voicing_tokens.py`: extracts aligned stop-window
  acoustic features from SBCSAE MFA TextGrids and WAVs.
- `prepare_yoruba_slr86_tone_tokens.py`: builds an OpenSLR 86 Yoruba manifest
  and parses orthographic H/M/L tone-bearing units. Its current tone-unit CSV is
  intentionally marked `pending_alignment`; it does not contain segment timings
  or acoustic F0 features yet.
