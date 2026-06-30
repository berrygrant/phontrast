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
- `prepare_aishell_mandarin_tone_tokens.py`: builds an AISHELL-1 manifest and
  pypinyin-derived Mandarin tone-unit inventory. Its tone-unit CSV is also
  marked `pending_alignment`.
- `prepare_aishell_mfa_pilot.py`: samples the AISHELL manifest and materializes
  an MFA-ready pilot corpus with matching WAV and transcript files.

## AISHELL MFA pilot

On the alignment machine, first make sure the AISHELL manifest exists. If it was
created on another machine with a different path convention, rerun the manifest
builder with the AISHELL root as seen by the alignment machine.

Example for a NAS-mounted AISHELL root:

```sh
python analysis/validation/prepare_aishell_mandarin_tone_tokens.py \
  --root /volume1/Corpus_Studies/Corpora/Mandarin/data_aishell \
  --out-dir analysis/validation/outputs/aishell_mandarin \
  --manifest-only
```

Then materialize a small dev-set corpus for MFA validation:

```sh
python analysis/validation/prepare_aishell_mfa_pilot.py \
  --root /volume1/Corpus_Studies/Corpora/Mandarin/data_aishell \
  --manifest analysis/validation/outputs/aishell_mandarin/aishell_manifest.csv \
  --out-dir analysis/validation/outputs/aishell_mfa_pilot \
  --splits dev \
  --max-speakers 5 \
  --max-utterances-per-speaker 20 \
  --mode copy \
  --overwrite
```

The MFA corpus is written to
`analysis/validation/outputs/aishell_mfa_pilot/corpus`. Use `--mode symlink`
when the NAS mount is fast and stable; use `--mode copy` when aligning from a
local SSD or scratch directory.

If MFA reports many Mandarin dictionary OOVs, rebuild the pilot with
`--transcript-format chars-spaced` and re-run MFA validation. That produces one
Chinese character token per transcript word position while preserving the same
audio sample.

After downloading the Mandarin MFA dictionary and acoustic models, confirm their
local names. The language model is not required for this forced-alignment pilot.

```sh
mfa model list acoustic | grep -i mandarin
mfa model list dictionary | grep -i mandarin
```

Then run validation before full alignment:

```sh
mfa validate \
  analysis/validation/outputs/aishell_mfa_pilot/corpus \
  <mandarin_dictionary_model> \
  <mandarin_acoustic_model> \
  --clean \
  --num_jobs 8
```

If validation looks clean enough for a pilot, align the same subset:

```sh
mfa align \
  analysis/validation/outputs/aishell_mfa_pilot/corpus \
  <mandarin_dictionary_model> \
  <mandarin_acoustic_model> \
  analysis/validation/outputs/aishell_mfa_pilot/aligned \
  --clean \
  --num_jobs 8
```

For the scratch-directory pilot used on the GPU host, the corresponding command
is:

```sh
mfa align \
  /tmp/phonjsd_aishell_mfa_pilot/corpus \
  <mandarin_dictionary_model> \
  <mandarin_acoustic_model> \
  /tmp/phonjsd_aishell_mfa_pilot/aligned \
  --clean \
  --num_jobs 8
```

Once alignment completes, extract the first-pass Mandarin tone features:

```sh
python analysis/validation/extract_aishell_mandarin_tone_features.py \
  --pilot-manifest /tmp/phonjsd_aishell_mfa_pilot/aishell_mfa_pilot_manifest.csv \
  --tone-units analysis/validation/outputs/aishell_mandarin/aishell_mandarin_tone_units_pending_alignment.csv \
  --aligned-dir /tmp/phonjsd_aishell_mfa_pilot/aligned \
  --out-dir analysis/validation/outputs/aishell_mandarin_tone_features \
  --pitch-method auto
```

Use `--pitch-method pyin` on the GPU host if you want to require
`librosa.pyin`; `--pitch-method auto` falls back to a lightweight
autocorrelation tracker if librosa is unavailable. The extractor needs lexical
tone labels. If the tone-unit CSV is not present, it can rebuild labels from the
pilot manifest, but that fallback requires `pypinyin` in the active Python
environment.

This pilot extractor uses MFA word intervals and equal character windows inside
each aligned word. That is sufficient for an initial JSD pass over F0 contour
features; a later extractor should replace this with syllable or vowel-nucleus
intervals from a phone-tier parse.

Run the common validator on the extracted tone table:

```sh
Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/aishell_mandarin_tone_features/aishell_mandarin_tone_features.csv \
  --out-dir analysis/validation/outputs/aishell_mandarin_tone_validation \
  --domain tone \
  --feature-sets analysis/validation/outputs/aishell_mandarin_tone_features/aishell_mandarin_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```
