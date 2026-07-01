# AISHELL Mandarin Tone Validation Workflow

This document records the reproducible AISHELL-1 Mandarin tone workflow used for
the phonJSD common-contrast validation. It is intended to support both reruns
and paper-methods reporting.

## Goal

Validate whether phonJSD distinguishes Mandarin lexical tone contrasts, and test
how that distinguishability changes across acoustic feature parameterizations.

The workflow produces token-level tone features from AISHELL-1 recordings,
aligns a sampled corpus with MFA, extracts F0/prosodic features, and runs the
shared validation runner.

## Data And Paths

Corpus:

- AISHELL-1 Mandarin speech corpus.
- GPU host AISHELL root:
  `/mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Mandarin/data_aishell`
- Earlier NAS-native path:
  `/volume1/Corpus_Studies/Corpora/Mandarin/data_aishell`
- Earlier Mac path:
  `/Volumes/Corpus_Studies/Corpora/Mandarin/data_aishell`

Use the path visible from the machine running the command. On the GPU host, use
`/mnt/LUV_LAB_NAS/...`, not `/volume1/...`.

AISHELL manifest audit from the available corpus:

| Metric | Value |
| --- | ---: |
| recordings_total | 141600 |
| speakers_total | 400 |
| recordings_train | 120098 |
| recordings_dev | 14326 |
| recordings_test | 7176 |
| recordings_with_audio | 141600 |
| recordings_missing_audio | 0 |

Tone-unit parse audit from `pypinyin`:

| Metric | Value |
| --- | ---: |
| tone_units_total | 2040219 |
| tone_units_recommended_exclude | 85325 |
| tone_units_T1 | 447236 |
| tone_units_T2 | 424620 |
| tone_units_T3 | 348672 |
| tone_units_T4 | 734366 |
| tone_units_T5 | 85318 |

Neutral tone T5 is excluded by default for validation because it is marked
`recommended_exclude`.

## Repository Branch

Current development branch:

```sh
git checkout validation-corpus-adapters
git pull
git rev-parse HEAD
```

Record the commit hash in lab notes and in any archived output bundle. The
analysis outputs under `analysis/validation/outputs/` are intentionally ignored
by git.

## Software Environments

### Python / MFA Environment

The extraction scripts need Python packages for audio and pitch tracking. On the
GPU host, the MFA environment used:

```sh
/opt/miniconda3/envs/mfa/bin/python -c "import sys; print(sys.executable)"
```

Install or verify dependencies:

```sh
mamba install -n mfa -c conda-forge numpy pysoundfile librosa pooch pypinyin

/opt/miniconda3/envs/mfa/bin/python -c \
  "import numpy, soundfile, librosa, pooch, pypinyin; print('python_deps_ok')"
```

Use `pysoundfile` for conda/mamba package installation. The Python import name
is still `soundfile`.

### R / phonJSD Environment

Use a conda R environment to avoid system-library permissions and source
compilation failures:

```sh
mamba create -n phonjsd-r -c conda-forge r-base r-ks r-dplyr r-purrr r-tibble r-rlang
mamba run -n phonjsd-r R CMD INSTALL .
mamba run -n phonjsd-r Rscript -e 'library(phonJSD); packageVersion("phonJSD")'
```

Record R session information:

```sh
mamba run -n phonjsd-r Rscript -e 'sessionInfo()'
```

## MFA Models

Mandarin MFA dictionary and acoustic models are required. A language model is
not required for this forced-alignment workflow.

Record the exact local MFA model names:

```sh
mfa model list acoustic | grep -i mandarin
mfa model list dictionary | grep -i mandarin
mfa version
```

Use those names in the `mfa validate` and `mfa align` commands below.

## Step 1: Build The AISHELL Manifest

Run this on the machine where the AISHELL root path is valid:

```sh
python analysis/validation/prepare_aishell_mandarin_tone_tokens.py \
  --root /mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Mandarin/data_aishell \
  --out-dir analysis/validation/outputs/aishell_mandarin \
  --manifest-only
```

Expected outputs:

- `analysis/validation/outputs/aishell_mandarin/aishell_manifest.csv`
- `analysis/validation/outputs/aishell_mandarin/aishell_parse_summary.csv`

If `pypinyin` is available and the full tone-unit parse is needed, omit
`--manifest-only`.

## Step 2: Materialize An MFA Corpus Sample

The sample builder copies or symlinks WAV files and writes matching transcript
files for MFA.

### Small Dev Pilot

This was used for the first 100-utterance sanity check:

```sh
python analysis/validation/prepare_aishell_mfa_pilot.py \
  --root /mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Mandarin/data_aishell \
  --manifest analysis/validation/outputs/aishell_mandarin/aishell_manifest.csv \
  --out-dir /tmp/phonjsd_aishell_mfa_pilot \
  --splits dev \
  --max-speakers 5 \
  --max-utterances-per-speaker 20 \
  --mode copy \
  --overwrite
```

Observed pilot summary:

| Metric | Value |
| --- | ---: |
| recordings_selected | 100 |
| speakers_selected | 5 |
| recordings_dev | 100 |
| audio_copied | 100 |
| transcript_written | 100 |
| recordings_missing_audio | 0 |

### Scaled Matched Run

This larger sample is intended to make `control_group = pinyin_base` contrasts
estimable:

```sh
python analysis/validation/prepare_aishell_mfa_pilot.py \
  --root /mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Mandarin/data_aishell \
  --manifest analysis/validation/outputs/aishell_mandarin/aishell_manifest.csv \
  --out-dir /tmp/phonjsd_aishell_mfa_scaled \
  --splits train \
  --max-speakers 50 \
  --max-utterances-per-speaker 50 \
  --mode copy \
  --overwrite
```

This selects up to 2500 train utterances. The exact sample is deterministic
given the script default seed `20260630`.

Check the materialized sample:

```sh
cat /tmp/phonjsd_aishell_mfa_scaled/aishell_mfa_pilot_summary.csv
```

Observed scaled sample summary:

| Metric | Value |
| --- | ---: |
| recordings_selected | 2500 |
| speakers_selected | 50 |
| recordings_train | 2500 |
| audio_copied | 2500 |
| transcript_written | 2500 |
| recordings_missing_audio | 0 |

## Step 3: Validate And Align With MFA

Run MFA validation first:

```sh
mfa validate \
  /tmp/phonjsd_aishell_mfa_scaled/corpus \
  <mandarin_dictionary_model> \
  <mandarin_acoustic_model> \
  --clean \
  --num_jobs 16
```

For the 100-utterance pilot, validation reported only 80 OOVs and completed
successfully. If OOVs are unexpectedly high, rebuild the sample with
`--transcript-format chars-spaced` and rerun validation.

Align:

```sh
mfa align \
  /tmp/phonjsd_aishell_mfa_scaled/corpus \
  <mandarin_dictionary_model> \
  <mandarin_acoustic_model> \
  /tmp/phonjsd_aishell_mfa_scaled/aligned \
  --clean \
  --num_jobs 16
```

MFA mainly benefits from CPU cores, RAM, and local/scratch disk throughput. The
GPU itself is not expected to accelerate MFA.

## Step 4: Extract Tone Features

For the scaled run:

```sh
/opt/miniconda3/envs/mfa/bin/python \
  analysis/validation/extract_aishell_mandarin_tone_features.py \
  --pilot-manifest /tmp/phonjsd_aishell_mfa_scaled/aishell_mfa_pilot_manifest.csv \
  --aligned-dir /tmp/phonjsd_aishell_mfa_scaled/aligned \
  --out-dir analysis/validation/outputs/aishell_mandarin_tone_features_scaled \
  --pitch-method pyin \
  --frame-ms 50
```

For the 100-utterance pilot, use `/tmp/phonjsd_aishell_mfa_pilot` and
`analysis/validation/outputs/aishell_mandarin_tone_features`.

Feature extraction details:

- Tone labels are generated by `pypinyin` using `Style.TONE3` with neutral tone
  encoded as tone 5.
- Default validation excludes T5 and parse-flagged tokens.
- MFA provides word intervals.
- The current pilot assigns equal-duration character windows inside each
  aligned word interval.
- F0 is tracked with `librosa.pyin` when `--pitch-method pyin` is used.
- Pitch frame length is 50 ms with 10 ms hop; this avoids low-F0 frame warnings
  observed with 40 ms frames at 16 kHz.
- F0 values are converted to semitone-like log units and centered by speaker for
  the primary level/contour features.
- Tokens shorter than 35 ms or with fewer than two voiced F0 frames are skipped.

Expected extractor outputs:

- `aishell_mandarin_tone_features.csv`
- `aishell_mandarin_tone_feature_summary.csv`
- `aishell_mandarin_tone_feature_sets.csv`

Check extraction health:

```sh
cat analysis/validation/outputs/aishell_mandarin_tone_features_scaled/aishell_mandarin_tone_feature_summary.csv
```

Observed 100-utterance pilot extraction:

| Metric | Value |
| --- | ---: |
| pilot_recordings | 100 |
| textgrids_found | 100 |
| tone_source | built_from_pilot_manifest_with_pypinyin |
| tone_units_loaded | 1378 |
| tokens_written | 1315 |
| tokens_skipped_short_token_interval | 7 |
| tokens_skipped_too_few_voiced_frames | 56 |
| tokens_T1 | 303 |
| tokens_T2 | 276 |
| tokens_T3 | 241 |
| tokens_T4 | 495 |

Observed 2,500-recording scaled extraction:

| Metric | Value |
| --- | ---: |
| pilot_recordings | 2500 |
| textgrids_found | 2500 |
| tone_source | built_from_pilot_manifest_with_pypinyin |
| tone_units_loaded | 34412 |
| tokens_written | 32850 |
| tokens_skipped_no_word_interval | 2 |
| tokens_skipped_short_token_interval | 89 |
| tokens_skipped_too_few_voiced_frames | 1471 |
| tokens_T1 | 7564 |
| tokens_T2 | 7240 |
| tokens_T3 | 5814 |
| tokens_T4 | 12232 |

## Step 5: Run Validation

### Pooled Sanity Check

Use pooled tone contrasts when the sample is too small for same-base matching:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/aishell_mandarin_tone_features/aishell_mandarin_tone_features.csv \
  --out-dir analysis/validation/outputs/aishell_mandarin_tone_validation_pooled \
  --domain tone \
  --control-col pooled_tone \
  --feature-sets analysis/validation/outputs/aishell_mandarin_tone_features/aishell_mandarin_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

`--control-col pooled_tone` intentionally names a non-existent column, causing
the runner to estimate global pairwise tone contrasts.

Observed 100-utterance pooled pilot:

| Metric | Value |
| --- | ---: |
| n_input_rows | 1315 |
| n_feature_sets | 4 |
| n_contrasts | 6 |
| n_attempted | 24 |
| n_successful | 24 |
| n_skipped | 0 |

All 24 rows used `metric_mode = all_metrics` after removing the redundant
`f0_delta_st` feature from the contour bundle.

### Same-Base Matched Validation

Use this for the scaled sample:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/aishell_mandarin_tone_features_scaled/aishell_mandarin_tone_features.csv \
  --out-dir analysis/validation/outputs/aishell_mandarin_tone_validation_scaled_matched \
  --domain tone \
  --feature-sets analysis/validation/outputs/aishell_mandarin_tone_features_scaled/aishell_mandarin_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

This uses the default `control_group` column, which is the pinyin base with tone
number removed. It therefore compares tones only within matched syllable base,
provided each base has at least `--min-per-category` tokens in both categories.

Check:

```sh
cat analysis/validation/outputs/aishell_mandarin_tone_validation_scaled_matched/validation_run_summary.csv
head -30 analysis/validation/outputs/aishell_mandarin_tone_validation_scaled_matched/validation_metrics.csv
head -30 analysis/validation/outputs/aishell_mandarin_tone_validation_scaled_matched/validation_skipped_contrasts.csv
```

Summarize the successful rows:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/summarize_validation_metrics.R \
  --validation-dir analysis/validation/outputs/aishell_mandarin_tone_validation_scaled_matched
```

Observed 2,500-recording matched validation summary:

| Metric | Value |
| --- | ---: |
| n_input_rows | 32850 |
| n_feature_sets | 4 |
| n_contrasts | 861 |
| n_attempted | 3444 |
| n_successful | 750 |
| n_skipped | 2694 |
| bw | scott.diag |
| eval_on | pooled_sample |
| eval_n | 300 |
| engine | fast_diag |
| min_per_category | 20 |
| min_tokens | 40 |
| cv_folds | 5 |

Most skipped rows are expected under strict same-pinyin-base matching: many
syllable bases do not have at least 20 complete tokens in both tone categories
after filtering.

## Feature Sets

The Mandarin extractor writes four feature sets:

| Feature set | Contents | Role |
| --- | --- | --- |
| `f0_static` | speaker-centered F0 level, F0 spread/range, voiced proportion | Baseline pitch-level cue set |
| `f0_contour` | speaker-centered start/mid/end F0 and F0 slope | Primary tonal contour cue set |
| `energy_duration` | token/word duration, character position, RMS level/change/slope | Non-F0 prosodic control set |
| `tone_prosody` | combined F0 static, F0 contour, energy, and duration features | Full cue bundle |

The original contour bundle included `f0_delta_st`, but this is exactly derived
from `f0_end_st - f0_start_st`; it was removed from feature-set definitions to
avoid rank-deficient classical tests. The column may still exist in the token
table, but it is not used in the generated feature sets.

## Observed Pooled Pilot Pattern

The 100-utterance pooled pilot showed:

| Contrast | `f0_contour` JSD / AUC | `tone_prosody` JSD / AUC |
| --- | ---: | ---: |
| T1-T2 | 0.155 / 0.774 | 0.549 / 0.788 |
| T1-T3 | 0.204 / 0.802 | 0.642 / 0.808 |
| T1-T4 | 0.061 / 0.777 | 0.396 / 0.770 |
| T2-T3 | 0.058 / 0.630 | 0.474 / 0.666 |
| T2-T4 | 0.102 / 0.780 | 0.459 / 0.801 |
| T3-T4 | 0.107 / 0.759 | 0.509 / 0.782 |

Interpretation:

- F0 contour features carry the primary contrast.
- Energy/duration features alone are weaker, which is expected for lexical tone.
- The combined `tone_prosody` bundle produces substantially larger JSD than the
  F0-only contour set, while classifier AUC improves more modestly. This is
  consistent with the broader validation theme: JSD is informative, but strongly
  parameterization-sensitive.
- T2-T3 is the weakest pair, plausibly because both are contour tones and are
  harder to separate acoustically in short, word-derived windows.

## Observed Scaled Matched Pattern

The 2,500-recording same-base matched run produced 750 successful metric rows
from 32,850 extracted tone tokens. Only 3 of 750 rows used
`jsd_overlap_fallback`; the remaining rows retained the full classical metric
bundle.

| Feature set | Rows | Contrasts | Control groups | All metrics | Fallback | Median JSD | Median overlap | Median AUC | Median balanced accuracy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `energy_duration` | 199 | 199 | 102 | 198 | 1 | 0.221 | 0.568 | 0.689 | 0.601 |
| `f0_contour` | 176 | 176 | 89 | 176 | 0 | 0.139 | 0.648 | 0.786 | 0.674 |
| `f0_static` | 199 | 199 | 102 | 199 | 0 | 0.087 | 0.739 | 0.680 | 0.567 |
| `tone_prosody` | 176 | 176 | 89 | 174 | 2 | 0.748 | 0.134 | 0.793 | 0.732 |

Interpretation:

- The matched run preserves the central validation result: phonJSD is strongly
  sensitive to feature parameterization.
- `f0_contour` provides the clearest F0-only classifier reference, with median
  AUC around 0.79.
- `f0_static` is weaker than contour features, as expected for Mandarin tones.
- `energy_duration` alone carries nontrivial information, but the lower balanced
  accuracy suggests it should be treated as a secondary prosodic cue rather than
  the primary tone contrast.
- `tone_prosody` yields much larger median JSD and lower overlap than the
  narrower feature sets, while classifier AUC improves only modestly over
  `f0_contour`. This is useful evidence that JSD is capturing distributional
  separation in the chosen feature space, not simply reproducing classifier
  separability.
- The low fallback count means rank-deficient classical tests are no longer a
  major blocker for the Mandarin matched analysis.

## Paper-Methods Notes

Report the following explicitly:

- Tone labels are lexical/provisional and derived from `pypinyin`, not manually
  annotated tone labels.
- Neutral tone T5 and parse-flagged tokens are excluded by default.
- MFA alignments are word-level for the current pilot extraction.
- Character-level tone windows are equal partitions of the aligned word
  interval; these are not phone-level nuclei.
- F0 is extracted with `librosa.pyin` using 50 ms frames and 10 ms hop.
- F0 features are speaker-centered before validation.
- Validation includes both JSD/overlap metrics and a held-out binomial GLM
  classifier reference.
- Pooled tone validation is a sanity check; same-pinyin-base matched validation
  is the target design for the paper.

Suggested wording:

> Mandarin tone validation used AISHELL-1 recordings aligned with Montreal
> Forced Aligner. Lexical tone labels were assigned from transcript characters
> using pypinyin tone-number output. In the current extraction, each MFA-aligned
> word interval was divided evenly across its constituent transcript
> characters, and F0 trajectories were extracted with librosa.pyin. F0 features
> were speaker-centered prior to validation. Pooled tone contrasts were used as
> an initial sanity check; the preregistered analysis target is matched
> same-pinyin-base tone contrasts when token counts permit.

## Limitations And Next Refinements

The current extractor is suitable for pipeline validation and preliminary JSD
patterns. It is not yet the final paper-quality Mandarin tone extractor.

Needed refinements:

- Replace equal character windows with syllable or vowel-nucleus intervals from
  an MFA phone-tier parse.
- Audit word-label mismatches and OOV-aligned intervals.
- Consider a larger or full-corpus matched run after nucleus-level extraction,
  if the paper needs tighter confidence around per-base tone-pair summaries.
- Consider per-speaker or mixed/grouped validation summaries after the matched
  analysis is stable.
- Archive exact MFA model names, package versions, command history, summary
  CSVs, and git commit hash for final reproducibility.

## Output Bundle To Archive

For each run, archive:

- `aishell_mfa_pilot_manifest.csv`
- `aishell_mfa_pilot_summary.csv`
- MFA validation logs
- MFA alignment output directory or a compressed copy of TextGrids
- `aishell_mandarin_tone_features.csv`
- `aishell_mandarin_tone_feature_summary.csv`
- `aishell_mandarin_tone_feature_sets.csv`
- `validation_metrics.csv`
- `validation_skipped_contrasts.csv`
- `validation_run_summary.csv`
- `git rev-parse HEAD`
- `mfa version`
- MFA acoustic/dictionary model names
- `python` and R session/package versions
