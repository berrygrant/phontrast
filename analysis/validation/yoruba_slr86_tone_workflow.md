# Yoruba OpenSLR 86 Tone Workflow

This document records the reproducible Yoruba OpenSLR 86 workflow for the
phonJSD common-contrast validation. It covers manifest construction,
orthographic tone-unit parsing, MFA corpus materialization, pilot alignment,
F0/prosodic feature extraction, validation, and the audit steps needed before
paper-scale reporting.

## Goal

Use OpenSLR 86 Yoruba as the first non-Mandarin lexical tone validation corpus.
The corpus provides tone-marked Yoruba transcripts and WAV recordings, allowing
H/M/L tone labels to be parsed directly from orthography before acoustic
alignment.

The current workflow produces:

- a recording-level manifest;
- an orthographic tone-bearing-unit table marked `pending_alignment`;
- parse summaries and count tables for audit;
- an MFA-ready Yoruba corpus sample or full-clean corpus;
- a generated grapheme-style dictionary for MFA;
- aligned phone-level tone-token acoustic features after MFA;
- pooled and vowel-quality matched validation outputs;
- a deterministic hand-audit sample for alignment/token quality control.

## Data And Paths

Corpus:

- OpenSLR 86 Yoruba.
- Local Mac root: `/Volumes/Corpus_Studies/Corpora/Yoruba`
- NAS root: `/volume1/Corpus_Studies/Corpora/Yoruba`

Expected source files:

- `line_index_female.tsv`
- `line_index_male.tsv`
- `yo_ng_female.zip`
- `yo_ng_male.zip`
- extracted audio directory `yo_ng_female/`
- extracted audio directory `yo_ng_male/`

Current completed NAS extraction status:

| Metric | Value |
| --- | ---: |
| recordings_with_extracted_audio | 3583 |
| recordings_with_zip_audio | 0 |
| recordings_with_available_audio | 3583 |
| recordings_missing_audio | 0 |
| speakers_total | 36 |
| tone_units_total | 51412 |

The earlier partial male extraction issue has been resolved. All 3,583 WAVs are
available as extracted audio in the current NAS corpus copy.

## Repository Branch

Current development branch:

```sh
git checkout validation-corpus-adapters
git pull
git rev-parse HEAD
```

Record the commit hash in lab notes and in any archived output bundle. Output
CSVs under `analysis/validation/outputs/` are ignored by git.

## Software Environment

The current Yoruba parser uses only the Python standard library. It does not
require `librosa`, `pypinyin`, MFA, R, or phonJSD.

Recommended command to verify the Python interpreter:

```sh
python --version
python analysis/validation/prepare_yoruba_slr86_tone_tokens.py --help
```

Later acoustic extraction will require an audio/F0 environment comparable to the
Mandarin extractor environment.

The MFA pilot materializer also uses only the Python standard library. MFA
itself is required only for the subsequent validation/alignment commands.

## Step 1: Parse The Corpus Manifest And Tone Units

Run on the NAS or a machine with fast local access to the extracted audio:

```sh
python analysis/validation/prepare_yoruba_slr86_tone_tokens.py \
  --root /volume1/Corpus_Studies/Corpora/Yoruba \
  --skip-audio-headers
```

On the GPU host where the NAS is mounted at `/mnt/LUV_LAB_NAS`, use:

```sh
python analysis/validation/prepare_yoruba_slr86_tone_tokens.py \
  --root /mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Yoruba \
  --skip-audio-headers
```

Use the Mac-mounted path if running locally:

```sh
python analysis/validation/prepare_yoruba_slr86_tone_tokens.py \
  --root /Volumes/Corpus_Studies/Corpora/Yoruba \
  --skip-audio-headers
```

`--skip-audio-headers` checks presence without opening each WAV file. It is
faster over mounted storage and sufficient for manifest/tone-label auditing.
For full audio duration/sample-rate audit, rerun without `--skip-audio-headers`
on data-local storage.

Expected outputs:

- `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv`
- `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv`
- `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_parse_summary.csv`
- `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_unit_counts.csv`

Observed command output:

```text
Wrote 3583 manifest rows to analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv
Wrote 51412 tone-unit rows to analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv
Wrote summaries to analysis/validation/outputs/yoruba_slr86/yoruba_slr86_parse_summary.csv and analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_unit_counts.csv
```

Check the summary:

```sh
cat analysis/validation/outputs/yoruba_slr86/yoruba_slr86_parse_summary.csv
```

## Tone-Label Parsing

Tone labels come directly from Yoruba orthography after Unicode normalization.
No acoustic values or timing values are inferred at this stage.

Rules:

- acute accent: high tone `H`;
- grave accent: low tone `L`;
- no acute/grave mark on a vowel-bearing unit: mid tone `M`;
- simultaneous acute and grave marks: `ambiguous`, flagged for exclusion;
- syllabic nasal `n` with acute/grave is treated as a tone-bearing unit;
- vowel quality is the `control_group` for first-pass validation:
  `a`, `e`, `e_dot`, `i`, `o`, `o_dot`, `u`, or `syllabic_n`.

The parser uses Unicode decomposition internally to detect combining tone marks
and dot-below marks, then writes normalized output.

## Speaker And Metadata Handling

Recording IDs are split into three fields. Speaker IDs are constructed from the
recording prefix plus the middle recording-ID field:

```text
yof_06136
yom_06136
```

The raw middle field overlaps across female and male files, so using it alone
would collapse distinct speakers. The manifest keeps both:

- `speaker`: prefix plus middle field;
- `speaker_code`: raw middle field;
- `speaker_source`: provenance note for the speaker assignment.

The parser also writes:

- `sex`
- `recording_prefix`
- `item_code`
- `source_line_index`
- `source_line_number`
- `transcript`
- `transcript_clean`
- `annotation_labels`
- `recommended_exclude`

## Annotation And Exclusion Rules

OpenSLR 86 transcripts include bracketed annotations such as:

- `[abrupt]`
- `[breath]`
- `[external]`
- `[hesitation]`
- `[snap]`

The parser strips bracketed annotations from `transcript_clean`, keeps the labels
in `annotation_labels`, and marks the recording as `recommended_exclude`.

Default downstream validation should exclude `recommended_exclude = TRUE` until
we have an explicit robustness check showing that annotations do not affect the
result.

## Current Output Schema

The manifest output is recording-level. Important columns include:

- `recording_id`
- `sex`
- `speaker`
- `speaker_code`
- `transcript`
- `transcript_clean`
- `annotation_labels`
- `recommended_exclude`
- `file`
- `audio_exists`
- `audio_available`
- `audio_location`

The tone-unit output is token-level but not yet acoustically aligned. Important
columns include:

- `token_id`
- `domain = tone`
- `language = Yoruba`
- `source_corpus = OpenSLR86`
- `speaker`
- `sex`
- `recording_id`
- `file`
- `word`
- `word_index`
- `category`: `H`, `M`, `L`, or `ambiguous`
- `control_group`: orthographic vowel quality
- `orthographic_unit`
- `orthographic_base`
- `vowel_quality`
- `tbu_type`
- `tone_unit_position`
- `word_tone_unit_count`
- `word_unit_index`
- `preceding_letter`
- `following_letter`
- `following_n_letter`
- `nasal_oral_status = not_assessed`
- `start`, `end`: blank until alignment
- `measurement_method = orthographic_tone_parse_pending_alignment`
- `quality_flag`
- `recommended_exclude`

Because `start` and `end` are blank and there are no acoustic feature columns,
this table should not be passed to `run_validation_metrics.R` yet.

## Step 2: Materialize An MFA Pilot

Create a small MFA-ready pilot corpus before attempting full alignment. The
pilot materializer samples speakers from the manifest, excludes annotation-
flagged recordings by default, writes cleaned tone-stripped `.lab` transcripts,
and creates a deterministic grapheme-style dictionary. Tone marks are removed
for alignment, while segmental Yoruba distinctions such as dot-below vowels are
preserved.

Example 10-speaker pilot on the Mac mount, rewriting NAS paths from the stored
manifest:

```sh
python analysis/validation/prepare_yoruba_mfa_pilot.py \
  --manifest analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv \
  --out-dir /tmp/phonjsd_yoruba_mfa_pilot_10 \
  --max-speakers 10 \
  --max-utterances-per-speaker 20 \
  --balance-sex \
  --source-root-rewrite /volume1/Corpus_Studies/Corpora/Yoruba=/Volumes/Corpus_Studies/Corpora/Yoruba \
  --mode copy \
  --overwrite
```

On the GPU/NAS host, use the path convention visible there, for example:

```sh
python analysis/validation/prepare_yoruba_mfa_pilot.py \
  --manifest analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv \
  --out-dir /tmp/phonjsd_yoruba_mfa_pilot_10 \
  --max-speakers 10 \
  --max-utterances-per-speaker 20 \
  --balance-sex \
  --source-root-rewrite /volume1/Corpus_Studies/Corpora/Yoruba=/mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Yoruba \
  --mode copy \
  --overwrite
```

Observed 10-speaker pilot summary:

| Metric | Value |
| --- | ---: |
| recordings_selected | 200 |
| speakers_selected | 10 |
| recordings_female | 100 |
| recordings_male | 100 |
| audio_copied | 200 |
| transcript_written | 200 |
| recordings_missing_audio | 0 |
| recordings_skipped_empty_source_audio | 1 |
| recordings_skipped_recommended_exclude | 1160 |
| dictionary_words | 599 |

The pilot outputs are:

- `/tmp/phonjsd_yoruba_mfa_pilot_10/corpus`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_manifest.csv`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_summary.csv`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_dictionary.txt`

For the next paper-oriented run, materialize the full clean corpus: all
recordings with available audio and no bracketed annotation flags. OpenSLR 86
has only 36 speakers, so a full clean run is preferable to an arbitrary larger
speaker sample.

```sh
python analysis/validation/prepare_yoruba_mfa_pilot.py \
  --manifest analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv \
  --out-dir /tmp/phonjsd_yoruba_mfa_full_clean \
  --max-speakers 0 \
  --max-utterances-per-speaker 0 \
  --balance-sex \
  --source-root-rewrite /volume1/Corpus_Studies/Corpora/Yoruba=/mnt/LUV_LAB_NAS/Corpus_Studies/Corpora/Yoruba \
  --mode copy \
  --overwrite
```

Check the materialized sample before alignment:

```sh
cat /tmp/phonjsd_yoruba_mfa_full_clean/yoruba_mfa_pilot_summary.csv
```

## Step 3: Alignment

No pretrained Yoruba MFA acoustic model was available, so the pilot used MFA
training from the generated grapheme-style dictionary:

```sh
mfa train \
  /tmp/phonjsd_yoruba_mfa_pilot_10/corpus \
  /tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_dictionary.txt \
  /tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_acoustic_model.zip \
  --output_directory /tmp/phonjsd_yoruba_mfa_pilot_10/aligned \
  --clean \
  --overwrite \
  --num_jobs 8
```

Observed pilot alignment output:

| Metric | Value |
| --- | ---: |
| TextGrids written | 200 |
| word tier | `words` |
| phone tier | `phones` |

The pilot is therefore sufficient for phone-level tone-token F0 extraction. A
larger or full-corpus run should reuse this training path, then audit a sample
of word/phone alignments before final paper analyses.

Full clean alignment:

```sh
mfa train \
  /tmp/phonjsd_yoruba_mfa_full_clean/corpus \
  /tmp/phonjsd_yoruba_mfa_full_clean/yoruba_mfa_pilot_dictionary.txt \
  /tmp/phonjsd_yoruba_mfa_full_clean/yoruba_mfa_full_clean_acoustic_model.zip \
  --output_directory /tmp/phonjsd_yoruba_mfa_full_clean/aligned \
  --clean \
  --overwrite \
  --num_jobs 16
```

## Step 4: Yoruba Tone Feature Extraction

The extractor mirrors the Mandarin tone feature contract, but with
Yoruba-specific labels and controls. It joins the parsed H/M/L tone-bearing
units to MFA `words` and `phones` tiers. Within each word, parsed tone-bearing
units are matched in sequence to aligned vowel or syllabic-nasal phone
intervals.

Run on the 10-speaker pilot:

```sh
/opt/miniconda3/envs/mfa/bin/python \
  analysis/validation/extract_yoruba_slr86_tone_features.py \
  --pilot-manifest /tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_manifest.csv \
  --tone-units analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv \
  --aligned-dir /tmp/phonjsd_yoruba_mfa_pilot_10/aligned \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10 \
  --pitch-method pyin \
  --frame-ms 50
```

If `librosa.pyin` is not available in the active environment, use
`--pitch-method autocorr` as a dependency-light smoke test, then rerun with
`pyin` before interpreting results.

Observed 10-speaker extraction:

| Metric | Value |
| --- | ---: |
| pilot_recordings | 200 |
| textgrids_found | 200 |
| tone_source | loaded_yoruba_tone_units_csv |
| tone_units_loaded | 2959 |
| tokens_written | 2622 |
| tokens_skipped_short_token_interval | 310 |
| tokens_skipped_too_few_voiced_frames | 27 |
| tokens_H | 984 |
| tokens_L | 687 |
| tokens_M | 951 |

Full clean extraction:

```sh
/opt/miniconda3/envs/mfa/bin/python \
  analysis/validation/extract_yoruba_slr86_tone_features.py \
  --pilot-manifest /tmp/phonjsd_yoruba_mfa_full_clean/yoruba_mfa_pilot_manifest.csv \
  --tone-units analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv \
  --aligned-dir /tmp/phonjsd_yoruba_mfa_full_clean/aligned \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_features_full_clean \
  --pitch-method pyin \
  --frame-ms 50 \
  --progress-every 100
```

Observed full-clean extraction:

| Metric | Value |
| --- | ---: |
| pilot_recordings | 2422 |
| textgrids_found | 2422 |
| tone_source | loaded_yoruba_tone_units_csv |
| tone_units_loaded | 33976 |
| tokens_written | 31343 |
| tokens_skipped_short_token_interval | 2090 |
| tokens_skipped_too_few_voiced_frames | 543 |
| tokens_H | 12051 |
| tokens_L | 8266 |
| tokens_M | 11026 |

Recommended token-feature columns:

- existing tone-unit metadata from
  `yoruba_slr86_tone_units_pending_alignment.csv`;
- aligned `start` and `end` for each tone-bearing unit;
- aligned interval source: word, phone, syllable, or vowel nucleus;
- `duration_ms`;
- `f0_total_frames`;
- `f0_voiced_frames`;
- `f0_voiced_prop`;
- speaker-centered F0 level features;
- time-normalized F0 contour points or start/mid/end F0;
- F0 slope and range;
- RMS/intensity controls;
- quality flags for low voicing, short intervals, annotation-bearing recordings,
  mismatched alignment, and nasal contexts.

Generated feature sets:

| Feature set | Contents | Role |
| --- | --- | --- |
| `f0_static` | speaker-centered F0 level, F0 range/spread, voiced proportion | Baseline tone cue set |
| `f0_contour` | start/mid/end F0 or 5-10 time-normalized F0 points plus slope | Primary Yoruba tone trajectory set |
| `energy_duration` | interval duration, position, RMS level/change/slope | Non-F0 prosodic controls |
| `tone_prosody` | combined F0, energy, and duration features | Full cue bundle |

## Step 5: Validation Commands

Once a Yoruba acoustic token-feature table exists, run pooled and matched
validation separately. For the 10-speaker pilot, the input path is
`analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10/yoruba_slr86_tone_features.csv`.

Pooled sanity check:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_pilot_10_pooled \
  --domain tone \
  --control-col pooled_tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10/yoruba_slr86_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

Matched vowel-quality validation:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_pilot_10_vowel_matched \
  --domain tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features_pilot_10/yoruba_slr86_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

The matched run uses the default `control_group` column, currently vowel
quality. We may later refine `control_group` to include vowel quality plus
position or local context if counts support it.

Full clean pooled validation:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/yoruba_slr86_tone_features_full_clean/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_full_clean_pooled \
  --domain tone \
  --control-col pooled_tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features_full_clean/yoruba_slr86_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

Full clean vowel-quality matched validation:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/yoruba_slr86_tone_features_full_clean/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_full_clean_vowel_matched \
  --domain tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features_full_clean/yoruba_slr86_tone_feature_sets.csv \
  --min-per-category 20 \
  --bw scott.diag \
  --eval-on pooled_sample \
  --eval-n 300 \
  --engine fast_diag \
  --cv-folds 5
```

Summarize each validation run:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/summarize_validation_metrics.R \
  --validation-dir analysis/validation/outputs/yoruba_slr86_tone_validation_full_clean_pooled

mamba run -n phonjsd-r Rscript analysis/validation/summarize_validation_metrics.R \
  --validation-dir analysis/validation/outputs/yoruba_slr86_tone_validation_full_clean_vowel_matched
```

## Observed 10-Speaker Pilot Pattern

The 10-speaker pilot produced a usable phone-level tone validation. All pooled
and vowel-matched rows used `metric_mode = all_metrics`; no fallback rows were
needed.

Pooled H/M/L validation:

| Feature set | Rows | Contrasts | All metrics | Median JSD | Median overlap | Median AUC | Median balanced accuracy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `energy_duration` | 3 | 3 | 3 | 0.083 | 0.741 | 0.650 | 0.586 |
| `f0_contour` | 3 | 3 | 3 | 0.133 | 0.637 | 0.753 | 0.699 |
| `f0_static` | 3 | 3 | 3 | 0.107 | 0.681 | 0.746 | 0.669 |
| `tone_prosody` | 3 | 3 | 3 | 0.295 | 0.476 | 0.766 | 0.705 |

Vowel-quality matched validation:

| Feature set | Rows | Contrasts | Control groups | All metrics | Median JSD | Median overlap | Median AUC | Median balanced accuracy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `energy_duration` | 21 | 21 | 7 | 21 | 0.194 | 0.586 | 0.664 | 0.602 |
| `f0_contour` | 21 | 21 | 7 | 21 | 0.132 | 0.646 | 0.764 | 0.701 |
| `f0_static` | 21 | 21 | 7 | 21 | 0.097 | 0.700 | 0.742 | 0.680 |
| `tone_prosody` | 21 | 21 | 7 | 21 | 0.563 | 0.261 | 0.780 | 0.703 |

Interpretation:

- The Yoruba pilot supports the Mandarin result that phonJSD can quantify
  lexical tone contrasts when the feature space encodes tone-relevant F0 cues.
- The result survives vowel-quality matching across seven vowel controls.
- `f0_static` and `f0_contour` both perform well, consistent with Yoruba as a
  register-tone system. `f0_contour` has the stronger classifier reference in
  this pilot, while `f0_static` remains a plausible compact register-tone
  parameterization.
- `tone_prosody` yields the largest JSD and lowest overlap, but it combines F0,
  energy, and duration cues and is therefore less clean as the primary
  phonological-tone specification.
- `energy_duration` is weaker by classifier AUC and should be treated as a
  control/sensitivity feature set rather than a primary tone cue.

## Observed Full-Clean Vowel-Matched Pattern

The full-clean Yoruba run confirms the pilot result on 2,422 recordings and
31,343 aligned tone tokens. The vowel-quality matched validation produced 88
metric rows: 22 contrasts for each of four feature sets across eight control
groups. Every row used `metric_mode = all_metrics`; no fallback rows were
needed.

Vowel-quality matched full-clean validation:

| Feature set | Rows | Contrasts | Control groups | All metrics | Median JSD | Median overlap | Median AUC | Median balanced accuracy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `energy_duration` | 22 | 22 | 8 | 22 | 0.083 | 0.726 | 0.667 | 0.606 |
| `f0_contour` | 22 | 22 | 8 | 22 | 0.140 | 0.644 | 0.761 | 0.691 |
| `f0_static` | 22 | 22 | 8 | 22 | 0.115 | 0.675 | 0.761 | 0.676 |
| `tone_prosody` | 22 | 22 | 8 | 22 | 0.282 | 0.512 | 0.812 | 0.733 |

Interpretation:

- The full-clean run supports the claim that phonJSD distinguishes Yoruba
  lexical tone when the feature space encodes F0.
- The result is robust to vowel-quality matching and is not driven by fallback
  metric mode.
- `f0_contour` and `f0_static` perform similarly by classifier AUC in the
  full-clean run. This is expected for a register-tone language where F0 level
  is a central cue; `f0_contour` remains the most portable primary tone
  specification across Mandarin and Yoruba.
- `tone_prosody` is the strongest cue bundle by JSD and classifier reference,
  but it combines F0, duration, and energy. It should be reported as a
  sensitivity/full-cue result rather than the primary phonological-tone
  parameterization.
- `energy_duration` remains weaker than F0-bearing feature sets and functions
  as a useful non-F0 control.

## Step 6: Alignment And Token Audit

Before treating a larger Yoruba run as paper-final, generate a deterministic
audit bundle from the extracted feature table:

```sh
/opt/miniconda3/envs/mfa/bin/python \
  analysis/validation/audit_yoruba_tone_alignment.py \
  --features analysis/validation/outputs/yoruba_slr86_tone_features_full_clean/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_alignment_audit_full_clean \
  --sample-size 120 \
  --max-per-category 40 \
  --max-per-category-control 6
```

The audit script writes:

- `yoruba_tone_alignment_audit_summary.csv`: corpus-wide counts, flag counts,
  and feature-quality medians;
- `yoruba_tone_alignment_audit_by_category.csv`: category-level duration,
  voicing, and F0 summaries;
- `yoruba_tone_alignment_audit_by_category_control.csv`: same summaries within
  tone-by-vowel-quality cells;
- `yoruba_tone_alignment_audit_flag_counts.csv`: extracted quality flags and
  derived audit flags;
- `yoruba_tone_alignment_audit_sample.csv`: stratified hand-audit sample with
  source WAV and TextGrid paths;
- `yoruba_tone_alignment_audit_issues.csv`: deterministic sample of rows with
  audit flags such as word-label mismatch, phone-label mismatch, low voiced
  proportion, or unexpectedly short duration.

For the hand audit, inspect the sampled `file` and `textgrid` pairs and verify:

- the TextGrid word interval corresponds to the token's orthographic word;
- the phone interval label matches `vowel_quality` or syllabic `n`;
- the interval is a plausible vowel or syllabic-nasal nucleus;
- voiced-frame coverage is adequate for F0 interpretation;
- any `word_label_mismatch` rows are rare and explainable.

Observed full-clean audit summary:

| Metric | Value |
| --- | ---: |
| tokens_total | 31343 |
| speakers_total | 36 |
| recordings_total | 2422 |
| control_groups_total | 8 |
| tokens_with_audit_flags | 502 |
| tokens_flagged_low_voiced_prop | 502 |
| word_label_mismatch | 0 |
| phone_label_mismatch | 0 |
| median_duration_ms | 100 |
| median_f0_voiced_prop | 1.0 |
| median_abs_f0_delta_st | 0.55 |

Observed full-clean audit by tone:

| Category | Tokens | Speakers | Recordings | Control groups | Low voiced-prop flags | Median duration ms | Median voiced prop | Median F0 mean ST |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| H | 12051 | 36 | 2411 | 8 | 81 | 100 | 1.0 | 90.14 |
| L | 8266 | 36 | 2236 | 8 | 344 | 100 | 1.0 | 87.46 |
| M | 11026 | 36 | 2384 | 7 | 77 | 90 | 1.0 | 88.48 |

The audit found no word-label or phone-label mismatches. The only derived audit
flag was low voiced-frame proportion in 502 tokens, or about 1.6% of retained
tokens. Because the extractor already skipped tokens with too few voiced frames,
these are lower-confidence retained tokens rather than alignment failures. The
issue sample should still be inspected before final paper submission.

## Step 7: Manual Audit Packet

Generate a manual review packet from the balanced sample and issue sample:

```sh
/opt/miniconda3/envs/mfa/bin/python \
  analysis/validation/prepare_yoruba_tone_manual_audit.py \
  --audit-dir analysis/validation/outputs/yoruba_slr86_tone_alignment_audit_full_clean \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_manual_audit_full_clean \
  --media-mode symlink
```

This writes:

- `yoruba_tone_manual_audit_items.csv`: one row per sampled token, with blank
  manual decision columns;
- `yoruba_tone_manual_audit_summary.csv`: counts by source, priority, tone, and
  speaker/recording coverage;
- `README.md`: coding instructions;
- optional `audio/` and `textgrid/` symlinks when `--media-mode symlink` is
  used.

Manual decisions should prioritize `review_priority = issue`, then inspect the
balanced sample rows. Fill:

- `manual_word_alignment_ok`: `yes`, `no`, or `unclear`;
- `manual_phone_alignment_ok`: `yes`, `no`, or `unclear`;
- `manual_interval_quality`: `good`, `minor_issue`, `bad`, or `unclear`;
- `manual_f0_usable`: `yes`, `no`, or `unclear`;
- `manual_decision`: `accept`, `exclude`, or `unclear`.

## Paper-Methods Notes

Report the following explicitly:

- Yoruba tone labels are orthographic labels derived from tone diacritics, not
  manually annotated acoustic tone labels.
- High tone is acute, low tone is grave, and mid tone is unmarked on
  vowel-bearing units.
- Syllabic nasal tone-bearing units are parsed explicitly.
- Recordings with bracketed annotations are flagged and should be excluded by
  default or reported separately.
- Speaker IDs are inferred from recording IDs as prefix plus middle field; this
  should be described as a corpus-internal speaker identifier.
- The parser-stage table is `pending_alignment` and contains no acoustic
  features; acoustic features are written only after MFA alignment.
- Paper-quality reporting should pair the full-clean aligned validation with a
  hand-audited alignment/token subset.

Suggested wording:

> Yoruba tone labels were derived from OpenSLR 86 orthographic transcripts after
> Unicode normalization. Acute-marked tone-bearing units were labeled high,
> grave-marked units low, and unmarked vowel-bearing units mid; syllabic nasal
> tone-bearing units were parsed separately. Bracketed recording annotations
> were retained as quality flags. A corpus-specific MFA acoustic model was
> trained from a grapheme-style dictionary, and orthographic tone labels were
> joined to aligned vowel or syllabic-nasal phone intervals before extracting
> speaker-normalized F0 trajectories. Vowel-quality matched validation was run
> on the full clean aligned subset, excluding annotation-flagged recordings.

## Limitations And Next Refinements

Completed:

- recording manifest;
- extracted-audio availability audit;
- speaker-code parsing;
- bracket-annotation flagging;
- orthographic H/M/L tone-unit parsing;
- vowel-quality controls;
- tone-unit count tables;
- 10-speaker MFA corpus materialization;
- corpus-specific MFA pilot training;
- phone-level F0 extraction;
- pooled and vowel-matched pilot JSD validation;
- full-clean scale-up command sequence;
- full-clean Yoruba alignment and feature extraction;
- full-clean vowel-matched JSD validation;
- deterministic alignment/token audit script and output contract;
- full-clean alignment/token audit summary;
- manual audit packet generation.

Not completed:

- manual inspection of the generated alignment/token audit sample;
- optional full-clean pooled validation summary, if wanted for supplement.

Needed refinements:

- validate the speaker-ID assumption against any OpenSLR metadata available;
- manually inspect a sample of parsed tone labels after Unicode normalization;
- audit a sample of generated dictionary entries and TextGrid phone alignments;
- decide whether nasal/oral status can be robustly derived from orthography;
- complete and archive manual decisions for the generated audit sample.

## Output Bundle To Archive

For the current parser stage, archive:

- `yoruba_slr86_manifest.csv`
- `yoruba_slr86_tone_units_pending_alignment.csv`
- `yoruba_slr86_parse_summary.csv`
- `yoruba_slr86_tone_unit_counts.csv`
- `git rev-parse HEAD`
- exact corpus root path used
- Python version
- command used to run the parser

For acoustic validation stages, also archive:

- alignment logs and TextGrids;
- aligned token-feature table;
- feature-set CSV;
- validation metrics/skips/summary CSVs;
- alignment/token audit summary and sample CSVs;
- alignment model and dictionary provenance;
- audio/F0 extraction package versions.

## Source

- OpenSLR 86: https://www.openslr.org/86/
