# Yoruba OpenSLR 86 Tone Workflow

This document records the reproducible Yoruba OpenSLR 86 workflow for the
phonJSD common-contrast validation. It covers the current completed stage:
manifest construction and orthographic tone-unit parsing. Acoustic alignment and
F0 extraction are the next development stage.

## Goal

Use OpenSLR 86 Yoruba as the first non-Mandarin lexical tone validation corpus.
The corpus provides tone-marked Yoruba transcripts and WAV recordings, allowing
H/M/L tone labels to be parsed directly from orthography before acoustic
alignment.

The current workflow produces:

- a recording-level manifest;
- an orthographic tone-bearing-unit table marked `pending_alignment`;
- parse summaries and count tables for audit;
- a documented plan for alignment, F0 extraction, and validation.

It does not yet produce acoustic tone features or JSD-ready token tables.

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

Observed local 10-speaker pilot summary:

| Metric | Value |
| --- | ---: |
| recordings_selected | 200 |
| speakers_selected | 10 |
| recordings_female | 100 |
| recordings_male | 100 |
| audio_copied | 200 |
| transcript_written | 200 |
| recordings_missing_audio | 0 |
| dictionary_words | 598 |

The pilot outputs are:

- `/tmp/phonjsd_yoruba_mfa_pilot_10/corpus`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_manifest.csv`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_summary.csv`
- `/tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_dictionary.txt`

## Step 3: Alignment Plan

The next stage is to create aligned intervals for tone-bearing units or vowel
nuclei. No final alignment method has been committed yet.

Preferred direction:

1. Build or identify a Yoruba pronunciation dictionary.
2. Align a small pilot subset with MFA or another phone-level aligner.
3. Audit word/phone alignment manually on a sample of files.
4. Map each parsed tone-bearing orthographic unit to a vowel or syllabic-nasal
   interval.
5. Extract F0 features over those intervals.

If a robust Yoruba acoustic model/dictionary is not available, a word-level
pilot can be used only as a smoke test. Paper-quality validation should use
syllable or vowel-nucleus intervals, not utterance-level or whole-word windows.

Validate the pilot corpus once an acoustic model choice is available:

```sh
mfa validate \
  /tmp/phonjsd_yoruba_mfa_pilot_10/corpus \
  /tmp/phonjsd_yoruba_mfa_pilot_10/yoruba_mfa_pilot_dictionary.txt \
  <yoruba_or_trainable_acoustic_model> \
  --clean \
  --num_jobs 8
```

If no usable pretrained Yoruba acoustic model is available, use the generated
dictionary to train a corpus-specific MFA model before extracting phone or vowel
intervals. The 10-speaker pilot is the smoke test for whether that path is
viable.

## Step 4: Proposed Yoruba Tone Feature Extraction

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

## Step 5: Proposed Validation Commands

Once a Yoruba acoustic token-feature table exists, run pooled and matched
validation separately.

Pooled sanity check:

```sh
mamba run -n phonjsd-r Rscript analysis/validation/run_validation_metrics.R \
  --input analysis/validation/outputs/yoruba_slr86_tone_features/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_pooled \
  --domain tone \
  --control-col pooled_tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features/yoruba_slr86_tone_feature_sets.csv \
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
  --input analysis/validation/outputs/yoruba_slr86_tone_features/yoruba_slr86_tone_features.csv \
  --out-dir analysis/validation/outputs/yoruba_slr86_tone_validation_vowel_matched \
  --domain tone \
  --feature-sets analysis/validation/outputs/yoruba_slr86_tone_features/yoruba_slr86_tone_feature_sets.csv \
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
- The current table is `pending_alignment` and contains no acoustic features.
- Paper-quality JSD validation requires aligned vowel or syllable intervals and
  speaker-normalized F0 features.

Suggested wording:

> Yoruba tone labels were derived from OpenSLR 86 orthographic transcripts after
> Unicode normalization. Acute-marked tone-bearing units were labeled high,
> grave-marked units low, and unmarked vowel-bearing units mid; syllabic nasal
> tone-bearing units were parsed separately. Bracketed recording annotations
> were retained as quality flags. The initial parser produces a pending-alignment
> token inventory; acoustic validation will join these labels to aligned
> vowel- or syllable-level intervals before extracting speaker-normalized F0
> trajectories.

## Limitations And Next Refinements

Completed:

- recording manifest;
- extracted-audio availability audit;
- speaker-code parsing;
- bracket-annotation flagging;
- orthographic H/M/L tone-unit parsing;
- vowel-quality controls;
- tone-unit count tables.

Not completed:

- forced alignment;
- phone/vowel nucleus mapping;
- acoustic F0 extraction;
- JSD validation on Yoruba acoustic features.

Needed refinements:

- validate the speaker-ID assumption against any OpenSLR metadata available;
- manually inspect a sample of parsed tone labels after Unicode normalization;
- define a Yoruba dictionary/alignment strategy;
- decide whether nasal/oral status can be robustly derived from orthography;
- create a Yoruba-specific acoustic feature extractor after alignment is stable;
- run pooled and vowel-matched JSD validation;
- add a small hand-audited alignment/token sample for paper quality control.

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

For the future acoustic validation stage, also archive:

- alignment logs and TextGrids;
- aligned token-feature table;
- feature-set CSV;
- validation metrics/skips/summary CSVs;
- alignment model and dictionary provenance;
- audio/F0 extraction package versions.

## Source

- OpenSLR 86: https://www.openslr.org/86/
