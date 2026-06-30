# Yoruba SLR86 development plan

## Corpus role

Use OpenSLR 86 as the active lexical-tone validation corpus while the larger
Mandarin corpus is deferred.

OpenSLR 86 is a crowdsourced Yoruba sentence corpus with manually checked WAV
recordings, separate female and male line-index files, and tone-marked Yoruba
orthographic transcriptions. The line-index files provide anonymized file IDs;
the file IDs appear to contain sex and speaker-like subfields, but that should
be verified before treating the middle code as `speaker`. The download footprint
is moderate: 462 MB for the female audio archive and 445 MB for the male audio
archive.

This replaces BibleTTS as the first Yoruba target. BibleTTS remains useful as a
later replication source, but SLR86 is better for this project because its
transcripts preserve Yoruba diacritics and speaker/file IDs directly.

## First-pass scope

### Tone

SLR86 is immediately useful for tone validation.

Token labels should be derived from tone marks on vowel-bearing units:

- acute accent: high tone;
- grave accent: low tone;
- unmarked vowel: mid tone;
- syllabic nasal/tone-bearing nasal should be handled explicitly rather than
  silently discarded.

The first validation pass should compare H/M/L tone categories within matched
segmental controls where possible:

- same vowel quality;
- same nasal/oral status if recoverable;
- same local consonant context when token counts allow;
- speaker- or sex-balanced subsets.

### Nasal/oral vowels

SLR86 is an exploratory nasal/oral source, not yet the main P3 corpus.

Yoruba has orthographic material that can expose nasalized vowels and syllabic
nasals, but we should validate the mapping before treating the corpus as a
nasal/oral vowel benchmark. The nasal-vowel pass should start as a small,
audited subset with hand-checkable tokens.

## Implementation steps

1. Download `line_index_female.tsv`, `line_index_male.tsv`,
   `annotation_info.txt`, `yo_ng_female.zip`, and `yo_ng_male.zip`.
2. Build `prepare_yoruba_slr86_tone_tokens.py` to:
   - parse line-index TSVs;
   - locate WAV files by file ID;
   - remove or flag bracketed annotations such as `[breath]`, `[snap]`, and
     `[external]`;
   - normalize Unicode to NFC;
   - emit a manifest with `recording_id`, provisional `speaker`, `sex`,
     `transcript`, and audio path.
3. Add a Yoruba orthography parser for tone-bearing units:
   - vowel base;
   - tone category;
   - nasal/oral flag where recoverable;
   - word position and neighboring characters.
4. Align text to audio on a pilot subset.
   - Preferred: MFA-style forced alignment if a Yoruba model/dictionary is
     available or can be derived.
   - Fallback: segment-level F0 validation on words or utterance regions only
     for smoke testing, not final JSD claims.
5. Extract tone features over aligned vowel-bearing units:
   - speaker-normalized F0 mean/min/max/range/slope;
   - time-normalized F0 trajectory points;
   - duration, intensity, and voiced-frame coverage.
6. Run the P0 validator with tone feature sets:
   - `f0_static`;
   - `f0_trajectory`;
   - `tone_prosody`.
7. Report pairwise H/M/L JSD against held-out classifier AUC/balanced accuracy,
   with skipped-token diagnostics for noisy annotations and poor F0 tracking.

## Current local status

`prepare_yoruba_slr86_tone_tokens.py` now builds the SLR86 manifest and
orthographic tone-unit inventory from the local corpus root:

- local root: `/Volumes/Corpus_Studies/Corpora/Yoruba`;
- manifest output:
  `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv`;
- tone-unit output:
  `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv`;
- summary output:
  `analysis/validation/outputs/yoruba_slr86/yoruba_slr86_parse_summary.csv`.

The current parse produced 3,583 manifest rows and 51,412 orthographic
tone-bearing units. Tone labels come only from the transcript diacritics: high
from acute marks, low from grave marks, and mid from unmarked vowel-bearing
units. No segment timings or acoustic values have been inferred.

Speaker IDs should use the recording prefix plus the middle file-ID field
(`yof_06136`, `yom_06136`, etc.). The raw middle field overlaps across sex, so
using it alone collapses distinct male and female speakers. The manifest keeps
both `speaker` and `speaker_code` so downstream analyses can audit this choice.

The female and male WAV directories are now fully extracted in the current NAS
copy. The manifest still distinguishes `audio_location = extracted_wav` from
`audio_location = zip_member` for auditability, but the current parse reports all
3,583 recordings as extracted WAV files with no missing audio.

## Quality checks

- Exclude or separately report lines with `[abrupt]`.
- Keep `[external]`, `[breath]`, `[hesitation]`, and `[snap]` rows in the
  manifest but flag them; default validation should exclude them until we know
  their impact.
- Verify that tone labels are not being inferred after Unicode decomposition
  has stripped diacritics.
- Audit a small sample of aligned tokens before running the full corpus.

## Source

- OpenSLR 86: https://www.openslr.org/86/
