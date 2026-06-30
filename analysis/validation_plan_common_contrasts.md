# JSD validation plan for additional phonological contrasts

This plan extends the existing vowel/MFCC validation architecture to three
contrast families: lexical tone, stop voicing, and nasal/oral vowels. The goal is
not to change the core JSD estimator, but to add reproducible validation modules
that prepare token-level acoustic feature tables and run the same metric
comparison logic across new phonological domains.

## Prioritization checklist

| Priority | Contrast | First target | Why this order | Main blocker | Done when |
| --- | --- | --- | --- | --- | --- |
| P0 | Shared validation scaffold | Domain-neutral runner | Prevents one-off scripts and lets each contrast reuse metric, ablation, and reporting code | Need to formalize a common token schema | A script can read any compliant token table, compute JSD/classical metrics, run ablations, and write skipped-contrast notes |
| P1 | Stop voicing beyond VOT | English /p b/, /t d/, /k g/ with place-controlled pairs | The cited VOT papers give a concrete cue inventory and event structure; it is the cleanest first test of consonantal contrasts | Need aligned stop-vowel windows and either measured or automatically estimated events | VOT-only, non-VOT cue-bundle, and full-cue JSD are compared against an independent separability reference |
| P2 | Yoruba lexical tone | OpenSLR 86 Yoruba H/M/L vowel-tone contrasts | Moderate download size, tone-marked transcripts, and manually checked WAV recordings make it the fastest tone path while larger Mandarin data are deferred | Need tone-bearing-unit parsing plus alignment | Yoruba tone contrasts run through the shared tone pipeline with speaker-normalized F0 trajectories and vowel/context controls where possible |
| P3 | Nasal vs oral vowels | French or Portuguese nasal/oral vowel pairs | Extends JSD to a well-known vowel contrast outside F1/F2-only oral vowel overlap | Nasality labels and oral controls are confounded by vowel quality and nasal-consonant context | Nasal/oral JSD is estimated from MFCC/formant/nasality-proxy features with oral tokens near nasal consonants excluded or modeled separately |
| P4 | Mandarin lexical tone replication | AISHELL-1 tone pairs, vowel/syllable nucleus features | Large multi-speaker corpus and tone labels derivable from pinyin/lexicon resources; deferred until a fast connection is available | Large corpus download plus syllable or vowel-nucleus alignment and robust F0 tracking | Mandarin tone-pair JSD matrix is produced with speaker-normalized F0 trajectories and same-rime/same-syllable controls where possible |

## Common validation contract

Each module should produce a token-level table with these required columns:

- `token_id`
- `domain`: `tone`, `stop_voicing`, or `nasal_vowels`
- `language`
- `source_corpus`
- `speaker`
- `category`: the phonological category being contrasted
- `control_group`: matched comparison bin, such as toneless syllable, stop place,
  vowel quality, word position, or speaker-normalized stratum
- numeric acoustic feature columns

Recommended optional columns:

- `file`, `start`, `end`
- `word`, `phone`, `syllable`, `rime`
- `following_phone`, `preceding_phone`
- `measurement_method`
- `quality_flag`

The validation runner should:

- enumerate binary contrasts within `domain` and optional `control_group`;
- enforce minimum token counts per category;
- compute `compare_overlap_metrics()` and `estimate_jsd()` with explicit KDE
  settings;
- compute at least one non-KDE reference separability measure, preferably
  held-out classifier AUC or balanced accuracy;
- run feature-set ablations;
- write both successful estimates and skipped-contrast diagnostics.

## Contrast-specific plans

### P1 stop voicing beyond VOT

Feature sets:

- `vot_only`: VOT, prevoicing duration/proportion if available.
- `event_durations`: closure duration, burst duration, aspiration duration,
  vowel-onset lag.
- `spectral_voicing`: low/high-frequency energy, spectral flatness/Wiener
  entropy, zero-crossing rate, autocorrelation or voicing probability, HNR.
- `vowel_onset`: onset F0, F1, intensity, and short trajectories after release.
- `full_cue_bundle`: all usable stop voicing cues.
- `mfcc_window`: MFCC means/trajectories over release plus vowel onset.

Adopter-facing parameter inventory:

- Time/event measures: stop interval duration, closure duration if recoverable,
  burst/release duration, aspiration duration, VOT, negative VOT/prevoicing
  duration, following-vowel lag, and following-vowel duration.
- Periodicity/voicing measures: autocorrelation voicing strength, voicing-frame
  proportion, HNR where available, F0 near vowel onset, and prevoicing
  proportion before release.
- Spectral measures: RMS intensity, zero-crossing rate, spectral centroid,
  spectral bandwidth, spectral skewness, spectral kurtosis, spectral flatness,
  low/mid/high-band energy proportions, and low/high energy ratios.
- Vowel-onset measures: onset F0, F1/F2/F3, early intensity trajectory, and
  release-to-vowel relative changes.
- General representations: MFCC means/trajectories and optional learned
  embeddings. These should be reported separately from interpretable phonetic
  cue sets because high-dimensional density estimation can over-separate at
  small sample sizes.

Validation design:

- Compare /p b/, /t d/, and /k g/ separately.
- Stratify by word position and speaker where data allow.
- Report whether the full cue bundle improves the relationship with held-out
  separability compared with VOT alone.
- Keep VOT-only results as a baseline, not as the full definition of voicing.

Source notes:

- Sonderegger & Keshet's AutoVOT work frames VOT as burst onset to following
  voicing onset and uses fine-grained spectral, autocorrelation, pitch, and
  voicing-detector features.
- Adi et al. extend automatic measurement to prevoicing, burst duration, and
  positive/negative VOT with frame labels for silence, prevoicing, burst, and
  vowel.

### P2 Yoruba lexical tone

Feature sets:

- `f0_static`: mean, min, max, range, and slope of speaker-normalized F0.
- `f0_trajectory`: 5-10 time-normalized F0 points over the voiced nucleus.
- `tone_prosody`: F0 trajectory plus duration, intensity, and voicing coverage.
- `ssl_optional`: learned frame or segment embeddings if available later.

Validation design:

- Start with Yoruba H/M/L pairwise contrasts.
- Prefer matched comparisons within the same vowel quality and local segmental
  context.
- Normalize F0 per speaker/session in semitones or z-scores.
- Flag tokens with low voiced-frame coverage or unreliable pitch tracking.
- Report both pooled and speaker/grouped JSD when token counts permit.

Data path:

- OpenSLR 86 is the preferred first tone corpus while Mandarin data are
  deferred. It provides manually checked WAV sentence recordings, separate
  female and male line indexes, and tone-marked Yoruba orthographic transcripts
  under CC BY-SA 4.0. The audio archives are moderate-sized: 462 MB for female
  speakers and 445 MB for male speakers. Its line indexes expose anonymized file
  IDs; verify the apparent speaker-like subfield before using it for
  speaker-normalized analyses.
- The first adapter should parse tone marks from Unicode-normalized Yoruba
  orthography, flag bracketed annotations such as `[breath]`, `[snap]`,
  `[external]`, `[hesitation]`, and `[abrupt]`, then align vowel-bearing units
  before extracting F0 trajectories.

### P3 nasal and non-nasal vowels

Feature sets:

- `formant_duration`: F1/F2/F3, bandwidths where available, and vowel duration.
- `mfcc_vowel`: MFCC means and optionally time-normalized trajectories over the
  vowel nucleus.
- `nasality_proxy`: spectral tilt, low-frequency nasal resonance proxies, HNR,
  and anti-formant-related measures where robust extraction is available.
- `full_vowel_bundle`: formants, duration, MFCCs, and nasality proxies.

Validation design:

- Use phonemic nasal vowels, not only coarticulatory nasalization.
- Compare language-specific nasal/oral pairs and document when oral and nasal
  categories are not quality-matched.
- Exclude or separately label oral vowels adjacent to nasal consonants.
- Stratify by speaker and vowel quality where data allow.

Candidate data:

- OpenSLR 86 Yoruba is an exploratory nasal/oral source after tone parsing is
  in place. It has WAV files and tone-marked Yoruba transcriptions, but it is
  not pre-aligned and nasal/oral vowel labels require a validated Yoruba
  orthography-to-vowel mapping. Treat this as a pilot, not the main P3 corpus,
  until labels and alignment are audited.
- UCLA Phonetics Lab Archive Portuguese and French are good calibration sources
  because the raw archive exposes targeted word lists, IPA transcriptions, and
  WAV files for nasal-vowel examples. Treat them as smoke-test/calibration data,
  not the main validation corpus, because they are short elicited recordings
  with small speaker counts, uneven metadata, and reel-tape provenance. Use WAV
  only and verify recording details per file. If the relevant language is
  available in VoxAngeles, prefer that processed release for the smoke test
  because it adds phone-level alignments and vowel measurements.
- MLS Portuguese is the preferred first corpus. It has a phonemic nasal/oral
  vowel system, open licensing, speaker structure, and a smaller footprint than
  MLS French. Use FLAC for final acoustic work and Opus only for a smoke test.
- FLEURS Brazilian Portuguese is a useful small pilot if we want to validate the
  ingestion and alignment pipeline before downloading MLS.
- VoxClamantis French or Portuguese can be a broad pre-aligned smoke-test source
  because it provides aligned phonetic resources and vowel measures, but its
  automatically estimated labels and non-commercial terms make it lower
  confidence for the main validation claim.
- French PFC is linguistically attractive if access is already available, but it
  is not the first target because new account creation is currently unavailable
  and usage is restricted to non-commercial scientific/pedagogical work.

### P4 Mandarin lexical tone replication

Feature sets and validation design should mirror the Yoruba tone module, but
with Mandarin T1/T2/T3/T4 categories and same-rime or same-syllable controls.

Data path:

- AISHELL-1 remains the preferred Mandarin corpus because it is large and
  multi-speaker. Its lexicon resource exposes pinyin entries with tone numbers,
  which can drive tone labels after alignment.
- Mandarin ingestion is deferred until the corpus can be downloaded over a fast
  connection. This should not block development of the shared tone feature
  extraction and validation runner.

## Development plan

1. Add `analysis/validation/` with a README, input schema, and reusable runner.
2. Implement metric computation on a precomputed token-feature table first.
3. Add feature-set ablation support and held-out classifier reference metrics.
4. Implement the stop-voicing feature table contract and run a small pilot.
5. Add Yoruba tone ingestion and F0 trajectory extraction.
6. Add nasal/oral vowel ingestion and feature extraction.
7. Add Mandarin tone ingestion as a replication once the large corpus is
   available.
8. Promote stable pieces into vignettes or package examples only after the
   analysis scripts have produced reproducible outputs.

## References and source leads

- Sonderegger & Keshet, AutoVOT PDF:
  https://people.linguistics.mcgill.ca/~morgan/interspeechVot.pdf
- Adi et al. 2016, DeepVOT / prevoicing:
  https://www.isca-archive.org/interspeech_2016/adi16_interspeech.html
- AISHELL-1:
  https://arxiv.org/abs/1709.05522
- Lexical tone quantization study using Mandarin and Yoruba:
  https://arxiv.org/abs/2604.07467
- BibleTTS:
  https://arxiv.org/abs/2207.03546
- VoxClamantis:
  https://arxiv.org/abs/2005.13962
- PFC project:
  https://www.projet-pfc.net/
