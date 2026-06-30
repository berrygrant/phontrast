# SBCSAE stop-voicing pilot

Date: 2026-06-30

## Inputs

- Audio: `/Volumes/Corpus_Studies/Corpora/English/SBCSAE/SBC001.wav` through
  `SBC060.wav`
- MFA TextGrids:
  `/Volumes/Corpus_Studies/Corpora/English/SBCSAE/SBCSAE_MFA_Aligned_Final`
- Adapter:
  `analysis/validation/prepare_sbcsae_stop_voicing_tokens.py`
- Validator:
  `analysis/validation/run_validation_metrics.R`

The aligned TextGrids contain `words` and `phones` tiers. They do not expose
speaker tiers, so the adapter currently uses `recording_id` as `speaker` and
sets `quality_flag = speaker_unavailable_recording_id_used`.

No acoustic values in this pilot are synthetic. Labels and time intervals are
read from the TextGrids; duration features are computed from those intervals;
and acoustic features are computed from the matching WAV samples in windows
anchored to those intervals. The only stochastic operation is balanced token
sampling with `--seed 2026`.

## Token inventory

Candidate stop tokens after requiring a following vowel within 200 ms:

| place | category | phone | n |
| --- | --- | --- | ---: |
| coronal | voiced | d | 23,721 |
| coronal | voiceless | t | 27,799 |
| dorsal | voiced | g | 7,618 |
| dorsal | voiceless | k | 12,948 |
| labial | voiced | b | 10,352 |
| labial | voiceless | p | 9,019 |

Pilot sample: 75 tokens per place/category, 450 tokens total.

## Feature sets

- `duration_context`: MFA stop/context/vowel interval timing.
- `spectral_stop`: pre-stop and stop-window energy, spectral-shape,
  low/high-ratio, and periodicity cues.
- `release_window`: 30 ms and 80 ms post-release/vowel-onset cues, including
  relative release/intensity changes.
- `f0_formants`: autocorrelation F0 and LPC F1/F2/F3 near vowel onset.
- `mfcc_context`: 13 MFCC means over the pre-stop through post-release context.
- `full_bundle`: all standardized pilot cues.

Audio features were extracted with `librosa` for this run. The current local
Anaconda `librosa` needed a temporary `pooch` dependency target plus
`NUMBA_DISABLE_JIT=1`; the adapter still has a `scipy` fallback.

Extracted acoustic parameters now include RMS intensity, zero-crossing rate,
spectral centroid/bandwidth/skewness/kurtosis, low/mid/high-band energy
proportions, low/high log energy ratio, spectral flatness, autocorrelation F0,
autocorrelation voicing strength, LPC formants, release-minus-stop deltas, and
MFCC means. True VOT and prevoicing onset are not yet estimated.

## Validation summary

All 18 place-by-feature-set contrasts completed; no skipped contrasts.
`full_bundle` rows used JSD/overlap fallback because Pillai fails from rank
deficiency in the 64-dimensional feature space. This is expected and now
explicitly recorded with `metric_mode = jsd_overlap_fallback`.

Mean by feature set:

| feature_set | mean JSD | mean classifier AUC | mean balanced accuracy | mean overlap |
| --- | ---: | ---: | ---: | ---: |
| duration_context | 0.066 | 0.636 | 0.602 | 0.794 |
| spectral_stop | 0.569 | 0.654 | 0.617 | 0.253 |
| release_window | 0.724 | 0.632 | 0.615 | 0.128 |
| f0_formants | 0.060 | 0.486 | 0.490 | 0.810 |
| mfcc_context | 0.886 | 0.568 | 0.551 | 0.041 |
| full_bundle | 1.000 | 0.562 | 0.556 | 0.000 |

## Interpretation

The pilot supports using SBCSAE for a stop-voicing validation module, but the
high-dimensional KDE estimates should not be interpreted at face value yet. The
`mfcc_context` and `full_bundle` JSD values are near-maximal while held-out
classifier AUC is only modest. That mismatch is consistent with sparse
high-dimensional KDE over-separation at this pilot sample size.

The most credible first-pass acoustic result remains `spectral_stop`: it has
the best mean classifier AUC and a plausible JSD/overlap profile. Expanded
`release_window` features now carry signal as well. `f0_formants` alone are weak
in this pilot, which is plausible for English stop voicing because onset F0 is a
secondary cue and the current F0/formant estimates are simple window-level
proxies.

## Next decision

Add true VOT/prevoicing estimation before treating this as a substantive
validation result. The current adapter uses MFA phone intervals and short
windows, so it does not identify burst onset, voicing onset, or prevoicing
onset. Those event measurements are the missing bridge to the AutoVOT/DeepVOT
literature and should be the next increment.
