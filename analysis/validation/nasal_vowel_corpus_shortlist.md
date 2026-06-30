# Nasal/oral vowel corpus shortlist

## Recommendation

Use Multilingual LibriSpeech Portuguese as the primary P3 nasal/oral vowel
validation corpus.

Rationale:

- Portuguese has a phonemic nasal/oral vowel contrast, so this is a real target
  contrast rather than English-style coarticulatory nasalization.
- MLS Portuguese is large enough for controlled vowel-quality contrasts, but
  much smaller than MLS French.
- MLS preserves speaker-level structure, which lets us use speaker-normalized
  formants, spectral measures, and held-out speaker checks.
- The corpus is openly released under CC BY 4.0 and has a manageable download
  size: 9.3 GB for FLAC, 2.5 GB for Opus.

Use FLAC for the final acoustic validation. Opus is acceptable only for a quick
pipeline smoke test because nasality proxies may be sensitive to compression.

## Alternatives

### OpenSLR 86 Yoruba

Best for a dual-use Yoruba tone corpus and an exploratory nasal/oral pilot.

- Pros: WAV files, moderate download size, manually checked sentence
  transcriptions, anonymized file IDs with apparent speaker-like subfields, and
  Yoruba diacritics preserved in the line indexes.
- Cons: not pre-aligned; nasal/oral vowel labels require a validated Yoruba
  orthography-to-vowel mapping; lexical tone and nasalization can interact, so
  controls must be stricter than in Portuguese.
- Role: use as the first Yoruba tone validation corpus. Treat nasal/oral vowels
  as an audited pilot only until the mapping and alignment are validated.

### UCLA Phonetics Lab Archive

Best for a small, hand-checkable calibration set.

- Pros: direct WAV downloads, IPA word-list transcriptions, targeted examples
  for Portuguese and French nasal vowels, and related Yoruba examples that are
  useful for checking tone/nasalization extraction logic.
- Cons: short elicited word lists, vintage reel-tape recordings, small speaker
  counts, uneven metadata, and no built-in token-level alignment in the raw
  archive pages.
- Quality rule: use WAV only. Verify each recording's details page before
  analysis; at least the Yoruba archive details report 44.1 kHz, 16-bit WAV
  digitization from reel tape, while MP3 copies are low bitrate.
- Role: use this before MLS as a sanity check for labeling, forced alignment,
  and feature extraction. Do not use it as the main evidence for population-
  level JSD behavior.

If the relevant language is present in VoxAngeles, prefer that processed release
for the smoke test because it adds audited transcriptions, phone-level
alignments, and vowel measurements for part of the UCLA archive. Its
non-commercial license still makes it lower priority than MLS for the main
validation path.

### FLEURS Brazilian Portuguese

Best for a small pilot before downloading MLS.

- Pros: easy Hugging Face access, CC BY 4.0, about 4.1k rows for `pt_br`.
- Cons: much smaller; metadata expose gender, but not enough speaker structure
  for the same normalization and speaker-held-out checks we can run on MLS.

### VoxClamantis French or Portuguese

Best for a pre-aligned smoke test, not for the main validation claim.

- Pros: already includes phone-level alignments and vowel measures; includes
  French and Portuguese readings with large vowel-token counts.
- Cons: labels are automatically estimated; data are for non-commercial use;
  readings are religious text; should be treated as lower-confidence evidence.

### PFC French

Best linguistically if access is already available.

- Pros: designed for contemporary French phonology; public subset covers
  multiple surveys and speakers.
- Cons: current project site says new accounts are not being created; usage is
  restricted to non-commercial scientific/pedagogical work; extraction would be
  more corpus-specific than MLS.

### MLS French

Good later replication target.

- Pros: large open French read-speech corpus with many speakers and clear nasal
  vowel contrasts.
- Cons: much larger than needed for the first pass: 61 GB FLAC, 16 GB Opus.

## Development plan for MLS Portuguese

1. Download and unpack MLS Portuguese, preferably the FLAC archive.
2. Build a small manifest reader for utterance paths, transcripts, speaker IDs,
   and splits.
3. Use UCLA Portuguese/French examples as a hand-checkable calibration set for
   phone labeling and nasality feature extraction.
4. Run forced alignment on a bounded MLS pilot subset first, then scale up.
5. Convert aligned phones to nasal/oral categories:
   - `category`: `nasal` or `oral`;
   - `control_group`: vowel quality or nasal/oral pair family;
   - exclude or flag oral vowels adjacent to nasal consonants.
6. Extract vowel-window features with the same librosa-preferred pattern used in
   the SBCSAE stop-voicing adapter:
   - formants and duration;
   - MFCC means and optional trajectories;
   - low-frequency energy, spectral tilt, flatness, HNR/voicing strength;
   - anti-formant proxy features where robust.
7. Run the P0 validator with feature-set ablations:
   - `formant_duration`;
   - `mfcc_vowel`;
   - `nasality_proxy`;
   - `full_vowel_bundle`.
8. Report JSD against independent classifier references and include speaker- or
   split-held-out checks.

## Notes

Do not use SBCSAE as the P3 corpus. It can only provide oral vowels in nasal
versus non-nasal contexts, which is useful for coarticulation but does not
validate phonemic nasal/oral vowel contrasts.
