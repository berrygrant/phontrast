#!/usr/bin/env python3
"""Materialize an OpenSLR 86 Yoruba corpus sample for MFA pilot alignment.

The Yoruba parser writes a recording manifest and orthographic tone-unit labels.
MFA expects a corpus directory containing WAV files and matching transcript
files, plus a pronunciation dictionary. This script samples speakers from the
manifest, writes an MFA-ready corpus, and creates a deterministic grapheme-style
dictionary from tone-stripped Yoruba orthography.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import re
import shutil
import unicodedata
import wave
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_MANIFEST = Path("analysis/validation/outputs/yoruba_slr86/yoruba_slr86_manifest.csv")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/yoruba_slr86_mfa_pilot")

ANNOTATION_RE = re.compile(r"\[([^\]]+)\]")
COMBINING_ACUTE = "\u0301"
COMBINING_GRAVE = "\u0300"
COMBINING_DOT_BELOW = "\u0323"
RIGHT_SINGLE_QUOTE = "\u2019"
LETTERLIKE_JOINERS = {"'", RIGHT_SINGLE_QUOTE}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Yoruba manifest produced by prepare_yoruba_slr86_tone_tokens.py.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-speakers", type=int, default=10, help="0 means no speaker cap.")
    parser.add_argument(
        "--max-utterances-per-speaker",
        type=int,
        default=20,
        help="0 means no utterance cap within selected speakers.",
    )
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument(
        "--balance-sex",
        action="store_true",
        help="Select speakers evenly across manifest sex labels when possible.",
    )
    parser.add_argument(
        "--include-recommended-exclude",
        action="store_true",
        help="Include annotation-flagged recordings. Default excludes them.",
    )
    parser.add_argument(
        "--skip-audio-validation",
        action="store_true",
        help="Do not pre-check WAV headers before sampling.",
    )
    parser.add_argument(
        "--mode",
        choices=("copy", "symlink", "manifest-only"),
        default="copy",
        help="Use copy for scratch alignment, symlink for stable local mounts.",
    )
    parser.add_argument(
        "--source-root-rewrite",
        default="",
        metavar="FROM=TO",
        help="Rewrite manifest file paths before reading audio, for cross-machine mounts.",
    )
    parser.add_argument(
        "--transcript-extension",
        default=".lab",
        help="Transcript suffix to write. Include the leading dot.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Select rows and write manifests without touching audio/transcript files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing files in the output corpus.",
    )
    return parser.parse_args()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_rewrite(value: str) -> tuple[str, str] | None:
    if not value:
        return None
    if "=" not in value:
        raise ValueError("--source-root-rewrite must have the form FROM=TO")
    old, new = value.split("=", 1)
    if not old:
        raise ValueError("--source-root-rewrite FROM cannot be empty")
    return old, new


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def strip_annotations(text: str) -> str:
    clean = ANNOTATION_RE.sub(" ", text)
    return re.sub(r"\s+", " ", clean).strip()


def strip_tone_marks(text: str) -> str:
    decomposed = unicodedata.normalize("NFD", text)
    kept = [ch for ch in decomposed if ch not in {COMBINING_ACUTE, COMBINING_GRAVE}]
    return unicodedata.normalize("NFC", "".join(kept))


def normalize_mfa_text(text: str) -> str:
    text = strip_tone_marks(strip_annotations(text)).lower().replace(RIGHT_SINGLE_QUOTE, "'")
    out: list[str] = []
    for ch in text:
        category = unicodedata.category(ch)
        if category.startswith("L") or category.startswith("M"):
            out.append(ch)
        elif ch in LETTERLIKE_JOINERS:
            continue
        else:
            out.append(" ")
    return re.sub(r"\s+", " ", "".join(out)).strip()


def phone_for_char(ch: str) -> str:
    decomposed = unicodedata.normalize("NFD", ch)
    base_chars = [part for part in decomposed if unicodedata.category(part).startswith("L")]
    if not base_chars:
        return ""
    base = base_chars[0].lower()
    has_dot = COMBINING_DOT_BELOW in decomposed
    if has_dot and base in {"e", "o", "s"}:
        return f"{base}_dot"
    return base


def word_to_phones(word: str) -> list[str]:
    phones: list[str] = []
    for ch in word:
        phone = phone_for_char(ch)
        if phone:
            phones.append(phone)

    collapsed: list[str] = []
    i = 0
    while i < len(phones):
        if i + 1 < len(phones) and phones[i] == "g" and phones[i + 1] == "b":
            collapsed.append("gb")
            i += 2
        elif i + 1 < len(phones) and phones[i] == "k" and phones[i + 1] == "p":
            collapsed.append("kp")
            i += 2
        else:
            collapsed.append(phones[i])
            i += 1
    return collapsed


def source_path(row: dict[str, str], rewrite: tuple[str, str] | None) -> Path:
    raw_path = row.get("file", "")
    if rewrite:
        old, new = rewrite
        if raw_path.startswith(old):
            raw_path = new + raw_path[len(old) :]
    return Path(raw_path)


def readable_wav_status(path: Path) -> str:
    if not path.exists():
        return "missing_source_audio"
    if path.stat().st_size <= 0:
        return "empty_source_audio"
    try:
        with wave.open(str(path), "rb") as handle:
            if handle.getframerate() <= 0 or handle.getnframes() <= 0:
                return "invalid_source_audio_header"
    except Exception as exc:
        return f"unreadable_source_audio:{type(exc).__name__}"
    return "readable"


def eligible_rows(
    rows: list[dict[str, str]],
    include_recommended_exclude: bool,
) -> tuple[list[dict[str, str]], Counter[str]]:
    out: list[dict[str, str]] = []
    skipped: Counter[str] = Counter()
    for row in rows:
        if not include_recommended_exclude and parse_bool(row.get("recommended_exclude", "")):
            skipped["recommended_exclude"] += 1
            continue
        if not parse_bool(row.get("audio_available", "")):
            skipped["manifest_audio_unavailable"] += 1
            continue
        transcript = normalize_mfa_text(row.get("transcript_clean") or row.get("transcript", ""))
        if not transcript:
            skipped["empty_mfa_transcript"] += 1
            continue
        if not row.get("speaker", ""):
            skipped["missing_speaker"] += 1
            continue
        out.append(row)
    return out, skipped


def balanced_speakers(by_speaker: dict[str, list[dict[str, str]]], max_speakers: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    by_sex: dict[str, list[str]] = defaultdict(list)
    for speaker, rows in by_speaker.items():
        sex = rows[0].get("sex", "unspecified") or "unspecified"
        by_sex[sex].append(speaker)

    for speakers in by_sex.values():
        speakers.sort()
        rng.shuffle(speakers)

    sexes = sorted(by_sex)
    selected: list[str] = []
    while len(selected) < max_speakers and any(by_sex.values()):
        for sex in sexes:
            if len(selected) >= max_speakers:
                break
            if by_sex[sex]:
                selected.append(by_sex[sex].pop(0))
    return selected


def select_rows(
    rows: list[dict[str, str]],
    max_speakers: int,
    max_utterances_per_speaker: int,
    seed: int,
    balance_sex: bool,
    rewrite: tuple[str, str] | None,
    validate_audio: bool,
) -> tuple[list[dict[str, str]], Counter[str]]:
    by_speaker: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_speaker[row["speaker"]].append(row)

    rng = random.Random(seed)
    speakers = sorted(by_speaker)
    if max_speakers > 0:
        if balance_sex:
            speakers = balanced_speakers(by_speaker, max_speakers, seed)
        else:
            rng.shuffle(speakers)
            speakers = speakers[:max_speakers]

    selected: list[dict[str, str]] = []
    skipped: Counter[str] = Counter()
    for speaker in speakers:
        speaker_rows = by_speaker[speaker][:]
        rng.shuffle(speaker_rows)
        speaker_selected = 0
        for row in speaker_rows:
            if max_utterances_per_speaker > 0 and speaker_selected >= max_utterances_per_speaker:
                break
            if validate_audio:
                audio_status = readable_wav_status(source_path(row, rewrite))
                if audio_status != "readable":
                    skipped[audio_status] += 1
                    continue
            selected.append(row)
            speaker_selected += 1

    return sorted(selected, key=lambda row: (row.get("sex", ""), row["speaker"], row["recording_id"])), skipped


def ensure_clean_target(path: Path, overwrite: bool) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if not overwrite:
        raise FileExistsError(f"{path} already exists; pass --overwrite to replace it")
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def materialize_audio(src: Path, dst: Path, mode: str, overwrite: bool) -> str:
    if mode == "manifest-only":
        return "manifest_only"
    if not src.exists():
        return "missing_source_audio"
    ensure_clean_target(dst, overwrite)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if mode == "copy":
        shutil.copy2(src, dst)
        return "copied"
    os.symlink(src, dst)
    return "symlinked"


def write_text(path: Path, text: str, overwrite: bool, dry_run: bool) -> str:
    if dry_run:
        return "dry_run"
    ensure_clean_target(path, overwrite)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n", encoding="utf-8")
    return "written"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_dictionary(path: Path, words: set[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8") as handle:
        for word in sorted(words):
            phones = word_to_phones(word)
            if not phones:
                continue
            handle.write(f"{word}\t{' '.join(phones)}\n")
            n += 1
    return n


def summary_rows(
    rows: list[dict[str, object]],
    dictionary_words: int,
    skipped_eligible: Counter[str],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        out.append({"metric": metric, "value": value})

    add("recordings_selected", len(rows))
    add("speakers_selected", len({row["speaker"] for row in rows}))
    for sex, n in sorted(Counter(str(row["sex"]) for row in rows).items()):
        add(f"recordings_{sex}", n)
    for status, n in sorted(Counter(str(row["audio_status"]) for row in rows).items()):
        add(f"audio_{status}", n)
    for status, n in sorted(Counter(str(row["transcript_status"]) for row in rows).items()):
        add(f"transcript_{status}", n)
    add("recordings_missing_audio", sum(1 for row in rows if row["audio_status"] == "missing_source_audio"))
    for status, n in sorted(skipped_eligible.items()):
        add(f"recordings_skipped_{status}", n)
    add("dictionary_words", dictionary_words)
    return out


def main() -> None:
    args = parse_args()
    rewrite = parse_rewrite(args.source_root_rewrite)
    if not args.transcript_extension.startswith("."):
        raise ValueError("--transcript-extension must include the leading dot")

    rows, skipped_eligible = eligible_rows(
        read_manifest(args.manifest),
        args.include_recommended_exclude,
    )
    selected, skipped_selected = select_rows(
        rows,
        max_speakers=args.max_speakers,
        max_utterances_per_speaker=args.max_utterances_per_speaker,
        seed=args.seed,
        balance_sex=args.balance_sex,
        rewrite=rewrite,
        validate_audio=not args.skip_audio_validation,
    )
    skipped_eligible.update(skipped_selected)

    corpus_dir = args.out_dir / "corpus"
    if args.overwrite and not args.dry_run and corpus_dir.exists():
        shutil.rmtree(corpus_dir)
    output_rows: list[dict[str, object]] = []
    dictionary_words: set[str] = set()
    for row in selected:
        src = source_path(row, rewrite)
        speaker_dir = corpus_dir / row["speaker"]
        dst_audio = speaker_dir / f"{row['recording_id']}{src.suffix or '.wav'}"
        dst_lab = speaker_dir / f"{row['recording_id']}{args.transcript_extension}"
        transcript = normalize_mfa_text(row.get("transcript_clean") or row.get("transcript", ""))
        dictionary_words.update(transcript.split())

        if args.dry_run:
            audio_status = "dry_run"
            transcript_status = "dry_run"
        elif args.mode == "manifest-only":
            audio_status = "manifest_only"
            transcript_status = "manifest_only"
        else:
            audio_status = materialize_audio(src, dst_audio, args.mode, args.overwrite)
            transcript_status = write_text(dst_lab, transcript, args.overwrite, args.dry_run)

        output_rows.append(
            {
                "recording_id": row["recording_id"],
                "sex": row.get("sex", ""),
                "speaker": row["speaker"],
                "source_audio_file": str(src),
                "mfa_audio_file": str(dst_audio),
                "mfa_transcript_file": str(dst_lab),
                "transcript": row.get("transcript", ""),
                "transcript_clean": row.get("transcript_clean", ""),
                "mfa_transcript": transcript,
                "audio_status": audio_status,
                "transcript_status": transcript_status,
            }
        )

    manifest_out = args.out_dir / "yoruba_mfa_pilot_manifest.csv"
    summary_out = args.out_dir / "yoruba_mfa_pilot_summary.csv"
    dictionary_out = args.out_dir / "yoruba_mfa_pilot_dictionary.txt"
    dictionary_count = write_dictionary(dictionary_out, dictionary_words)

    write_csv(
        manifest_out,
        output_rows,
        [
            "recording_id",
            "sex",
            "speaker",
            "source_audio_file",
            "mfa_audio_file",
            "mfa_transcript_file",
            "transcript",
            "transcript_clean",
            "mfa_transcript",
            "audio_status",
            "transcript_status",
        ],
    )
    write_csv(summary_out, summary_rows(output_rows, dictionary_count, skipped_eligible), ["metric", "value"])

    print(f"Wrote {len(output_rows)} pilot rows to {manifest_out}")
    print(f"Wrote dictionary with {dictionary_count} words to {dictionary_out}")
    print(f"Wrote summary to {summary_out}")
    print(f"MFA corpus directory: {corpus_dir}")


if __name__ == "__main__":
    main()
