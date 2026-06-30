#!/usr/bin/env python3
"""Prepare AISHELL-1 Mandarin tone manifests and pinyin-derived tone units.

This script does not create acoustic tone measurements. It reads the AISHELL
transcript file, maps recording IDs to split/speaker paths, derives provisional
lexical tone labels with pypinyin, and writes audit tables that can later be
joined to forced-alignment intervals.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_ROOT = Path("/Volumes/Corpus_Studies/Corpora/Mandarin/data_aishell")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/aishell_mandarin")

RECORDING_RE = re.compile(r"^(?P<prefix>[A-Z0-9]+)(?P<speaker>S\d+)(?P<item>W\d+)$")
PINYIN_TONE_RE = re.compile(r"^(?P<base>.+?)(?P<tone>[1-5])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--transcript", type=Path, default=None)
    parser.add_argument("--wav-root", type=Path, default=None)
    parser.add_argument(
        "--file-root-override",
        type=Path,
        default=None,
        help=(
            "Root to use when writing output file paths. Use this when parsing "
            "locally but preparing CSVs for NAS-side audio extraction."
        ),
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--manifest-out", type=Path, default=None)
    parser.add_argument("--tone-units-out", type=Path, default=None)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--tone-counts-out", type=Path, default=None)
    parser.add_argument("--max-recordings", type=int, default=0)
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Write only the utterance manifest, not pinyin-derived tone units.",
    )
    parser.add_argument(
        "--scan-audio",
        action="store_true",
        help="Inventory WAV basenames under split/speaker directories. Faster on data-local storage.",
    )
    return parser.parse_args()


def load_pypinyin():
    try:
        from pypinyin import Style, lazy_pinyin  # type: ignore

        return lazy_pinyin, Style
    except Exception as exc:  # pragma: no cover - depends on local env
        return None, exc


def parse_recording_id(recording_id: str) -> dict[str, str]:
    match = RECORDING_RE.match(recording_id)
    if not match:
        return {"recording_prefix": "", "speaker": "", "item_code": ""}
    groups = match.groupdict()
    return {
        "recording_prefix": groups["prefix"],
        "speaker": groups["speaker"],
        "item_code": groups["item"],
    }


def split_speaker_dirs(wav_root: Path) -> dict[str, str]:
    speaker_to_split: dict[str, str] = {}
    for split in ("train", "dev", "test"):
        split_dir = wav_root / split
        if not split_dir.exists():
            continue
        with os.scandir(split_dir) as iterator:
            for entry in iterator:
                if entry.is_dir():
                    speaker_to_split[entry.name] = split
    return speaker_to_split


def audio_inventory(wav_root: Path, speaker_to_split: dict[str, str]) -> set[str]:
    ids: set[str] = set()
    for speaker, split in speaker_to_split.items():
        speaker_dir = wav_root / split / speaker
        if not speaker_dir.exists():
            continue
        with os.scandir(speaker_dir) as iterator:
            for entry in iterator:
                if entry.is_file() and entry.name.lower().endswith(".wav"):
                    ids.add(Path(entry.name).stem)
    return ids


def load_manifest_rows(
    transcript_path: Path,
    wav_root: Path,
    file_wav_root: Path,
    speaker_to_split: dict[str, str],
    wav_ids: set[str] | None,
    max_recordings: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with transcript_path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.strip().split()
            if not fields:
                continue
            recording_id = fields[0]
            words = fields[1:]
            parsed = parse_recording_id(recording_id)
            speaker = parsed["speaker"]
            split = speaker_to_split.get(speaker, "")
            expected_path = file_wav_root / split / speaker / f"{recording_id}.wav" if split else Path("")
            audio_exists = recording_id in wav_ids if wav_ids is not None else ""
            rows.append(
                {
                    "recording_id": recording_id,
                    "source_line_number": line_number,
                    "split": split,
                    "speaker": speaker,
                    "recording_prefix": parsed["recording_prefix"],
                    "item_code": parsed["item_code"],
                    "transcript": " ".join(words),
                    "word_count": len(words),
                    "character_count": sum(len(word) for word in words),
                    "file": str(expected_path),
                    "audio_inventory_checked": wav_ids is not None,
                    "audio_exists": audio_exists,
                    "source_transcript": str(transcript_path),
                }
            )
            if max_recordings > 0 and len(rows) >= max_recordings:
                break
    return rows


def pinyin_for_word(word: str, lazy_pinyin, Style) -> list[str]:
    try:
        syllables = lazy_pinyin(
            word,
            style=Style.TONE3,
            neutral_tone_with_five=True,
            errors="default",
        )
    except Exception:
        syllables = []
    if len(syllables) == len(word):
        return syllables
    out: list[str] = []
    for char in word:
        try:
            values = lazy_pinyin(
                char,
                style=Style.TONE3,
                neutral_tone_with_five=True,
                errors="default",
            )
        except Exception:
            values = []
        out.append(values[0] if values else char)
    return out


def parse_pinyin_tone(pinyin_value: str) -> tuple[str, str, str]:
    match = PINYIN_TONE_RE.match(pinyin_value)
    if not match:
        return pinyin_value, "unknown", "missing_tone_digit"
    tone = match.group("tone")
    return match.group("base"), f"T{tone}", ""


def build_tone_units(manifest_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], str]:
    lazy_pinyin, style_or_error = load_pypinyin()
    if lazy_pinyin is None:
        return [], f"pypinyin_unavailable:{style_or_error}"
    Style = style_or_error
    rows: list[dict[str, object]] = []
    for manifest in manifest_rows:
        recording_unit_index = 0
        for word_index, word in enumerate(str(manifest["transcript"]).split(), start=1):
            syllables = pinyin_for_word(word, lazy_pinyin, Style)
            for char_index, char in enumerate(word, start=1):
                recording_unit_index += 1
                pinyin_value = syllables[char_index - 1] if char_index - 1 < len(syllables) else char
                pinyin_base, tone_category, flag = parse_pinyin_tone(pinyin_value)
                token_id = f"{manifest['recording_id']}_w{word_index:03d}_c{char_index:02d}"
                rows.append(
                    {
                        "token_id": token_id,
                        "domain": "tone",
                        "language": "Mandarin",
                        "source_corpus": "AISHELL-1",
                        "speaker": manifest["speaker"],
                        "split": manifest["split"],
                        "recording_id": manifest["recording_id"],
                        "file": manifest["file"],
                        "word": word,
                        "word_index": word_index,
                        "character": char,
                        "character_index": char_index,
                        "category": tone_category,
                        "control_group": pinyin_base,
                        "pinyin_tone3": pinyin_value,
                        "pinyin_base": pinyin_base,
                        "recording_unit_index": recording_unit_index,
                        "start": "",
                        "end": "",
                        "measurement_method": "pypinyin_lexical_tone_parse_pending_alignment",
                        "quality_flag": flag,
                        "recommended_exclude": bool(flag) or tone_category == "T5",
                    }
                )
    return rows, ""


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summary_rows(
    manifest_rows: list[dict[str, object]],
    tone_rows: list[dict[str, object]],
    pinyin_status: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        rows.append({"metric": metric, "value": value})

    add("recordings_total", len(manifest_rows))
    add("speakers_total", len({row["speaker"] for row in manifest_rows if row["speaker"]}))
    split_counts = Counter(str(row["split"]) for row in manifest_rows)
    for split, n in sorted(split_counts.items()):
        add(f"recordings_{split or 'unknown_split'}", n)
    checked = [row for row in manifest_rows if row["audio_inventory_checked"]]
    add("audio_inventory_checked", bool(checked))
    if checked:
        add("recordings_with_audio", sum(1 for row in checked if row["audio_exists"]))
        add("recordings_missing_audio", sum(1 for row in checked if not row["audio_exists"]))
    add("tone_units_total", len(tone_rows))
    add("tone_units_recommended_exclude", sum(1 for row in tone_rows if row["recommended_exclude"]))
    tone_counts = Counter(str(row["category"]) for row in tone_rows)
    for tone, n in sorted(tone_counts.items()):
        add(f"tone_units_{tone}", n)
    if pinyin_status:
        add("pinyin_status", pinyin_status)
    return rows


def tone_count_rows(tone_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for row in tone_rows:
        grouped[(str(row["split"]), str(row["speaker"]), str(row["control_group"]), str(row["category"]))] += 1
    return [
        {
            "split": split,
            "speaker": speaker,
            "control_group": control_group,
            "category": category,
            "n": n,
        }
        for (split, speaker, control_group, category), n in sorted(grouped.items())
    ]


def main() -> None:
    args = parse_args()
    root = args.root
    transcript_path = args.transcript or root / "transcript" / "aishell_transcript_v0.8.txt"
    wav_root = args.wav_root or root / "wav"
    file_wav_root = (args.file_root_override / "wav") if args.file_root_override else wav_root
    out_dir = args.out_dir
    manifest_out = args.manifest_out or out_dir / "aishell_manifest.csv"
    tone_units_out = args.tone_units_out or out_dir / "aishell_mandarin_tone_units_pending_alignment.csv"
    summary_out = args.summary_out or out_dir / "aishell_parse_summary.csv"
    tone_counts_out = args.tone_counts_out or out_dir / "aishell_mandarin_tone_unit_counts.csv"

    speaker_to_split = split_speaker_dirs(wav_root)
    wav_ids = audio_inventory(wav_root, speaker_to_split) if args.scan_audio else None
    manifest_rows = load_manifest_rows(
        transcript_path,
        wav_root,
        file_wav_root,
        speaker_to_split,
        wav_ids,
        args.max_recordings,
    )
    tone_rows: list[dict[str, object]] = []
    pinyin_status = ""
    if not args.manifest_only:
        tone_rows, pinyin_status = build_tone_units(manifest_rows)

    manifest_fields = [
        "recording_id",
        "source_line_number",
        "split",
        "speaker",
        "recording_prefix",
        "item_code",
        "transcript",
        "word_count",
        "character_count",
        "file",
        "audio_inventory_checked",
        "audio_exists",
        "source_transcript",
    ]
    tone_fields = [
        "token_id",
        "domain",
        "language",
        "source_corpus",
        "speaker",
        "split",
        "recording_id",
        "file",
        "word",
        "word_index",
        "character",
        "character_index",
        "category",
        "control_group",
        "pinyin_tone3",
        "pinyin_base",
        "recording_unit_index",
        "start",
        "end",
        "measurement_method",
        "quality_flag",
        "recommended_exclude",
    ]
    write_csv(manifest_out, manifest_rows, manifest_fields)
    if not args.manifest_only:
        write_csv(tone_units_out, tone_rows, tone_fields)
        write_csv(tone_counts_out, tone_count_rows(tone_rows), ["split", "speaker", "control_group", "category", "n"])
    write_csv(summary_out, summary_rows(manifest_rows, tone_rows, pinyin_status), ["metric", "value"])

    print(f"Wrote {len(manifest_rows)} manifest rows to {manifest_out}")
    if not args.manifest_only:
        print(f"Wrote {len(tone_rows)} tone-unit rows to {tone_units_out}")
        print(f"Wrote tone counts to {tone_counts_out}")
    print(f"Wrote summary to {summary_out}")


if __name__ == "__main__":
    main()
