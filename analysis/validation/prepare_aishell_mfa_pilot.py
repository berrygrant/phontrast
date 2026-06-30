#!/usr/bin/env python3
"""Materialize a small AISHELL-1 Mandarin corpus for MFA pilot alignment.

The AISHELL tone-unit parser writes labels and paths, but MFA expects a corpus
directory containing audio files and matching transcript files. This script
samples utterances from the AISHELL manifest and creates that MFA-ready
directory by copying or symlinking WAVs and writing matching `.lab` files.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_ROOT = Path("/Volumes/Corpus_Studies/Corpora/Mandarin/data_aishell")
DEFAULT_MANIFEST = Path("analysis/validation/outputs/aishell_mandarin/aishell_manifest.csv")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/aishell_mfa_pilot")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="AISHELL root as seen from the machine running this script.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="AISHELL manifest produced by prepare_aishell_mandarin_tone_tokens.py.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--splits",
        default="dev",
        help="Comma-separated AISHELL splits to sample from, e.g. dev or dev,test.",
    )
    parser.add_argument("--max-speakers", type=int, default=5, help="0 means no speaker cap.")
    parser.add_argument(
        "--max-utterances-per-speaker",
        type=int,
        default=20,
        help="0 means no utterance cap within selected speakers.",
    )
    parser.add_argument("--seed", type=int, default=20260630)
    parser.add_argument(
        "--mode",
        choices=("copy", "symlink", "manifest-only"),
        default="copy",
        help="Use copy for local scratch alignment, symlink for data-local mounts.",
    )
    parser.add_argument(
        "--source",
        choices=("root", "manifest"),
        default="root",
        help="Use --root-derived WAV paths or the manifest file column.",
    )
    parser.add_argument(
        "--source-root-rewrite",
        default="",
        metavar="FROM=TO",
        help="Rewrite manifest file paths before reading audio, for cross-machine mounts.",
    )
    parser.add_argument(
        "--transcript-format",
        choices=("original", "nospace", "chars-spaced"),
        default="original",
        help="Controls transcript tokenization written to MFA .lab files.",
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


def parse_csv_list(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


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


def select_rows(
    rows: list[dict[str, str]],
    splits: set[str],
    max_speakers: int,
    max_utterances_per_speaker: int,
    seed: int,
) -> list[dict[str, str]]:
    eligible = [
        row
        for row in rows
        if row.get("split", "") in splits
        and row.get("speaker", "")
        and row.get("recording_id", "")
        and row.get("transcript", "").strip()
    ]
    by_speaker: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        by_speaker[row["speaker"]].append(row)

    rng = random.Random(seed)
    speakers = sorted(by_speaker)
    rng.shuffle(speakers)
    if max_speakers > 0:
        speakers = speakers[:max_speakers]

    selected: list[dict[str, str]] = []
    for speaker in speakers:
        speaker_rows = by_speaker[speaker][:]
        rng.shuffle(speaker_rows)
        if max_utterances_per_speaker > 0:
            speaker_rows = speaker_rows[:max_utterances_per_speaker]
        selected.extend(speaker_rows)

    return sorted(selected, key=lambda row: (row["split"], row["speaker"], row["recording_id"]))


def source_path(row: dict[str, str], root: Path, source: str, rewrite: tuple[str, str] | None) -> Path:
    if source == "manifest":
        raw_path = row.get("file", "")
        if rewrite:
            old, new = rewrite
            if raw_path.startswith(old):
                raw_path = new + raw_path[len(old) :]
        return Path(raw_path)
    return root / "wav" / row["split"] / row["speaker"] / f"{row['recording_id']}.wav"


def format_transcript(text: str, transcript_format: str) -> str:
    words = text.split()
    if transcript_format == "original":
        return " ".join(words)
    joined = "".join(words)
    if transcript_format == "nospace":
        return joined
    return " ".join(joined)


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


def summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        out.append({"metric": metric, "value": value})

    add("recordings_selected", len(rows))
    add("speakers_selected", len({row["speaker"] for row in rows}))
    for split, n in sorted(Counter(str(row["split"]) for row in rows).items()):
        add(f"recordings_{split}", n)
    for status, n in sorted(Counter(str(row["audio_status"]) for row in rows).items()):
        add(f"audio_{status}", n)
    for status, n in sorted(Counter(str(row["transcript_status"]) for row in rows).items()):
        add(f"transcript_{status}", n)
    add("recordings_missing_audio", sum(1 for row in rows if row["audio_status"] == "missing_source_audio"))
    return out


def main() -> None:
    args = parse_args()
    rewrite = parse_rewrite(args.source_root_rewrite)
    splits = parse_csv_list(args.splits)
    if not splits:
        raise ValueError("--splits selected no splits")
    if not args.transcript_extension.startswith("."):
        raise ValueError("--transcript-extension must include the leading dot")

    manifest_rows = read_manifest(args.manifest)
    selected = select_rows(
        manifest_rows,
        splits=splits,
        max_speakers=args.max_speakers,
        max_utterances_per_speaker=args.max_utterances_per_speaker,
        seed=args.seed,
    )

    corpus_dir = args.out_dir / "corpus"
    output_rows: list[dict[str, object]] = []
    for row in selected:
        src = source_path(row, args.root, args.source, rewrite)
        speaker_dir = corpus_dir / row["speaker"]
        dst_audio = speaker_dir / f"{row['recording_id']}{src.suffix or '.wav'}"
        dst_lab = speaker_dir / f"{row['recording_id']}{args.transcript_extension}"
        transcript = format_transcript(row["transcript"], args.transcript_format)

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
                "split": row["split"],
                "speaker": row["speaker"],
                "source_audio_file": str(src),
                "mfa_audio_file": str(dst_audio),
                "mfa_transcript_file": str(dst_lab),
                "transcript": row["transcript"],
                "mfa_transcript": transcript,
                "audio_status": audio_status,
                "transcript_status": transcript_status,
            }
        )

    manifest_out = args.out_dir / "aishell_mfa_pilot_manifest.csv"
    summary_out = args.out_dir / "aishell_mfa_pilot_summary.csv"
    write_csv(
        manifest_out,
        output_rows,
        [
            "recording_id",
            "split",
            "speaker",
            "source_audio_file",
            "mfa_audio_file",
            "mfa_transcript_file",
            "transcript",
            "mfa_transcript",
            "audio_status",
            "transcript_status",
        ],
    )
    write_csv(summary_out, summary_rows(output_rows), ["metric", "value"])

    print(f"Wrote {len(output_rows)} pilot rows to {manifest_out}")
    print(f"Wrote summary to {summary_out}")
    print(f"MFA corpus directory: {corpus_dir}")


if __name__ == "__main__":
    main()
