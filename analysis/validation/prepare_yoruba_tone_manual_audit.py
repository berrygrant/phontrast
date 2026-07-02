#!/usr/bin/env python3
"""Prepare a manual review packet for Yoruba tone alignment audit rows.

The automated audit creates a balanced sample and an issue-focused sample. This
script combines those rows, adds blank manual decision columns, and optionally
materializes media symlinks/copies so the review can be done in Praat or any
TextGrid/audio viewer.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import shlex
from collections import Counter
from pathlib import Path


DEFAULT_AUDIT_DIR = Path("analysis/validation/outputs/yoruba_slr86_tone_alignment_audit_full_clean")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/yoruba_slr86_tone_manual_audit_full_clean")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, default=DEFAULT_AUDIT_DIR)
    parser.add_argument("--sample", type=Path, default=None)
    parser.add_argument("--issues", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--media-mode",
        choices=("none", "symlink", "copy"),
        default="none",
        help="Whether to materialize referenced WAV/TextGrid files in the packet.",
    )
    parser.add_argument(
        "--overwrite-media",
        action="store_true",
        help="Replace existing symlinks/copies in the packet media directories.",
    )
    return parser.parse_args()


def read_csv(path: Path, source_label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Audit CSV does not exist: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["_audit_source"] = source_label
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sanitize_filename(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return out.strip("._") or "item"


def token_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    token_id = row.get("token_id", "")
    if token_id:
        return ("token_id", token_id, "", "")
    return (
        row.get("recording_id", ""),
        row.get("category", ""),
        row.get("start", ""),
        row.get("end", ""),
    )


def combine_rows(sample_rows: list[dict[str, str]], issue_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key: dict[tuple[str, str, str, str], dict[str, str]] = {}
    source_order: dict[tuple[str, str, str, str], list[str]] = {}
    for row in sample_rows + issue_rows:
        key = token_key(row)
        if key not in by_key:
            by_key[key] = dict(row)
            source_order[key] = []
        source = row.get("_audit_source", "")
        if source and source not in source_order[key]:
            source_order[key].append(source)
        if source == "issue":
            by_key[key].update(row)
    for key, row in by_key.items():
        row["_audit_source"] = ";".join(source_order[key])
    return sorted(
        by_key.values(),
        key=lambda row: (
            "0" if "issue" in row.get("_audit_source", "") else "1",
            row.get("category", ""),
            row.get("control_group", ""),
            row.get("recording_id", ""),
            row.get("start", ""),
        ),
    )


def materialize_file(src_value: str, dst_dir: Path, media_mode: str, overwrite: bool) -> str:
    src = Path(src_value)
    if media_mode == "none" or not src_value:
        return ""
    if not src.exists():
        return ""
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    if dst.exists() or dst.is_symlink():
        if not overwrite:
            return str(dst)
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    if media_mode == "copy":
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src)
    return str(dst)


def praat_command(audio_path: str, textgrid_path: str) -> str:
    if not audio_path or not textgrid_path:
        return ""
    return "praat --open " + " ".join([shlex.quote(audio_path), shlex.quote(textgrid_path)])


def packet_row(
    row: dict[str, str],
    index: int,
    local_audio: str,
    local_textgrid: str,
) -> dict[str, object]:
    source = row.get("_audit_source", "")
    recording_id = row.get("recording_id", "")
    audit_item_id = f"yo_audit_{index:04d}_{sanitize_filename(recording_id)}"
    audio_for_open = local_audio or row.get("file", "")
    textgrid_for_open = local_textgrid or row.get("textgrid", "")
    return {
        "audit_item_id": audit_item_id,
        "audit_source": source,
        "review_priority": "issue" if "issue" in source else "sample",
        "token_id": row.get("token_id", ""),
        "recording_id": recording_id,
        "speaker": row.get("speaker", ""),
        "sex": row.get("sex", ""),
        "category": row.get("category", ""),
        "control_group": row.get("control_group", ""),
        "vowel_quality": row.get("vowel_quality", ""),
        "word": row.get("word", ""),
        "word_index": row.get("word_index", ""),
        "orthographic_unit": row.get("orthographic_unit", ""),
        "orthographic_base": row.get("orthographic_base", ""),
        "word_unit_index": row.get("word_unit_index", ""),
        "word_interval_label": row.get("word_interval_label", ""),
        "phone_interval_label": row.get("phone_interval_label", ""),
        "expected_phone_label": row.get("expected_phone_label", ""),
        "start": row.get("start", ""),
        "end": row.get("end", ""),
        "duration_ms": row.get("duration_ms", ""),
        "f0_voiced_frames": row.get("f0_voiced_frames", ""),
        "f0_voiced_prop": row.get("f0_voiced_prop", ""),
        "f0_mean_st": row.get("f0_mean_st", ""),
        "f0_start_st": row.get("f0_start_st", ""),
        "f0_mid_st": row.get("f0_mid_st", ""),
        "f0_end_st": row.get("f0_end_st", ""),
        "quality_flag": row.get("quality_flag", ""),
        "audit_flags": row.get("audit_flags", ""),
        "source_audio": row.get("file", ""),
        "source_textgrid": row.get("textgrid", ""),
        "local_audio": local_audio,
        "local_textgrid": local_textgrid,
        "praat_open_command": praat_command(audio_for_open, textgrid_for_open),
        "manual_word_alignment_ok": "",
        "manual_phone_alignment_ok": "",
        "manual_interval_quality": "",
        "manual_f0_usable": "",
        "manual_decision": "",
        "manual_notes": "",
        "reviewer": "",
        "review_date": "",
    }


def summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        out.append({"metric": metric, "value": value})

    add("manual_audit_items", len(rows))
    for source, n in sorted(Counter(str(row.get("audit_source", "")) for row in rows).items()):
        add(f"items_source_{source or 'unknown'}", n)
    for priority, n in sorted(Counter(str(row.get("review_priority", "")) for row in rows).items()):
        add(f"items_priority_{priority or 'unknown'}", n)
    for tone, n in sorted(Counter(str(row.get("category", "")) for row in rows).items()):
        add(f"items_{tone or 'missing_category'}", n)
    add("recordings_total", len({str(row.get("recording_id", "")) for row in rows}))
    add("speakers_total", len({str(row.get("speaker", "")) for row in rows}))
    return out


def write_readme(path: Path, items_csv: Path, summary_csv: Path, media_mode: str) -> None:
    text = f"""# Yoruba Tone Manual Audit Packet

Review table:

- `{items_csv.name}`
- `{summary_csv.name}`

Media mode: `{media_mode}`

Recommended manual coding values:

- `manual_word_alignment_ok`: `yes`, `no`, or `unclear`
- `manual_phone_alignment_ok`: `yes`, `no`, or `unclear`
- `manual_interval_quality`: `good`, `minor_issue`, `bad`, or `unclear`
- `manual_f0_usable`: `yes`, `no`, or `unclear`
- `manual_decision`: `accept`, `exclude`, or `unclear`

For each row, open the WAV/TextGrid pair with `praat_open_command`, navigate to
`start`-`end`, and check whether the word and phone intervals match the token.
Prioritize rows with `review_priority = issue`, then inspect the balanced
sample rows.
"""
    path.write_text(text, encoding="utf-8")


ITEM_FIELDS = [
    "audit_item_id",
    "audit_source",
    "review_priority",
    "token_id",
    "recording_id",
    "speaker",
    "sex",
    "category",
    "control_group",
    "vowel_quality",
    "word",
    "word_index",
    "orthographic_unit",
    "orthographic_base",
    "word_unit_index",
    "word_interval_label",
    "phone_interval_label",
    "expected_phone_label",
    "start",
    "end",
    "duration_ms",
    "f0_voiced_frames",
    "f0_voiced_prop",
    "f0_mean_st",
    "f0_start_st",
    "f0_mid_st",
    "f0_end_st",
    "quality_flag",
    "audit_flags",
    "source_audio",
    "source_textgrid",
    "local_audio",
    "local_textgrid",
    "praat_open_command",
    "manual_word_alignment_ok",
    "manual_phone_alignment_ok",
    "manual_interval_quality",
    "manual_f0_usable",
    "manual_decision",
    "manual_notes",
    "reviewer",
    "review_date",
]


def main() -> None:
    args = parse_args()
    sample_path = args.sample or args.audit_dir / "yoruba_tone_alignment_audit_sample.csv"
    issues_path = args.issues or args.audit_dir / "yoruba_tone_alignment_audit_issues.csv"
    sample_rows = read_csv(sample_path, "sample")
    issue_rows = read_csv(issues_path, "issue")
    combined = combine_rows(sample_rows, issue_rows)

    audio_dir = args.out_dir / "audio"
    textgrid_dir = args.out_dir / "textgrid"
    packet_rows: list[dict[str, object]] = []
    for index, row in enumerate(combined, start=1):
        local_audio = materialize_file(row.get("file", ""), audio_dir, args.media_mode, args.overwrite_media)
        local_textgrid = materialize_file(row.get("textgrid", ""), textgrid_dir, args.media_mode, args.overwrite_media)
        packet_rows.append(packet_row(row, index, local_audio, local_textgrid))

    items_out = args.out_dir / "yoruba_tone_manual_audit_items.csv"
    summary_out = args.out_dir / "yoruba_tone_manual_audit_summary.csv"
    readme_out = args.out_dir / "README.md"
    write_csv(items_out, packet_rows, ITEM_FIELDS)
    write_csv(summary_out, summary_rows(packet_rows), ["metric", "value"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_readme(readme_out, items_out, summary_out, args.media_mode)

    print(f"Read {len(sample_rows)} sample rows from {sample_path}")
    print(f"Read {len(issue_rows)} issue rows from {issues_path}")
    print(f"Wrote {len(packet_rows)} manual audit rows to {items_out}")
    print(f"Wrote summary to {summary_out}")
    print(f"Wrote instructions to {readme_out}")


if __name__ == "__main__":
    main()
