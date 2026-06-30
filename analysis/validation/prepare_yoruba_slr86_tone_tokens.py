#!/usr/bin/env python3
"""Prepare OpenSLR 86 Yoruba tone manifests and orthographic tone units.

This script does not create acoustic tone measurements. It reads the SLR86 line
indexes and WAV headers, parses tone-bearing units from Yoruba orthography, and
writes audit tables that can later be joined to forced-alignment intervals.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import unicodedata
import wave
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


DEFAULT_ROOT = Path("/Volumes/Corpus_Studies/Corpora/Yoruba")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/yoruba_slr86")

ANNOTATION_RE = re.compile(r"\[([^\]]+)\]")

COMBINING_ACUTE = "\u0301"
COMBINING_GRAVE = "\u0300"
COMBINING_DOT_BELOW = "\u0323"
RIGHT_SINGLE_QUOTE = "\u2019"

LETTERLIKE_JOINERS = {"'", RIGHT_SINGLE_QUOTE}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--female-index", type=Path, default=None)
    parser.add_argument("--male-index", type=Path, default=None)
    parser.add_argument("--female-audio-dir", type=Path, default=None)
    parser.add_argument("--male-audio-dir", type=Path, default=None)
    parser.add_argument("--female-zip", type=Path, default=None)
    parser.add_argument("--male-zip", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--manifest-out", type=Path, default=None)
    parser.add_argument("--tone-units-out", type=Path, default=None)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--tone-counts-out", type=Path, default=None)
    parser.add_argument(
        "--max-recordings",
        type=int,
        default=0,
        help="Limit rows per sex for quick parser audits. Zero means all rows.",
    )
    parser.add_argument(
        "--skip-missing-audio",
        action="store_true",
        help="Drop manifest rows whose expected WAV is neither extracted nor present in the source zip.",
    )
    parser.add_argument(
        "--no-zip-fallback",
        action="store_true",
        help="Do not use source zip inventories to mark unextracted WAV files as available.",
    )
    parser.add_argument(
        "--skip-audio-headers",
        action="store_true",
        help="Only check audio presence; do not open extracted WAVs for sample rate or duration.",
    )
    return parser.parse_args()


def normalize_text(text: str) -> str:
    return unicodedata.normalize("NFC", text.strip())


def strip_annotations(text: str) -> tuple[str, list[str]]:
    labels = [match.group(1).strip().lower() for match in ANNOTATION_RE.finditer(text)]
    clean = ANNOTATION_RE.sub(" ", text)
    clean = re.sub(r"\s+", " ", clean).strip()
    return clean, labels


def split_recording_id(recording_id: str) -> tuple[str, str, str]:
    parts = recording_id.split("_", 2)
    if len(parts) != 3:
        return "", "", ""
    return parts[0], parts[1], parts[2]


def read_audio_info(
    path: Path,
    read_headers: bool = True,
    exists_override: bool | None = None,
) -> dict[str, object]:
    exists = path.exists() if exists_override is None else exists_override
    info: dict[str, object] = {
        "audio_exists": exists,
        "audio_in_zip": False,
        "audio_available": exists,
        "audio_location": "extracted_wav" if exists else "missing",
        "zip_path": "",
        "zip_member": "",
        "sample_rate": "",
        "channels": "",
        "sample_width_bytes": "",
        "n_frames": "",
        "duration_s": "",
        "audio_error": "",
    }
    if not exists or not read_headers:
        return info
    try:
        with wave.open(str(path), "rb") as handle:
            sample_rate = handle.getframerate()
            n_frames = handle.getnframes()
            info.update(
                {
                    "sample_rate": sample_rate,
                    "channels": handle.getnchannels(),
                    "sample_width_bytes": handle.getsampwidth(),
                    "n_frames": n_frames,
                    "duration_s": n_frames / sample_rate if sample_rate else math.nan,
                }
            )
    except Exception as exc:  # pragma: no cover - depends on local files
        info["audio_error"] = str(exc)
    return info


def zip_inventory(path: Path | None) -> set[str]:
    if path is None or not path.exists():
        return set()
    with zipfile.ZipFile(path) as archive:
        return {name for name in archive.namelist() if name.lower().endswith(".wav")}


def extracted_wav_inventory(audio_dir: Path) -> set[str]:
    if not audio_dir.exists():
        return set()
    names: set[str] = set()
    with os.scandir(audio_dir) as iterator:
        for entry in iterator:
            if entry.is_file() and entry.name.lower().endswith(".wav"):
                names.add(Path(entry.name).stem)
    return names


def mark_zip_fallback(
    info: dict[str, object],
    zip_path: Path | None,
    zip_members: set[str],
    recording_id: str,
) -> None:
    if info["audio_exists"] or zip_path is None:
        return
    wav_name = f"{recording_id}.wav"
    candidates = (wav_name, f"yo_ng/{wav_name}", f"yo_ng_female/{wav_name}", f"yo_ng_male/{wav_name}")
    for member in candidates:
        if member in zip_members:
            info["audio_in_zip"] = True
            info["audio_available"] = True
            info["audio_location"] = "zip_member"
            info["zip_path"] = str(zip_path)
            info["zip_member"] = member
            return
    suffix = f"/{wav_name}"
    for member in zip_members:
        if member.endswith(suffix):
            info["audio_in_zip"] = True
            info["audio_available"] = True
            info["audio_location"] = "zip_member"
            info["zip_path"] = str(zip_path)
            info["zip_member"] = member
            return


def load_index_rows(
    index_path: Path,
    audio_dir: Path,
    extracted_ids: set[str],
    zip_path: Path | None,
    zip_members: set[str],
    sex: str,
    max_recordings: int,
    skip_missing_audio: bool,
    read_audio_headers: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with index_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw_index, fields in enumerate(reader, start=1):
            if not fields:
                continue
            if len(fields) < 2:
                raise ValueError(f"Expected two TSV columns in {index_path}:{raw_index}")
            recording_id = fields[0].strip()
            transcript = normalize_text(fields[1])
            prefix, speaker_code, item_code = split_recording_id(recording_id)
            clean_transcript, annotation_labels = strip_annotations(transcript)
            audio_path = audio_dir / f"{recording_id}.wav"
            audio_info = read_audio_info(
                audio_path,
                read_headers=read_audio_headers,
                exists_override=recording_id in extracted_ids,
            )
            mark_zip_fallback(audio_info, zip_path, zip_members, recording_id)
            if skip_missing_audio and not audio_info["audio_available"]:
                continue
            speaker = f"{prefix}_{speaker_code}" if prefix and speaker_code else speaker_code
            row = {
                "recording_id": recording_id,
                "sex": sex,
                "speaker": speaker,
                "speaker_code": speaker_code,
                "speaker_source": "recording_prefix_plus_middle_recording_id_field",
                "recording_prefix": prefix,
                "item_code": item_code,
                "source_line_index": str(index_path),
                "source_line_number": raw_index,
                "transcript": transcript,
                "transcript_clean": clean_transcript,
                "annotation_labels": ";".join(annotation_labels),
                "has_annotation": bool(annotation_labels),
                "recommended_exclude": bool(annotation_labels),
                "file": str(audio_path),
                **audio_info,
            }
            rows.append(row)
            if max_recordings > 0 and len(rows) >= max_recordings:
                break
    return rows


def is_word_char(ch: str) -> bool:
    category = unicodedata.category(ch)
    return category.startswith("L") or category.startswith("M") or ch in LETTERLIKE_JOINERS


def iter_words(text: str) -> Iterable[dict[str, object]]:
    current: list[str] = []
    start = 0
    word_index = 0
    for idx, ch in enumerate(text):
        if is_word_char(ch):
            if not current:
                start = idx
            current.append(ch)
            continue
        if current:
            word = normalize_text("".join(current)).strip("".join(LETTERLIKE_JOINERS))
            if word:
                word_index += 1
                yield {
                    "word": word,
                    "word_index": word_index,
                    "word_start_char": start,
                    "word_end_char": idx,
                }
            current = []
    if current:
        word = normalize_text("".join(current)).strip("".join(LETTERLIKE_JOINERS))
        if word:
            word_index += 1
            yield {
                "word": word,
                "word_index": word_index,
                "word_start_char": start,
                "word_end_char": len(text),
            }


def grapheme_clusters(word: str) -> list[dict[str, object]]:
    clusters: list[dict[str, object]] = []
    current: list[str] = []
    start = 0
    for idx, ch in enumerate(word):
        if unicodedata.category(ch).startswith("M"):
            if current:
                current.append(ch)
            continue
        if current:
            clusters.append({"text": normalize_text("".join(current)), "start": start, "end": idx})
        current = [ch]
        start = idx
    if current:
        clusters.append({"text": normalize_text("".join(current)), "start": start, "end": len(word)})
    return clusters


def cluster_base_and_marks(cluster_text: str) -> tuple[str, set[str]]:
    decomposed = unicodedata.normalize("NFD", cluster_text.lower())
    base = ""
    marks: set[str] = set()
    for ch in decomposed:
        if unicodedata.category(ch).startswith("M"):
            marks.add(ch)
        elif not base:
            base = ch
    return base, marks


def vowel_quality(base: str, marks: set[str]) -> str | None:
    if base not in {"a", "e", "i", "o", "u"}:
        return None
    if base == "e" and COMBINING_DOT_BELOW in marks:
        return "e_dot"
    if base == "o" and COMBINING_DOT_BELOW in marks:
        return "o_dot"
    return base


def tone_category(marks: set[str]) -> tuple[str, str]:
    has_high = COMBINING_ACUTE in marks
    has_low = COMBINING_GRAVE in marks
    if has_high and has_low:
        return "ambiguous", "both_acute_and_grave"
    if has_high:
        return "H", ""
    if has_low:
        return "L", ""
    return "M", ""


def parse_tone_units_for_word(word_info: dict[str, object]) -> list[dict[str, object]]:
    word = str(word_info["word"])
    clusters = grapheme_clusters(word)
    units: list[dict[str, object]] = []
    cluster_meta = []
    for idx, cluster in enumerate(clusters):
        base, marks = cluster_base_and_marks(str(cluster["text"]))
        quality = vowel_quality(base, marks)
        tbu_type = "vowel"
        if quality is None and base == "n" and (COMBINING_ACUTE in marks or COMBINING_GRAVE in marks):
            quality = "syllabic_n"
            tbu_type = "syllabic_nasal"
        cluster_meta.append(
            {
                "cluster_index": idx,
                "text": cluster["text"],
                "start": cluster["start"],
                "end": cluster["end"],
                "base": base,
                "marks": marks,
                "quality": quality,
                "tbu_type": tbu_type,
            }
        )

    for meta in cluster_meta:
        if meta["quality"] is None:
            continue
        tone, tone_flag = tone_category(meta["marks"])
        preceding = previous_letterlike_cluster(cluster_meta, int(meta["cluster_index"]))
        following = next_letterlike_cluster(cluster_meta, int(meta["cluster_index"]))
        units.append(
            {
                **word_info,
                "orthographic_unit": meta["text"],
                "orthographic_base": meta["base"],
                "vowel_quality": meta["quality"],
                "tbu_type": meta["tbu_type"],
                "tone_category": tone,
                "tone_parse_flag": tone_flag,
                "word_unit_index": len(units) + 1,
                "cluster_index": int(meta["cluster_index"]) + 1,
                "unit_start_char_in_word": meta["start"],
                "unit_end_char_in_word": meta["end"],
                "preceding_letter": preceding["text"] if preceding else "",
                "following_letter": following["text"] if following else "",
                "following_n_letter": bool(following and following["base"] == "n"),
                "nasal_oral_status": "not_assessed",
            }
        )

    unit_count = len(units)
    for unit in units:
        if unit_count == 1:
            position = "monosyllabic"
        elif unit["word_unit_index"] == 1:
            position = "initial"
        elif unit["word_unit_index"] == unit_count:
            position = "final"
        else:
            position = "medial"
        unit["word_tone_unit_count"] = unit_count
        unit["tone_unit_position"] = position
    return units


def previous_letterlike_cluster(
    clusters: list[dict[str, object]],
    current_index: int,
) -> dict[str, object] | None:
    for idx in range(current_index - 1, -1, -1):
        base = clusters[idx]["base"]
        if base and base not in LETTERLIKE_JOINERS:
            return clusters[idx]
    return None


def next_letterlike_cluster(
    clusters: list[dict[str, object]],
    current_index: int,
) -> dict[str, object] | None:
    for idx in range(current_index + 1, len(clusters)):
        base = clusters[idx]["base"]
        if base and base not in LETTERLIKE_JOINERS:
            return clusters[idx]
    return None


def build_tone_units(manifest_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    tone_rows: list[dict[str, object]] = []
    for row in manifest_rows:
        recording_unit_index = 0
        for word_info in iter_words(str(row["transcript_clean"])):
            word_units = parse_tone_units_for_word(word_info)
            for unit in word_units:
                recording_unit_index += 1
                quality_flag = []
                if row["has_annotation"]:
                    quality_flag.append("recording_has_annotation")
                if unit["tone_parse_flag"]:
                    quality_flag.append(str(unit["tone_parse_flag"]))
                token_id = (
                    f"{row['recording_id']}_w{int(unit['word_index']):03d}"
                    f"_u{int(unit['word_unit_index']):02d}"
                )
                tone_rows.append(
                    {
                        "token_id": token_id,
                        "domain": "tone",
                        "language": "Yoruba",
                        "source_corpus": "OpenSLR86",
                        "speaker": row["speaker"],
                        "speaker_code": row["speaker_code"],
                        "speaker_source": row["speaker_source"],
                        "sex": row["sex"],
                        "recording_id": row["recording_id"],
                        "file": row["file"],
                        "word": unit["word"],
                        "word_index": unit["word_index"],
                        "category": unit["tone_category"],
                        "control_group": unit["vowel_quality"],
                        "orthographic_unit": unit["orthographic_unit"],
                        "orthographic_base": unit["orthographic_base"],
                        "vowel_quality": unit["vowel_quality"],
                        "tbu_type": unit["tbu_type"],
                        "tone_unit_position": unit["tone_unit_position"],
                        "word_tone_unit_count": unit["word_tone_unit_count"],
                        "word_unit_index": unit["word_unit_index"],
                        "recording_unit_index": recording_unit_index,
                        "preceding_letter": unit["preceding_letter"],
                        "following_letter": unit["following_letter"],
                        "following_n_letter": unit["following_n_letter"],
                        "nasal_oral_status": unit["nasal_oral_status"],
                        "start": "",
                        "end": "",
                        "measurement_method": "orthographic_tone_parse_pending_alignment",
                        "quality_flag": ";".join(quality_flag),
                        "recommended_exclude": bool(row["recommended_exclude"]) or unit["tone_category"] == "ambiguous",
                    }
                )
    return tone_rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summary_rows(manifest_rows: list[dict[str, object]], tone_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        rows.append({"metric": metric, "value": value})

    add("recordings_total", len(manifest_rows))
    add("recordings_with_extracted_audio", sum(1 for row in manifest_rows if row["audio_exists"]))
    add("recordings_with_zip_audio", sum(1 for row in manifest_rows if row["audio_in_zip"]))
    add("recordings_with_available_audio", sum(1 for row in manifest_rows if row["audio_available"]))
    add("recordings_missing_audio", sum(1 for row in manifest_rows if not row["audio_available"]))
    add("recordings_with_annotations", sum(1 for row in manifest_rows if row["has_annotation"]))
    add("speakers_total", len({row["speaker"] for row in manifest_rows}))
    for sex in ("female", "male"):
        sex_rows = [row for row in manifest_rows if row["sex"] == sex]
        add(f"{sex}_recordings", len(sex_rows))
        add(f"{sex}_speakers", len({row["speaker"] for row in sex_rows}))
    add("tone_units_total", len(tone_rows))
    add("tone_units_recommended_exclude", sum(1 for row in tone_rows if row["recommended_exclude"]))
    tone_counts = Counter(str(row["category"]) for row in tone_rows)
    for tone in ("H", "M", "L", "ambiguous"):
        add(f"tone_units_{tone}", tone_counts.get(tone, 0))
    audio_rates = Counter(str(row["sample_rate"]) for row in manifest_rows if row["sample_rate"] != "")
    for sample_rate, count in sorted(audio_rates.items()):
        add(f"sample_rate_{sample_rate}_recordings", count)
    annotation_counts = Counter()
    for row in manifest_rows:
        for label in str(row["annotation_labels"]).split(";"):
            if label:
                annotation_counts[label] += 1
    for label, count in sorted(annotation_counts.items()):
        add(f"annotation_{label}", count)
    return rows


def tone_count_rows(tone_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for row in tone_rows:
        key = (
            str(row["sex"]),
            str(row["speaker"]),
            str(row["control_group"]),
            str(row["category"]),
        )
        grouped[key] += 1
    return [
        {
            "sex": sex,
            "speaker": speaker,
            "control_group": control_group,
            "category": category,
            "n": n,
        }
        for (sex, speaker, control_group, category), n in sorted(grouped.items())
    ]


def main() -> None:
    args = parse_args()
    root = args.root
    female_index = args.female_index or root / "line_index_female.tsv"
    male_index = args.male_index or root / "line_index_male.tsv"
    female_audio_dir = args.female_audio_dir or root / "yo_ng_female"
    male_audio_dir = args.male_audio_dir or root / "yo_ng_male"
    female_zip = None if args.no_zip_fallback else (args.female_zip or root / "yo_ng_female.zip")
    male_zip = None if args.no_zip_fallback else (args.male_zip or root / "yo_ng_male.zip")
    female_extracted_ids = extracted_wav_inventory(female_audio_dir)
    male_extracted_ids = extracted_wav_inventory(male_audio_dir)
    female_zip_members = zip_inventory(female_zip)
    male_zip_members = zip_inventory(male_zip)

    out_dir = args.out_dir
    manifest_out = args.manifest_out or out_dir / "yoruba_slr86_manifest.csv"
    tone_units_out = args.tone_units_out or out_dir / "yoruba_slr86_tone_units_pending_alignment.csv"
    summary_out = args.summary_out or out_dir / "yoruba_slr86_parse_summary.csv"
    tone_counts_out = args.tone_counts_out or out_dir / "yoruba_slr86_tone_unit_counts.csv"

    manifest_rows = []
    manifest_rows.extend(
        load_index_rows(
            female_index,
            female_audio_dir,
            female_extracted_ids,
            female_zip,
            female_zip_members,
            "female",
            args.max_recordings,
            args.skip_missing_audio,
            not args.skip_audio_headers,
        )
    )
    manifest_rows.extend(
        load_index_rows(
            male_index,
            male_audio_dir,
            male_extracted_ids,
            male_zip,
            male_zip_members,
            "male",
            args.max_recordings,
            args.skip_missing_audio,
            not args.skip_audio_headers,
        )
    )
    tone_rows = build_tone_units(manifest_rows)

    manifest_fields = [
        "recording_id",
        "sex",
        "speaker",
        "speaker_code",
        "speaker_source",
        "recording_prefix",
        "item_code",
        "source_line_index",
        "source_line_number",
        "transcript",
        "transcript_clean",
        "annotation_labels",
        "has_annotation",
        "recommended_exclude",
        "file",
        "audio_exists",
        "audio_in_zip",
        "audio_available",
        "audio_location",
        "zip_path",
        "zip_member",
        "sample_rate",
        "channels",
        "sample_width_bytes",
        "n_frames",
        "duration_s",
        "audio_error",
    ]
    tone_fields = [
        "token_id",
        "domain",
        "language",
        "source_corpus",
        "speaker",
        "speaker_code",
        "speaker_source",
        "sex",
        "recording_id",
        "file",
        "word",
        "word_index",
        "category",
        "control_group",
        "orthographic_unit",
        "orthographic_base",
        "vowel_quality",
        "tbu_type",
        "tone_unit_position",
        "word_tone_unit_count",
        "word_unit_index",
        "recording_unit_index",
        "preceding_letter",
        "following_letter",
        "following_n_letter",
        "nasal_oral_status",
        "start",
        "end",
        "measurement_method",
        "quality_flag",
        "recommended_exclude",
    ]
    write_csv(manifest_out, manifest_rows, manifest_fields)
    write_csv(tone_units_out, tone_rows, tone_fields)
    write_csv(summary_out, summary_rows(manifest_rows, tone_rows), ["metric", "value"])
    write_csv(tone_counts_out, tone_count_rows(tone_rows), ["sex", "speaker", "control_group", "category", "n"])

    print(f"Wrote {len(manifest_rows)} manifest rows to {manifest_out}")
    print(f"Wrote {len(tone_rows)} tone-unit rows to {tone_units_out}")
    print(f"Wrote summaries to {summary_out} and {tone_counts_out}")


if __name__ == "__main__":
    main()
