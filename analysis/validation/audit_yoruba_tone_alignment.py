#!/usr/bin/env python3
"""Audit aligned Yoruba tone-token feature tables.

The Yoruba extractor writes one row per aligned tone-bearing unit. This script
summarizes alignment/feature quality and writes a deterministic hand-audit
sample that can be checked against the source WAV and TextGrid files before
paper-scale validation is treated as final.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_FEATURES = Path("analysis/validation/outputs/yoruba_slr86_tone_features/yoruba_slr86_tone_features.csv")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/yoruba_slr86_tone_alignment_audit")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--by-category-out", type=Path, default=None)
    parser.add_argument("--by-category-control-out", type=Path, default=None)
    parser.add_argument("--flag-counts-out", type=Path, default=None)
    parser.add_argument("--sample-out", type=Path, default=None)
    parser.add_argument("--issues-out", type=Path, default=None)
    parser.add_argument("--sample-size", type=int, default=120)
    parser.add_argument("--max-per-category", type=int, default=40)
    parser.add_argument("--max-per-category-control", type=int, default=6)
    parser.add_argument("--issue-sample-size", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument("--low-voiced-prop", type=float, default=0.5)
    parser.add_argument("--short-duration-ms", type=float, default=35.0)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    fieldnames.append(key)
                    seen.add(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def finite_float(value: object) -> float:
    try:
        out = float(value)
    except Exception:
        return math.nan
    return out if math.isfinite(out) else math.nan


def median(values: list[float]) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return (clean[mid - 1] + clean[mid]) / 2.0


def split_quality_flags(row: dict[str, str]) -> list[str]:
    return [flag for flag in str(row.get("quality_flag", "")).split(";") if flag]


def expected_phone_label(row: dict[str, str]) -> str:
    quality = str(row.get("vowel_quality") or row.get("control_group") or "")
    if quality == "syllabic_n":
        return "n"
    return quality


def audit_flags(row: dict[str, str], low_voiced_prop: float, short_duration_ms: float) -> list[str]:
    flags: list[str] = []
    quality_flags = set(split_quality_flags(row))
    if "word_label_mismatch" in quality_flags:
        flags.append("word_label_mismatch")
    expected = expected_phone_label(row)
    observed = str(row.get("phone_interval_label", "")).strip().lower()
    if expected and observed and expected != observed:
        flags.append("phone_label_mismatch")
    if parse_bool(row.get("recommended_exclude", "")):
        flags.append("recommended_exclude")
    voiced_prop = finite_float(row.get("f0_voiced_prop", ""))
    if math.isfinite(voiced_prop) and voiced_prop < low_voiced_prop:
        flags.append("low_voiced_prop")
    duration_ms = finite_float(row.get("duration_ms", ""))
    if math.isfinite(duration_ms) and duration_ms < short_duration_ms:
        flags.append("short_duration_ms")
    if not str(row.get("file", "")).strip():
        flags.append("missing_audio_path")
    if not str(row.get("textgrid", "")).strip():
        flags.append("missing_textgrid_path")
    return flags


def enrich_rows(rows: list[dict[str, str]], low_voiced_prop: float, short_duration_ms: float) -> list[dict[str, object]]:
    enriched: list[dict[str, object]] = []
    for row in rows:
        out: dict[str, object] = dict(row)
        expected = expected_phone_label(row)
        observed = str(row.get("phone_interval_label", "")).strip().lower()
        out["expected_phone_label"] = expected
        out["phone_label_match"] = bool(expected and observed and expected == observed)
        out["audit_flags"] = ";".join(audit_flags(row, low_voiced_prop, short_duration_ms))
        enriched.append(out)
    return enriched


def count_distinct(rows: list[dict[str, object]], col: str) -> int:
    return len({str(row.get(col, "")) for row in rows if str(row.get(col, "")).strip()})


def metric_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    flag_counts = Counter()
    quality_flag_counts = Counter()
    for row in rows:
        flag_counts.update(flag for flag in str(row.get("audit_flags", "")).split(";") if flag)
        quality_flag_counts.update(split_quality_flags({k: str(v) for k, v in row.items()}))

    out: list[dict[str, object]] = []

    def add(metric: str, value: object) -> None:
        out.append({"metric": metric, "value": value})

    add("tokens_total", len(rows))
    add("speakers_total", count_distinct(rows, "speaker"))
    add("recordings_total", count_distinct(rows, "recording_id"))
    add("categories_total", count_distinct(rows, "category"))
    add("control_groups_total", count_distinct(rows, "control_group"))
    add("textgrids_total", count_distinct(rows, "textgrid"))
    add("audio_files_total", count_distinct(rows, "file"))
    add("tokens_with_audit_flags", sum(1 for row in rows if str(row.get("audit_flags", "")).strip()))
    for flag, n in sorted(flag_counts.items()):
        add(f"tokens_flagged_{flag}", n)
    for flag, n in sorted(quality_flag_counts.items()):
        add(f"tokens_quality_flag_{flag}", n)
    for category, n in sorted(Counter(str(row.get("category", "")) for row in rows).items()):
        add(f"tokens_{category or 'missing_category'}", n)
    add("median_duration_ms", median([finite_float(row.get("duration_ms", "")) for row in rows]))
    add("median_f0_voiced_prop", median([finite_float(row.get("f0_voiced_prop", "")) for row in rows]))
    add("median_f0_mean_st", median([finite_float(row.get("f0_mean_st", "")) for row in rows]))
    add(
        "median_abs_f0_delta_st",
        median([abs(finite_float(row.get("f0_delta_st", ""))) for row in rows]),
    )
    return out


def summarize_group(rows: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(str(row.get(key, "")) for key in keys)].append(row)

    out: list[dict[str, object]] = []
    for group_key, group_rows in sorted(grouped.items()):
        flag_counts = Counter(
            flag
            for row in group_rows
            for flag in str(row.get("audit_flags", "")).split(";")
            if flag
        )
        row = {key: value for key, value in zip(keys, group_key)}
        row.update(
            {
                "n_tokens": len(group_rows),
                "n_speakers": count_distinct(group_rows, "speaker"),
                "n_recordings": count_distinct(group_rows, "recording_id"),
                "n_control_groups": count_distinct(group_rows, "control_group"),
                "n_word_label_mismatch": flag_counts.get("word_label_mismatch", 0),
                "n_phone_label_mismatch": flag_counts.get("phone_label_mismatch", 0),
                "n_low_voiced_prop": flag_counts.get("low_voiced_prop", 0),
                "median_duration_ms": median([finite_float(item.get("duration_ms", "")) for item in group_rows]),
                "median_f0_voiced_prop": median([finite_float(item.get("f0_voiced_prop", "")) for item in group_rows]),
                "median_f0_mean_st": median([finite_float(item.get("f0_mean_st", "")) for item in group_rows]),
                "median_f0_start_st": median([finite_float(item.get("f0_start_st", "")) for item in group_rows]),
                "median_f0_mid_st": median([finite_float(item.get("f0_mid_st", "")) for item in group_rows]),
                "median_f0_end_st": median([finite_float(item.get("f0_end_st", "")) for item in group_rows]),
            }
        )
        out.append(row)
    return out


def flag_count_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    counts = Counter()
    for row in rows:
        counts.update(flag for flag in str(row.get("audit_flags", "")).split(";") if flag)
        counts.update(f"quality_flag:{flag}" for flag in split_quality_flags({k: str(v) for k, v in row.items()}))
    return [{"flag": flag, "n": n} for flag, n in sorted(counts.items())]


def stratified_sample(
    rows: list[dict[str, object]],
    sample_size: int,
    max_per_category: int,
    max_per_category_control: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    by_group: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_group[(str(row.get("category", "")), str(row.get("control_group", "")))].append(row)
    for group_rows in by_group.values():
        rng.shuffle(group_rows)

    selected: list[dict[str, object]] = []
    category_counts: Counter[str] = Counter()
    group_counts: Counter[tuple[str, str]] = Counter()
    groups = sorted(by_group)
    made_progress = True
    while len(selected) < sample_size and made_progress:
        made_progress = False
        for group in groups:
            if len(selected) >= sample_size:
                break
            category = group[0]
            if category_counts[category] >= max_per_category:
                continue
            if group_counts[group] >= max_per_category_control:
                continue
            group_rows = by_group[group]
            if not group_rows:
                continue
            selected.append(group_rows.pop())
            category_counts[category] += 1
            group_counts[group] += 1
            made_progress = True
    return selected


AUDIT_SAMPLE_FIELDS = [
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
    "phone_interval_label",
    "expected_phone_label",
    "phone_label_match",
    "start",
    "end",
    "duration_ms",
    "f0_voiced_frames",
    "f0_voiced_prop",
    "f0_mean_st",
    "f0_start_st",
    "f0_mid_st",
    "f0_end_st",
    "f0_slope_st_per_s",
    "word_interval_label",
    "word_start",
    "word_end",
    "quality_flag",
    "audit_flags",
    "file",
    "textgrid",
]


def main() -> None:
    args = parse_args()
    rows = enrich_rows(read_csv(args.features), args.low_voiced_prop, args.short_duration_ms)
    if not rows:
        raise RuntimeError(f"No rows found in {args.features}")

    summary_out = args.summary_out or args.out_dir / "yoruba_tone_alignment_audit_summary.csv"
    by_category_out = args.by_category_out or args.out_dir / "yoruba_tone_alignment_audit_by_category.csv"
    by_category_control_out = (
        args.by_category_control_out
        or args.out_dir / "yoruba_tone_alignment_audit_by_category_control.csv"
    )
    flag_counts_out = args.flag_counts_out or args.out_dir / "yoruba_tone_alignment_audit_flag_counts.csv"
    sample_out = args.sample_out or args.out_dir / "yoruba_tone_alignment_audit_sample.csv"
    issues_out = args.issues_out or args.out_dir / "yoruba_tone_alignment_audit_issues.csv"

    sample_rows = stratified_sample(
        rows,
        sample_size=args.sample_size,
        max_per_category=args.max_per_category,
        max_per_category_control=args.max_per_category_control,
        seed=args.seed,
    )
    issue_rows = [row for row in rows if str(row.get("audit_flags", "")).strip()]
    random.Random(args.seed).shuffle(issue_rows)
    issue_rows = issue_rows[: args.issue_sample_size]

    write_csv(summary_out, metric_rows(rows), ["metric", "value"])
    write_csv(by_category_out, summarize_group(rows, ["category"]))
    write_csv(by_category_control_out, summarize_group(rows, ["category", "control_group"]))
    write_csv(flag_counts_out, flag_count_rows(rows), ["flag", "n"])
    write_csv(sample_out, sample_rows, AUDIT_SAMPLE_FIELDS)
    write_csv(issues_out, issue_rows, AUDIT_SAMPLE_FIELDS)

    print(f"Read {len(rows)} aligned tone-token rows from {args.features}")
    print(f"Wrote summary to {summary_out}")
    print(f"Wrote category summary to {by_category_out}")
    print(f"Wrote category/control summary to {by_category_control_out}")
    print(f"Wrote flag counts to {flag_counts_out}")
    print(f"Wrote hand-audit sample to {sample_out}")
    print(f"Wrote issue sample to {issues_out}")


if __name__ == "__main__":
    main()
