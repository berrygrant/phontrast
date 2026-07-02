#!/usr/bin/env python3
"""Extract Yoruba tone-token features from OpenSLR 86 MFA alignments.

This script joins orthographic Yoruba H/M/L tone labels to MFA word/phone
TextGrids and extracts F0/energy features from aligned vowel or syllabic-nasal
intervals. It expects the grapheme-style dictionary produced by
prepare_yoruba_mfa_pilot.py, where tone marks are stripped for alignment while
dot-below vowel distinctions are preserved.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
import time
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import soundfile as sf


DEFAULT_PILOT_MANIFEST = Path("analysis/validation/outputs/yoruba_slr86_mfa_pilot/yoruba_mfa_pilot_manifest.csv")
DEFAULT_TONE_UNITS = Path("analysis/validation/outputs/yoruba_slr86/yoruba_slr86_tone_units_pending_alignment.csv")
DEFAULT_ALIGNED_DIR = Path("analysis/validation/outputs/yoruba_slr86_mfa_pilot/aligned")
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/yoruba_slr86_tone_features")

NON_SPEECH = {"", "sil", "sp", "spn", "<eps>", "<unk>", "<oov>"}
COMBINING_ACUTE = "\u0301"
COMBINING_GRAVE = "\u0300"
RIGHT_SINGLE_QUOTE = "\u2019"
LETTERLIKE_JOINERS = {"'", RIGHT_SINGLE_QUOTE}
VOWEL_PHONE_LABELS = {"a", "e", "e_dot", "i", "o", "o_dot", "u"}


@dataclass(frozen=True)
class Interval:
    xmin: float
    xmax: float
    text: str

    @property
    def duration(self) -> float:
        return self.xmax - self.xmin


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-manifest", type=Path, default=DEFAULT_PILOT_MANIFEST)
    parser.add_argument("--tone-units", type=Path, default=DEFAULT_TONE_UNITS)
    parser.add_argument("--aligned-dir", type=Path, default=DEFAULT_ALIGNED_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--tokens-out", type=Path, default=None)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--feature-sets-out", type=Path, default=None)
    parser.add_argument("--max-recordings", type=int, default=0)
    parser.add_argument("--include-recommended-exclude", action="store_true")
    parser.add_argument("--min-token-ms", type=float, default=35.0)
    parser.add_argument("--min-voiced-frames", type=int, default=2)
    parser.add_argument("--frame-ms", type=float, default=50.0)
    parser.add_argument("--hop-ms", type=float, default=10.0)
    parser.add_argument("--fmin-hz", type=float, default=50.0)
    parser.add_argument("--fmax-hz", type=float, default=500.0)
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print extraction progress every N recordings. Use 0 to disable.",
    )
    parser.add_argument(
        "--pitch-method",
        choices=("auto", "pyin", "yin", "autocorr"),
        default="auto",
        help="Use librosa.pyin when available; autocorr is the dependency-light fallback.",
    )
    return parser.parse_args()


_LIBROSA = None
_LIBROSA_ERROR: Exception | None = None


def load_librosa(required: bool = False):
    global _LIBROSA, _LIBROSA_ERROR
    if _LIBROSA is not None:
        return _LIBROSA
    if _LIBROSA_ERROR is not None:
        if required:
            raise RuntimeError(f"librosa could not be imported: {_LIBROSA_ERROR}") from _LIBROSA_ERROR
        return None

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    for name, value in (("complex", complex), ("float", float), ("int", int)):
        if not hasattr(np, name):
            setattr(np, name, value)

    try:
        import librosa  # type: ignore

        _LIBROSA = librosa
        return _LIBROSA
    except Exception as exc:  # pragma: no cover - depends on local env
        _LIBROSA_ERROR = exc
        if required:
            raise RuntimeError(f"librosa could not be imported: {exc}") from exc
        return None


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_int(value: object, default: int = 0) -> int:
    try:
        return int(str(value))
    except Exception:
        return default


def finite_float(value: object) -> float:
    try:
        out = float(value)
    except Exception:
        return np.nan
    return out if np.isfinite(out) else np.nan


def parse_textgrid(path: Path) -> dict[str, list[Interval]]:
    tiers: dict[str, list[Interval]] = {}
    current_tier: str | None = None
    current_xmin: float | None = None
    current_xmax: float | None = None

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("name = "):
                current_tier = line.split("=", 1)[1].strip().strip('"')
                tiers.setdefault(current_tier, [])
                current_xmin = None
                current_xmax = None
            elif current_tier is not None and line.startswith("xmin = "):
                current_xmin = float(line.split("=", 1)[1].strip())
            elif current_tier is not None and line.startswith("xmax = "):
                current_xmax = float(line.split("=", 1)[1].strip())
            elif current_tier is not None and line.startswith("text = "):
                label = line.split("=", 1)[1].strip().strip('"')
                if current_xmin is not None and current_xmax is not None:
                    tiers[current_tier].append(Interval(current_xmin, current_xmax, label))
                current_xmin = None
                current_xmax = None
    return tiers


def choose_word_tier(tiers: dict[str, list[Interval]]) -> list[Interval]:
    lowered = {name.lower(): intervals for name, intervals in tiers.items()}
    if "words" in lowered:
        return lowered["words"]
    if "word" in lowered:
        return lowered["word"]
    for name, intervals in lowered.items():
        if "word" in name:
            return intervals
    raise KeyError(f"No word tier found. Available tiers: {', '.join(sorted(tiers))}")


def choose_phone_tier(tiers: dict[str, list[Interval]]) -> list[Interval]:
    lowered = {name.lower(): intervals for name, intervals in tiers.items()}
    if "phones" in lowered:
        return lowered["phones"]
    if "phone" in lowered:
        return lowered["phone"]
    for name, intervals in lowered.items():
        if "phone" in name:
            return intervals
    raise KeyError(f"No phone tier found. Available tiers: {', '.join(sorted(tiers))}")


def speech_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    return [interval for interval in intervals if interval.text.strip().lower() not in NON_SPEECH]


def normalize_label(value: str) -> str:
    return "".join(value.split()).lower()


def strip_tone_marks(text: str) -> str:
    decomposed = unicodedata.normalize("NFD", text)
    kept = [ch for ch in decomposed if ch not in {COMBINING_ACUTE, COMBINING_GRAVE}]
    return unicodedata.normalize("NFC", "".join(kept))


def normalize_mfa_text(text: str) -> str:
    text = strip_tone_marks(text).lower().replace(RIGHT_SINGLE_QUOTE, "'")
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


def textgrid_index(aligned_dir: Path) -> dict[str, Path]:
    paths = sorted(aligned_dir.rglob("*.TextGrid"))
    out: dict[str, Path] = {}
    for path in paths:
        out.setdefault(path.stem, path)
    return out


def load_pilot_manifest(path: Path, max_recordings: int = 0) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    if max_recordings > 0:
        rows = rows[:max_recordings]
    return {row["recording_id"]: row for row in rows if row.get("recording_id")}


def load_tone_units(
    tone_units_path: Path,
    pilot_rows: dict[str, dict[str, str]],
    include_recommended_exclude: bool,
) -> tuple[list[dict[str, object]], str]:
    recording_ids = set(pilot_rows)
    if not tone_units_path.exists():
        raise FileNotFoundError(f"Yoruba tone-unit CSV does not exist: {tone_units_path}")
    rows: list[dict[str, object]] = []
    with tone_units_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("recording_id") not in recording_ids:
                continue
            if parse_bool(row.get("recommended_exclude", "")) and not include_recommended_exclude:
                continue
            row["word_index"] = parse_int(row.get("word_index"))
            row["word_unit_index"] = parse_int(row.get("word_unit_index"))
            row["recording_unit_index"] = parse_int(row.get("recording_unit_index"))
            row["word_tone_unit_count"] = parse_int(row.get("word_tone_unit_count"), 0)
            row["recommended_exclude"] = parse_bool(row.get("recommended_exclude", ""))
            rows.append(row)
    return rows, "loaded_yoruba_tone_units_csv"


def audio_path_for_recording(manifest: dict[str, str]) -> Path:
    for col in ("mfa_audio_file", "source_audio_file", "file"):
        value = manifest.get(col, "")
        if value:
            return Path(value)
    return Path("")


def compute_track_features(
    audio_path: Path,
    frame_ms: float,
    hop_ms: float,
    fmin_hz: float,
    fmax_hz: float,
    pitch_method: str,
) -> dict[str, object]:
    if pitch_method == "autocorr":
        return compute_track_features_autocorr(audio_path, frame_ms, hop_ms, fmin_hz, fmax_hz)

    librosa = load_librosa(required=pitch_method in {"pyin", "yin"})
    if librosa is None:
        return compute_track_features_autocorr(audio_path, frame_ms, hop_ms, fmin_hz, fmax_hz)

    try:
        y, sr = librosa.load(str(audio_path), sr=None, mono=True)
    except Exception:
        if pitch_method in {"pyin", "yin"}:
            raise
        return compute_track_features_autocorr(audio_path, frame_ms, hop_ms, fmin_hz, fmax_hz)
    y = np.asarray(y, dtype=float)
    frame_length = max(64, int(round(frame_ms * sr / 1000.0)))
    hop_length = max(16, int(round(hop_ms * sr / 1000.0)))

    method_used = "pyin"
    if pitch_method in {"auto", "pyin"} and hasattr(librosa, "pyin"):
        try:
            f0, voiced_flag, voiced_prob = librosa.pyin(
                y,
                fmin=fmin_hz,
                fmax=fmax_hz,
                sr=sr,
                frame_length=frame_length,
                hop_length=hop_length,
            )
        except Exception:
            if pitch_method == "pyin":
                raise
            f0, voiced_flag, voiced_prob = None, None, None
    else:
        f0, voiced_flag, voiced_prob = None, None, None

    if f0 is None:
        method_used = "yin"
        f0 = librosa.yin(
            y,
            fmin=fmin_hz,
            fmax=fmax_hz,
            sr=sr,
            frame_length=frame_length,
            hop_length=hop_length,
        )
        voiced_flag = np.isfinite(f0)
        voiced_prob = np.full_like(f0, np.nan, dtype=float)

    rms = librosa.feature.rms(y=y, frame_length=frame_length, hop_length=hop_length, center=True)[0]
    rms_db = 20.0 * np.log10(np.maximum(rms, np.finfo(float).eps))
    times = librosa.frames_to_time(np.arange(len(f0)), sr=sr, hop_length=hop_length)
    usable = min(len(times), len(f0), len(rms_db))

    return {
        "sample_rate_hz": sr,
        "pitch_method": method_used,
        "frame_ms": frame_ms,
        "hop_ms": hop_ms,
        "times": times[:usable],
        "f0": np.asarray(f0[:usable], dtype=float),
        "voiced_flag": np.asarray(voiced_flag[:usable], dtype=bool),
        "voiced_prob": np.asarray(voiced_prob[:usable], dtype=float),
        "rms_db": np.asarray(rms_db[:usable], dtype=float),
    }


def frame_audio(y: np.ndarray, frame_length: int, hop_length: int) -> tuple[np.ndarray, np.ndarray]:
    if y.size == 0:
        return np.empty((0, frame_length), dtype=float), np.empty(0, dtype=int)
    if y.size < frame_length:
        y = np.pad(y, (0, frame_length - y.size))
    starts = np.arange(0, max(1, y.size - frame_length + 1), hop_length)
    if starts.size == 0:
        starts = np.asarray([0])
    frames = []
    for start in starts:
        frame = y[start : start + frame_length]
        if frame.size < frame_length:
            frame = np.pad(frame, (0, frame_length - frame.size))
        frames.append(frame)
    return np.vstack(frames), starts


def autocorr_f0(y: np.ndarray, sr: int, fmin_hz: float, fmax_hz: float) -> tuple[float, float]:
    if y.size < max(16, int(sr / fmin_hz)):
        return np.nan, np.nan
    y = y.astype(float)
    y = y - np.mean(y)
    energy = float(np.dot(y, y))
    if energy <= np.finfo(float).eps:
        return np.nan, 0.0
    y = y * np.hanning(len(y))
    ac = np.correlate(y, y, mode="full")[len(y) - 1 :]
    if ac[0] <= np.finfo(float).eps:
        return np.nan, 0.0
    min_lag = max(1, int(math.floor(sr / fmax_hz)))
    max_lag = min(len(ac) - 1, int(math.ceil(sr / fmin_hz)))
    if max_lag <= min_lag:
        return np.nan, np.nan
    segment = ac[min_lag : max_lag + 1] / ac[0]
    peak_offset = int(np.argmax(segment))
    peak = float(segment[peak_offset])
    lag = min_lag + peak_offset
    if peak <= 0:
        return np.nan, peak
    return float(sr / lag), peak


def compute_track_features_autocorr(
    audio_path: Path,
    frame_ms: float,
    hop_ms: float,
    fmin_hz: float,
    fmax_hz: float,
) -> dict[str, object]:
    y, sr = sf.read(audio_path, dtype="float32", always_2d=True)
    if y.shape[1] > 1:
        y = y.mean(axis=1)
    else:
        y = y[:, 0]
    y = np.asarray(y, dtype=float)
    frame_length = max(64, int(round(frame_ms * sr / 1000.0)))
    hop_length = max(16, int(round(hop_ms * sr / 1000.0)))
    frames, starts = frame_audio(y, frame_length, hop_length)
    times = (starts + frame_length / 2.0) / sr
    f0 = np.full(frames.shape[0], np.nan, dtype=float)
    voiced_prob = np.full(frames.shape[0], np.nan, dtype=float)
    rms_db = np.full(frames.shape[0], np.nan, dtype=float)
    for i, frame in enumerate(frames):
        f0_value, strength = autocorr_f0(frame, sr, fmin_hz, fmax_hz)
        f0[i] = f0_value
        voiced_prob[i] = strength
        rms = math.sqrt(float(np.mean(frame * frame)))
        rms_db[i] = 20.0 * math.log10(max(rms, np.finfo(float).eps))
    voiced_flag = np.isfinite(f0) & (voiced_prob >= 0.25)
    return {
        "sample_rate_hz": sr,
        "pitch_method": "autocorr",
        "frame_ms": frame_ms,
        "hop_ms": hop_ms,
        "times": times,
        "f0": f0,
        "voiced_flag": voiced_flag,
        "voiced_prob": voiced_prob,
        "rms_db": rms_db,
    }


def median_or_nan(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.median(values))


def mean_or_nan(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.mean(values))


def sd_or_nan(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan
    return float(np.std(values, ddof=1))


def slope_or_nan(times: np.ndarray, values: np.ndarray) -> float:
    ok = np.isfinite(times) & np.isfinite(values)
    if int(ok.sum()) < 2:
        return np.nan
    rel = times[ok] - float(times[ok][0])
    if np.nanmax(rel) <= 0:
        return np.nan
    slope, _ = np.polyfit(rel, values[ok], 1)
    return float(slope)


def segment_median(times: np.ndarray, values: np.ndarray, start: float, end: float, part: int) -> float:
    width = (end - start) / 3.0
    lo = start + part * width
    hi = start + (part + 1) * width
    if part == 2:
        mask = (times >= lo) & (times <= hi)
    else:
        mask = (times >= lo) & (times < hi)
    return median_or_nan(values[mask])


def token_features(
    track: dict[str, object],
    start: float,
    end: float,
    min_voiced_frames: int,
) -> tuple[dict[str, float], str]:
    times = track["times"]
    f0 = track["f0"]
    voiced_flag = track["voiced_flag"]
    voiced_prob = track["voiced_prob"]
    rms_db = track["rms_db"]
    assert isinstance(times, np.ndarray)
    assert isinstance(f0, np.ndarray)
    assert isinstance(voiced_flag, np.ndarray)
    assert isinstance(voiced_prob, np.ndarray)
    assert isinstance(rms_db, np.ndarray)

    mask = (times >= start) & (times <= end)
    total_frames = int(mask.sum())
    if total_frames == 0:
        return {}, "no_pitch_frames"

    token_f0 = f0[mask]
    token_voiced = voiced_flag[mask] & np.isfinite(token_f0)
    voiced_frames = int(token_voiced.sum())
    if voiced_frames < min_voiced_frames:
        return {}, "too_few_voiced_frames"

    token_times = times[mask]
    token_st = 12.0 * np.log2(token_f0[token_voiced])
    voiced_times = token_times[token_voiced]
    token_rms = rms_db[mask]

    features = {
        "f0_total_frames": float(total_frames),
        "f0_voiced_frames": float(voiced_frames),
        "f0_voiced_prop": float(voiced_frames / total_frames),
        "f0_voicing_prob_mean": mean_or_nan(voiced_prob[mask]),
        "f0_mean_hz": mean_or_nan(token_f0[token_voiced]),
        "f0_median_hz": median_or_nan(token_f0[token_voiced]),
        "f0_sd_hz": sd_or_nan(token_f0[token_voiced]),
        "f0_min_hz": float(np.nanmin(token_f0[token_voiced])),
        "f0_max_hz": float(np.nanmax(token_f0[token_voiced])),
        "f0_range_hz": float(np.nanmax(token_f0[token_voiced]) - np.nanmin(token_f0[token_voiced])),
        "f0_mean_st": mean_or_nan(token_st),
        "f0_median_st": median_or_nan(token_st),
        "f0_sd_st": sd_or_nan(token_st),
        "f0_range_st": float(np.nanmax(token_st) - np.nanmin(token_st)),
        "f0_start_st": segment_median(voiced_times, token_st, start, end, 0),
        "f0_mid_st": segment_median(voiced_times, token_st, start, end, 1),
        "f0_end_st": segment_median(voiced_times, token_st, start, end, 2),
        "f0_slope_st_per_s": slope_or_nan(voiced_times, token_st),
        "rms_db_mean": mean_or_nan(token_rms),
        "rms_db_sd": sd_or_nan(token_rms),
        "rms_db_start": segment_median(token_times, token_rms, start, end, 0),
        "rms_db_mid": segment_median(token_times, token_rms, start, end, 1),
        "rms_db_end": segment_median(token_times, token_rms, start, end, 2),
        "rms_db_slope_per_s": slope_or_nan(token_times, token_rms),
    }
    features["f0_delta_st"] = features["f0_end_st"] - features["f0_start_st"]
    features["rms_db_delta"] = features["rms_db_end"] - features["rms_db_start"]
    return features, ""


def word_unit_count(row: dict[str, object]) -> int:
    explicit = parse_int(row.get("word_tone_unit_count", 0), 0)
    if explicit > 0:
        return explicit
    return max(1, len(str(row.get("word", ""))))


def phone_label_for_quality(quality: str) -> str:
    if quality == "syllabic_n":
        return "n"
    return quality


def phone_intervals_for_word(phone_intervals: list[Interval], word_interval: Interval) -> list[Interval]:
    tol = 0.005
    return [
        phone
        for phone in phone_intervals
        if phone.xmax > word_interval.xmin + tol and phone.xmin < word_interval.xmax - tol
    ]


def align_tone_units_to_phones(
    tone_rows: list[dict[str, object]],
    phone_intervals: list[Interval],
) -> tuple[dict[str, Interval], Counter[str]]:
    out: dict[str, Interval] = {}
    skipped: Counter[str] = Counter()
    cursor = 0
    for row in sorted(tone_rows, key=lambda item: parse_int(item.get("word_unit_index"))):
        target = phone_label_for_quality(str(row.get("vowel_quality", "")))
        if not target:
            skipped["missing_vowel_quality"] += 1
            continue
        matched: Interval | None = None
        while cursor < len(phone_intervals):
            phone = phone_intervals[cursor]
            cursor += 1
            label = phone.text.strip().lower()
            if label == target:
                matched = phone
                break
        if matched is None:
            skipped["phone_sequence_mismatch"] += 1
            continue
        out[str(row.get("token_id", ""))] = matched
    return out, skipped


def add_centered_and_zscores(rows: list[dict[str, object]]) -> list[str]:
    st_cols = [
        "f0_mean_st",
        "f0_median_st",
        "f0_start_st",
        "f0_mid_st",
        "f0_end_st",
    ]
    rms_cols = ["rms_db_mean", "rms_db_start", "rms_db_mid", "rms_db_end"]
    by_speaker: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_speaker[str(row.get("speaker", ""))].append(row)

    for speaker_rows in by_speaker.values():
        f0_center = median_or_nan(np.asarray([finite_float(row.get("f0_mean_st")) for row in speaker_rows]))
        rms_center = median_or_nan(np.asarray([finite_float(row.get("rms_db_mean")) for row in speaker_rows]))
        for row in speaker_rows:
            for col in st_cols:
                row[f"{col}_speaker_centered"] = finite_float(row.get(col)) - f0_center
            for col in rms_cols:
                row[f"{col}_speaker_centered"] = finite_float(row.get(col)) - rms_center

    metadata_numeric = {
        "word_index",
        "word_unit_index",
        "recording_unit_index",
        "word_tone_unit_count",
        "start",
        "end",
        "word_start",
        "word_end",
        "phone_start",
        "phone_end",
        "sample_rate_hz",
        "frame_ms",
        "hop_ms",
    }
    numeric_cols: list[str] = []
    for col in rows[0]:
        if col in metadata_numeric:
            continue
        values = [finite_float(row.get(col)) for row in rows]
        if any(np.isfinite(value) for value in values):
            numeric_cols.append(col)

    for col in numeric_cols:
        values = np.asarray([finite_float(row.get(col)) for row in rows], dtype=float)
        mean = float(np.nanmean(values))
        sd = float(np.nanstd(values, ddof=1)) if np.isfinite(values).sum() > 1 else np.nan
        z_col = f"z_{col}"
        for row, value in zip(rows, values):
            row[z_col] = (value - mean) / sd if np.isfinite(value) and np.isfinite(sd) and sd > 0 else np.nan
    return [f"z_{col}" for col in numeric_cols]


def extract_features(args: argparse.Namespace) -> tuple[list[dict[str, object]], list[dict[str, object]], str]:
    pilot_rows = load_pilot_manifest(args.pilot_manifest, args.max_recordings)
    tg_by_id = textgrid_index(args.aligned_dir)
    tone_rows, tone_source = load_tone_units(args.tone_units, pilot_rows, args.include_recommended_exclude)
    tone_by_recording: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in tone_rows:
        if row.get("recording_id") in pilot_rows:
            tone_by_recording[str(row["recording_id"])].append(row)
    for rows in tone_by_recording.values():
        rows.sort(key=lambda row: (parse_int(row.get("word_index")), parse_int(row.get("word_unit_index"))))

    out_rows: list[dict[str, object]] = []
    skipped = Counter()
    librosa = load_librosa()
    _ = librosa
    start_time = time.monotonic()
    recording_items = sorted(pilot_rows.items())

    def report_progress(processed: int) -> None:
        if args.progress_every <= 0:
            return
        if processed != 1 and processed % args.progress_every != 0 and processed != len(recording_items):
            return
        top_skips = ", ".join(f"{reason}={n}" for reason, n in skipped.most_common(5)) or "none"
        elapsed_s = time.monotonic() - start_time
        print(
            "[yoruba-extract] "
            f"processed_recordings={processed}/{len(recording_items)} "
            f"tokens_written={len(out_rows)} "
            f"tokens_skipped={sum(skipped.values())} "
            f"top_skips={top_skips} "
            f"elapsed_s={elapsed_s:.1f}",
            file=sys.stderr,
            flush=True,
        )

    for processed, (recording_id, manifest) in enumerate(recording_items, start=1):
        try:
            tg_path = tg_by_id.get(recording_id)
            if tg_path is None:
                skipped["no_textgrid"] += len(tone_by_recording.get(recording_id, []))
                continue
            audio_path = audio_path_for_recording(manifest)
            if not audio_path.exists():
                skipped["missing_audio"] += len(tone_by_recording.get(recording_id, []))
                continue
            try:
                tiers = parse_textgrid(tg_path)
                word_intervals = speech_intervals(choose_word_tier(tiers))
                phone_intervals = speech_intervals(choose_phone_tier(tiers))
            except Exception:
                skipped["missing_textgrid_tier"] += len(tone_by_recording.get(recording_id, []))
                continue
            track = compute_track_features(
                audio_path,
                frame_ms=args.frame_ms,
                hop_ms=args.hop_ms,
                fmin_hz=args.fmin_hz,
                fmax_hz=args.fmax_hz,
                pitch_method=args.pitch_method,
            )
            rows_by_word: dict[int, list[dict[str, object]]] = defaultdict(list)
            for row in tone_by_recording.get(recording_id, []):
                rows_by_word[parse_int(row.get("word_index"))].append(row)

            token_phone_intervals: dict[str, Interval] = {}
            word_interval_by_index: dict[int, Interval] = {}
            for word_index, word_rows in sorted(rows_by_word.items()):
                if word_index <= 0 or word_index > len(word_intervals):
                    skipped["no_word_interval"] += len(word_rows)
                    continue
                word_interval = word_intervals[word_index - 1]
                word_interval_by_index[word_index] = word_interval
                word_phones = phone_intervals_for_word(phone_intervals, word_interval)
                if not word_phones:
                    skipped["no_phone_intervals_in_word"] += len(word_rows)
                    continue
                aligned, phone_skips = align_tone_units_to_phones(word_rows, word_phones)
                token_phone_intervals.update(aligned)
                skipped.update(phone_skips)

            for tone_row in tone_by_recording.get(recording_id, []):
                word_index = parse_int(tone_row.get("word_index"))
                word_unit_index = parse_int(tone_row.get("word_unit_index"))
                word_interval = word_interval_by_index.get(word_index)
                phone_interval = token_phone_intervals.get(str(tone_row.get("token_id", "")))
                if word_interval is None:
                    continue
                if phone_interval is None:
                    continue
                unit_count = word_unit_count(tone_row)
                if word_unit_index <= 0 or word_unit_index > unit_count:
                    skipped["bad_word_unit_index"] += 1
                    continue
                unit_start = phone_interval.xmin
                unit_end = phone_interval.xmax
                duration_ms = (unit_end - unit_start) * 1000.0
                quality_flags = ["mfa_phone_interval"]
                expected_word = normalize_mfa_text(str(tone_row.get("word", "")))
                if normalize_label(word_interval.text) != normalize_label(expected_word):
                    quality_flags.append("word_label_mismatch")
                if parse_bool(tone_row.get("recommended_exclude", "")):
                    quality_flags.append("recommended_exclude")
                if duration_ms < args.min_token_ms:
                    skipped["short_token_interval"] += 1
                    continue

                features, skip_reason = token_features(track, unit_start, unit_end, args.min_voiced_frames)
                if skip_reason:
                    skipped[skip_reason] += 1
                    continue

                out = {
                    "token_id": tone_row.get("token_id", ""),
                    "domain": "tone",
                    "language": "Yoruba",
                    "source_corpus": "OpenSLR86",
                    "speaker": manifest.get("speaker", tone_row.get("speaker", "")),
                    "speaker_code": tone_row.get("speaker_code", ""),
                    "speaker_source": tone_row.get("speaker_source", ""),
                    "sex": manifest.get("sex", tone_row.get("sex", "")),
                    "recording_id": recording_id,
                    "file": str(audio_path),
                    "textgrid": str(tg_path),
                    "word": tone_row.get("word", ""),
                    "word_index": word_index,
                    "orthographic_unit": tone_row.get("orthographic_unit", ""),
                    "orthographic_base": tone_row.get("orthographic_base", ""),
                    "vowel_quality": tone_row.get("vowel_quality", ""),
                    "tbu_type": tone_row.get("tbu_type", ""),
                    "word_unit_index": word_unit_index,
                    "category": tone_row.get("category", ""),
                    "control_group": tone_row.get("control_group", tone_row.get("vowel_quality", "")),
                    "tone_unit_position": tone_row.get("tone_unit_position", ""),
                    "vowel_quality_control": tone_row.get("vowel_quality", ""),
                    "nasal_oral_status": tone_row.get("nasal_oral_status", ""),
                    "preceding_letter": tone_row.get("preceding_letter", ""),
                    "following_letter": tone_row.get("following_letter", ""),
                    "following_n_letter": tone_row.get("following_n_letter", ""),
                    "recording_unit_index": parse_int(tone_row.get("recording_unit_index")),
                    "word_tone_unit_count": unit_count,
                    "word_interval_label": word_interval.text,
                    "phone_interval_label": phone_interval.text,
                    "start": unit_start,
                    "end": unit_end,
                    "word_start": word_interval.xmin,
                    "word_end": word_interval.xmax,
                    "phone_start": phone_interval.xmin,
                    "phone_end": phone_interval.xmax,
                    "duration_ms": duration_ms,
                    "word_duration_ms": word_interval.duration * 1000.0,
                    "tone_unit_relative_position": (word_unit_index - 0.5) / unit_count,
                    "sample_rate_hz": track["sample_rate_hz"],
                    "frame_ms": track["frame_ms"],
                    "hop_ms": track["hop_ms"],
                    "pitch_method": track["pitch_method"],
                    "measurement_method": "mfa_phone_interval_yoruba_tone_unit_librosa_f0",
                    "quality_flag": ";".join(quality_flags),
                    "recommended_exclude": parse_bool(tone_row.get("recommended_exclude", "")),
                }
                out.update(features)
                out_rows.append(out)
        finally:
            report_progress(processed)

    summary = [
        {"metric": "pilot_recordings", "value": len(pilot_rows)},
        {"metric": "textgrids_found", "value": len(tg_by_id)},
        {"metric": "tone_source", "value": tone_source},
        {"metric": "tone_units_loaded", "value": len(tone_rows)},
        {"metric": "tokens_written", "value": len(out_rows)},
    ]
    for reason, n in sorted(skipped.items()):
        summary.append({"metric": f"tokens_skipped_{reason}", "value": n})
    for tone, n in sorted(Counter(str(row.get("category", "")) for row in out_rows).items()):
        summary.append({"metric": f"tokens_{tone}", "value": n})
    return out_rows, summary, tone_source


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


def write_feature_sets(path: Path) -> None:
    f0_static = [
        "z_f0_mean_st_speaker_centered",
        "z_f0_sd_st",
        "z_f0_range_st",
        "z_f0_voiced_prop",
    ]
    f0_contour = [
        "z_f0_start_st_speaker_centered",
        "z_f0_mid_st_speaker_centered",
        "z_f0_end_st_speaker_centered",
        "z_f0_slope_st_per_s",
    ]
    energy_duration = [
        "z_duration_ms",
        "z_word_duration_ms",
        "z_tone_unit_relative_position",
        "z_rms_db_mean_speaker_centered",
        "z_rms_db_delta",
        "z_rms_db_slope_per_s",
    ]
    rows = [
        ("f0_static", f0_static, "Speaker-centered F0 level, spread, and voiced-frame coverage."),
        ("f0_contour", f0_contour, "Speaker-centered F0 thirds plus token-level F0 movement."),
        ("energy_duration", energy_duration, "Duration, position, and intensity-envelope controls."),
        ("tone_prosody", f0_static + f0_contour + energy_duration, "Combined pilot Yoruba tone bundle."),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature_set", "features", "note"])
        for name, features, note in rows:
            writer.writerow([name, "; ".join(features), note])


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir
    tokens_out = args.tokens_out or out_dir / "yoruba_slr86_tone_features.csv"
    summary_out = args.summary_out or out_dir / "yoruba_slr86_tone_feature_summary.csv"
    feature_sets_out = args.feature_sets_out or out_dir / "yoruba_slr86_tone_feature_sets.csv"

    rows, summary, _ = extract_features(args)
    if not rows:
        write_csv(summary_out, summary, ["metric", "value"])
        raise RuntimeError(f"No tone feature rows were extracted. Wrote summary to {summary_out}")
    add_centered_and_zscores(rows)
    write_csv(tokens_out, rows)
    write_csv(summary_out, summary, ["metric", "value"])
    write_feature_sets(feature_sets_out)

    print(f"tokens_written={len(rows)}")
    print(f"tokens_out={tokens_out}")
    print(f"summary_out={summary_out}")
    print(f"feature_sets_out={feature_sets_out}")


if __name__ == "__main__":
    main()
