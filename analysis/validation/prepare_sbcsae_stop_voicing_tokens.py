#!/usr/bin/env python3
"""Prepare SBCSAE stop-voicing token features for the common validator.

The script reads MFA TextGrids with `words` and `phones` tiers, matches them to
SBCSAE WAV files by basename, extracts stop tokens for /p b/, /t d/, and /k g/,
and writes a P0-compatible token-feature CSV.

This is a pilot adapter: it uses aligned phone intervals and short acoustic
windows, but it does not estimate true burst onset, VOT, or prevoicing onset.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

# Avoid noisy CPU feature probes from pyarrow when pandas imports optional
# Arrow support under the managed sandbox.
os.environ.setdefault("ARROW_USER_SIMD_LEVEL", "NONE")

import numpy as np
import pandas as pd
from scipy.fftpack import dct
import soundfile as sf


VOLUME_ROOT = Path("/Volumes/Corpus_Studies/Corpora/English/SBCSAE")
DEFAULT_TEXTGRID_DIR = VOLUME_ROOT / "SBCSAE_MFA_Aligned_Final"
DEFAULT_OUT_DIR = Path("analysis/validation/outputs/sbcsae_stop_voicing")

STOP_MAP = {
    "p": ("labial", "voiceless", "/p/~/b/", {"p", "pʰ", "pʲ"}),
    "b": ("labial", "voiced", "/p/~/b/", {"b", "bʲ"}),
    "t": ("coronal", "voiceless", "/t/~/d/", {"t", "tʰ", "tʲ", "t̪", "tʷ"}),
    "d": ("coronal", "voiced", "/t/~/d/", {"d", "dʲ", "d̪"}),
    "k": ("dorsal", "voiceless", "/k/~/g/", {"k", "kʰ", "kʷ"}),
    "g": ("dorsal", "voiced", "/k/~/g/", {"ɡ", "ɡʷ"}),
}

LABEL_TO_STOP = {}
for canonical, (_, _, _, labels) in STOP_MAP.items():
    for label in labels:
        LABEL_TO_STOP[label] = canonical

VOWEL_CHARS = set("aeiouæɑɒɔɛɪʊəɐɚɝʉɜ")
NON_SPEECH = {"", "spn", "sil", "<eps>"}


@dataclass(frozen=True)
class Interval:
    xmin: float
    xmax: float
    text: str

    @property
    def duration(self) -> float:
        return self.xmax - self.xmin

    @property
    def midpoint(self) -> float:
        return (self.xmin + self.xmax) / 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--textgrid-dir", type=Path, default=DEFAULT_TEXTGRID_DIR)
    parser.add_argument("--audio-dir", type=Path, default=VOLUME_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--tokens-out", type=Path, default=None)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--feature-sets-out", type=Path, default=None)
    parser.add_argument("--max-per-place-category", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--min-stop-ms", type=float, default=10.0)
    parser.add_argument("--max-following-vowel-lag-ms", type=float, default=200.0)
    parser.add_argument("--allow-no-following-vowel", action="store_true")
    parser.add_argument("--pre-ms", type=float, default=30.0)
    parser.add_argument("--post-ms", type=float, default=80.0)
    parser.add_argument("--mfcc-count", type=int, default=13)
    parser.add_argument(
        "--feature-engine",
        choices=("auto", "librosa", "scipy"),
        default="auto",
        help="Prefer librosa when available; scipy is a dependency-light fallback.",
    )
    return parser.parse_args()


_LIBROSA = None
_LIBROSA_ERROR: Exception | None = None


def load_librosa(required: bool = False):
    """Load librosa with compatibility shims for the local Anaconda version."""
    global _LIBROSA, _LIBROSA_ERROR
    if _LIBROSA is not None:
        return _LIBROSA
    if _LIBROSA_ERROR is not None:
        if required:
            raise RuntimeError(f"librosa could not be imported: {_LIBROSA_ERROR}") from _LIBROSA_ERROR
        return None

    # The local librosa is old enough to use NumPy aliases removed in newer
    # NumPy, and numba caching is not available in the managed sandbox.
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


def nonempty_index(intervals: list[Interval]) -> list[int]:
    return [i for i, interval in enumerate(intervals) if interval.text not in NON_SPEECH]


def is_vowel(label: str) -> bool:
    if label in NON_SPEECH:
        return False
    return any(ch in VOWEL_CHARS for ch in label)


def surrounding_nonempty(intervals: list[Interval], index: int) -> tuple[Interval | None, Interval | None]:
    prev_interval = None
    for i in range(index - 1, -1, -1):
        if intervals[i].text not in NON_SPEECH:
            prev_interval = intervals[i]
            break
    next_interval = None
    for i in range(index + 1, len(intervals)):
        if intervals[i].text not in NON_SPEECH:
            next_interval = intervals[i]
            break
    return prev_interval, next_interval


def first_following_vowel(
    intervals: list[Interval],
    index: int,
    max_lag_s: float,
) -> Interval | None:
    stop_end = intervals[index].xmax
    for i in range(index + 1, len(intervals)):
        interval = intervals[i]
        if interval.text in NON_SPEECH:
            continue
        if interval.xmin - stop_end > max_lag_s:
            return None
        if is_vowel(interval.text):
            return interval
    return None


def word_at_time(words: list[Interval], time_s: float) -> str:
    # TextGrids are time-ordered, but linear scan is fine at this scale.
    for interval in words:
        if interval.xmin <= time_s < interval.xmax and interval.text:
            return interval.text
    return ""


def find_audio(audio_dir: Path, file_id: str) -> Path | None:
    for ext in (".wav", ".WAV", ".flac", ".FLAC"):
        path = audio_dir / f"{file_id}{ext}"
        if path.exists():
            return path
    return None


def collect_candidates(args: argparse.Namespace) -> pd.DataFrame:
    textgrid_paths = sorted(args.textgrid_dir.glob("*.TextGrid"))
    if not textgrid_paths:
        raise FileNotFoundError(f"No TextGrid files found in {args.textgrid_dir}")

    rows = []
    max_lag_s = args.max_following_vowel_lag_ms / 1000.0
    min_stop_s = args.min_stop_ms / 1000.0

    for tg_path in textgrid_paths:
        file_id = tg_path.stem
        audio_path = find_audio(args.audio_dir, file_id)
        if audio_path is None:
            continue
        tiers = parse_textgrid(tg_path)
        phones = tiers.get("phones")
        words = tiers.get("words", [])
        if not phones:
            continue
        for i, phone in enumerate(phones):
            canonical = LABEL_TO_STOP.get(phone.text)
            if canonical is None or phone.duration < min_stop_s:
                continue
            place, category, contrast_pair, _ = STOP_MAP[canonical]
            prev_phone, next_phone = surrounding_nonempty(phones, i)
            following_vowel = first_following_vowel(phones, i, max_lag_s)
            if following_vowel is None and not args.allow_no_following_vowel:
                continue
            rows.append(
                {
                    "token_id": f"{file_id}_{i:06d}_{phone.text}",
                    "domain": "stop_voicing",
                    "language": "English",
                    "source_corpus": "SBCSAE",
                    # Speaker tiers are not present in these MFA TextGrids.
                    "speaker": file_id,
                    "recording_id": file_id,
                    "file": str(audio_path),
                    "word": word_at_time(words, phone.midpoint),
                    "phone": canonical,
                    "phone_variant": phone.text,
                    "category": category,
                    "control_group": place,
                    "contrast_pair": contrast_pair,
                    "start": phone.xmin,
                    "end": phone.xmax,
                    "stop_duration_ms": phone.duration * 1000.0,
                    "preceding_phone": prev_phone.text if prev_phone else "",
                    "following_phone": next_phone.text if next_phone else "",
                    "prev_phone_duration_ms": prev_phone.duration * 1000.0 if prev_phone else np.nan,
                    "next_phone_duration_ms": next_phone.duration * 1000.0 if next_phone else np.nan,
                    "following_vowel_phone": following_vowel.text if following_vowel else "",
                    "following_vowel_lag_ms": (following_vowel.xmin - phone.xmax) * 1000.0
                    if following_vowel
                    else np.nan,
                    "following_vowel_duration_ms": following_vowel.duration * 1000.0
                    if following_vowel
                    else np.nan,
                    "measurement_method": "mfa_phone_interval_plus_audio_windows",
                    "quality_flag": "speaker_unavailable_recording_id_used",
                }
            )
    if not rows:
        raise RuntimeError("No stop tokens were found after filtering.")
    return pd.DataFrame(rows)


def sample_candidates(df: pd.DataFrame, max_per_place_category: int, seed: int) -> pd.DataFrame:
    if max_per_place_category <= 0:
        return df.sort_values(["recording_id", "start", "phone_variant"]).reset_index(drop=True)
    pieces = []
    for _, group in df.groupby(["control_group", "category"], sort=True):
        pieces.append(
            group.sample(
                n=min(len(group), max_per_place_category),
                random_state=seed,
            )
        )
    sampled = pd.concat(pieces, ignore_index=True)
    sampled = sampled.sort_values(["recording_id", "start", "phone_variant"]).reset_index(drop=True)
    return sampled


def read_window(handle: sf.SoundFile, start_s: float, end_s: float) -> np.ndarray:
    sr = handle.samplerate
    start_frame = max(0, int(math.floor(start_s * sr)))
    end_frame = min(len(handle), int(math.ceil(end_s * sr)))
    frames = max(0, end_frame - start_frame)
    if frames <= 0:
        return np.asarray([], dtype=np.float32)
    handle.seek(start_frame)
    audio = handle.read(frames, dtype="float32", always_2d=True)
    if audio.shape[1] > 1:
        audio = audio.mean(axis=1)
    else:
        audio = audio[:, 0]
    return np.asarray(audio, dtype=np.float32)


def librosa_frame_params(y: np.ndarray, sr: int) -> tuple[int, int]:
    target = max(64, min(int(round(0.025 * sr)), max(1, len(y))))
    n_fft = 1 << int(math.floor(math.log2(target)))
    n_fft = max(64, min(2048, n_fft, max(64, len(y))))
    hop_length = max(16, n_fft // 4)
    return n_fft, hop_length


def audio_summary_features(
    y: np.ndarray,
    sr: int,
    prefix: str,
    librosa_mod=None,
) -> dict[str, float]:
    out = {
        f"{prefix}_rms_db": np.nan,
        f"{prefix}_zcr": np.nan,
        f"{prefix}_spectral_centroid_hz": np.nan,
        f"{prefix}_spectral_bandwidth_hz": np.nan,
        f"{prefix}_spectral_skewness": np.nan,
        f"{prefix}_spectral_kurtosis": np.nan,
        f"{prefix}_low_energy_prop": np.nan,
        f"{prefix}_mid_energy_prop": np.nan,
        f"{prefix}_high_energy_prop": np.nan,
        f"{prefix}_low_high_log_ratio": np.nan,
        f"{prefix}_spectral_flatness": np.nan,
        f"{prefix}_f0_hz": np.nan,
        f"{prefix}_voicing_strength": np.nan,
    }
    if y.size < 4:
        return out
    eps = np.finfo(float).eps
    y = y.astype(float)
    out[f"{prefix}_rms_db"] = 10.0 * math.log10(float(np.mean(y * y)) + eps)
    windowed = y * np.hanning(len(y))
    power = np.abs(np.fft.rfft(windowed)) ** 2
    freqs = np.fft.rfftfreq(len(windowed), d=1.0 / sr)
    total = float(power.sum())
    if total <= eps:
        return out
    out[f"{prefix}_low_energy_prop"] = float(power[(freqs >= 50) & (freqs < 1000)].sum() / total)
    out[f"{prefix}_mid_energy_prop"] = float(power[(freqs >= 1000) & (freqs < 3000)].sum() / total)
    out[f"{prefix}_high_energy_prop"] = float(power[freqs >= 3000].sum() / total)
    out[f"{prefix}_low_high_log_ratio"] = float(
        math.log((out[f"{prefix}_low_energy_prop"] + eps) / (out[f"{prefix}_high_energy_prop"] + eps))
    )
    centroid_direct = float(np.sum(freqs * power) / total)
    spread = np.sqrt(np.sum(((freqs - centroid_direct) ** 2) * power) / total)
    if spread > eps:
        standardized = (freqs - centroid_direct) / spread
        out[f"{prefix}_spectral_skewness"] = float(np.sum((standardized**3) * power) / total)
        out[f"{prefix}_spectral_kurtosis"] = float(np.sum((standardized**4) * power) / total)

    f0_hz, voicing_strength = autocorr_f0(y, sr)
    out[f"{prefix}_f0_hz"] = f0_hz
    out[f"{prefix}_voicing_strength"] = voicing_strength

    if librosa_mod is not None:
        n_fft, hop_length = librosa_frame_params(y, sr)
        try:
            out[f"{prefix}_zcr"] = float(
                np.nanmean(
                    librosa_mod.feature.zero_crossing_rate(
                        y,
                        frame_length=n_fft,
                        hop_length=hop_length,
                        center=True,
                    )
                )
            )
            out[f"{prefix}_spectral_centroid_hz"] = float(
                np.nanmean(
                    librosa_mod.feature.spectral_centroid(
                        y=y,
                        sr=sr,
                        n_fft=n_fft,
                        hop_length=hop_length,
                    )
                )
            )
            out[f"{prefix}_spectral_bandwidth_hz"] = float(
                np.nanmean(
                    librosa_mod.feature.spectral_bandwidth(
                        y=y,
                        sr=sr,
                        n_fft=n_fft,
                        hop_length=hop_length,
                    )
                )
            )
            out[f"{prefix}_spectral_flatness"] = float(
                np.nanmean(
                    librosa_mod.feature.spectral_flatness(
                        y=y,
                        n_fft=n_fft,
                        hop_length=hop_length,
                    )
                )
            )
            return out
        except Exception:
            pass

    out[f"{prefix}_zcr"] = float(np.mean(np.abs(np.diff(np.signbit(y)))))
    centroid = float(np.sum(freqs * power) / total)
    bandwidth = float(np.sqrt(np.sum(((freqs - centroid) ** 2) * power) / total))
    out[f"{prefix}_spectral_centroid_hz"] = centroid
    out[f"{prefix}_spectral_bandwidth_hz"] = bandwidth
    out[f"{prefix}_spectral_flatness"] = float(np.exp(np.mean(np.log(power + eps))) / (np.mean(power) + eps))
    return out


def autocorr_f0(y: np.ndarray, sr: int, fmin: float = 50.0, fmax: float = 400.0) -> tuple[float, float]:
    if y.size < max(16, int(sr / fmin)):
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
    min_lag = max(1, int(math.floor(sr / fmax)))
    max_lag = min(len(ac) - 1, int(math.ceil(sr / fmin)))
    if max_lag <= min_lag:
        return np.nan, np.nan
    segment = ac[min_lag : max_lag + 1] / ac[0]
    peak_offset = int(np.argmax(segment))
    peak = float(segment[peak_offset])
    lag = min_lag + peak_offset
    if peak <= 0:
        return np.nan, peak
    return float(sr / lag), peak


def formant_features(y: np.ndarray, sr: int, librosa_mod=None, prefix: str = "post80") -> dict[str, float]:
    out = {
        f"{prefix}_f1_hz": np.nan,
        f"{prefix}_f2_hz": np.nan,
        f"{prefix}_f3_hz": np.nan,
    }
    if librosa_mod is None or y.size < int(0.035 * sr):
        return out
    try:
        x = y.astype(float)
        x = x - np.mean(x)
        x = np.append(x[0], x[1:] - 0.97 * x[:-1])
        x = x * np.hamming(len(x))
        order = min(16, max(8, int(sr / 1000) + 2))
        coeffs = librosa_mod.lpc(x, order)
        roots = np.roots(coeffs)
        roots = roots[np.imag(roots) >= 0.01]
        freqs = np.sort(np.angle(roots) * (sr / (2 * np.pi)))
        freqs = freqs[(freqs > 90) & (freqs < 5500)]
        for i, freq in enumerate(freqs[:3], start=1):
            out[f"{prefix}_f{i}_hz"] = float(freq)
    except Exception:
        pass
    return out


def hz_to_mel(hz: np.ndarray) -> np.ndarray:
    return 2595.0 * np.log10(1.0 + hz / 700.0)


def mel_to_hz(mel: np.ndarray) -> np.ndarray:
    return 700.0 * (10.0 ** (mel / 2595.0) - 1.0)


def mel_filterbank(sr: int, n_fft: int, n_mels: int = 26) -> np.ndarray:
    low_mel = hz_to_mel(np.asarray([0.0]))[0]
    high_mel = hz_to_mel(np.asarray([sr / 2.0]))[0]
    mel_points = np.linspace(low_mel, high_mel, n_mels + 2)
    hz_points = mel_to_hz(mel_points)
    bins = np.floor((n_fft + 1) * hz_points / sr).astype(int)
    filters = np.zeros((n_mels, n_fft // 2 + 1), dtype=float)
    for m in range(1, n_mels + 1):
        left, center, right = bins[m - 1], bins[m], bins[m + 1]
        if center <= left:
            center = left + 1
        if right <= center:
            right = center + 1
        for k in range(left, min(center, filters.shape[1])):
            filters[m - 1, k] = (k - left) / (center - left)
        for k in range(center, min(right, filters.shape[1])):
            filters[m - 1, k] = (right - k) / (right - center)
    return filters


def frame_audio(y: np.ndarray, frame_length: int, hop_length: int) -> np.ndarray:
    if y.size < frame_length:
        y = np.pad(y, (0, frame_length - y.size))
    starts = np.arange(0, max(1, y.size - frame_length + 1), hop_length)
    if starts.size == 0:
        starts = np.asarray([0])
    frames = np.vstack([y[start : start + frame_length] for start in starts])
    if frames.shape[1] < frame_length:
        frames = np.pad(frames, ((0, 0), (0, frame_length - frames.shape[1])))
    return frames


def mfcc_means(y: np.ndarray, sr: int, n_mfcc: int, librosa_mod=None) -> np.ndarray:
    if y.size < 8:
        return np.full(n_mfcc, np.nan)
    y = y.astype(float)
    if librosa_mod is not None:
        n_fft, hop_length = librosa_frame_params(y, sr)
        try:
            coeffs = librosa_mod.feature.mfcc(
                y=y,
                sr=sr,
                n_mfcc=n_mfcc,
                n_fft=n_fft,
                hop_length=hop_length,
            )
            return np.asarray(coeffs.mean(axis=1), dtype=float)
        except Exception:
            pass

    frame_length = max(64, int(round(0.025 * sr)))
    hop_length = max(32, int(round(0.010 * sr)))
    n_fft = 1 << int(math.ceil(math.log2(frame_length)))
    frames = frame_audio(y, frame_length, hop_length)
    frames = frames * np.hamming(frame_length)
    power = np.abs(np.fft.rfft(frames, n=n_fft, axis=1)) ** 2
    filters = mel_filterbank(sr, n_fft)
    mel_energy = np.maximum(power @ filters.T, np.finfo(float).eps)
    coeffs = dct(np.log(mel_energy), type=2, axis=1, norm="ortho")[:, :n_mfcc]
    return coeffs.mean(axis=0)


def extract_audio_features(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    pre_s = args.pre_ms / 1000.0
    post_s = args.post_ms / 1000.0
    librosa_mod = None
    backend = "scipy"
    if args.feature_engine in {"auto", "librosa"}:
        librosa_mod = load_librosa(required=args.feature_engine == "librosa")
        if librosa_mod is not None:
            backend = "librosa"
        elif args.feature_engine == "auto":
            print(f"librosa_unavailable_falling_back_to_scipy={_LIBROSA_ERROR}")

    for file_path, group in df.groupby("file", sort=True):
        with sf.SoundFile(file_path) as handle:
            sr = handle.samplerate
            for row in group.itertuples(index=False):
                start = float(row.start)
                end = float(row.end)
                pre = read_window(handle, start - pre_s, start)
                stop = read_window(handle, start, end)
                post30 = read_window(handle, end, end + 0.030)
                post = read_window(handle, end, end + post_s)
                whole = read_window(handle, start - pre_s, end + post_s)
                features = {
                    "sample_rate_hz": sr,
                    "pre_window_ms": args.pre_ms,
                    "post_window_ms": args.post_ms,
                    "audio_feature_backend": backend,
                }
                features.update(audio_summary_features(pre, sr, "pre30", librosa_mod=librosa_mod))
                features.update(audio_summary_features(stop, sr, "stop", librosa_mod=librosa_mod))
                features.update(audio_summary_features(post30, sr, "post30", librosa_mod=librosa_mod))
                features.update(audio_summary_features(post, sr, "post80", librosa_mod=librosa_mod))
                features.update(formant_features(post, sr, librosa_mod=librosa_mod, prefix="post80"))
                features["post30_minus_stop_rms_db"] = features["post30_rms_db"] - features["stop_rms_db"]
                features["post30_minus_pre30_rms_db"] = features["post30_rms_db"] - features["pre30_rms_db"]
                features["post30_minus_stop_high_energy_prop"] = (
                    features["post30_high_energy_prop"] - features["stop_high_energy_prop"]
                )
                features["post80_minus_stop_low_high_log_ratio"] = (
                    features["post80_low_high_log_ratio"] - features["stop_low_high_log_ratio"]
                )
                features["post80_minus_stop_voicing_strength"] = (
                    features["post80_voicing_strength"] - features["stop_voicing_strength"]
                )
                mfcc = mfcc_means(whole, sr, args.mfcc_count, librosa_mod=librosa_mod)
                for i, value in enumerate(mfcc, start=1):
                    features[f"mfcc{i}_mean"] = float(value)
                rows.append(features)

    audio_df = pd.DataFrame(rows)
    if len(audio_df) != len(df):
        raise RuntimeError("Audio feature row count does not match token row count.")
    return pd.concat([df.reset_index(drop=True), audio_df.reset_index(drop=True)], axis=1)


def add_zscore_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    metadata = {
        "start",
        "end",
        "sample_rate_hz",
        "pre_window_ms",
        "post_window_ms",
    }
    numeric_cols = [
        col
        for col in df.columns
        if col not in metadata and pd.api.types.is_numeric_dtype(df[col])
    ]
    z_cols = []
    for col in numeric_cols:
        values = pd.to_numeric(df[col], errors="coerce")
        mean = values.mean(skipna=True)
        sd = values.std(skipna=True)
        z_col = f"z_{col}"
        if not np.isfinite(sd) or sd == 0:
            df[z_col] = np.nan
        else:
            df[z_col] = (values - mean) / sd
        z_cols.append(z_col)
    return df, z_cols


def write_feature_sets(path: Path, mfcc_count: int) -> None:
    duration = [
        "z_stop_duration_ms",
        "z_prev_phone_duration_ms",
        "z_next_phone_duration_ms",
        "z_following_vowel_lag_ms",
        "z_following_vowel_duration_ms",
    ]
    spectral_stop = [
        "z_pre30_rms_db",
        "z_pre30_zcr",
        "z_pre30_low_energy_prop",
        "z_pre30_high_energy_prop",
        "z_pre30_low_high_log_ratio",
        "z_pre30_voicing_strength",
        "z_stop_rms_db",
        "z_stop_zcr",
        "z_stop_spectral_centroid_hz",
        "z_stop_spectral_bandwidth_hz",
        "z_stop_spectral_skewness",
        "z_stop_spectral_kurtosis",
        "z_stop_low_energy_prop",
        "z_stop_mid_energy_prop",
        "z_stop_high_energy_prop",
        "z_stop_low_high_log_ratio",
        "z_stop_spectral_flatness",
        "z_stop_voicing_strength",
    ]
    release = [
        "z_post30_rms_db",
        "z_post30_zcr",
        "z_post30_spectral_centroid_hz",
        "z_post30_high_energy_prop",
        "z_post30_low_high_log_ratio",
        "z_post30_voicing_strength",
        "z_post80_rms_db",
        "z_post80_zcr",
        "z_post80_spectral_centroid_hz",
        "z_post80_spectral_bandwidth_hz",
        "z_post80_spectral_skewness",
        "z_post80_spectral_kurtosis",
        "z_post80_low_energy_prop",
        "z_post80_mid_energy_prop",
        "z_post80_high_energy_prop",
        "z_post80_low_high_log_ratio",
        "z_post80_spectral_flatness",
        "z_post80_voicing_strength",
        "z_post30_minus_stop_rms_db",
        "z_post30_minus_pre30_rms_db",
        "z_post30_minus_stop_high_energy_prop",
        "z_post80_minus_stop_low_high_log_ratio",
        "z_post80_minus_stop_voicing_strength",
    ]
    f0_formants = [
        "z_post30_f0_hz",
        "z_post80_f0_hz",
        "z_post80_f1_hz",
        "z_post80_f2_hz",
        "z_post80_f3_hz",
    ]
    mfcc = [f"z_mfcc{i}_mean" for i in range(1, mfcc_count + 1)]
    rows = [
        ("duration_context", duration, "MFA interval durations and following-vowel timing."),
        ("spectral_stop", spectral_stop, "Pre-stop and stop-window spectral, energy, and periodicity cues."),
        ("release_window", release, "Short post-release, vowel-onset, and relative-change acoustic cues."),
        ("f0_formants", f0_formants, "Autocorrelation F0 and LPC formants near vowel onset."),
        ("mfcc_context", mfcc, "MFCC means over pre-stop through post-release context."),
        ("full_bundle", duration + spectral_stop + release + f0_formants + mfcc, "All standardized pilot cues."),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature_set", "features", "note"])
        for name, features, note in rows:
            writer.writerow([name, "; ".join(features), note])


def write_summary(path: Path, candidates: pd.DataFrame, sampled: pd.DataFrame) -> None:
    rows = []
    for label, frame in (("candidates", candidates), ("sampled", sampled)):
        counts = (
            frame.groupby(["control_group", "category", "phone"], dropna=False)
            .size()
            .reset_index(name="n")
        )
        counts.insert(0, "stage", label)
        rows.append(counts)
    summary = pd.concat(rows, ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path, index=False)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tokens_out = args.tokens_out or args.out_dir / "sbcsae_stop_voicing_tokens.csv"
    summary_out = args.summary_out or args.out_dir / "sbcsae_stop_voicing_summary.csv"
    feature_sets_out = args.feature_sets_out or args.out_dir / "sbcsae_stop_voicing_feature_sets.csv"

    candidates = collect_candidates(args)
    sampled = sample_candidates(candidates, args.max_per_place_category, args.seed)
    featured = extract_audio_features(sampled, args)
    featured, _ = add_zscore_columns(featured)

    featured.to_csv(tokens_out, index=False)
    write_summary(summary_out, candidates, sampled)
    write_feature_sets(feature_sets_out, args.mfcc_count)

    print(f"candidate_tokens={len(candidates)}")
    print(f"sampled_tokens={len(sampled)}")
    print(f"tokens_out={tokens_out}")
    print(f"summary_out={summary_out}")
    print(f"feature_sets_out={feature_sets_out}")


if __name__ == "__main__":
    main()
