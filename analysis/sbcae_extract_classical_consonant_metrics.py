import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import soundfile as sf


STOP_MAP = {
    "b": "b",
    "bʲ": "b",
    "p": "p",
    "pʰ": "p",
    "pʲ": "p",
    "d": "d",
    "d̪": "d",
    "dʲ": "d",
    "t": "t",
    "tʰ": "t",
    "tʲ": "t",
    "t̪": "t",
    "ɡ": "g",
    "g": "g",
    "k": "k",
    "kʰ": "k",
    "kʷ": "k",
}

FRICATIVE_MAP = {
    "s": "s",
    "ʃ": "sh",
    "z": "z",
    "ʒ": "zh",
    "f": "f",
    "fʲ": "f",
    "θ": "th",
}

TARGET_MAP = {}
TARGET_MAP.update(STOP_MAP)
TARGET_MAP.update(FRICATIVE_MAP)


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract classical consonant metrics (VOT + fricative spectral measures)."
    )
    p.add_argument(
        "--align-dir",
        default="analysis/sbcae_mfa_segments_aligned",
        help="Directory containing aligned TextGrid files.",
    )
    p.add_argument(
        "--wav-dir",
        default="analysis/sbcae_mfa_segments",
        help="Directory containing segmented WAV files.",
    )
    p.add_argument(
        "--out-csv",
        default="analysis/sbcae_classical_consonant_metrics.csv",
        help="Output CSV path.",
    )
    p.add_argument(
        "--file-limit",
        type=int,
        default=0,
        help="If > 0, only process this many TextGrids.",
    )
    p.add_argument(
        "--min-dur",
        type=float,
        default=0.015,
        help="Minimum segment duration in seconds.",
    )
    return p.parse_args()


def parse_phones(textgrid_path: Path):
    lines = textgrid_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    phones = []
    in_phones = False
    for i, line in enumerate(lines):
        if 'name = "phones"' in line:
            in_phones = True
        if in_phones and line.strip().startswith("intervals ["):
            try:
                xmin = float(lines[i + 1].split("=")[1].strip())
                xmax = float(lines[i + 2].split("=")[1].strip())
                text = lines[i + 3].split("=")[1].strip().strip('"')
            except Exception:
                continue
            phones.append((xmin, xmax, text))
    return phones


def canonical_label(label):
    if label is None:
        return None
    base = re.sub(r"\d", "", str(label)).strip()
    return TARGET_MAP.get(base), base


def frame_rms(x, frame_len, hop):
    if len(x) < frame_len:
        return np.array([math.sqrt(float(np.mean(x * x)))]) if len(x) else np.array([0.0])
    starts = np.arange(0, len(x) - frame_len + 1, hop)
    out = np.empty(len(starts), dtype=float)
    for i, s in enumerate(starts):
        w = x[s : s + frame_len]
        out[i] = math.sqrt(float(np.mean(w * w)))
    return out


def estimate_vot_ms(y, sr, seg_start, seg_end):
    seg = y[seg_start:seg_end]
    if len(seg) < int(0.008 * sr):
        return np.nan

    frame_len = max(16, int(0.001 * sr))
    hop = max(8, int(0.0005 * sr))

    hp = np.concatenate(([seg[0]], seg[1:] - 0.97 * seg[:-1]))
    e = frame_rms(hp, frame_len, hop)
    if len(e) < 4:
        return np.nan

    loge = np.log(e + 1e-12)
    de = np.diff(loge)
    search_start = max(1, int(0.3 * len(de)))
    rel_idx = search_start + int(np.argmax(de[search_start:]))
    release_sample = seg_start + rel_idx * hop + frame_len // 2

    pre = y[max(seg_start, release_sample - int(0.03 * sr)) : release_sample]
    pre_rms = math.sqrt(float(np.mean(pre * pre))) if len(pre) else 0.0
    min_rms = max(0.002, 2.0 * pre_rms)

    post_end = min(len(y), release_sample + int(0.12 * sr))
    post = y[release_sample:post_end]
    if len(post) < int(0.006 * sr):
        return np.nan

    v_frame = max(24, int(0.005 * sr))
    v_hop = max(12, int(0.001 * sr))
    starts = np.arange(0, max(1, len(post) - v_frame + 1), v_hop)
    voiced_flags = []
    onset_idx = None

    lag_min = max(1, int(sr / 400))
    lag_max = max(lag_min + 1, int(sr / 60))

    for s in starts:
        w = post[s : s + v_frame]
        if len(w) < v_frame:
            break
        rms = math.sqrt(float(np.mean(w * w)))
        if rms < min_rms:
            voiced_flags.append(False)
            continue

        zcr = float(np.mean(w[:-1] * w[1:] < 0))
        w0 = w - np.mean(w)
        ac = np.correlate(w0, w0, mode="full")[len(w0) - 1 :]
        if ac[0] <= 0:
            voiced_flags.append(False)
            continue
        ac = ac / (ac[0] + 1e-12)
        periodicity = float(np.max(ac[lag_min:min(lag_max, len(ac))]))
        voiced = periodicity > 0.30 and zcr < 0.24
        voiced_flags.append(voiced)

    for i in range(0, max(0, len(voiced_flags) - 2)):
        if voiced_flags[i] and voiced_flags[i + 1] and voiced_flags[i + 2]:
            onset_idx = starts[i]
            break

    if onset_idx is None:
        return np.nan

    vot_ms = 1000.0 * (onset_idx / sr)
    if vot_ms < -20 or vot_ms > 150:
        return np.nan
    return vot_ms


def spectral_metrics_fricative(seg, sr):
    if len(seg) < int(0.008 * sr):
        return {
            "cog_hz": np.nan,
            "spec_sd_hz": np.nan,
            "spec_skew": np.nan,
            "spec_kurt": np.nan,
            "peak_hz": np.nan,
            "spec_slope_db_per_khz": np.nan,
            "band_ratio_hi_lo_db": np.nan,
            "intensity_db": np.nan,
        }

    win_len = int(0.04 * sr)
    if len(seg) > win_len:
        s = (len(seg) - win_len) // 2
        x = seg[s : s + win_len]
    else:
        x = seg

    x = np.concatenate(([x[0]], x[1:] - 0.97 * x[:-1]))
    x = x * np.hamming(len(x))

    nfft = 1
    while nfft < len(x):
        nfft *= 2
    nfft = max(1024, nfft)

    spec = np.fft.rfft(x, n=nfft)
    power = np.abs(spec) ** 2
    freqs = np.fft.rfftfreq(nfft, d=1.0 / sr)

    fmax = min(11000.0, float(sr) / 2.0 - 10.0)
    mask = (freqs >= 500.0) & (freqs <= fmax)
    f = freqs[mask]
    p = power[mask]
    if len(f) < 4 or float(np.sum(p)) <= 0:
        return {
            "cog_hz": np.nan,
            "spec_sd_hz": np.nan,
            "spec_skew": np.nan,
            "spec_kurt": np.nan,
            "peak_hz": np.nan,
            "spec_slope_db_per_khz": np.nan,
            "band_ratio_hi_lo_db": np.nan,
            "intensity_db": np.nan,
        }

    w = p / (np.sum(p) + 1e-12)
    cog = float(np.sum(f * w))
    sd = float(np.sqrt(np.sum(((f - cog) ** 2) * w)))
    skew = float(np.sum(((f - cog) ** 3) * w) / (sd**3 + 1e-12))
    kurt = float(np.sum(((f - cog) ** 4) * w) / (sd**4 + 1e-12))
    peak_hz = float(f[np.argmax(w)])

    xk = f / 1000.0
    yk = 10.0 * np.log10(p + 1e-12)
    try:
        slope = float(np.polyfit(xk, yk, 1)[0])
    except Exception:
        slope = np.nan

    lo = float(np.sum(p[(f >= 500) & (f < 3000)]))
    hi = float(np.sum(p[(f >= 3000) & (f <= fmax)]))
    band_ratio = float(10.0 * np.log10((hi + 1e-12) / (lo + 1e-12)))

    rms = math.sqrt(float(np.mean(seg * seg)))
    intensity = float(20.0 * np.log10(rms + 1e-12))

    return {
        "cog_hz": cog,
        "spec_sd_hz": sd,
        "spec_skew": skew,
        "spec_kurt": kurt,
        "peak_hz": peak_hz,
        "spec_slope_db_per_khz": slope,
        "band_ratio_hi_lo_db": band_ratio,
        "intensity_db": intensity,
    }


def main():
    args = parse_args()
    align_dir = Path(args.align_dir)
    wav_dir = Path(args.wav_dir)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(align_dir.glob("*.TextGrid"))
    if args.file_limit and args.file_limit > 0:
        files = files[: args.file_limit]

    rows_written = 0
    counts = Counter()
    skipped_audio = 0

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "file",
                "speaker",
                "segment",
                "raw_label",
                "start",
                "end",
                "duration_ms",
                "vot_ms",
                "cog_hz",
                "spec_sd_hz",
                "spec_skew",
                "spec_kurt",
                "peak_hz",
                "spec_slope_db_per_khz",
                "band_ratio_hi_lo_db",
                "intensity_db",
            ]
        )

        for tg_path in files:
            base = tg_path.stem
            speaker = base.split("_")[0]
            wav_path = wav_dir / f"{base}.wav"
            if not wav_path.exists():
                skipped_audio += 1
                continue

            try:
                y, sr = sf.read(str(wav_path))
            except Exception:
                skipped_audio += 1
                continue

            if y.ndim > 1:
                y = y.mean(axis=1)

            phones = parse_phones(tg_path)
            if not phones:
                continue

            for xmin, xmax, label in phones:
                canon, raw = canonical_label(label)
                if canon is None:
                    continue

                dur = xmax - xmin
                if dur < args.min_dur:
                    continue

                start = int(max(0, xmin) * sr)
                end = int(min(xmax, len(y) / sr) * sr)
                if end <= start or end > len(y):
                    continue

                seg = y[start:end]
                if len(seg) < 64:
                    continue

                vot_ms = np.nan
                if canon in {"b", "p", "d", "t", "g", "k"}:
                    vot_ms = estimate_vot_ms(y, sr, start, end)

                fric_metrics = {
                    "cog_hz": np.nan,
                    "spec_sd_hz": np.nan,
                    "spec_skew": np.nan,
                    "spec_kurt": np.nan,
                    "peak_hz": np.nan,
                    "spec_slope_db_per_khz": np.nan,
                    "band_ratio_hi_lo_db": np.nan,
                    "intensity_db": np.nan,
                }
                if canon in {"s", "sh", "z", "zh", "f", "th"}:
                    fric_metrics = spectral_metrics_fricative(seg, sr)

                writer.writerow(
                    [
                        base,
                        speaker,
                        canon,
                        raw,
                        xmin,
                        xmax,
                        1000.0 * dur,
                        vot_ms,
                        fric_metrics["cog_hz"],
                        fric_metrics["spec_sd_hz"],
                        fric_metrics["spec_skew"],
                        fric_metrics["spec_kurt"],
                        fric_metrics["peak_hz"],
                        fric_metrics["spec_slope_db_per_khz"],
                        fric_metrics["band_ratio_hi_lo_db"],
                        fric_metrics["intensity_db"],
                    ]
                )
                counts[canon] += 1
                rows_written += 1

    print(f"Wrote: {out_csv}")
    print(f"Rows: {rows_written}")
    print(f"TextGrids scanned: {len(files)}")
    print(f"Missing/bad audio files skipped: {skipped_audio}")
    print("Counts by canonical segment:")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")


if __name__ == "__main__":
    main()
