#!/usr/bin/env python3

"""Summarise flexible-tail tip-displacement histories."""

import argparse
from pathlib import Path
import sys

import numpy as np


def read_history(path):
    samples = {}
    if path.is_dir():
        paths = sorted(path.glob("*/solidPointDisplacement_pointDisp.dat"))
    else:
        paths = [path]
    if not paths:
        raise ValueError(f"{path}: no tip-displacement histories found")

    for history in paths:
        data = np.loadtxt(history, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 3:
            raise ValueError(f"{history}: expected at least three columns")

        # A restarted run may repeat its first time. Retain the newest sample.
        for row in data:
            samples[row[0]] = row[2]

    times = np.asarray(sorted(samples))
    displacement = np.asarray([samples[time] for time in times])
    return times, displacement


def upward_crossing_frequency(times, values):
    trend = np.polyval(np.polyfit(times, values, 1), times)
    centred = values - trend
    indices = np.flatnonzero((centred[:-1] <= 0) & (centred[1:] > 0))
    if len(indices) < 2:
        return np.nan, len(indices)

    crossings = []
    for index in indices:
        fraction = -centred[index]/(centred[index + 1] - centred[index])
        crossings.append(
            times[index] + fraction*(times[index + 1] - times[index])
        )

    periods = np.diff(crossings)
    return 1.0/np.mean(periods), len(periods)


def fft_frequency(times, values):
    delta_t = np.median(np.diff(times))
    uniform_times = np.arange(times[0], times[-1] + 0.5*delta_t, delta_t)
    uniform_values = np.interp(uniform_times, times, values)
    trend = np.polyval(
        np.polyfit(uniform_times, uniform_values, 1),
        uniform_times,
    )
    signal = (uniform_values - trend)*np.hanning(len(uniform_values))
    amplitudes = np.abs(np.fft.rfft(signal))
    frequencies = np.fft.rfftfreq(len(signal), delta_t)
    if len(amplitudes) < 3:
        return np.nan

    peak = np.argmax(amplitudes[1:]) + 1
    if peak == len(amplitudes) - 1:
        return frequencies[peak]

    # Quadratic interpolation reduces the coarse FFT-bin bias in short windows.
    neighbours = np.maximum(amplitudes[peak - 1:peak + 2], np.finfo(float).tiny)
    logarithms = np.log(neighbours)
    denominator = logarithms[0] - 2*logarithms[1] + logarithms[2]
    if denominator == 0:
        offset = 0.0
    else:
        offset = 0.5*(logarithms[0] - logarithms[2])/denominator
    return frequencies[peak] + offset*(frequencies[1] - frequencies[0])


def analyse(label, path, start, end):
    times, displacement = read_history(path)
    selected = (times >= start) & (times <= end)
    if np.count_nonzero(selected) < 4:
        duration = end - start
        fallback_start = max(times[0], times[-1] - duration)
        selected = times >= fallback_start
        print(
            f"{path}: requested window [{start:g}, {end:g}] s is unavailable; "
            f"using the final available {duration:g} s window",
            file=sys.stderr,
        )
    times = times[selected]
    displacement = displacement[selected]
    if len(times) < 4:
        raise ValueError(f"{path}: fewer than four samples in [{start}, {end}] s")

    crossing_frequency, periods = upward_crossing_frequency(times, displacement)
    return {
        "label": label,
        "start": times[0],
        "end": times[-1],
        "half_range": 0.5*(np.max(displacement) - np.min(displacement)),
        "peak": np.max(np.abs(displacement)),
        "crossing_frequency": crossing_frequency,
        "periods": periods,
        "fft_frequency": fft_frequency(times, displacement),
    }


def markdown(results):
    lines = [
        "| Case | Window (s) | Half peak-to-peak (m) | Peak magnitude (m) | "
        "Zero-crossing frequency (Hz) | FFT frequency (Hz) |",
        "| :-- | :--: | --: | --: | --: | --: |",
    ]
    for result in results:
        crossing = result["crossing_frequency"]
        if np.isnan(crossing):
            crossing_text = "-"
        else:
            crossing_text = f"{crossing:.3f} ({result['periods']} periods)"
        fft = result["fft_frequency"]
        fft_text = "-" if np.isnan(fft) else f"{fft:.3f}"
        lines.append(
            f"| {result['label']} | {result['start']:g}--{result['end']:g} | "
            f"{result['half_range']:.6g} | {result['peak']:.6g} | "
            f"{crossing_text} | {fft_text} |"
        )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        nargs=4,
        metavar=("LABEL", "FILE", "START", "END"),
        required=True,
        help="case label, displacement file, and inclusive analysis window",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="also write the Markdown table here",
    )
    args = parser.parse_args()

    results = []
    for label, filename, start, end in args.case:
        results.append(analyse(label, Path(filename), float(start), float(end)))

    table = markdown(results)
    print(table, end="")
    if args.output:
        args.output.write_text(table)


if __name__ == "__main__":
    main()
