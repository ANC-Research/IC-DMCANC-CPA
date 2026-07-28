"""MATLAB-compatible signal construction and display-only smoothing."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import firwin, lfilter

from .config import DTYPE, NUM_NODES, require_integer_samples


@dataclass(frozen=True)
class ReferenceSignals:
    reference_clean: np.ndarray
    reference_input: np.ndarray
    awgn_noise: np.ndarray
    fir_coefficients: np.ndarray | None


def design_fir1_bandpass(
    order: int,
    low_hz: float,
    high_hz: float,
    fs: int,
) -> np.ndarray:
    """Match ``fir1(order, [2*low/Fs, 2*high/Fs])`` for this case.

    MATLAB ``fir1`` interprets its first argument as filter order, so the
    returned filter contains ``order + 1`` Hamming-window taps. SciPy's
    ``firwin`` uses the tap count directly and applies the same passband
    scaling convention for this band-pass design.
    """

    if order < 1:
        raise ValueError("FIR order must be at least one.")
    if not 0 < low_hz < high_hz < fs / 2:
        raise ValueError("Band edges must satisfy 0 < low < high < fs/2.")
    taps = firwin(
        order + 1,
        [low_hz, high_hz],
        pass_zero=False,
        window="hamming",
        scale=True,
        fs=fs,
    )
    return np.ascontiguousarray(taps, dtype=DTYPE)


def add_measured_awgn(
    signal: np.ndarray,
    snr_db: float,
    *,
    seed: int,
    noise: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Implement ``awgn(x, snr_db, 'measured')`` without de-meaning ``x``."""

    x = np.ascontiguousarray(np.asarray(signal, dtype=DTYPE).reshape(-1))
    if not np.isfinite(x).all():
        raise ValueError("Signal contains NaN or Inf.")
    signal_power = float(np.mean(x * x))
    if signal_power <= 0:
        raise ValueError("Measured-SNR AWGN requires nonzero signal power.")
    noise_power = signal_power / (10.0 ** (float(snr_db) / 10.0))
    if noise is None:
        rng = np.random.default_rng(seed)
        standard_noise = rng.standard_normal(x.size)
    else:
        standard_noise = np.asarray(noise, dtype=DTYPE).reshape(-1)
        if standard_noise.shape != x.shape:
            raise ValueError(
                f"Provided AWGN vector has shape {standard_noise.shape}; "
                f"expected {x.shape}."
            )
    standard_power = float(np.mean(standard_noise * standard_noise))
    if standard_power <= 0:
        raise ValueError("Provided AWGN vector has zero power.")
    # Scaling by its measured power makes the finite-record SNR exact.
    awgn_noise = standard_noise * np.sqrt(noise_power / standard_power)
    return (
        np.ascontiguousarray(x + awgn_noise, dtype=DTYPE),
        np.ascontiguousarray(awgn_noise, dtype=DTYPE),
    )


def generate_synthetic_reference(
    *,
    fs: int,
    duration_seconds: float,
    low_hz: float,
    high_hz: float,
    fir_order: int,
    random_seed: int,
    awgn_seed: int,
    snr_db: float,
) -> ReferenceSignals:
    """Generate the case1/3/4/5 band-limited reference deterministically."""

    num_samples = (
        require_integer_samples(
            fs, duration_seconds, name="duration_seconds"
        )
        + 1
    )
    rng = np.random.default_rng(random_seed)
    white = rng.standard_normal(num_samples)
    taps = design_fir1_bandpass(fir_order, low_hz, high_hz, fs)
    clean = np.ascontiguousarray(lfilter(taps, [1.0], white), dtype=DTYPE)
    noisy, awgn_noise = add_measured_awgn(
        clean, snr_db, seed=awgn_seed
    )
    return ReferenceSignals(
        reference_clean=clean,
        reference_input=noisy,
        awgn_noise=awgn_noise,
        fir_coefficients=taps,
    )


def prepare_recorded_reference(
    signal: np.ndarray,
    *,
    fs: int,
    duration_seconds: float,
    awgn_seed: int,
    snr_db: float,
) -> ReferenceSignals:
    """Select ``yr(1:T*Fs)`` and add measured-SNR AWGN as in case 2."""

    sample_count = require_integer_samples(
        fs, duration_seconds, name="duration_seconds"
    )
    source = np.asarray(signal, dtype=DTYPE).reshape(-1)
    if source.size < sample_count:
        raise ValueError(
            f"Recording has {source.size} samples; {sample_count} required."
        )
    clean = np.ascontiguousarray(source[:sample_count])
    noisy, awgn_noise = add_measured_awgn(
        clean, snr_db, seed=awgn_seed
    )
    return ReferenceSignals(
        reference_clean=clean,
        reference_input=noisy,
        awgn_noise=awgn_noise,
        fir_coefficients=None,
    )


def filter_primary_paths(
    primary_path: np.ndarray,
    reference_clean: np.ndarray,
) -> np.ndarray:
    """Create ``Dis(m,:) = filter(PrimaryPath(m,:), 1, Ref)``."""

    primary = np.asarray(primary_path, dtype=DTYPE)
    reference = np.asarray(reference_clean, dtype=DTYPE).reshape(-1)
    if primary.ndim != 2 or primary.shape[0] != NUM_NODES:
        raise ValueError(
            f"primary_path must have shape (6, p_len); got {primary.shape}."
        )
    disturbance = np.empty((NUM_NODES, reference.size), dtype=DTYPE)
    for node in range(NUM_NODES):
        disturbance[node] = lfilter(primary[node], [1.0], reference)
    if not np.isfinite(disturbance).all():
        raise FloatingPointError("Generated disturbance contains NaN or Inf.")
    return disturbance


def _smooth_vector(values: np.ndarray, span: int) -> np.ndarray:
    data = np.asarray(values, dtype=DTYPE).reshape(-1)
    if span <= 0:
        raise ValueError("Smoothing span must be positive.")
    if data.size == 0 or span == 1:
        return data.copy()
    # MATLAB smooth reduces an even moving-average span by one.
    effective = min(span, data.size)
    if effective % 2 == 0:
        effective -= 1
    if effective <= 1:
        return data.copy()
    half = effective // 2
    cumsum = np.empty(data.size + 1, dtype=DTYPE)
    cumsum[0] = 0.0
    np.cumsum(data, out=cumsum[1:])
    indices = np.arange(data.size)
    starts = np.maximum(0, indices - half)
    stops = np.minimum(data.size, indices + half + 1)
    return (cumsum[stops] - cumsum[starts]) / (stops - starts)


def matlab_smooth(values: np.ndarray, span: int) -> np.ndarray:
    """MATLAB-style centered moving average along the final axis.

    This helper is display/metric-only. It is never called by an adaptive
    controller update.
    """

    array = np.asarray(values, dtype=DTYPE)
    if array.ndim == 1:
        return _smooth_vector(array, span)
    if array.ndim != 2:
        raise ValueError("matlab_smooth supports one- or two-dimensional data.")
    output = np.empty_like(array)
    for row in range(array.shape[0]):
        output[row] = _smooth_vector(array[row], span)
    return output
