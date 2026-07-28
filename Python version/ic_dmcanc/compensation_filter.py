"""Compensation-filter system identification and MWD tail filtering."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.signal import lfilter

from .config import DTYPE, NUM_NODES

try:  # Optional acceleration; the Python function remains the authority.
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - environment dependent
    njit = None
    NUMBA_AVAILABLE = False


@dataclass(frozen=True)
class CompensationIdentificationResult:
    """Results for all 30 off-diagonal secondary-path pairs."""

    compensation_filters: np.ndarray
    identified_coefficients: np.ndarray
    final_mse: np.ndarray
    convergence_sample_indices: np.ndarray
    convergence_error_squared: np.ndarray
    white_noise: np.ndarray
    full_error: np.ndarray | None
    engine: str
    seed: int | None
    output_transform: str = "negative_then_flip"
    shared_white_noise: bool = True
    skipped_diagonal: bool = True


def _identify_pair_core(
    white_noise: np.ndarray,
    desired: np.ndarray,
    diagonal_path: np.ndarray,
    c_len: int,
    step_size: float,
    convergence_stride: int,
    save_full_error: bool,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray, np.ndarray]:
    """Sample-by-sample Filtered-x LMS system identification.

    The adaptive FIR output is passed through the actual diagonal secondary
    path before computing ``error = desired - secondary_output``. The same
    diagonal path filters the reference for the coefficient update:

    ``identified[n+1] = identified[n] + mu * filtered_x[n] * error[n]``.

    This function deliberately returns the raw identified coefficient order.
    The MATLAB repository stores ``-flip(Coefficients)`` afterwards.
    """

    num_samples = white_noise.size
    s_len = diagonal_path.size
    coefficients = np.zeros(c_len, dtype=np.float64)
    input_buffer = np.zeros(c_len, dtype=np.float64)
    controller_output_buffer = np.zeros(s_len, dtype=np.float64)
    path_input_buffer = np.zeros(s_len, dtype=np.float64)
    filtered_input_buffer = np.zeros(c_len, dtype=np.float64)
    full_error = np.empty(
        num_samples if save_full_error else 0, dtype=np.float64
    )

    point_count = (num_samples - 1) // convergence_stride + 1
    if (num_samples - 1) % convergence_stride != 0:
        point_count += 1
    convergence_indices = np.empty(point_count, dtype=np.int64)
    convergence_error_squared = np.empty(point_count, dtype=np.float64)
    convergence_position = 0

    final_window = max(1, num_samples // 10)
    final_squared_sum = 0.0

    for sample in range(num_samples):
        for tap in range(c_len - 1, 0, -1):
            input_buffer[tap] = input_buffer[tap - 1]
        input_buffer[0] = white_noise[sample]

        adaptive_output = 0.0
        for tap in range(c_len):
            adaptive_output += coefficients[tap] * input_buffer[tap]

        for tap in range(s_len - 1, 0, -1):
            controller_output_buffer[tap] = controller_output_buffer[tap - 1]
        controller_output_buffer[0] = adaptive_output

        secondary_output = 0.0
        for tap in range(s_len):
            secondary_output += (
                diagonal_path[tap] * controller_output_buffer[tap]
            )
        error = desired[sample] - secondary_output

        for tap in range(s_len - 1, 0, -1):
            path_input_buffer[tap] = path_input_buffer[tap - 1]
        path_input_buffer[0] = white_noise[sample]

        filtered_input_sample = 0.0
        for tap in range(s_len):
            filtered_input_sample += (
                diagonal_path[tap] * path_input_buffer[tap]
            )
        for tap in range(c_len - 1, 0, -1):
            filtered_input_buffer[tap] = filtered_input_buffer[tap - 1]
        filtered_input_buffer[0] = filtered_input_sample

        for tap in range(c_len):
            coefficients[tap] += (
                step_size * filtered_input_buffer[tap] * error
            )

        if save_full_error:
            full_error[sample] = error
        if sample >= num_samples - final_window:
            final_squared_sum += error * error
        if sample % convergence_stride == 0:
            convergence_indices[convergence_position] = sample + 1
            convergence_error_squared[convergence_position] = error * error
            convergence_position += 1

    if (num_samples - 1) % convergence_stride != 0:
        convergence_indices[convergence_position] = num_samples
        convergence_error_squared[convergence_position] = error * error

    return (
        coefficients,
        final_squared_sum / final_window,
        convergence_indices,
        convergence_error_squared,
        full_error,
    )


if NUMBA_AVAILABLE:  # pragma: no branch - selected at import time
    _identify_pair_numba = njit(cache=True)(_identify_pair_core)
else:  # pragma: no cover - trivial alias
    _identify_pair_numba = None


def identify_single_pair(
    white_noise: np.ndarray,
    desired: np.ndarray,
    diagonal_path: np.ndarray,
    *,
    c_len: int,
    step_size: float,
    convergence_stride: int = 1_000,
    save_full_error: bool = False,
    use_numba: bool = True,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray, np.ndarray | None, str]:
    """Identify one pair and expose both raw and diagnostic state."""

    white = np.ascontiguousarray(
        np.asarray(white_noise, dtype=DTYPE).reshape(-1)
    )
    target = np.ascontiguousarray(
        np.asarray(desired, dtype=DTYPE).reshape(-1)
    )
    diagonal = np.ascontiguousarray(
        np.asarray(diagonal_path, dtype=DTYPE).reshape(-1)
    )
    if white.size == 0 or target.shape != white.shape:
        raise ValueError("white_noise and desired must be equal nonempty vectors.")
    if diagonal.size == 0:
        raise ValueError("diagonal_path must not be empty.")
    if c_len <= 0 or convergence_stride <= 0:
        raise ValueError("c_len and convergence_stride must be positive.")
    if not np.isfinite(step_size) or step_size < 0:
        raise ValueError("step_size must be finite and nonnegative.")
    if not (
        np.isfinite(white).all()
        and np.isfinite(target).all()
        and np.isfinite(diagonal).all()
    ):
        raise ValueError("System-identification inputs contain NaN or Inf.")

    if use_numba and NUMBA_AVAILABLE:
        runner = _identify_pair_numba
        engine = "numba"
    else:
        runner = _identify_pair_core
        engine = "numpy_reference"
    raw, final_mse, indices, convergence, full = runner(
        white,
        target,
        diagonal,
        int(c_len),
        float(step_size),
        int(convergence_stride),
        bool(save_full_error),
    )
    return (
        np.asarray(raw, dtype=DTYPE),
        float(final_mse),
        np.asarray(indices, dtype=np.int64),
        np.asarray(convergence, dtype=DTYPE),
        np.asarray(full, dtype=DTYPE) if save_full_error else None,
        engine,
    )


def identify_compensation_filters(
    secondary_path: np.ndarray,
    *,
    c_len: int,
    step_size: float,
    num_samples: int = 200_000,
    seed: int = 0,
    white_noise: np.ndarray | None = None,
    convergence_stride: int = 1_000,
    save_full_error: bool = False,
    use_numba: bool = True,
) -> CompensationIdentificationResult:
    """Reproduce both MATLAB ``CompensateSP`` methods.

    A single white-noise vector is shared by all off-diagonal pairs. Each
    pair starts with independent zero adaptive-filter state. Diagonal entries
    are skipped. The stored compensation filter is exactly:

    ``compensation_filter = -identified_coefficients[::-1]``.
    """

    secondary = np.ascontiguousarray(
        np.asarray(secondary_path, dtype=DTYPE)
    )
    if (
        secondary.ndim != 3
        or secondary.shape[0] != NUM_NODES
        or secondary.shape[1] != NUM_NODES
    ):
        raise ValueError(
            "secondary_path must have shape (6, 6, s_len); "
            f"got {secondary.shape}."
        )
    if not np.isfinite(secondary).all():
        raise ValueError("secondary_path contains NaN or Inf.")
    if num_samples <= 0:
        raise ValueError("num_samples must be positive.")
    if white_noise is None:
        rng = np.random.default_rng(seed)
        white = np.ascontiguousarray(
            rng.standard_normal(num_samples), dtype=DTYPE
        )
        stored_seed: int | None = int(seed)
    else:
        white = np.ascontiguousarray(
            np.asarray(white_noise, dtype=DTYPE).reshape(-1)
        )
        if white.size != num_samples:
            raise ValueError(
                f"white_noise has {white.size} samples; expected {num_samples}."
            )
        stored_seed = None
    if not np.isfinite(white).all():
        raise ValueError("white_noise contains NaN or Inf.")

    if use_numba and not NUMBA_AVAILABLE:
        warnings.warn(
            "Numba is unavailable; compensation identification is using the "
            "slower NumPy reference kernel.",
            RuntimeWarning,
            stacklevel=2,
        )

    coefficients = np.zeros((NUM_NODES, NUM_NODES, c_len), dtype=DTYPE)
    compensation = np.zeros_like(coefficients)
    # MATLAB preallocates the skipped diagonal error entries as zeros.
    final_mse = np.zeros((NUM_NODES, NUM_NODES), dtype=DTYPE)
    convergence_indices: np.ndarray | None = None
    convergence = None
    full_error = (
        np.zeros((NUM_NODES, NUM_NODES, num_samples), dtype=DTYPE)
        if save_full_error
        else None
    )
    engine = "numpy_reference"

    for error_node in range(NUM_NODES):
        diagonal = secondary[error_node, error_node]
        for source_node in range(NUM_NODES):
            if error_node == source_node:
                continue
            desired = np.ascontiguousarray(
                lfilter(
                    secondary[error_node, source_node],
                    [1.0],
                    white,
                ),
                dtype=DTYPE,
            )
            (
                raw,
                pair_mse,
                pair_indices,
                pair_convergence,
                pair_full_error,
                pair_engine,
            ) = identify_single_pair(
                white,
                desired,
                diagonal,
                c_len=c_len,
                step_size=step_size,
                convergence_stride=convergence_stride,
                save_full_error=save_full_error,
                use_numba=use_numba,
            )
            if convergence_indices is None:
                convergence_indices = pair_indices
                convergence = np.zeros(
                    (
                        NUM_NODES,
                        NUM_NODES,
                        pair_convergence.size,
                    ),
                    dtype=DTYPE,
                )
            elif not np.array_equal(convergence_indices, pair_indices):
                raise AssertionError("Pair convergence indices unexpectedly differ.")
            coefficients[error_node, source_node] = raw
            compensation[error_node, source_node] = -raw[::-1]
            final_mse[error_node, source_node] = pair_mse
            convergence[error_node, source_node] = pair_convergence
            if save_full_error and pair_full_error is not None:
                full_error[error_node, source_node] = pair_full_error
            engine = pair_engine

    if convergence_indices is None or convergence is None:
        raise AssertionError("No off-diagonal path pairs were identified.")
    if not (
        np.isfinite(compensation).all()
        and np.isfinite(coefficients).all()
        and np.isfinite(final_mse).all()
    ):
        raise FloatingPointError("Compensation identification produced NaN/Inf.")

    return CompensationIdentificationResult(
        compensation_filters=compensation,
        identified_coefficients=coefficients,
        final_mse=final_mse,
        convergence_sample_indices=convergence_indices,
        convergence_error_squared=convergence,
        white_noise=white,
        full_error=full_error,
        engine=engine,
        seed=stored_seed,
    )


def mwd_tail_lfilter(
    compensation_filter: np.ndarray,
    weight_difference: np.ndarray,
    w_len: int,
) -> np.ndarray:
    """Literal MATLAB mapping: ``filter(C,1,Nabla); a(end-wLen+1:end)``."""

    coefficients = np.asarray(compensation_filter, dtype=DTYPE).reshape(-1)
    difference = np.asarray(weight_difference, dtype=DTYPE).reshape(-1)
    expected = w_len + coefficients.size - 1
    if difference.size != expected:
        raise ValueError(
            f"len(weight_difference) must be w_len + c_len - 1 = "
            f"{expected}; got {difference.size}."
        )
    filtered = lfilter(coefficients, [1.0], difference)
    return np.ascontiguousarray(filtered[-w_len:], dtype=DTYPE)


def mwd_tail_direct(
    compensation_filter: np.ndarray,
    weight_difference: np.ndarray,
    w_len: int,
) -> np.ndarray:
    """Allocation-light expression exactly equal to ``mwd_tail_lfilter``.

    This direct form is used inside the optional Numba kernels. It is *not*
    ``convolve(..., mode='same')`` and it selects the final ``w_len`` samples.
    """

    coefficients = np.asarray(compensation_filter, dtype=DTYPE).reshape(-1)
    difference = np.asarray(weight_difference, dtype=DTYPE).reshape(-1)
    c_len = coefficients.size
    if difference.size != w_len + c_len - 1:
        raise ValueError(
            "weight_difference length must equal w_len + c_len - 1."
        )
    tail = np.empty(w_len, dtype=DTYPE)
    for output_tap in range(w_len):
        signal_index = c_len - 1 + output_tap
        total = 0.0
        for filter_tap in range(c_len):
            total += (
                coefficients[filter_tap]
                * difference[signal_index - filter_tap]
            )
        tail[output_tap] = total
    return tail
