"""Centralized 1×6×6 FxLMS mapping of ``McFxLMS_SIMO_166``."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np

from .config import DTYPE, NUM_NODES, SimulationConfig

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - environment dependent
    njit = None
    NUMBA_AVAILABLE = False


@dataclass(frozen=True)
class CentralizedFxLMSResult:
    residual_error: np.ndarray
    control_output: np.ndarray | None
    final_weights: np.ndarray
    diagnostic_sample_indices: np.ndarray
    controller_norm: np.ndarray
    engine: str


def _centralized_core(
    reference: np.ndarray,
    disturbance: np.ndarray,
    secondary_path: np.ndarray,
    w_len: int,
    mu_w: float,
    save_control_output: bool,
    diagnostic_stride: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    num_samples = reference.size
    s_len = secondary_path.shape[2]
    residual_error = np.empty((NUM_NODES, num_samples), dtype=np.float64)
    control_output = np.empty(
        (NUM_NODES, num_samples)
        if save_control_output
        else (0, 0),
        dtype=np.float64,
    )
    weights = np.zeros((NUM_NODES, w_len), dtype=np.float64)
    control_reference_buffer = np.zeros(w_len, dtype=np.float64)
    path_reference_buffer = np.zeros(s_len, dtype=np.float64)
    control_path_buffers = np.zeros(
        (NUM_NODES, s_len), dtype=np.float64
    )
    filtered_reference = np.zeros(
        (NUM_NODES, NUM_NODES, w_len), dtype=np.float64
    )
    current_control = np.zeros(NUM_NODES, dtype=np.float64)

    diagnostic_count = (num_samples - 1) // diagnostic_stride + 1
    if (num_samples - 1) % diagnostic_stride != 0:
        diagnostic_count += 1
    diagnostic_indices = np.empty(diagnostic_count, dtype=np.int64)
    controller_norm = np.empty(
        (NUM_NODES, diagnostic_count), dtype=np.float64
    )
    diagnostic_position = 0

    for sample in range(num_samples):
        for tap in range(w_len - 1, 0, -1):
            control_reference_buffer[tap] = (
                control_reference_buffer[tap - 1]
            )
        control_reference_buffer[0] = reference[sample]

        for source in range(NUM_NODES):
            value = 0.0
            for tap in range(w_len):
                value += weights[source, tap] * control_reference_buffer[tap]
            current_control[source] = value
            if save_control_output:
                control_output[source, sample] = value

        for source in range(NUM_NODES):
            for tap in range(s_len - 1, 0, -1):
                control_path_buffers[source, tap] = (
                    control_path_buffers[source, tap - 1]
                )
            control_path_buffers[source, 0] = current_control[source]

        for error_node in range(NUM_NODES):
            secondary_output = 0.0
            for source in range(NUM_NODES):
                for tap in range(s_len):
                    secondary_output += (
                        secondary_path[error_node, source, tap]
                        * control_path_buffers[source, tap]
                    )
            residual_error[error_node, sample] = (
                disturbance[error_node, sample] - secondary_output
            )

        for tap in range(s_len - 1, 0, -1):
            path_reference_buffer[tap] = path_reference_buffer[tap - 1]
        path_reference_buffer[0] = reference[sample]

        # MATLAB xf(source,error,:) uses SecP(error,source,:).
        for source in range(NUM_NODES):
            for error_node in range(NUM_NODES):
                filtered_sample = 0.0
                for tap in range(s_len):
                    filtered_sample += (
                        secondary_path[error_node, source, tap]
                        * path_reference_buffer[tap]
                    )
                for tap in range(w_len - 1, 0, -1):
                    filtered_reference[source, error_node, tap] = (
                        filtered_reference[source, error_node, tap - 1]
                    )
                filtered_reference[source, error_node, 0] = filtered_sample

        for source in range(NUM_NODES):
            for tap in range(w_len):
                gradient = 0.0
                for error_node in range(NUM_NODES):
                    gradient += (
                        filtered_reference[source, error_node, tap]
                        * residual_error[error_node, sample]
                    )
                weights[source, tap] += mu_w * gradient

        if (
            sample % diagnostic_stride == 0
            or sample == num_samples - 1
        ):
            diagnostic_indices[diagnostic_position] = sample + 1
            for source in range(NUM_NODES):
                squared = 0.0
                for tap in range(w_len):
                    squared += weights[source, tap] * weights[source, tap]
                controller_norm[source, diagnostic_position] = np.sqrt(
                    squared
                )
            diagnostic_position += 1

    return (
        residual_error,
        control_output,
        weights,
        diagnostic_indices,
        controller_norm,
    )


if NUMBA_AVAILABLE:  # pragma: no branch
    _centralized_numba = njit(cache=True)(_centralized_core)
else:  # pragma: no cover
    _centralized_numba = None


def run_centralized_fxlms(
    reference: np.ndarray,
    disturbance: np.ndarray,
    secondary_path: np.ndarray,
    config: SimulationConfig,
    *,
    save_control_output: bool = True,
) -> CentralizedFxLMSResult:
    """Run the case1/case2 1×6×6 centralized baseline."""

    config = config.validate()
    reference_array = np.ascontiguousarray(
        np.asarray(reference, dtype=DTYPE).reshape(-1)
    )
    disturbance_array = np.ascontiguousarray(
        np.asarray(disturbance, dtype=DTYPE)
    )
    secondary = np.ascontiguousarray(
        np.asarray(secondary_path, dtype=DTYPE)
    )
    expected_disturbance = (NUM_NODES, reference_array.size)
    if disturbance_array.shape != expected_disturbance:
        raise ValueError(
            f"disturbance must have shape {expected_disturbance}; "
            f"got {disturbance_array.shape}."
        )
    if secondary.shape != (NUM_NODES, NUM_NODES, config.s_len):
        raise ValueError("secondary_path shape does not match the configuration.")
    if not all(
        np.isfinite(array).all()
        for array in (reference_array, disturbance_array, secondary)
    ):
        raise ValueError("Centralized FxLMS input contains NaN or Inf.")

    if config.use_numba and NUMBA_AVAILABLE:
        runner = _centralized_numba
        engine = "numba"
    else:
        if config.use_numba and not NUMBA_AVAILABLE:
            warnings.warn(
                "Numba is unavailable; centralized FxLMS is using the slower "
                "NumPy reference kernel.",
                RuntimeWarning,
                stacklevel=2,
            )
        runner = _centralized_core
        engine = "numpy_reference"
    residual, control, weights, indices, norms = runner(
        reference_array,
        disturbance_array,
        secondary,
        int(config.w_len),
        float(config.mu_w),
        bool(save_control_output),
        int(config.diagnostic_stride),
    )
    arrays = [residual, weights, norms]
    if save_control_output:
        arrays.append(control)
    if not all(np.isfinite(array).all() for array in arrays):
        raise FloatingPointError("Centralized FxLMS produced NaN or Inf.")
    return CentralizedFxLMSResult(
        residual_error=np.asarray(residual, dtype=DTYPE),
        control_output=(
            np.asarray(control, dtype=DTYPE)
            if save_control_output
            else None
        ),
        final_weights=np.asarray(weights, dtype=DTYPE),
        diagnostic_sample_indices=np.asarray(indices, dtype=np.int64),
        controller_norm=np.asarray(norms, dtype=DTYPE),
        engine=engine,
    )
