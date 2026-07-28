"""Strict six-node FedDMCANC/WCFxLMS and communication schedules."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np

from .compensation_filter import mwd_tail_lfilter
from .config import (
    DTYPE,
    NUM_NODES,
    CommunicationMode,
    SimulationConfig,
)

try:  # Optional acceleration; both paths execute the same core function.
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - environment dependent
    njit = None
    NUMBA_AVAILABLE = False


@dataclass(frozen=True)
class CommunicationSchedule:
    """One-based communication events matching MATLAB loop indexing."""

    mode: CommunicationMode
    sample_indices: np.ndarray
    node_indices: np.ndarray
    interval_samples: np.ndarray
    per_node_counts: np.ndarray
    total_node_uploads: int

    @property
    def total_event_rows(self) -> int:
        return int(self.sample_indices.size)


@dataclass(frozen=True)
class FedMCANCResult:
    residual_error: np.ndarray
    control_output: np.ndarray | None
    secondary_output: np.ndarray | None
    final_local_weights: np.ndarray
    final_center_weights: np.ndarray
    final_weight_difference: np.ndarray
    diagnostic_sample_indices: np.ndarray
    local_weight_norm: np.ndarray
    center_weight_norm: np.ndarray
    communication: CommunicationSchedule
    engine: str


def build_communication_schedule(
    num_samples: int,
    mode: CommunicationMode,
    interval_samples: np.ndarray,
) -> CommunicationSchedule:
    """Build deterministic events before simulation for logging and tests."""

    if num_samples <= 0:
        raise ValueError("num_samples must be positive.")
    mode = CommunicationMode(mode)
    intervals = np.asarray(interval_samples, dtype=np.int64).reshape(-1)
    if np.any(intervals <= 0):
        raise ValueError("All communication intervals must be positive.")

    if mode is CommunicationMode.IDEAL:
        sample_indices = np.arange(1, num_samples + 1, dtype=np.int64)
        node_indices = np.full(sample_indices.size, -1, dtype=np.int64)
        counts = np.full(NUM_NODES, sample_indices.size, dtype=np.int64)
    elif mode in {
        CommunicationMode.FIXED_10_SECONDS,
        CommunicationMode.COMMON_PERIODIC,
    }:
        if intervals.size != 1:
            raise ValueError("Common communication requires one interval.")
        sample_indices = np.arange(
            intervals[0], num_samples + 1, intervals[0], dtype=np.int64
        )
        node_indices = np.full(sample_indices.size, -1, dtype=np.int64)
        counts = np.full(NUM_NODES, sample_indices.size, dtype=np.int64)
    else:
        if intervals.size != NUM_NODES:
            raise ValueError("Individual communication requires six intervals.")
        samples: list[int] = []
        nodes: list[int] = []
        for node in range(NUM_NODES):
            triggered = np.arange(
                intervals[node],
                num_samples + 1,
                intervals[node],
                dtype=np.int64,
            )
            samples.extend(int(value) for value in triggered)
            nodes.extend([node] * triggered.size)
        if samples:
            sample_array = np.asarray(samples, dtype=np.int64)
            node_array = np.asarray(nodes, dtype=np.int64)
            order = np.lexsort((node_array, sample_array))
            sample_indices = sample_array[order]
            node_indices = node_array[order]
        else:
            sample_indices = np.empty(0, dtype=np.int64)
            node_indices = np.empty(0, dtype=np.int64)
        counts = np.bincount(
            node_indices, minlength=NUM_NODES
        ).astype(np.int64)

    return CommunicationSchedule(
        mode=mode,
        sample_indices=sample_indices,
        node_indices=node_indices,
        interval_samples=intervals,
        per_node_counts=counts,
        total_node_uploads=int(np.sum(counts)),
    )


def apply_common_communication(
    local_weights: np.ndarray,
    center_weights: np.ndarray,
    compensation_filters: np.ndarray,
    *,
    w_len: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Literal common/ideal MWD aggregation using SciPy ``lfilter``."""

    local = np.asarray(local_weights, dtype=DTYPE)
    center = np.asarray(center_weights, dtype=DTYPE)
    compensation = np.asarray(compensation_filters, dtype=DTYPE)
    if local.shape != center.shape:
        raise ValueError("local_weights and center_weights shapes differ.")
    if local.shape[0] != NUM_NODES or compensation.shape[:2] != (
        NUM_NODES,
        NUM_NODES,
    ):
        raise ValueError("Common communication requires six nodes.")
    c_len = compensation.shape[2]
    if local.shape[1] != w_len + c_len - 1:
        raise ValueError("Extended weight length must equal w_len+c_len-1.")
    difference = local - center
    updated_center = center.copy()
    for target in range(NUM_NODES):
        delta = difference[target, :w_len].copy()
        for other in range(NUM_NODES):
            if other == target:
                continue
            delta += mwd_tail_lfilter(
                compensation[other, target],
                difference[other],
                w_len,
            )
        updated_center[target, :w_len] += delta
    updated_local = updated_center.copy()
    return updated_local, updated_center, difference


def apply_individual_communication(
    local_weights: np.ndarray,
    center_weights: np.ndarray,
    stale_weight_difference: np.ndarray,
    compensation_filters: np.ndarray,
    *,
    target_node: int,
    w_len: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply one individual event and preserve all other stale differences."""

    if target_node not in range(NUM_NODES):
        raise ValueError("target_node must be in 0..5.")
    local = np.asarray(local_weights, dtype=DTYPE).copy()
    center = np.asarray(center_weights, dtype=DTYPE).copy()
    difference = np.asarray(stale_weight_difference, dtype=DTYPE).copy()
    compensation = np.asarray(compensation_filters, dtype=DTYPE)
    c_len = compensation.shape[2]
    expected = (NUM_NODES, w_len + c_len - 1)
    if local.shape != expected or center.shape != expected:
        raise ValueError(f"Weight states must have shape {expected}.")
    if difference.shape != expected:
        raise ValueError(f"stale_weight_difference must have shape {expected}.")

    difference[target_node] = local[target_node] - center[target_node]
    delta = difference[target_node, :w_len].copy()
    for other in range(NUM_NODES):
        if other == target_node:
            continue
        delta += mwd_tail_lfilter(
            compensation[other, target_node],
            difference[other],
            w_len,
        )
    center[target_node, :w_len] += delta
    local[target_node] = center[target_node]
    return local, center, difference


def _mode_code(mode: CommunicationMode) -> int:
    mode = CommunicationMode(mode)
    if mode is CommunicationMode.IDEAL:
        return 0
    if mode in {
        CommunicationMode.FIXED_10_SECONDS,
        CommunicationMode.COMMON_PERIODIC,
    }:
        return 1
    return 2


def _diagnostic_point_count(num_samples: int, stride: int) -> int:
    count = (num_samples - 1) // stride + 1
    if (num_samples - 1) % stride != 0:
        count += 1
    return count


def _fed_mcanc_core(
    reference: np.ndarray,
    disturbance: np.ndarray,
    secondary_path: np.ndarray,
    compensation_filters: np.ndarray,
    initial_center_weights: np.ndarray,
    w_len: int,
    mu_w: float,
    alpha: float,
    mode_code: int,
    interval_samples: np.ndarray,
    save_control_output: bool,
    save_secondary_output: bool,
    diagnostic_stride: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Single source of truth used as Python and optional Numba kernel."""

    num_samples = reference.size
    s_len = secondary_path.shape[2]
    c_len = compensation_filters.shape[2]
    extended_len = w_len + c_len - 1

    residual_error = np.empty((NUM_NODES, num_samples), dtype=np.float64)
    control_output = np.empty(
        (NUM_NODES, num_samples)
        if save_control_output
        else (0, 0),
        dtype=np.float64,
    )
    secondary_output_store = np.empty(
        (NUM_NODES, num_samples)
        if save_secondary_output
        else (0, 0),
        dtype=np.float64,
    )
    local_weights = initial_center_weights.copy()
    center_weights = initial_center_weights.copy()
    weight_difference = np.zeros(
        (NUM_NODES, extended_len), dtype=np.float64
    )

    reference_control_buffer = np.zeros(w_len, dtype=np.float64)
    control_path_buffers = np.zeros(
        (NUM_NODES, s_len), dtype=np.float64
    )
    reference_path_buffer = np.zeros(s_len, dtype=np.float64)
    filtered_reference_buffers = np.zeros(
        (NUM_NODES, extended_len), dtype=np.float64
    )
    current_control = np.zeros(NUM_NODES, dtype=np.float64)
    current_secondary = np.zeros(NUM_NODES, dtype=np.float64)

    diagnostic_count = (num_samples - 1) // diagnostic_stride + 1
    if (num_samples - 1) % diagnostic_stride != 0:
        diagnostic_count += 1
    diagnostic_indices = np.empty(diagnostic_count, dtype=np.int64)
    local_norm = np.empty(
        (NUM_NODES, diagnostic_count), dtype=np.float64
    )
    center_norm = np.empty(
        (NUM_NODES, diagnostic_count), dtype=np.float64
    )
    diagnostic_position = 0

    for sample in range(num_samples):
        # xc = [Ref(i), xc(1:end-1)]
        for tap in range(w_len - 1, 0, -1):
            reference_control_buffer[tap] = reference_control_buffer[tap - 1]
        reference_control_buffer[0] = reference[sample]

        # MATLAB only uses the first w_len entries for acoustic output.
        for source in range(NUM_NODES):
            value = 0.0
            for tap in range(w_len):
                value += (
                    local_weights[source, tap]
                    * reference_control_buffer[tap]
                )
            current_control[source] = value
            if save_control_output:
                control_output[source, sample] = value

        # ys(k,:) = [yc(k,i), ys(k,1:end-1)]
        for source in range(NUM_NODES):
            for tap in range(s_len - 1, 0, -1):
                control_path_buffers[source, tap] = (
                    control_path_buffers[source, tap - 1]
                )
            control_path_buffers[source, 0] = current_control[source]

        # Full 6x6 physical secondary-path coupling.
        for error_node in range(NUM_NODES):
            value = 0.0
            for source in range(NUM_NODES):
                for tap in range(s_len):
                    value += (
                        secondary_path[error_node, source, tap]
                        * control_path_buffers[source, tap]
                    )
            current_secondary[error_node] = value
            residual_error[error_node, sample] = (
                disturbance[error_node, sample] - value
            )
            if save_secondary_output:
                secondary_output_store[error_node, sample] = value

        # xs = [Ref(i), xs(1:end-1)]
        for tap in range(s_len - 1, 0, -1):
            reference_path_buffer[tap] = reference_path_buffer[tap - 1]
        reference_path_buffer[0] = reference[sample]

        # Six diagonal filtered references and local WCFxLMS updates.
        for node in range(NUM_NODES):
            filtered_sample = 0.0
            for tap in range(s_len):
                filtered_sample += (
                    secondary_path[node, node, tap]
                    * reference_path_buffer[tap]
                )
            for tap in range(extended_len - 1, 0, -1):
                filtered_reference_buffers[node, tap] = (
                    filtered_reference_buffers[node, tap - 1]
                )
            filtered_reference_buffers[node, 0] = filtered_sample
            error = residual_error[node, sample]
            for tap in range(extended_len):
                old_local = local_weights[node, tap]
                local_weights[node, tap] = (
                    old_local
                    + mu_w
                    * filtered_reference_buffers[node, tap]
                    * error
                    + alpha
                    * mu_w
                    * (center_weights[node, tap] - old_local)
                )

        matlab_sample_index = sample + 1
        common_trigger = (
            mode_code == 0
            or (
                mode_code == 1
                and matlab_sample_index % interval_samples[0] == 0
            )
        )

        if common_trigger:
            # Snapshot all six differences before updating any center row.
            for node in range(NUM_NODES):
                for tap in range(extended_len):
                    weight_difference[node, tap] = (
                        local_weights[node, tap] - center_weights[node, tap]
                    )

            for target in range(NUM_NODES):
                for output_tap in range(w_len):
                    delta = weight_difference[target, output_tap]
                    signal_index = c_len - 1 + output_tap
                    for other in range(NUM_NODES):
                        if other == target:
                            continue
                        for filter_tap in range(c_len):
                            delta += (
                                compensation_filters[
                                    other, target, filter_tap
                                ]
                                * weight_difference[
                                    other, signal_index - filter_tap
                                ]
                            )
                    center_weights[target, output_tap] += delta

            # MATLAB: obj.Wc = obj.Wcsubopt
            for node in range(NUM_NODES):
                for tap in range(extended_len):
                    local_weights[node, tap] = center_weights[node, tap]

        elif mode_code == 2:
            # The loop order is semantically observable when triggers coincide.
            for target in range(NUM_NODES):
                if matlab_sample_index % interval_samples[target] != 0:
                    continue

                # Refresh only the triggering node; all others remain stale.
                for tap in range(extended_len):
                    weight_difference[target, tap] = (
                        local_weights[target, tap]
                        - center_weights[target, tap]
                    )

                for output_tap in range(w_len):
                    delta = weight_difference[target, output_tap]
                    signal_index = c_len - 1 + output_tap
                    for other in range(NUM_NODES):
                        if other == target:
                            continue
                        for filter_tap in range(c_len):
                            delta += (
                                compensation_filters[
                                    other, target, filter_tap
                                ]
                                * weight_difference[
                                    other, signal_index - filter_tap
                                ]
                            )
                    center_weights[target, output_tap] += delta

                # Reset only the triggering local controller.
                for tap in range(extended_len):
                    local_weights[target, tap] = center_weights[target, tap]

        if (
            sample % diagnostic_stride == 0
            or sample == num_samples - 1
        ):
            diagnostic_indices[diagnostic_position] = matlab_sample_index
            for node in range(NUM_NODES):
                local_squared = 0.0
                center_squared = 0.0
                for tap in range(extended_len):
                    local_squared += (
                        local_weights[node, tap] * local_weights[node, tap]
                    )
                    center_squared += (
                        center_weights[node, tap] * center_weights[node, tap]
                    )
                local_norm[node, diagnostic_position] = np.sqrt(local_squared)
                center_norm[node, diagnostic_position] = np.sqrt(
                    center_squared
                )
            diagnostic_position += 1

    return (
        residual_error,
        control_output,
        secondary_output_store,
        local_weights,
        center_weights,
        weight_difference,
        diagnostic_indices,
        local_norm,
        center_norm,
    )


if NUMBA_AVAILABLE:  # pragma: no branch - selected at import time
    _fed_mcanc_numba = njit(cache=True)(_fed_mcanc_core)
else:  # pragma: no cover - trivial alias
    _fed_mcanc_numba = None


def run_fed_mcanc(
    reference: np.ndarray,
    disturbance: np.ndarray,
    secondary_path: np.ndarray,
    compensation_filters: np.ndarray,
    config: SimulationConfig,
    *,
    initial_center_weights: np.ndarray | None = None,
    save_control_output: bool = True,
    save_secondary_output: bool = False,
) -> FedMCANCResult:
    """Run one independently initialized FedDMCANC experiment."""

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
    compensation = np.ascontiguousarray(
        np.asarray(compensation_filters, dtype=DTYPE)
    )
    num_samples = reference_array.size
    if num_samples == 0:
        raise ValueError("reference must not be empty.")
    if disturbance_array.shape != (NUM_NODES, num_samples):
        raise ValueError(
            f"disturbance must have shape {(NUM_NODES, num_samples)}; "
            f"got {disturbance_array.shape}."
        )
    if secondary.shape != (NUM_NODES, NUM_NODES, config.s_len):
        raise ValueError(
            "secondary_path shape does not match (6,6,s_len): "
            f"{secondary.shape}."
        )
    if compensation.shape != (
        NUM_NODES,
        NUM_NODES,
        config.c_len,
    ):
        raise ValueError(
            "compensation_filters shape does not match (6,6,c_len): "
            f"{compensation.shape}."
        )
    extended_len = config.extended_weight_len
    if initial_center_weights is None:
        initial_center = np.zeros(
            (NUM_NODES, extended_len), dtype=DTYPE
        )
    else:
        initial_center = np.ascontiguousarray(
            np.asarray(initial_center_weights, dtype=DTYPE)
        )
        if initial_center.shape != (NUM_NODES, extended_len):
            raise ValueError(
                "initial_center_weights must have shape "
                f"{(NUM_NODES, extended_len)}; got {initial_center.shape}."
            )
    if not all(
        np.isfinite(array).all()
        for array in (
            reference_array,
            disturbance_array,
            secondary,
            compensation,
            initial_center,
        )
    ):
        raise ValueError("FedDMCANC input contains NaN or Inf.")

    intervals = config.communication_intervals()
    schedule = build_communication_schedule(
        num_samples,
        CommunicationMode(config.communication_mode),
        intervals,
    )
    if config.use_numba and NUMBA_AVAILABLE:
        runner = _fed_mcanc_numba
        engine = "numba"
    else:
        if config.use_numba and not NUMBA_AVAILABLE:
            warnings.warn(
                "Numba is unavailable; FedDMCANC is using the slower NumPy "
                "reference kernel.",
                RuntimeWarning,
                stacklevel=2,
            )
        runner = _fed_mcanc_core
        engine = "numpy_reference"

    (
        residual_error,
        control_output,
        secondary_output_store,
        local_weights,
        center_weights,
        weight_difference,
        diagnostic_indices,
        local_norm,
        center_norm,
    ) = runner(
        reference_array,
        disturbance_array,
        secondary,
        compensation,
        initial_center,
        int(config.w_len),
        float(config.mu_w),
        float(config.alpha),
        _mode_code(config.communication_mode),
        intervals,
        bool(save_control_output),
        bool(save_secondary_output),
        int(config.diagnostic_stride),
    )

    arrays_to_check = [
        residual_error,
        local_weights,
        center_weights,
        weight_difference,
        local_norm,
        center_norm,
    ]
    if save_control_output:
        arrays_to_check.append(control_output)
    if save_secondary_output:
        arrays_to_check.append(secondary_output_store)
    if not all(np.isfinite(array).all() for array in arrays_to_check):
        raise FloatingPointError("FedDMCANC produced NaN or Inf.")

    return FedMCANCResult(
        residual_error=np.asarray(residual_error, dtype=DTYPE),
        control_output=(
            np.asarray(control_output, dtype=DTYPE)
            if save_control_output
            else None
        ),
        secondary_output=(
            np.asarray(secondary_output_store, dtype=DTYPE)
            if save_secondary_output
            else None
        ),
        final_local_weights=np.asarray(local_weights, dtype=DTYPE),
        final_center_weights=np.asarray(center_weights, dtype=DTYPE),
        final_weight_difference=np.asarray(weight_difference, dtype=DTYPE),
        diagnostic_sample_indices=np.asarray(
            diagnostic_indices, dtype=np.int64
        ),
        local_weight_norm=np.asarray(local_norm, dtype=DTYPE),
        center_weight_norm=np.asarray(center_norm, dtype=DTYPE),
        communication=schedule,
        engine=engine,
    )
