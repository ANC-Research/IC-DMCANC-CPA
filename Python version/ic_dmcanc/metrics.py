"""Raw and display-smoothed residual metrics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import DTYPE, NUM_NODES
from .signal_utils import matlab_smooth


@dataclass(frozen=True)
class MetricBundle:
    raw_relative_level_db: np.ndarray
    smoothed_relative_level_db: np.ndarray
    raw_mean_relative_level_db: np.ndarray
    smoothed_mean_relative_level_db: np.ndarray
    final_mse_per_node: np.ndarray
    disturbance_mse_per_node: np.ndarray
    final_noise_reduction_db_per_node: np.ndarray
    final_mean_noise_reduction_db: float
    evaluation_sample_count: int


def _safe_db_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    tiny = np.finfo(DTYPE).tiny
    return 10.0 * np.log10(
        np.maximum(numerator, tiny) / np.maximum(denominator, tiny)
    )


def compute_metrics(
    residual_error: np.ndarray,
    disturbance: np.ndarray,
    *,
    evaluation_sample_count: int | None = None,
    smoothing_span: int = 2_000,
    final_fraction: float = 0.1,
) -> MetricBundle:
    """Compute MATLAB-meaningful curves and unsmoothed counterparts."""

    residual = np.asarray(residual_error, dtype=DTYPE)
    primary = np.asarray(disturbance, dtype=DTYPE)
    if residual.ndim != 2 or residual.shape[0] != NUM_NODES:
        raise ValueError("residual_error must have shape (6, num_samples).")
    if primary.shape != residual.shape:
        raise ValueError("disturbance and residual_error shapes must match.")
    if evaluation_sample_count is None:
        count = residual.shape[1]
    else:
        count = int(evaluation_sample_count)
        if count <= 0 or count > residual.shape[1]:
            raise ValueError("evaluation_sample_count is outside the signal.")
    if not 0 < final_fraction <= 1:
        raise ValueError("final_fraction must be in (0, 1].")

    residual = residual[:, :count]
    primary = primary[:, :count]
    residual_squared = residual * residual
    disturbance_squared = primary * primary
    raw = _safe_db_ratio(residual_squared, disturbance_squared)

    smooth_residual = matlab_smooth(residual_squared, smoothing_span)
    smooth_disturbance = matlab_smooth(
        disturbance_squared, smoothing_span
    )
    smoothed = _safe_db_ratio(smooth_residual, smooth_disturbance)

    final_count = max(1, int(np.ceil(count * final_fraction)))
    final_residual_mse = np.mean(
        residual_squared[:, -final_count:], axis=1
    )
    final_disturbance_mse = np.mean(
        disturbance_squared[:, -final_count:], axis=1
    )
    final_nr = _safe_db_ratio(
        final_disturbance_mse, final_residual_mse
    )
    return MetricBundle(
        raw_relative_level_db=raw,
        smoothed_relative_level_db=smoothed,
        raw_mean_relative_level_db=np.mean(raw, axis=0),
        smoothed_mean_relative_level_db=np.mean(smoothed, axis=0),
        final_mse_per_node=final_residual_mse,
        disturbance_mse_per_node=final_disturbance_mse,
        final_noise_reduction_db_per_node=final_nr,
        final_mean_noise_reduction_db=float(np.mean(final_nr)),
        evaluation_sample_count=count,
    )


def metric_summary(bundle: MetricBundle) -> dict[str, object]:
    return {
        "evaluation_sample_count": bundle.evaluation_sample_count,
        "final_mse_per_node": bundle.final_mse_per_node.tolist(),
        "disturbance_mse_per_node": (
            bundle.disturbance_mse_per_node.tolist()
        ),
        "final_noise_reduction_db_per_node": (
            bundle.final_noise_reduction_db_per_node.tolist()
        ),
        "final_mean_noise_reduction_db": (
            bundle.final_mean_noise_reduction_db
        ),
        "all_finite": bool(
            np.isfinite(bundle.raw_relative_level_db).all()
            and np.isfinite(bundle.smoothed_relative_level_db).all()
            and np.isfinite(bundle.final_noise_reduction_db_per_node).all()
        ),
    }
