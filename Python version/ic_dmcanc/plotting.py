"""Headless plots matching the meanings of the five MATLAB cases."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ic_dmcanc_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .config import NUM_NODES
from .metrics import MetricBundle
from .signal_utils import matlab_smooth


def _display_curve(
    curve: np.ndarray,
    *,
    second_smoothing_span: int = 5_000,
) -> tuple[np.ndarray, int]:
    """Apply the scripts' display-only trim and second ``smooth`` call."""

    values = np.asarray(curve, dtype=np.float64).reshape(-1)
    start = 99 if values.size > 99 else 0
    stop = values.size - 1_000 if values.size > start + 1_000 else values.size
    trimmed = values[start:stop]
    if trimmed.size == 0:
        trimmed = values
        start = 0
    effective_span = min(second_smoothing_span, trimmed.size)
    if effective_span <= 1:
        return trimmed.copy(), start
    return matlab_smooth(trimmed, effective_span), start


def plot_metric_comparison(
    metrics: dict[str, MetricBundle],
    *,
    fs: int,
    output_directory: Path,
    title: str,
) -> list[Path]:
    """Plot six node curves and their node-wise average."""

    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []

    figure, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    for node, axis in enumerate(axes.flat):
        for label, bundle in metrics.items():
            curve, start = _display_curve(
                bundle.smoothed_relative_level_db[node]
            )
            time = (np.arange(curve.size) + start) / fs
            axis.plot(time, curve, label=label, linewidth=1.1)
        axis.set_title(f"Node {node + 1}")
        axis.set_ylabel("Residual / disturbance (dB)")
        axis.grid(True, alpha=0.3)
    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.suptitle(title, y=0.995)
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.96),
        ncol=max(1, len(labels)),
    )
    figure.tight_layout(rect=(0, 0, 1, 0.90))
    per_node_path = destination / "residual_level_per_node.png"
    figure.savefig(per_node_path, dpi=160)
    plt.close(figure)
    created.append(per_node_path)

    figure, axis = plt.subplots(figsize=(10, 4.8))
    for label, bundle in metrics.items():
        curve, start = _display_curve(
            bundle.smoothed_mean_relative_level_db
        )
        time = (np.arange(curve.size) + start) / fs
        axis.plot(time, curve, label=label, linewidth=1.2)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Six-node mean residual / disturbance (dB)")
    axis.set_title(f"{title} — six-node mean")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    average_path = destination / "residual_level_mean.png"
    figure.savefig(average_path, dpi=160)
    plt.close(figure)
    created.append(average_path)
    return created


def plot_controller_norms(
    diagnostics: dict[str, tuple[np.ndarray, np.ndarray]],
    *,
    fs: int,
    output_directory: Path,
) -> Path:
    """Plot the RMS-across-node controller norm diagnostics."""

    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(10, 4.8))
    for label, (sample_indices, per_node_norm) in diagnostics.items():
        norm = np.sqrt(np.mean(np.asarray(per_node_norm) ** 2, axis=0))
        axis.plot(np.asarray(sample_indices) / fs, norm, label=label)
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("RMS node controller L2 norm")
    axis.set_title("Controller norm diagnostics")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    path = destination / "controller_norms.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path


def plot_compensation_convergence(
    sample_indices: np.ndarray,
    error_squared: np.ndarray,
    *,
    output_directory: Path,
) -> Path:
    """Plot all 30 off-diagonal system-ID squared-error traces."""

    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    for error_node, axis in enumerate(axes.flat):
        for source in range(NUM_NODES):
            if source == error_node:
                continue
            axis.semilogy(
                sample_indices,
                np.maximum(
                    error_squared[error_node, source],
                    np.finfo(np.float64).tiny,
                ),
                label=f"{error_node + 1},{source + 1}",
                linewidth=0.9,
            )
        axis.set_title(f"Error node {error_node + 1}")
        axis.grid(True, which="both", alpha=0.25)
    axes[-1, 0].set_xlabel("Identification sample")
    axes[-1, 1].set_xlabel("Identification sample")
    figure.supylabel("System-ID error squared")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.suptitle("Compensation-filter identification", y=0.995)
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.96),
        ncol=5,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.90))
    path = destination / "compensation_identification.png"
    figure.savefig(path, dpi=160)
    plt.close(figure)
    return path
