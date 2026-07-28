"""Case orchestration, diagnostics, plots, and reproducible output files."""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import json
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np

from .centralized_fxlms import (
    CentralizedFxLMSResult,
    run_centralized_fxlms,
)
from .compensation_filter import (
    CompensationIdentificationResult,
    identify_compensation_filters,
)
from .config import (
    NUM_NODES,
    CaseConfig,
    CommunicationMode,
    OutputConfig,
    json_ready,
)
from .data_io import (
    AcousticPaths,
    load_acoustic_paths,
    load_compressor_recording,
)
from .fed_mcanc import FedMCANCResult, run_fed_mcanc
from .metrics import MetricBundle, compute_metrics, metric_summary
from .mgdfxlms import MGDFxLMSResult, run_mgdfxlms
from .plotting import (
    plot_compensation_convergence,
    plot_controller_norms,
    plot_metric_comparison,
)
from .signal_utils import (
    ReferenceSignals,
    filter_primary_paths,
    generate_synthetic_reference,
    prepare_recorded_reference,
)

AlgorithmResult = FedMCANCResult | CentralizedFxLMSResult | MGDFxLMSResult


@dataclass(frozen=True)
class CaseRunResult:
    output_directory: Path
    results_file: Path
    metrics_file: Path
    config_file: Path
    summary_file: Path
    summary: dict[str, Any]


def _safe_label(value: float) -> str:
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p").replace("+", "")


def _prepare_signals(
    case: CaseConfig,
    repository_root: Path | None,
) -> tuple[AcousticPaths, ReferenceSignals, np.ndarray, int]:
    simulation = case.simulation
    paths = load_acoustic_paths(
        repository_root, expected_s_len=simulation.s_len
    )
    if case.reference_kind == "synthetic":
        signals = generate_synthetic_reference(
            fs=simulation.fs,
            duration_seconds=simulation.duration_seconds,
            low_hz=simulation.synthetic_low_hz,
            high_hz=simulation.synthetic_high_hz,
            fir_order=simulation.synthetic_fir_order,
            random_seed=simulation.random_seed,
            awgn_seed=simulation.awgn_seed,
            snr_db=simulation.awgn_snr_db,
        )
        # MATLAB simulates T*Fs+1 points but plots/evaluates 1:T*Fs.
        evaluation_count = signals.reference_clean.size - 1
    else:
        recording = load_compressor_recording(repository_root)
        signals = prepare_recorded_reference(
            recording.signal,
            fs=simulation.fs,
            duration_seconds=simulation.duration_seconds,
            awgn_seed=simulation.awgn_seed,
            snr_db=simulation.awgn_snr_db,
        )
        evaluation_count = signals.reference_clean.size
    disturbance = filter_primary_paths(
        paths.primary_path, signals.reference_clean
    )
    return paths, signals, disturbance, evaluation_count


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _save_npz(path: Path, arrays: dict[str, np.ndarray], compress: bool) -> None:
    if compress:
        np.savez_compressed(path, **arrays)
    else:
        np.savez(path, **arrays)


def _result_arrays(
    label: str,
    result: AlgorithmResult,
    *,
    save_full_error: bool,
    save_full_control: bool,
    save_full_secondary: bool,
) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    prefix = f"{label}__"
    if save_full_error:
        arrays[prefix + "residual_error"] = result.residual_error
    if save_full_control and result.control_output is not None:
        arrays[prefix + "control_output"] = result.control_output

    if isinstance(result, FedMCANCResult):
        if save_full_secondary and result.secondary_output is not None:
            arrays[prefix + "secondary_output"] = result.secondary_output
        arrays[prefix + "final_local_weights"] = result.final_local_weights
        arrays[prefix + "final_center_weights"] = result.final_center_weights
        arrays[prefix + "final_weight_difference"] = (
            result.final_weight_difference
        )
        arrays[prefix + "diagnostic_sample_indices"] = (
            result.diagnostic_sample_indices
        )
        arrays[prefix + "local_weight_norm"] = result.local_weight_norm
        arrays[prefix + "center_weight_norm"] = result.center_weight_norm
        arrays[prefix + "communication_sample_indices"] = (
            result.communication.sample_indices
        )
        arrays[prefix + "communication_node_indices"] = (
            result.communication.node_indices
        )
        arrays[prefix + "communication_counts"] = (
            result.communication.per_node_counts
        )
    elif isinstance(result, CentralizedFxLMSResult):
        arrays[prefix + "final_weights"] = result.final_weights
        arrays[prefix + "diagnostic_sample_indices"] = (
            result.diagnostic_sample_indices
        )
        arrays[prefix + "controller_norm"] = result.controller_norm
    else:
        arrays[prefix + "final_weights"] = result.final_weights
        arrays[prefix + "final_gradient"] = result.final_gradient
        arrays[prefix + "diagnostic_sample_indices"] = (
            result.diagnostic_sample_indices
        )
        arrays[prefix + "controller_norm"] = result.controller_norm
    return arrays


def _metric_arrays(
    label: str, bundle: MetricBundle
) -> dict[str, np.ndarray]:
    prefix = f"{label}__"
    return {
        prefix + "raw_relative_level_db": bundle.raw_relative_level_db,
        prefix + "smoothed_relative_level_db": (
            bundle.smoothed_relative_level_db
        ),
        prefix + "raw_mean_relative_level_db": (
            bundle.raw_mean_relative_level_db
        ),
        prefix + "smoothed_mean_relative_level_db": (
            bundle.smoothed_mean_relative_level_db
        ),
        prefix + "final_mse_per_node": bundle.final_mse_per_node,
        prefix + "disturbance_mse_per_node": (
            bundle.disturbance_mse_per_node
        ),
        prefix + "final_noise_reduction_db_per_node": (
            bundle.final_noise_reduction_db_per_node
        ),
    }


def _write_communication_csv(
    path: Path,
    result: FedMCANCResult,
) -> None:
    schedule = result.communication
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "matlab_sample_index",
                "python_sample_index",
                "matlab_node_index",
                "python_node_index",
                "mode",
                "interval_samples",
            ]
        )
        for sample, node in zip(
            schedule.sample_indices, schedule.node_indices, strict=True
        ):
            if node < 0:
                matlab_node: int | str = "all"
                python_node: int | str = "all"
                interval = int(schedule.interval_samples[0])
            else:
                matlab_node = int(node) + 1
                python_node = int(node)
                interval = int(schedule.interval_samples[node])
            writer.writerow(
                [
                    int(sample),
                    int(sample) - 1,
                    matlab_node,
                    python_node,
                    schedule.mode.value,
                    interval,
                ]
            )


def _algorithm_summary(
    result: AlgorithmResult,
    metrics: MetricBundle,
    runtime_seconds: float,
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "engine": result.engine,
        "runtime_seconds": runtime_seconds,
        "residual_shape": list(result.residual_error.shape),
        "control_output_saved": result.control_output is not None,
        "all_finite": bool(
            np.isfinite(result.residual_error).all()
            and (
                result.control_output is None
                or np.isfinite(result.control_output).all()
            )
        ),
        "metrics": metric_summary(metrics),
    }
    if isinstance(result, FedMCANCResult):
        schedule = result.communication
        summary.update(
            {
                "final_local_weight_norm": np.linalg.norm(
                    result.final_local_weights, axis=1
                ).tolist(),
                "final_center_weight_norm": np.linalg.norm(
                    result.final_center_weights, axis=1
                ).tolist(),
                "communication": {
                    "mode": schedule.mode.value,
                    "interval_samples": schedule.interval_samples.tolist(),
                    "per_node_counts": schedule.per_node_counts.tolist(),
                    "total_node_uploads": schedule.total_node_uploads,
                    "total_event_rows": schedule.total_event_rows,
                    "first_matlab_sample_index": (
                        int(schedule.sample_indices[0])
                        if schedule.sample_indices.size
                        else None
                    ),
                    "last_matlab_sample_index": (
                        int(schedule.sample_indices[-1])
                        if schedule.sample_indices.size
                        else None
                    ),
                },
            }
        )
    else:
        summary["final_controller_norm"] = np.linalg.norm(
            result.final_weights, axis=1
        ).tolist()
    return summary


def _memory_bytes(arrays: dict[str, np.ndarray]) -> int:
    return int(sum(array.nbytes for array in arrays.values()))


def run_case(
    case: CaseConfig,
    output: OutputConfig = OutputConfig(),
    *,
    repository_root: Path | None = None,
) -> CaseRunResult:
    """Run one case with independent state for every sweep point."""

    case = case.validate()
    simulation = case.simulation
    case_directory = Path(output.output_directory) / f"case{case.case_id}"
    case_directory.mkdir(parents=True, exist_ok=True)

    total_start = perf_counter()
    paths, signals, disturbance, evaluation_count = _prepare_signals(
        case, repository_root
    )

    compensation_start = perf_counter()
    mgd_compensation: CompensationIdentificationResult | None = None
    mgd_compensation_runtime = 0.0
    if case.case_id in {1, 2} and case.include_baselines:
        mgd_start = perf_counter()
        mgd_compensation = identify_compensation_filters(
            paths.secondary_path,
            c_len=simulation.c_len,
            step_size=simulation.mu_c,
            num_samples=simulation.comp_id_num_samples,
            seed=simulation.comp_id_seed,
            convergence_stride=simulation.comp_id_convergence_stride,
            save_full_error=output.save_full_comp_id_error,
            use_numba=simulation.use_numba,
        )
        mgd_compensation_runtime = perf_counter() - mgd_start
        proposed_compensation_seed = simulation.comp_id_seed + 1
    else:
        proposed_compensation_seed = simulation.comp_id_seed
    proposed_start = perf_counter()
    compensation = identify_compensation_filters(
        paths.secondary_path,
        c_len=simulation.c_len,
        step_size=simulation.mu_c,
        num_samples=simulation.comp_id_num_samples,
        seed=proposed_compensation_seed,
        convergence_stride=simulation.comp_id_convergence_stride,
        save_full_error=output.save_full_comp_id_error,
        use_numba=simulation.use_numba,
    )
    compensation_runtime = perf_counter() - proposed_start
    compensation_total_runtime = perf_counter() - compensation_start

    algorithms: dict[str, AlgorithmResult] = {}
    runtimes: dict[str, float] = {}

    def run_timed(label: str, callback: Any) -> None:
        start = perf_counter()
        algorithms[label] = callback()
        runtimes[label] = perf_counter() - start

    if case.case_id in {1, 2}:
        if case.include_baselines:
            run_timed(
                "centralized",
                lambda: run_centralized_fxlms(
                    signals.reference_input,
                    disturbance,
                    paths.secondary_path,
                    simulation,
                    save_control_output=output.save_full_control_output,
                ),
            )
            run_timed(
                "mgdfxlms",
                lambda: run_mgdfxlms(
                    signals.reference_input,
                    disturbance,
                    paths.secondary_path,
                    mgd_compensation.compensation_filters,
                    simulation,
                    save_control_output=output.save_full_control_output,
                ),
            )
        ideal_config = simulation.with_updates(
            communication_mode=CommunicationMode.IDEAL,
            alpha=case.ideal_alpha,
        )
        run_timed(
            "feddmcanc_ideal",
            lambda: run_fed_mcanc(
                signals.reference_input,
                disturbance,
                paths.secondary_path,
                compensation.compensation_filters,
                ideal_config,
                save_control_output=output.save_full_control_output,
                save_secondary_output=output.save_full_secondary_output,
            ),
        )
    elif case.case_id == 3:
        ideal_config = simulation.with_updates(
            communication_mode=CommunicationMode.IDEAL,
            alpha=case.ideal_alpha,
        )
        run_timed(
            "feddmcanc_ideal",
            lambda: run_fed_mcanc(
                signals.reference_input,
                disturbance,
                paths.secondary_path,
                compensation.compensation_filters,
                ideal_config,
                save_control_output=output.save_full_control_output,
                save_secondary_output=output.save_full_secondary_output,
            ),
        )
        for tc in case.tc_sweep_seconds:
            periodic_config = simulation.with_updates(
                communication_mode=CommunicationMode.COMMON_PERIODIC,
                tc_seconds=float(tc),
                alpha=case.ideal_alpha,
            )
            label = f"feddmcanc_tc_{_safe_label(tc)}s"
            run_timed(
                label,
                lambda config=periodic_config: run_fed_mcanc(
                    signals.reference_input,
                    disturbance,
                    paths.secondary_path,
                    compensation.compensation_filters,
                    config,
                    save_control_output=output.save_full_control_output,
                    save_secondary_output=output.save_full_secondary_output,
                ),
            )
    elif case.case_id == 4:
        ideal_config = simulation.with_updates(
            communication_mode=CommunicationMode.IDEAL,
            alpha=case.ideal_alpha,
        )
        run_timed(
            "feddmcanc_ideal_alpha_1000",
            lambda: run_fed_mcanc(
                signals.reference_input,
                disturbance,
                paths.secondary_path,
                compensation.compensation_filters,
                ideal_config,
                save_control_output=output.save_full_control_output,
                save_secondary_output=output.save_full_secondary_output,
            ),
        )
        tc = case.tc_sweep_seconds[0]
        for alpha in case.alpha_sweep:
            periodic_config = simulation.with_updates(
                communication_mode=CommunicationMode.COMMON_PERIODIC,
                tc_seconds=float(tc),
                alpha=float(alpha),
            )
            label = f"feddmcanc_alpha_{_safe_label(alpha)}"
            run_timed(
                label,
                lambda config=periodic_config: run_fed_mcanc(
                    signals.reference_input,
                    disturbance,
                    paths.secondary_path,
                    compensation.compensation_filters,
                    config,
                    save_control_output=output.save_full_control_output,
                    save_secondary_output=output.save_full_secondary_output,
                ),
            )
    else:
        ideal_config = simulation.with_updates(
            communication_mode=CommunicationMode.IDEAL,
            alpha=case.ideal_alpha,
        )
        run_timed(
            "feddmcanc_ideal",
            lambda: run_fed_mcanc(
                signals.reference_input,
                disturbance,
                paths.secondary_path,
                compensation.compensation_filters,
                ideal_config,
                save_control_output=output.save_full_control_output,
                save_secondary_output=output.save_full_secondary_output,
            ),
        )
        individual_config = simulation.with_updates(
            communication_mode=CommunicationMode.INDIVIDUAL_PERIODIC,
            tc_seconds=tuple(case.individual_tc_seconds),
            alpha=case.ideal_alpha,
        )
        run_timed(
            "feddmcanc_individual",
            lambda: run_fed_mcanc(
                signals.reference_input,
                disturbance,
                paths.secondary_path,
                compensation.compensation_filters,
                individual_config,
                save_control_output=output.save_full_control_output,
                save_secondary_output=output.save_full_secondary_output,
            ),
        )

    metrics = {
        label: compute_metrics(
            result.residual_error,
            disturbance,
            evaluation_sample_count=evaluation_count,
        )
        for label, result in algorithms.items()
    }

    result_arrays: dict[str, np.ndarray] = {
        "reference_clean": signals.reference_clean,
        "reference_input": signals.reference_input,
        "awgn_noise": signals.awgn_noise,
        "disturbance": disturbance,
        "primary_path": paths.primary_path,
        "secondary_path": paths.secondary_path,
        "compensation_filters": compensation.compensation_filters,
        "compensation_identified_coefficients": (
            compensation.identified_coefficients
        ),
        "compensation_final_mse": compensation.final_mse,
        "compensation_convergence_sample_indices": (
            compensation.convergence_sample_indices
        ),
        "compensation_convergence_error_squared": (
            compensation.convergence_error_squared
        ),
        "compensation_identification_white_noise": compensation.white_noise,
    }
    if mgd_compensation is not None:
        result_arrays.update(
            {
                "mgdfxlms_compensation_filters": (
                    mgd_compensation.compensation_filters
                ),
                "mgdfxlms_compensation_identified_coefficients": (
                    mgd_compensation.identified_coefficients
                ),
                "mgdfxlms_compensation_final_mse": (
                    mgd_compensation.final_mse
                ),
                "mgdfxlms_compensation_convergence_sample_indices": (
                    mgd_compensation.convergence_sample_indices
                ),
                "mgdfxlms_compensation_convergence_error_squared": (
                    mgd_compensation.convergence_error_squared
                ),
                "mgdfxlms_compensation_identification_white_noise": (
                    mgd_compensation.white_noise
                ),
            }
        )
    if signals.fir_coefficients is not None:
        result_arrays["synthetic_fir_coefficients"] = (
            signals.fir_coefficients
        )
    if compensation.full_error is not None:
        result_arrays["compensation_full_error"] = compensation.full_error
    if (
        mgd_compensation is not None
        and mgd_compensation.full_error is not None
    ):
        result_arrays["mgdfxlms_compensation_full_error"] = (
            mgd_compensation.full_error
        )
    for label, result in algorithms.items():
        result_arrays.update(
            _result_arrays(
                label,
                result,
                save_full_error=output.save_full_error,
                save_full_control=output.save_full_control_output,
                save_full_secondary=output.save_full_secondary_output,
            )
        )

    metric_arrays: dict[str, np.ndarray] = {}
    for label, bundle in metrics.items():
        metric_arrays.update(_metric_arrays(label, bundle))

    results_path = case_directory / "results.npz"
    metrics_path = case_directory / "metrics.npz"
    _save_npz(results_path, result_arrays, output.compress_npz)
    _save_npz(metrics_path, metric_arrays, output.compress_npz)

    if output.save_communication_events:
        for label, result in algorithms.items():
            if isinstance(result, FedMCANCResult):
                _write_communication_csv(
                    case_directory / f"communication_events__{label}.csv",
                    result,
                )

    plot_files: list[str] = []
    if output.save_plots:
        plot_files.extend(
            str(path.name)
            for path in plot_metric_comparison(
                metrics,
                fs=simulation.fs,
                output_directory=case_directory,
                title=f"IC-DMCANC-CPA Python case {case.case_id}",
            )
        )
        norm_diagnostics: dict[
            str, tuple[np.ndarray, np.ndarray]
        ] = {}
        for label, result in algorithms.items():
            if isinstance(result, FedMCANCResult):
                norm_diagnostics[label + " local"] = (
                    result.diagnostic_sample_indices,
                    result.local_weight_norm,
                )
                norm_diagnostics[label + " center"] = (
                    result.diagnostic_sample_indices,
                    result.center_weight_norm,
                )
            else:
                norm_diagnostics[label] = (
                    result.diagnostic_sample_indices,
                    result.controller_norm,
                )
        plot_files.append(
            plot_controller_norms(
                norm_diagnostics,
                fs=simulation.fs,
                output_directory=case_directory,
            ).name
        )
        plot_files.append(
            plot_compensation_convergence(
                compensation.convergence_sample_indices,
                compensation.convergence_error_squared,
                output_directory=case_directory,
            ).name
        )

    total_runtime = perf_counter() - total_start
    summary: dict[str, Any] = {
        "case_id": case.case_id,
        "case_name": case.name,
        "reference_kind": case.reference_kind,
        "matlab_authority": True,
        "num_nodes": NUM_NODES,
        "sample_count": int(signals.reference_input.size),
        "evaluation_sample_count": int(evaluation_count),
        "input_shapes": {
            "primary_path": list(paths.primary_path.shape),
            "secondary_path": list(paths.secondary_path.shape),
            "disturbance": list(disturbance.shape),
            "reference_clean": list(signals.reference_clean.shape),
            "reference_input": list(signals.reference_input.shape),
        },
        "seeds": {
            "random_seed": simulation.random_seed,
            "awgn_seed": simulation.awgn_seed,
            "comp_id_seed_mgdfxlms": (
                simulation.comp_id_seed
                if mgd_compensation is not None
                else None
            ),
            "comp_id_seed_feddmcanc": proposed_compensation_seed,
        },
        "compensation_identification": {
            "engine": compensation.engine,
            "runtime_seconds": compensation_runtime,
            "num_samples": simulation.comp_id_num_samples,
            "seed": proposed_compensation_seed,
            "shared_white_noise": compensation.shared_white_noise,
            "skipped_diagonal": compensation.skipped_diagonal,
            "output_transform": compensation.output_transform,
            "all_finite": bool(
                np.isfinite(compensation.compensation_filters).all()
            ),
            "final_mse": compensation.final_mse.tolist(),
        },
        "compensation_total_runtime_seconds": compensation_total_runtime,
        "algorithms": {
            label: _algorithm_summary(result, metrics[label], runtimes[label])
            for label, result in algorithms.items()
        },
        "sweep_state_mode": "independent",
        "external_baselines": {
            "DFxLMS": "not implemented; optional external interface only",
            "ADFxLMS": "not implemented; optional external interface only",
        },
        "output": {
            "result_array_bytes_uncompressed": _memory_bytes(result_arrays),
            "metric_array_bytes_uncompressed": _memory_bytes(metric_arrays),
            "plot_files": plot_files,
        },
        "total_runtime_seconds": total_runtime,
    }
    if mgd_compensation is not None:
        summary["mgdfxlms_compensation_identification"] = {
            "engine": mgd_compensation.engine,
            "runtime_seconds": mgd_compensation_runtime,
            "num_samples": simulation.comp_id_num_samples,
            "seed": simulation.comp_id_seed,
            "shared_white_noise": mgd_compensation.shared_white_noise,
            "skipped_diagonal": mgd_compensation.skipped_diagonal,
            "output_transform": mgd_compensation.output_transform,
            "all_finite": bool(
                np.isfinite(
                    mgd_compensation.compensation_filters
                ).all()
            ),
            "final_mse": mgd_compensation.final_mse.tolist(),
        }

    config_path = case_directory / "config.json"
    summary_path = case_directory / "summary.json"
    _write_json(
        config_path,
        {
            "case": case,
            "output": output,
            "resolved_communication_intervals": {
                label: (
                    result.communication.interval_samples.tolist()
                    if isinstance(result, FedMCANCResult)
                    else None
                )
                for label, result in algorithms.items()
            },
        },
    )
    _write_json(summary_path, summary)
    return CaseRunResult(
        output_directory=case_directory,
        results_file=results_path,
        metrics_file=metrics_path,
        config_file=config_path,
        summary_file=summary_path,
        summary=summary,
    )


def with_output_directory(output: OutputConfig, path: Path) -> OutputConfig:
    return replace(output, output_directory=Path(path))
