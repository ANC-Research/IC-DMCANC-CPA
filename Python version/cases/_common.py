"""Shared command-line interface for all five case entry points."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from ic_dmcanc.config import (
    CaseConfig,
    OutputConfig,
    case_defaults,
    normalize_float_sequence,
)
from ic_dmcanc.runner import run_case


def _parser(case_id: int) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            f"Run the six-node IC-DMCANC-CPA Python reference case {case_id}."
        )
    )
    parser.add_argument("--fs", type=int, help="Sampling rate in Hz.")
    parser.add_argument("--duration", type=float, help="Duration in seconds.")
    parser.add_argument("--w-len", type=int, help="Controller FIR length.")
    parser.add_argument(
        "--s-len",
        type=int,
        help="Secondary-path length (must match the MAT asset).",
    )
    parser.add_argument("--c-len", type=int, help="Compensation FIR length.")
    parser.add_argument("--mu-w", type=float, help="Controller step size.")
    parser.add_argument("--mu-c", type=float, help="System-ID step size.")
    parser.add_argument("--alpha", type=float, help="Center-attraction value.")
    parser.add_argument(
        "--tc",
        type=float,
        nargs="+",
        help=(
            "Communication period(s) in seconds: case3 sweep, one value for "
            "case4, or six node periods for case5."
        ),
    )
    parser.add_argument(
        "--alpha-values",
        type=float,
        nargs="+",
        help="Case4 alpha sweep values.",
    )
    parser.add_argument("--random-seed", type=int)
    parser.add_argument("--awgn-seed", type=int)
    parser.add_argument("--comp-id-seed", type=int)
    parser.add_argument("--awgn-snr-db", type=float)
    parser.add_argument("--comp-id-samples", type=int)
    parser.add_argument("--comp-id-convergence-stride", type=int)
    parser.add_argument("--diagnostic-stride", type=int)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Parent directory; a caseN subdirectory is created.",
    )
    parser.add_argument(
        "--no-numba",
        action="store_true",
        help="Force the readable NumPy reference kernels.",
    )
    parser.add_argument(
        "--no-baselines",
        action="store_true",
        help="Skip centralized and MGDFxLMS in cases 1 and 2.",
    )
    parser.add_argument(
        "--no-save-error",
        action="store_true",
        help="Do not persist full residual arrays (metrics are still saved).",
    )
    parser.add_argument(
        "--no-save-control",
        action="store_true",
        help="Do not retain or persist full control-output arrays.",
    )
    parser.add_argument(
        "--save-secondary-output",
        action="store_true",
        help="Persist full physical secondary output for FedDMCANC.",
    )
    parser.add_argument(
        "--save-full-comp-id-error",
        action="store_true",
        help="Persist the 6x6xN compensation-ID error tensor.",
    )
    parser.add_argument("--no-events", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--uncompressed", action="store_true")
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Short smoke-test settings: 0.01 s, w_len=16, c_len=5, "
            "128 compensation-ID samples, and reference kernels."
        ),
    )
    return parser


def _apply_arguments(case: CaseConfig, args: argparse.Namespace) -> CaseConfig:
    updates = {}
    mapping = {
        "fs": "fs",
        "duration": "duration_seconds",
        "w_len": "w_len",
        "s_len": "s_len",
        "c_len": "c_len",
        "mu_w": "mu_w",
        "mu_c": "mu_c",
        "alpha": "alpha",
        "random_seed": "random_seed",
        "awgn_seed": "awgn_seed",
        "comp_id_seed": "comp_id_seed",
        "awgn_snr_db": "awgn_snr_db",
        "comp_id_samples": "comp_id_num_samples",
        "comp_id_convergence_stride": "comp_id_convergence_stride",
        "diagnostic_stride": "diagnostic_stride",
    }
    for argument_name, config_name in mapping.items():
        value = getattr(args, argument_name)
        if value is not None:
            updates[config_name] = value
    if args.no_numba:
        updates["use_numba"] = False
    if args.quick:
        updates.update(
            {
                "duration_seconds": 0.01,
                "w_len": 16,
                "c_len": 5,
                "comp_id_num_samples": 128,
                "comp_id_convergence_stride": 16,
                "diagnostic_stride": 16,
                "use_numba": False,
            }
        )
    updated = case.with_simulation_updates(**updates) if updates else case
    if args.alpha is not None:
        updated = replace(updated, ideal_alpha=float(args.alpha))
    if args.no_baselines:
        updated = replace(updated, include_baselines=False)
    if args.tc is not None:
        periods = normalize_float_sequence(args.tc)
        if updated.case_id == 3:
            updated = replace(updated, tc_sweep_seconds=periods)
        elif updated.case_id == 4:
            if len(periods) != 1:
                raise ValueError("Case 4 --tc requires exactly one value.")
            updated = replace(updated, tc_sweep_seconds=periods)
        elif updated.case_id == 5:
            if len(periods) != 6:
                raise ValueError("Case 5 --tc requires exactly six values.")
            updated = replace(updated, individual_tc_seconds=periods)
        else:
            raise ValueError("--tc applies to cases 3, 4, and 5.")
    if args.alpha_values is not None:
        if updated.case_id != 4:
            raise ValueError("--alpha-values applies only to case 4.")
        updated = replace(
            updated,
            alpha_sweep=normalize_float_sequence(args.alpha_values),
        )
    return updated.validate()


def main_for_case(
    case_id: int,
    argv: Sequence[str] | None = None,
) -> int:
    args = _parser(case_id).parse_args(argv)
    case = _apply_arguments(case_defaults(case_id), args)
    output = OutputConfig(
        output_directory=args.output_dir,
        save_full_error=not args.no_save_error,
        save_full_control_output=not args.no_save_control,
        save_full_secondary_output=args.save_secondary_output,
        save_full_comp_id_error=args.save_full_comp_id_error,
        save_communication_events=not args.no_events,
        save_plots=not args.no_plots,
        compress_npz=not args.uncompressed,
    )
    result = run_case(case, output)
    print(f"Case {case_id} completed.")
    print(f"Summary: {result.summary_file.resolve()}")
    print(f"Results: {result.results_file.resolve()}")
    return 0
