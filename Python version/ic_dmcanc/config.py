"""Central configuration and validation for the fixed six-node reference."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import numpy as np

NUM_NODES = 6
DTYPE = np.float64
COMM_INTERVAL_ROUNDING = "require_integer"
SWEEP_STATE_MODE = "independent"


class CommunicationMode(str, Enum):
    """Communication modes implemented by ``FedMCANC.m``."""

    IDEAL = "ideal"
    FIXED_10_SECONDS = "fixed_10_seconds"
    COMMON_PERIODIC = "common_periodic"
    INDIVIDUAL_PERIODIC = "individual_periodic"


def require_integer_samples(fs: int, seconds: float, *, name: str) -> int:
    """Convert seconds to samples without floor, ceil, round, or truncation.

    ``Fraction(str(...))`` treats user-facing decimal values such as ``0.1``
    exactly, while still rejecting a period whose sample count is not an
    integer. This is the Python equivalent of requiring a valid integer
    divisor in the MATLAB ``mod(i, Fs * tc)`` expression.
    """

    if fs <= 0:
        raise ValueError(f"fs must be positive; got {fs}.")
    if not np.isfinite(seconds) or seconds <= 0:
        raise ValueError(f"{name} must be finite and positive; got {seconds}.")
    raw = Fraction(fs) * Fraction(str(seconds))
    if raw.denominator != 1:
        raise ValueError(
            f"{name}={seconds} s gives {float(raw):.17g} samples at fs={fs}; "
            "COMM_INTERVAL_ROUNDING='require_integer' forbids implicit "
            "floor/ceil/round/truncation."
        )
    samples = raw.numerator
    if samples <= 0:
        raise ValueError(f"{name} must resolve to more than zero samples.")
    return samples


@dataclass(frozen=True)
class SimulationConfig:
    """Algorithm and experiment settings shared by the five cases."""

    fs: int = 16_000
    duration_seconds: float = 90.0
    num_nodes: int = NUM_NODES
    w_len: int = 512
    s_len: int = 256
    c_len: int = 33
    mu_w: float = 1e-6
    mu_c: float = 1e-5
    alpha: float = 1_000.0
    communication_mode: CommunicationMode = CommunicationMode.IDEAL
    tc_seconds: float | tuple[float, ...] = 0.5
    random_seed: int = 0
    awgn_seed: int = 1
    comp_id_seed: int = 0
    awgn_snr_db: float = 40.0
    synthetic_low_hz: float = 100.0
    synthetic_high_hz: float = 1_000.0
    synthetic_fir_order: int = 63
    comp_id_num_samples: int = 200_000
    comp_id_convergence_stride: int = 1_000
    diagnostic_stride: int = 1_600
    use_numba: bool = True

    def validate(self) -> "SimulationConfig":
        if self.num_nodes != NUM_NODES:
            raise ValueError(
                f"This reference is fixed at NUM_NODES={NUM_NODES}; "
                f"got {self.num_nodes}."
            )
        for name in ("fs", "w_len", "s_len", "c_len"):
            value = int(getattr(self, name))
            if value <= 0:
                raise ValueError(f"{name} must be positive; got {value}.")
        for name in ("mu_w", "mu_c"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative.")
        if not np.isfinite(self.alpha) or self.alpha < 0:
            raise ValueError("alpha must be finite and nonnegative.")
        if self.comp_id_num_samples <= 0:
            raise ValueError("comp_id_num_samples must be positive.")
        if self.comp_id_convergence_stride <= 0:
            raise ValueError("comp_id_convergence_stride must be positive.")
        if self.diagnostic_stride <= 0:
            raise ValueError("diagnostic_stride must be positive.")
        if not 0 < self.synthetic_low_hz < self.synthetic_high_hz < self.fs / 2:
            raise ValueError(
                "Synthetic band must satisfy "
                "0 < low_hz < high_hz < fs/2."
            )
        if self.synthetic_fir_order < 1:
            raise ValueError("synthetic_fir_order must be at least one.")
        require_integer_samples(
            self.fs, self.duration_seconds, name="duration_seconds"
        )
        self.communication_intervals()
        return self

    @property
    def synthetic_num_samples(self) -> int:
        """Match MATLAB ``0:1/Fs:T``: include both endpoints."""

        return (
            require_integer_samples(
                self.fs, self.duration_seconds, name="duration_seconds"
            )
            + 1
        )

    @property
    def extended_weight_len(self) -> int:
        return self.w_len + self.c_len - 1

    def communication_intervals(self) -> np.ndarray:
        """Return one or six integer periods, depending on the mode."""

        mode = CommunicationMode(self.communication_mode)
        if mode is CommunicationMode.IDEAL:
            return np.ones(1, dtype=np.int64)
        if mode is CommunicationMode.FIXED_10_SECONDS:
            return np.asarray(
                [require_integer_samples(self.fs, 10.0, name="fixed interval")],
                dtype=np.int64,
            )
        if mode is CommunicationMode.COMMON_PERIODIC:
            if isinstance(self.tc_seconds, tuple):
                if len(self.tc_seconds) != 1:
                    raise ValueError(
                        "common_periodic requires exactly one tc_seconds value."
                    )
                tc = self.tc_seconds[0]
            else:
                tc = self.tc_seconds
            return np.asarray(
                [require_integer_samples(self.fs, tc, name="tc_seconds")],
                dtype=np.int64,
            )
        if not isinstance(self.tc_seconds, tuple):
            raise ValueError(
                "individual_periodic requires a six-value tc_seconds tuple."
            )
        if len(self.tc_seconds) != NUM_NODES:
            raise ValueError(
                "individual_periodic requires exactly six tc_seconds values."
            )
        return np.asarray(
            [
                require_integer_samples(
                    self.fs, tc, name=f"tc_seconds[{index}]"
                )
                for index, tc in enumerate(self.tc_seconds)
            ],
            dtype=np.int64,
        )

    def with_updates(self, **updates: Any) -> "SimulationConfig":
        return replace(self, **updates).validate()


@dataclass(frozen=True)
class OutputConfig:
    """Controls diagnostic persistence without changing the algorithm."""

    output_directory: Path = Path("outputs")
    save_full_error: bool = True
    save_full_control_output: bool = True
    save_full_secondary_output: bool = False
    save_full_comp_id_error: bool = False
    save_communication_events: bool = True
    save_controller_history: bool = False
    save_plots: bool = True
    compress_npz: bool = True


@dataclass(frozen=True)
class CaseConfig:
    """Case-specific definitions copied from the MATLAB scripts."""

    case_id: int
    name: str
    reference_kind: str
    simulation: SimulationConfig
    include_baselines: bool = False
    tc_sweep_seconds: tuple[float, ...] = field(default_factory=tuple)
    alpha_sweep: tuple[float, ...] = field(default_factory=tuple)
    individual_tc_seconds: tuple[float, ...] = field(default_factory=tuple)
    ideal_alpha: float = 1_000.0

    def validate(self) -> "CaseConfig":
        if self.case_id not in range(1, 6):
            raise ValueError("case_id must be in 1..5.")
        if self.reference_kind not in {"synthetic", "compressor"}:
            raise ValueError("reference_kind must be synthetic or compressor.")
        self.simulation.validate()
        if self.case_id == 3 and not self.tc_sweep_seconds:
            raise ValueError("Case 3 requires tc_sweep_seconds.")
        if self.case_id == 4 and not self.alpha_sweep:
            raise ValueError("Case 4 requires alpha_sweep.")
        if self.case_id == 5 and len(self.individual_tc_seconds) != NUM_NODES:
            raise ValueError("Case 5 requires six individual periods.")
        return self

    def with_simulation_updates(self, **updates: Any) -> "CaseConfig":
        return replace(
            self, simulation=self.simulation.with_updates(**updates)
        ).validate()


def case_defaults(case_id: int) -> CaseConfig:
    """Return the parameters currently present in ``FedDMCANC_case*.m``."""

    base = SimulationConfig()
    if case_id == 1:
        return CaseConfig(
            case_id=1,
            name="synthetic_ideal_baseline_comparison",
            reference_kind="synthetic",
            simulation=base,
            include_baselines=True,
        ).validate()
    if case_id == 2:
        return CaseConfig(
            case_id=2,
            name="compressor_ideal_baseline_comparison",
            reference_kind="compressor",
            simulation=base.with_updates(
                duration_seconds=14.0,
                mu_w=3e-6,
            ),
            include_baselines=True,
        ).validate()
    if case_id == 3:
        return CaseConfig(
            case_id=3,
            name="common_communication_interval_sweep",
            reference_kind="synthetic",
            simulation=base,
            tc_sweep_seconds=(0.1, 0.5, 1.0, 3.0),
        ).validate()
    if case_id == 4:
        return CaseConfig(
            case_id=4,
            name="center_attraction_alpha_sweep",
            reference_kind="synthetic",
            simulation=base,
            alpha_sweep=(300.0, 600.0, 1_000.0, 2_000.0, 5_000.0, 2.1e6),
            tc_sweep_seconds=(0.5,),
        ).validate()
    if case_id == 5:
        return CaseConfig(
            case_id=5,
            name="individual_communication_periods",
            reference_kind="synthetic",
            simulation=base,
            individual_tc_seconds=(0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
        ).validate()
    raise ValueError(f"Unknown case_id={case_id}; expected 1..5.")


def json_ready(value: Any) -> Any:
    """Recursively convert dataclass configuration to JSON-safe values."""

    if hasattr(value, "__dataclass_fields__"):
        return json_ready(asdict(value))
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_ready(item) for item in value]
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def normalize_float_sequence(values: Sequence[float]) -> tuple[float, ...]:
    return tuple(float(value) for value in values)
