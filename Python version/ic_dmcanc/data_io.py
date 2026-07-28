"""Strict MATLAB asset loading with explicit variable and shape checks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from .config import DTYPE, NUM_NODES


@dataclass(frozen=True)
class AcousticPaths:
    primary_path: np.ndarray
    secondary_path: np.ndarray
    primary_variable: str = "Primary_path"
    secondary_variable: str = "Secondary_path"


@dataclass(frozen=True)
class CompressorRecording:
    time_seconds: np.ndarray
    signal: np.ndarray
    time_variable: str = "tr"
    signal_variable: str = "yr"


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_required_variable(path: Path, variable: str) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"Required MATLAB asset not found: {path}")
    contents = loadmat(path, squeeze_me=False, struct_as_record=False)
    public_names = sorted(name for name in contents if not name.startswith("__"))
    if variable not in contents:
        raise KeyError(
            f"{path} does not contain {variable!r}; found {public_names}."
        )
    array = np.asarray(contents[variable], dtype=DTYPE)
    if not np.isfinite(array).all():
        raise ValueError(f"{path}:{variable} contains NaN or Inf.")
    return np.ascontiguousarray(array)


def load_acoustic_paths(
    root: Path | None = None,
    *,
    expected_s_len: int = 256,
) -> AcousticPaths:
    """Load the committed 1×6 primary and 6×6 secondary paths."""

    base = repository_root() if root is None else Path(root)
    primary = _load_required_variable(
        base / "simulation path" / "PrimaryPath_1x6.mat", "Primary_path"
    )
    secondary = _load_required_variable(
        base / "simulation path" / "SecondaryPath_6x6.mat", "Secondary_path"
    )
    if primary.ndim != 2 or primary.shape[0] != NUM_NODES:
        raise ValueError(
            "Primary_path must have shape (6, p_len); "
            f"got {primary.shape}."
        )
    expected_secondary = (NUM_NODES, NUM_NODES, expected_s_len)
    if secondary.shape != expected_secondary:
        raise ValueError(
            f"Secondary_path must have shape {expected_secondary}; "
            f"got {secondary.shape}."
        )
    return AcousticPaths(primary_path=primary, secondary_path=secondary)


def load_compressor_recording(
    root: Path | None = None,
) -> CompressorRecording:
    """Load ``tr`` and ``yr`` without unconditional ``squeeze``."""

    base = repository_root() if root is None else Path(root)
    path = base / "compressor_16kHz.mat"
    time = _load_required_variable(path, "tr")
    signal = _load_required_variable(path, "yr")
    for name, array in (("tr", time), ("yr", signal)):
        if array.ndim != 2 or array.shape[1] != 1:
            raise ValueError(
                f"{name} must be an explicit MATLAB column vector; "
                f"got {array.shape}."
            )
    if time.shape != signal.shape:
        raise ValueError(
            f"tr and yr shapes differ: {time.shape} vs {signal.shape}."
        )
    return CompressorRecording(
        time_seconds=np.ascontiguousarray(time[:, 0]),
        signal=np.ascontiguousarray(signal[:, 0]),
    )


def load_vector_from_mat(path: Path, variable: str) -> np.ndarray:
    """Load a fixed cross-language test vector as one contiguous row."""

    array = _load_required_variable(Path(path), variable)
    if array.ndim != 2 or 1 not in array.shape:
        raise ValueError(
            f"{path}:{variable} must be a MATLAB row or column vector; "
            f"got {array.shape}."
        )
    return np.ascontiguousarray(array.reshape(-1, order="F"), dtype=DTYPE)
