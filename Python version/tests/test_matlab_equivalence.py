from __future__ import annotations

import os
from pathlib import Path
import unittest

import numpy as np
from scipy.io import loadmat
from scipy.signal import lfilter

from ic_dmcanc.config import CommunicationMode, SimulationConfig
from ic_dmcanc.fed_mcanc import run_fed_mcanc


def _matlab_style_reference(
    reference: np.ndarray,
    disturbance: np.ndarray,
    secondary: np.ndarray,
    compensation: np.ndarray,
    config: SimulationConfig,
    initial_center: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Independent literal translation using concatenation and lfilter."""

    w_len = config.w_len
    c_len = config.c_len
    extended_len = w_len + c_len - 1
    num_samples = reference.size
    local = initial_center.copy()
    center = initial_center.copy()
    nabla = np.zeros_like(local)
    xc = np.zeros(w_len)
    xs = np.zeros(config.s_len)
    ys = np.zeros((6, config.s_len))
    xf = np.zeros((6, extended_len))
    residual = np.zeros((6, num_samples))
    control = np.zeros((6, num_samples))
    intervals = config.communication_intervals()

    for sample in range(num_samples):
        xc = np.r_[reference[sample], xc[:-1]]
        control[:, sample] = np.sum(local[:, :w_len] * xc, axis=1)
        for source in range(6):
            ys[source] = np.r_[control[source, sample], ys[source, :-1]]
        for error_node in range(6):
            secondary_output = 0.0
            for source in range(6):
                secondary_output += np.sum(
                    secondary[error_node, source] * ys[source]
                )
            residual[error_node, sample] = (
                disturbance[error_node, sample] - secondary_output
            )

        xs = np.r_[reference[sample], xs[:-1]]
        for node in range(6):
            filtered_sample = np.sum(secondary[node, node] * xs)
            xf[node] = np.r_[filtered_sample, xf[node, :-1]]
            local[node] = (
                local[node]
                + config.mu_w * xf[node] * residual[node, sample]
                + config.alpha * config.mu_w * (center[node] - local[node])
            )

        matlab_index = sample + 1
        common = config.communication_mode is CommunicationMode.IDEAL or (
            config.communication_mode
            in {
                CommunicationMode.COMMON_PERIODIC,
                CommunicationMode.FIXED_10_SECONDS,
            }
            and matlab_index % intervals[0] == 0
        )
        if common:
            nabla = local - center
            updated = center.copy()
            for target in range(6):
                delta = nabla[target, :w_len].copy()
                for other in range(6):
                    if other == target:
                        continue
                    filtered = lfilter(
                        compensation[other, target], [1.0], nabla[other]
                    )
                    delta += filtered[-w_len:]
                updated[target, :w_len] += delta
            center = updated
            local = center.copy()
        elif config.communication_mode is CommunicationMode.INDIVIDUAL_PERIODIC:
            for target in range(6):
                if matlab_index % intervals[target] != 0:
                    continue
                nabla[target] = local[target] - center[target]
                delta = nabla[target, :w_len].copy()
                for other in range(6):
                    if other == target:
                        continue
                    filtered = lfilter(
                        compensation[other, target], [1.0], nabla[other]
                    )
                    delta += filtered[-w_len:]
                center[target, :w_len] += delta
                local[target] = center[target]
    return residual, control, local, center, nabla


class MatlabStyleEquivalenceTests(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(12)
        self.reference = rng.standard_normal(8) * 0.2
        self.disturbance = rng.standard_normal((6, 8)) * 0.1
        self.secondary = rng.standard_normal((6, 6, 3)) * 0.03
        for node in range(6):
            self.secondary[node, node, 0] += 0.3
        self.compensation = rng.standard_normal((6, 6, 2)) * 0.02
        for node in range(6):
            self.compensation[node, node] = 0.0
        self.initial = rng.standard_normal((6, 4)) * 0.01

    def _compare(
        self,
        mode: CommunicationMode,
        tc: float | tuple[float, ...],
    ) -> None:
        config = SimulationConfig(
            fs=10,
            duration_seconds=1.0,
            w_len=3,
            s_len=3,
            c_len=2,
            mu_w=0.01,
            mu_c=0.01,
            alpha=0.2,
            synthetic_low_hz=1.0,
            synthetic_high_hz=4.0,
            communication_mode=mode,
            tc_seconds=tc,
            comp_id_num_samples=8,
            comp_id_convergence_stride=1,
            diagnostic_stride=1,
            use_numba=False,
        )
        expected = _matlab_style_reference(
            self.reference,
            self.disturbance,
            self.secondary,
            self.compensation,
            config,
            self.initial,
        )
        actual = run_fed_mcanc(
            self.reference,
            self.disturbance,
            self.secondary,
            self.compensation,
            config,
            initial_center_weights=self.initial,
        )
        np.testing.assert_allclose(
            actual.residual_error, expected[0], rtol=1e-13, atol=1e-14
        )
        np.testing.assert_allclose(
            actual.control_output, expected[1], rtol=1e-13, atol=1e-14
        )
        np.testing.assert_allclose(
            actual.final_local_weights, expected[2], rtol=1e-13, atol=1e-14
        )
        np.testing.assert_allclose(
            actual.final_center_weights, expected[3], rtol=1e-13, atol=1e-14
        )
        np.testing.assert_allclose(
            actual.final_weight_difference,
            expected[4],
            rtol=1e-13,
            atol=1e-14,
        )

    def test_ideal_trace(self) -> None:
        self._compare(CommunicationMode.IDEAL, 0.2)

    def test_common_periodic_trace(self) -> None:
        self._compare(CommunicationMode.COMMON_PERIODIC, 0.2)

    def test_individual_periodic_trace(self) -> None:
        self._compare(
            CommunicationMode.INDIVIDUAL_PERIODIC,
            (0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
        )


class OptionalExportedMatlabFixtureTests(unittest.TestCase):
    def test_exported_ideal_fixture_when_supplied(self) -> None:
        fixture_name = os.environ.get("IC_DMCANC_MATLAB_FIXTURE")
        if not fixture_name:
            self.skipTest(
                "Set IC_DMCANC_MATLAB_FIXTURE to a fixture produced by "
                "matlab/export_equivalence_fixture.m."
            )
        fixture = loadmat(Path(fixture_name), squeeze_me=False)
        reference = fixture["reference"].reshape(-1, order="F")
        disturbance = fixture["disturbance"]
        secondary = fixture["secondary_path"]
        compensation = fixture["compensation_filters"]
        initial = fixture["initial_center"]
        config = SimulationConfig(
            fs=16_000,
            duration_seconds=1.0,
            w_len=int(fixture["w_len"].item()),
            s_len=int(fixture["s_len"].item()),
            c_len=int(fixture["c_len"].item()),
            mu_w=float(fixture["mu_w"].item()),
            mu_c=0.0,
            alpha=float(fixture["alpha"].item()),
            communication_mode=CommunicationMode.IDEAL,
            comp_id_num_samples=1,
            diagnostic_stride=1,
            use_numba=False,
        )
        actual = run_fed_mcanc(
            reference,
            disturbance,
            secondary,
            compensation,
            config,
            initial_center_weights=initial,
        )
        np.testing.assert_allclose(
            actual.residual_error,
            fixture["residual_error"],
            rtol=1e-11,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            actual.control_output,
            fixture["control_output"],
            rtol=1e-11,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            actual.final_local_weights,
            fixture["final_local_weights"],
            rtol=1e-11,
            atol=1e-12,
        )


if __name__ == "__main__":
    unittest.main()
