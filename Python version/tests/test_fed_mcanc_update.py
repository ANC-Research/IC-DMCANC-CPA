from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.config import CommunicationMode, SimulationConfig
from ic_dmcanc.fed_mcanc import run_fed_mcanc


def _config(
    *,
    mode: CommunicationMode = CommunicationMode.COMMON_PERIODIC,
    tc: float | tuple[float, ...] = 0.2,
    c_len: int = 1,
    alpha: float = 0.25,
) -> SimulationConfig:
    return SimulationConfig(
        fs=10,
        duration_seconds=1.0,
        w_len=2,
        s_len=1,
        c_len=c_len,
        mu_w=0.1,
        mu_c=0.01,
        alpha=alpha,
        synthetic_low_hz=1.0,
        synthetic_high_hz=4.0,
        communication_mode=mode,
        tc_seconds=tc,
        comp_id_num_samples=4,
        comp_id_convergence_stride=1,
        diagnostic_stride=1,
        use_numba=False,
    )


class FedMCANCUpdateTests(unittest.TestCase):
    def test_residual_full_coupling_and_positive_local_update(self) -> None:
        config = _config()
        reference = np.asarray([2.0])
        disturbance = np.arange(1.0, 7.0).reshape(6, 1)
        secondary = np.zeros((6, 6, 1))
        secondary[:, :, 0] = np.arange(1.0, 37.0).reshape(6, 6) / 100.0
        compensation = np.zeros((6, 6, 1))
        initial = np.zeros((6, 2))
        initial[:, 0] = np.linspace(0.1, 0.6, 6)

        result = run_fed_mcanc(
            reference,
            disturbance,
            secondary,
            compensation,
            config,
            initial_center_weights=initial,
            save_secondary_output=True,
        )
        control = initial[:, 0] * 2.0
        expected_secondary = secondary[:, :, 0] @ control
        expected_error = disturbance[:, 0] - expected_secondary
        np.testing.assert_allclose(result.control_output[:, 0], control)
        np.testing.assert_allclose(
            result.secondary_output[:, 0], expected_secondary
        )
        np.testing.assert_allclose(
            result.residual_error[:, 0], expected_error
        )

        diagonal_filtered_reference = 2.0 * np.diag(
            secondary[:, :, 0]
        )
        expected_first_tap = (
            initial[:, 0]
            + config.mu_w
            * diagonal_filtered_reference
            * expected_error
        )
        np.testing.assert_allclose(
            result.final_local_weights[:, 0], expected_first_tap
        )
        # No center attraction on sample one because local==center beforehand.
        np.testing.assert_allclose(
            result.final_center_weights, initial
        )

    def test_ideal_aggregation_occurs_after_local_update(self) -> None:
        config = _config(
            mode=CommunicationMode.IDEAL,
            c_len=2,
            alpha=0.0,
        )
        reference = np.asarray([1.0, 0.5])
        disturbance = np.ones((6, 2))
        secondary = np.zeros((6, 6, 1))
        for node in range(6):
            secondary[node, node, 0] = 1.0
        compensation = np.zeros((6, 6, 2))
        result = run_fed_mcanc(
            reference,
            disturbance,
            secondary,
            compensation,
            config,
        )
        np.testing.assert_allclose(
            result.final_local_weights, result.final_center_weights
        )
        # Center only accepts the first w_len taps; the extended tail resets.
        np.testing.assert_array_equal(
            result.final_center_weights[:, config.w_len :], 0.0
        )
        np.testing.assert_array_equal(
            result.communication.sample_indices, [1, 2]
        )

    def test_two_independent_runs_do_not_chain_state(self) -> None:
        config = _config(
            mode=CommunicationMode.COMMON_PERIODIC,
            tc=0.1,
        )
        reference = np.asarray([0.2, -0.4, 0.7])
        disturbance = np.ones((6, 3))
        secondary = np.zeros((6, 6, 1))
        for node in range(6):
            secondary[node, node, 0] = 0.5
        compensation = np.zeros((6, 6, 1))
        first = run_fed_mcanc(
            reference, disturbance, secondary, compensation, config
        )
        second = run_fed_mcanc(
            reference, disturbance, secondary, compensation, config
        )
        np.testing.assert_array_equal(
            first.residual_error, second.residual_error
        )
        np.testing.assert_array_equal(
            first.final_local_weights, second.final_local_weights
        )


if __name__ == "__main__":
    unittest.main()
