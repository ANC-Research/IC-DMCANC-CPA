from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.config import CommunicationMode, SimulationConfig
from ic_dmcanc.fed_mcanc import NUMBA_AVAILABLE, run_fed_mcanc


@unittest.skipUnless(NUMBA_AVAILABLE, "Numba optional dependency is absent.")
class NumbaEquivalenceTests(unittest.TestCase):
    def test_numba_and_reference_kernels_match(self) -> None:
        rng = np.random.default_rng(23)
        reference = rng.standard_normal(10) * 0.1
        disturbance = rng.standard_normal((6, 10)) * 0.1
        secondary = rng.standard_normal((6, 6, 2)) * 0.02
        compensation = rng.standard_normal((6, 6, 2)) * 0.01
        for node in range(6):
            compensation[node, node] = 0.0
        common = dict(
            fs=10,
            duration_seconds=1.0,
            w_len=3,
            s_len=2,
            c_len=2,
            mu_w=0.01,
            mu_c=0.01,
            alpha=0.3,
            synthetic_low_hz=1.0,
            synthetic_high_hz=4.0,
            communication_mode=CommunicationMode.COMMON_PERIODIC,
            tc_seconds=0.2,
            comp_id_num_samples=4,
            diagnostic_stride=1,
        )
        reference_result = run_fed_mcanc(
            reference,
            disturbance,
            secondary,
            compensation,
            SimulationConfig(**common, use_numba=False),
        )
        numba_result = run_fed_mcanc(
            reference,
            disturbance,
            secondary,
            compensation,
            SimulationConfig(**common, use_numba=True),
        )
        np.testing.assert_allclose(
            numba_result.residual_error,
            reference_result.residual_error,
            rtol=0,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            numba_result.final_local_weights,
            reference_result.final_local_weights,
            rtol=0,
            atol=1e-14,
        )


if __name__ == "__main__":
    unittest.main()
