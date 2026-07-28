from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.centralized_fxlms import run_centralized_fxlms
from ic_dmcanc.config import SimulationConfig
from ic_dmcanc.mgdfxlms import run_mgdfxlms


class BaselineSmokeTests(unittest.TestCase):
    def test_centralized_and_mgd_short_runs_are_finite(self) -> None:
        rng = np.random.default_rng(8)
        reference = rng.standard_normal(12) * 0.1
        disturbance = rng.standard_normal((6, 12)) * 0.1
        secondary = rng.standard_normal((6, 6, 2)) * 0.01
        for node in range(6):
            secondary[node, node, 0] += 0.4
        compensation = rng.standard_normal((6, 6, 2)) * 0.01
        for node in range(6):
            compensation[node, node] = 0.0
        config = SimulationConfig(
            fs=10,
            duration_seconds=1.0,
            w_len=3,
            s_len=2,
            c_len=2,
            mu_w=0.01,
            mu_c=0.01,
            synthetic_low_hz=1.0,
            synthetic_high_hz=4.0,
            comp_id_num_samples=4,
            diagnostic_stride=2,
            use_numba=False,
        )
        centralized = run_centralized_fxlms(
            reference, disturbance, secondary, config
        )
        mgd = run_mgdfxlms(
            reference,
            disturbance,
            secondary,
            compensation,
            config,
        )
        for result in (centralized, mgd):
            self.assertEqual(result.residual_error.shape, (6, 12))
            self.assertEqual(result.control_output.shape, (6, 12))
            self.assertTrue(np.isfinite(result.residual_error).all())
            self.assertTrue(np.isfinite(result.final_weights).all())


if __name__ == "__main__":
    unittest.main()
