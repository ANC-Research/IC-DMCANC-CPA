from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.compensation_filter import (
    identify_compensation_filters,
    identify_single_pair,
    mwd_tail_direct,
    mwd_tail_lfilter,
)


class CompensationFilterTests(unittest.TestCase):
    def test_single_sample_system_id_update_order(self) -> None:
        raw, mse, indices, convergence, error, engine = identify_single_pair(
            np.asarray([2.0]),
            np.asarray([1.0]),
            np.asarray([0.5]),
            c_len=2,
            step_size=0.1,
            convergence_stride=1,
            save_full_error=True,
            use_numba=False,
        )
        np.testing.assert_allclose(raw, [0.1, 0.0])
        self.assertAlmostEqual(mse, 1.0)
        np.testing.assert_array_equal(indices, [1])
        np.testing.assert_allclose(convergence, [1.0])
        np.testing.assert_allclose(error, [1.0])
        self.assertEqual(engine, "numpy_reference")

    def test_all_pairs_share_noise_skip_diagonal_and_negative_flip(self) -> None:
        secondary = np.full((6, 6, 1), 0.15, dtype=np.float64)
        for node in range(6):
            secondary[node, node, 0] = 1.0
        result = identify_compensation_filters(
            secondary,
            c_len=2,
            step_size=0.02,
            num_samples=500,
            seed=42,
            convergence_stride=50,
            use_numba=False,
        )
        self.assertEqual(result.white_noise.shape, (500,))
        np.testing.assert_allclose(
            result.compensation_filters,
            -result.identified_coefficients[:, :, ::-1],
        )
        for node in range(6):
            np.testing.assert_array_equal(
                result.compensation_filters[node, node], 0.0
            )
            self.assertEqual(result.final_mse[node, node], 0.0)
        self.assertTrue(np.isfinite(result.compensation_filters).all())

    def test_mwd_direct_is_exact_lfilter_tail(self) -> None:
        rng = np.random.default_rng(4)
        for w_len, c_len in ((1, 1), (4, 3), (17, 5), (32, 9)):
            compensation = rng.standard_normal(c_len)
            difference = rng.standard_normal(w_len + c_len - 1)
            expected = mwd_tail_lfilter(
                compensation, difference, w_len
            )
            actual = mwd_tail_direct(compensation, difference, w_len)
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-15)

    def test_mwd_rejects_same_or_center_cropping_lengths(self) -> None:
        with self.assertRaises(ValueError):
            mwd_tail_lfilter(
                np.ones(3), np.ones(4), w_len=4
            )


if __name__ == "__main__":
    unittest.main()
