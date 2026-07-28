from __future__ import annotations

import unittest

import numpy as np
from scipy.signal import freqz

from ic_dmcanc.signal_utils import (
    add_measured_awgn,
    design_fir1_bandpass,
    matlab_smooth,
)


class SignalUtilityTests(unittest.TestCase):
    def test_fir1_mapping_uses_order_plus_one_taps(self) -> None:
        taps = design_fir1_bandpass(63, 100.0, 1_000.0, 16_000)
        self.assertEqual(taps.shape, (64,))
        frequencies, response = freqz(taps, worN=16_384, fs=16_000)
        magnitude = np.abs(response)
        passband = magnitude[
            (frequencies >= 200.0) & (frequencies <= 900.0)
        ]
        upper_stop = magnitude[frequencies >= 2_000.0]
        self.assertGreater(float(np.max(passband)), 0.9)
        self.assertLess(float(np.max(upper_stop)), 0.05)

    def test_measured_awgn_hits_requested_finite_record_snr(self) -> None:
        signal = np.linspace(-0.8, 1.2, 10_001)
        noisy, noise = add_measured_awgn(signal, 40.0, seed=123)
        measured = 10.0 * np.log10(
            np.mean(signal * signal) / np.mean(noise * noise)
        )
        self.assertAlmostEqual(measured, 40.0, places=12)
        np.testing.assert_allclose(noisy, signal + noise)

    def test_measured_awgn_does_not_remove_mean(self) -> None:
        signal = np.full(1_000, 2.0)
        _, noise = add_measured_awgn(signal, 20.0, seed=7)
        expected_noise_power = np.mean(signal * signal) / 100.0
        self.assertAlmostEqual(
            float(np.mean(noise * noise)), expected_noise_power, places=13
        )

    def test_matlab_smooth_center_and_endpoints(self) -> None:
        values = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0])
        expected = np.asarray([1.5, 2.0, 3.0, 4.0, 4.5])
        np.testing.assert_allclose(matlab_smooth(values, 3), expected)
        # MATLAB reduces an even moving-average span by one.
        np.testing.assert_allclose(matlab_smooth(values, 4), expected)


if __name__ == "__main__":
    unittest.main()
