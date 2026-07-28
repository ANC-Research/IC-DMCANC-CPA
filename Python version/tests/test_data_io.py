from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.data_io import (
    load_acoustic_paths,
    load_compressor_recording,
)


class DataIOTests(unittest.TestCase):
    def test_committed_mat_shapes_and_variables(self) -> None:
        paths = load_acoustic_paths()
        self.assertEqual(paths.primary_path.shape, (6, 512))
        self.assertEqual(paths.secondary_path.shape, (6, 6, 256))
        self.assertEqual(paths.primary_path.dtype, np.float64)
        self.assertEqual(paths.secondary_path.dtype, np.float64)
        self.assertTrue(np.isfinite(paths.primary_path).all())
        self.assertTrue(np.isfinite(paths.secondary_path).all())

    def test_compressor_columns_are_explicitly_validated(self) -> None:
        recording = load_compressor_recording()
        self.assertEqual(recording.time_seconds.shape, (279_378,))
        self.assertEqual(recording.signal.shape, (279_378,))
        self.assertTrue(np.isfinite(recording.signal).all())


if __name__ == "__main__":
    unittest.main()
