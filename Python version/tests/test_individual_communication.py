from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.fed_mcanc import apply_individual_communication


class IndividualCommunicationTests(unittest.TestCase):
    def test_only_target_resets_and_other_differences_stay_stale(self) -> None:
        w_len = 2
        c_len = 2
        shape = (6, w_len + c_len - 1)
        local = np.arange(np.prod(shape), dtype=float).reshape(shape)
        center = np.zeros(shape)
        stale = np.full(shape, -7.0)
        compensation = np.zeros((6, 6, c_len))
        updated_local, updated_center, updated_stale = (
            apply_individual_communication(
                local,
                center,
                stale,
                compensation,
                target_node=3,
                w_len=w_len,
            )
        )
        np.testing.assert_array_equal(
            updated_stale[3], local[3] - center[3]
        )
        np.testing.assert_array_equal(
            updated_stale[[0, 1, 2, 4, 5]],
            stale[[0, 1, 2, 4, 5]],
        )
        np.testing.assert_array_equal(
            updated_local[[0, 1, 2, 4, 5]],
            local[[0, 1, 2, 4, 5]],
        )
        np.testing.assert_array_equal(
            updated_center[[0, 1, 2, 4, 5]],
            center[[0, 1, 2, 4, 5]],
        )
        np.testing.assert_array_equal(
            updated_local[3], updated_center[3]
        )

    def test_same_sample_sequence_uses_earlier_refreshed_difference(self) -> None:
        w_len = 1
        c_len = 2
        local = np.zeros((6, 2))
        center = np.zeros((6, 2))
        stale = np.zeros((6, 2))
        compensation = np.zeros((6, 6, c_len))
        local[0] = [2.0, 3.0]
        local[1] = [1.0, 0.0]
        compensation[0, 1] = [1.0, 0.0]

        local, center, stale = apply_individual_communication(
            local,
            center,
            stale,
            compensation,
            target_node=0,
            w_len=w_len,
        )
        np.testing.assert_array_equal(stale[0], [2.0, 3.0])

        local, center, stale = apply_individual_communication(
            local,
            center,
            stale,
            compensation,
            target_node=1,
            w_len=w_len,
        )
        # Own difference 1 plus lfilter([1,0], Nabla0)[-1] == 3.
        self.assertAlmostEqual(center[1, 0], 4.0)


if __name__ == "__main__":
    unittest.main()
