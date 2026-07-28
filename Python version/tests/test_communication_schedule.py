from __future__ import annotations

import unittest

import numpy as np

from ic_dmcanc.config import (
    CommunicationMode,
    SimulationConfig,
    require_integer_samples,
)
from ic_dmcanc.fed_mcanc import build_communication_schedule


class CommunicationScheduleTests(unittest.TestCase):
    def test_common_schedule_matches_matlab_one_based_mod(self) -> None:
        config = SimulationConfig(
            communication_mode=CommunicationMode.COMMON_PERIODIC,
            tc_seconds=0.1,
        )
        intervals = config.communication_intervals()
        np.testing.assert_array_equal(intervals, [1_600])
        schedule = build_communication_schedule(
            5_000, config.communication_mode, intervals
        )
        np.testing.assert_array_equal(
            schedule.sample_indices, [1_600, 3_200, 4_800]
        )
        self.assertNotIn(0, schedule.sample_indices)
        np.testing.assert_array_equal(schedule.per_node_counts, [3] * 6)
        self.assertEqual(schedule.total_node_uploads, 18)

    def test_fixed_interval_replaces_only_the_hard_coded_fs(self) -> None:
        config = SimulationConfig(
            communication_mode=CommunicationMode.FIXED_10_SECONDS
        )
        np.testing.assert_array_equal(
            config.communication_intervals(), [160_000]
        )

    def test_individual_simultaneous_events_are_node_ordered(self) -> None:
        intervals = np.asarray([2, 3, 4, 5, 6, 7], dtype=np.int64)
        schedule = build_communication_schedule(
            12, CommunicationMode.INDIVIDUAL_PERIODIC, intervals
        )
        at_six = schedule.node_indices[schedule.sample_indices == 6]
        np.testing.assert_array_equal(at_six, [0, 1, 4])
        np.testing.assert_array_equal(
            schedule.per_node_counts, [6, 4, 3, 2, 2, 1]
        )
        self.assertEqual(schedule.total_node_uploads, 18)

    def test_noninteger_interval_is_rejected_not_rounded(self) -> None:
        with self.assertRaisesRegex(ValueError, "require_integer"):
            require_integer_samples(16_000, 0.0001, name="tc")
        with self.assertRaises(ValueError):
            require_integer_samples(16_000, 0.0, name="tc")


if __name__ == "__main__":
    unittest.main()
