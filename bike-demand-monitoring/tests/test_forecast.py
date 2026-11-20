import unittest

import numpy as np
import pandas as pd

from forecast import baseline_predictions, make_features, validate_counts


class ForecastTests(unittest.TestCase):
    def test_perturbing_the_current_and_future_targets_does_not_change_features(self):
        timestamps = pd.date_range("2012-01-01", periods=200, freq="h")
        history = pd.Series(np.arange(len(timestamps)), index=timestamps)
        perturbed = history.copy()
        perturbed.iloc[190:] = 1_000_000
        target = timestamps[[190]]
        pd.testing.assert_frame_equal(make_features(history, target), make_features(perturbed, target))

    def test_missing_lag_is_not_filled_from_a_nonadjacent_time(self):
        timestamps = pd.to_datetime(["2012-01-01 08:00", "2012-01-01 10:00"])
        features = make_features(pd.Series([8, 20], index=timestamps), timestamps[1:])
        self.assertTrue(pd.isna(features["lag_1h"].item()))
        self.assertEqual(features["lag_2h"].item(), 8)

    def test_new_weekday_hour_defaults_to_training_median_for_baseline(self):
        training = pd.Series([2, 8], index=pd.date_range("2012-01-01", periods=2, freq="h"))
        target = pd.DatetimeIndex(["2012-07-01 12:00"])
        features = make_features(training, target)
        self.assertEqual(baseline_predictions(features, training)["seasonal_naive"].item(), 5)

    def test_duplicate_hour_is_rejected(self):
        history = pd.Series([3, 5], index=pd.DatetimeIndex(["2012-01-01", "2012-01-01"]))
        with self.assertRaisesRegex(ValueError, "unique"):
            validate_counts(history)


if __name__ == "__main__":
    unittest.main()
