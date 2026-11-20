import unittest

import numpy as np
import pandas as pd

from backtest import coverage_tests, forecasts


class BacktestTests(unittest.TestCase):
    def test_forecasts_are_unchanged_when_current_and_future_losses_change(self):
        random = np.random.default_rng(14)
        original = pd.Series(random.normal(size=180))
        perturbed = original.copy()
        perturbed.iloc[120:] *= 100
        prior = forecasts(original, window=35).drop(columns="loss")
        modified = forecasts(perturbed, window=35).drop(columns="loss")
        pd.testing.assert_frame_equal(prior.loc[:120], modified.loc[:120])
        self.assertFalse(prior.loc[121:].equals(modified.loc[121:]))

    def test_first_historical_forecast_uses_only_previous_window(self):
        losses = pd.Series(np.arange(60))
        self.assertAlmostEqual(forecasts(losses, window=5).loc[5, "historical"], 3.96)

    def test_matching_exception_frequency_has_perfect_coverage_score(self):
        hits = np.zeros(500, dtype=bool)
        hits[::100] = True
        p_coverage, _ = coverage_tests(hits)
        self.assertAlmostEqual(p_coverage, 1)

    def test_independence_flags_same_count_when_exceptions_cluster(self):
        dispersed = np.zeros(1000, dtype=bool)
        clustered = dispersed.copy()
        dispersed[::100] = True
        clustered[5:15] = True
        _, p_dispersed = coverage_tests(dispersed)
        _, p_clustered = coverage_tests(clustered)
        self.assertGreater(p_dispersed, 0.05)
        self.assertLess(p_clustered, 0.001)


if __name__ == "__main__":
    unittest.main()
