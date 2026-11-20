import unittest

import numpy as np

from simulation import SCENARIOS, analytic_rmst, run_simulation, true_difference
from survival import kaplan_meier, rmst_difference


class KaplanMeierTests(unittest.TestCase):
    def test_uncensored_rmst_matches_mean_truncated_event_time(self):
        time = np.array([1, 2, 3, 9])
        clipped = np.minimum(time, 5)
        fit = kaplan_meier(time, np.ones(len(time)))
        estimate = fit.rmst(5)
        self.assertAlmostEqual(estimate.days, clipped.mean())
        self.assertAlmostEqual(estimate.variance, clipped.var() / len(clipped))

    def test_tied_failure_and_censoring_use_full_risk_set(self):
        fit = kaplan_meier(np.array([1, 2, 2, 4]), np.array([1, 1, 0, 1]))
        np.testing.assert_allclose(fit.risk, [4, 3, 1])
        np.testing.assert_allclose(fit.survival, [0.75, 0.50, 0])
        self.assertAlmostEqual(fit.rmst(5).days, 2.75)
        self.assertAlmostEqual(fit.rmst(5).variance, 1.75**2 / 12 + 1 / 6)

    def test_event_at_horizon_does_not_change_area(self):
        fit = kaplan_meier(np.array([2, 3, 6]), np.array([1, 1, 0]))
        self.assertAlmostEqual(fit.rmst(3).days, 2 + 2 / 3)

    def test_support_ends_at_censored_tail(self):
        fit = kaplan_meier(np.array([2, 4]), np.array([1, 0]))
        with self.assertRaisesRegex(ValueError, "beyond observed"):
            fit.rmst(5)

    def test_arm_contrast_has_signed_difference_and_pooled_variance(self):
        standard = kaplan_meier(np.array([1, 3]), np.array([1, 1]))
        test = kaplan_meier(np.array([2, 4]), np.array([1, 1]))
        difference, se, low, high = rmst_difference(standard, test, 4)
        self.assertAlmostEqual(difference, 1)
        self.assertAlmostEqual(se, 1)
        self.assertLess(low, difference)
        self.assertGreater(high, difference)


class SimulationTests(unittest.TestCase):
    def test_piecewise_integral_when_hazards_match(self):
        hazard = 0.01
        for horizon in (20, 60, 180):
            self.assertAlmostEqual(
                analytic_rmst(hazard, hazard, horizon), -np.expm1(-hazard * horizon) / hazard
            )

    def test_informative_censoring_does_not_change_event_time_truth(self):
        self.assertAlmostEqual(true_difference(SCENARIOS[-1]), 0.0)

    def test_monte_carlo_is_reproducible(self):
        first = run_simulation(replicates=12, per_arm=40, seed=42)
        second = run_simulation(replicates=12, per_arm=40, seed=42)
        self.assertTrue(first.equals(second))
        self.assertTrue(first["coverage"].between(0, 1).all())


if __name__ == "__main__":
    unittest.main()
