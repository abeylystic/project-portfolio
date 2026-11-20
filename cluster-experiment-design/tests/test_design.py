import unittest

import numpy as np

from design import cluster_means, cluster_residual_variance, paired_power, paired_stats


class DesignTests(unittest.TestCase):
    def test_noncentral_t_under_null_has_nominal_size(self):
        self.assertAlmostEqual(paired_power(20, 0.2, 0), 0.05)

    def test_variance_increases_when_icc_increases_at_fixed_cluster_size(self):
        self.assertGreater(cluster_residual_variance(0.35, 30), cluster_residual_variance(0.02, 30))

    def test_pair_stats_uses_variance_between_pair_differences(self):
        estimate, se, critical = paired_stats(np.array([[0, 2, 4]]))
        self.assertAlmostEqual(estimate.item(), 2)
        self.assertAlmostEqual(se.item(), 2 / np.sqrt(3))
        self.assertGreater(critical, 4)

    def test_generate_one_treated_site_per_pair(self):
        baseline, followup, assignment = cluster_means(8, 30, 0.15, 100, np.random.default_rng(1))
        self.assertEqual(baseline.shape, (100, 8, 2))
        self.assertEqual(followup.shape, baseline.shape)
        self.assertEqual(assignment.shape, (100, 8, 1))
        self.assertEqual(set(assignment.ravel()), {0, 1})


if __name__ == "__main__":
    unittest.main()
