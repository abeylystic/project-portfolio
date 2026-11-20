"""Power and test calibration for a matched-pair cluster randomized trial."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import nct, t

ALPHA = 0.05
REPLICATES = 4000
WITHIN_SITE = 30
PERIOD_CORRELATION = 0.65
PAIR_SHARE = 0.30
PLANNING_EFFECT = 0.30
ICCS = (0.02, 0.15, 0.35)
PAIR_GRID = range(6, 65)


def cluster_residual_variance(icc: float, per_site: int) -> float:
    """Cluster mean variance after the shared pair component is removed."""
    return (1 - PAIR_SHARE) * icc + (1 - icc) / per_site


def paired_power(pairs: int, variance: float, effect: float) -> float:
    df = pairs - 1
    critical = t.ppf(1 - ALPHA / 2, df)
    noncentrality = effect / np.sqrt(variance / pairs)
    return float(nct.sf(critical, df, noncentrality) + nct.cdf(-critical, df, noncentrality))


def cluster_means(
    pairs: int, per_site: int, icc: float, size: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = (size, pairs, 2)
    site_sd = np.sqrt(cluster_residual_variance(icc, per_site))
    pair_sd = np.sqrt(PAIR_SHARE * icc)
    baseline = rng.normal(0, site_sd, shape) + rng.normal(0, pair_sd, (size, pairs, 1))
    innovation = rng.normal(0, site_sd, shape) + rng.normal(0, pair_sd, (size, pairs, 1))
    followup = (
        PERIOD_CORRELATION * baseline + np.sqrt(1 - PERIOD_CORRELATION**2) * innovation
    )
    assignment = rng.integers(0, 2, size=(size, pairs, 1))
    return baseline, followup, assignment


def paired_stats(differences: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    pairs = differences.shape[1]
    return (
        differences.mean(axis=1),
        differences.std(axis=1, ddof=1) / np.sqrt(pairs),
        float(t.ppf(1 - ALPHA / 2, pairs - 1)),
    )


def naive_stats(
    treated: np.ndarray, control: np.ndarray, per_site: int,
    icc: float, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Welch SE computed with independent-Gaussian within-site sums of squares."""
    simulations, pairs = treated.shape
    n = pairs * per_site
    arm_variances = []
    for means in (treated, control):
        within_ss = (1 - icc) * rng.chisquare(pairs * (per_site - 1), size=simulations)
        between_ss = per_site * np.square(means - means.mean(axis=1, keepdims=True)).sum(axis=1)
        arm_variances.append((within_ss + between_ss) / (n - 1))
    first, second = np.asarray(arm_variances) / n
    variance = first + second
    df = variance**2 / ((first**2 + second**2) / (n - 1))
    return treated.mean(axis=1) - control.mean(axis=1), np.sqrt(variance), t.ppf(0.975, df)


def assess(icc: float, effect: float, rng: np.random.Generator) -> list[dict]:
    pairs = 20
    baseline, followup, assignment = cluster_means(pairs, WITHIN_SITE, icc, REPLICATES, rng)
    treated = np.take_along_axis(followup, assignment, axis=2)[..., 0] + effect
    control = np.take_along_axis(followup, 1 - assignment, axis=2)[..., 0]
    baseline_diff = (
        np.take_along_axis(baseline, assignment, axis=2)
        - np.take_along_axis(baseline, 1 - assignment, axis=2)
    )[..., 0]
    followup_diff = treated - control
    methods = {
        "Naive participant Welch": naive_stats(treated, control, WITHIN_SITE, icc, rng),
        "Paired follow-up": paired_stats(followup_diff),
        "Paired baseline adjustment": paired_stats(followup_diff - PERIOD_CORRELATION * baseline_diff),
    }
    rows = []
    for model, (estimate, se, critical) in methods.items():
        half_width = critical * se
        reject = np.mean(np.abs(estimate) > half_width)
        rows.append({
            "icc": icc, "true_effect_sd": effect, "pairs": pairs, "participants_per_site": WITHIN_SITE,
            "replicates": REPLICATES, "method": model,
            "bias_sd": (estimate - effect).mean(),
            "empirical_sd": estimate.std(ddof=1), "mean_se": se.mean(),
            "rejection_rate": reject,
            "rejection_mcse": np.sqrt(reject * (1 - reject) / REPLICATES),
            "coverage": np.mean(np.abs(estimate - effect) <= half_width),
        })
    return rows


def power_grid() -> pd.DataFrame:
    rows = []
    for icc in ICCS:
        variance = 2 * cluster_residual_variance(icc, WITHIN_SITE)
        adjusted_variance = variance * (1 - PERIOD_CORRELATION**2)
        for count in PAIR_GRID:
            rows.append({
                "icc": icc, "pairs": count, "sites": 2 * count,
                "participants": 2 * count * WITHIN_SITE,
                "effect_sd": PLANNING_EFFECT,
                "power_followup": paired_power(count, variance, PLANNING_EFFECT),
                "power_adjusted": paired_power(count, adjusted_variance, PLANNING_EFFECT),
            })
    return pd.DataFrame(rows)


def plot_power(power: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True, layout="constrained")
    for (icc, values), axis in zip(power.groupby("icc"), axes, strict=True):
        axis.plot(values["pairs"], 100 * values["power_adjusted"], color="#126B68", label="Baseline-adjusted")
        axis.plot(values["pairs"], 100 * values["power_followup"], color="#9C3C32", label="Follow-up only")
        axis.axhline(80, color="#637083", linestyle="--", linewidth=0.9)
        axis.set(title=f"ICC {icc:.2f}", xlabel="Matched pairs", ylim=(0, 103))
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Power (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    fig.suptitle("Power to detect a 0.30 SD effect")
    fig.savefig(path, dpi=145)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--seed", type=int, default=20261120)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = [row for icc in ICCS for effect in (0.0, PLANNING_EFFECT) for row in assess(icc, effect, rng)]
    power = power_grid()
    args.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output / "simulation.csv", index=False, float_format="%.8f")
    power.to_csv(args.output / "power.csv", index=False, float_format="%.8f")
    plot_power(power, args.output / "power.png")


if __name__ == "__main__":
    main()
