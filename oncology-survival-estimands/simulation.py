"""Monte Carlo calibration against analytical restricted survival means."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from survival import UnsupportedHorizon, kaplan_meier, rmst_difference

SWITCH_DAY = 60
BASE_HAZARD = 0.008
HORIZON = 180


@dataclass(frozen=True)
class Scenario:
    name: str
    early_hr: float
    late_hr: float
    informative_censoring: bool = False


SCENARIOS = (
    Scenario("Null, independent censoring", 1.0, 1.0),
    Scenario("Proportional benefit", 0.70, 0.70),
    Scenario("Early harm, late benefit", 1.65, 0.42),
    Scenario("Null, informative censoring", 1.0, 1.0, informative_censoring=True),
)


def analytic_rmst(early_hazard: float, late_hazard: float, horizon: int) -> float:
    first = min(horizon, SWITCH_DAY)
    second = max(horizon - SWITCH_DAY, 0)
    first_area = -np.expm1(-early_hazard * first) / early_hazard
    second_area = -np.expm1(-late_hazard * second) / late_hazard
    return float(first_area + np.exp(-early_hazard * first) * second_area)


def true_difference(scenario: Scenario) -> float:
    # Both arms have an identical latent-risk mixture under the informative null.
    return analytic_rmst(
        BASE_HAZARD * scenario.early_hr, BASE_HAZARD * scenario.late_hr, HORIZON
    ) - analytic_rmst(BASE_HAZARD, BASE_HAZARD, HORIZON)


def draw_arm(
    rng: np.random.Generator, scenario: Scenario, size: int, test: bool
) -> tuple[np.ndarray, np.ndarray]:
    early_hr, late_hr = (scenario.early_hr, scenario.late_hr) if test else (1.0, 1.0)
    early_hazard = np.full(size, early_hr * BASE_HAZARD)
    late_hazard = np.full(size, late_hr * BASE_HAZARD)
    censor_mean = np.full(size, 360.0)
    if scenario.informative_censoring:
        high_risk = rng.random(size) < 0.45
        early_hazard[high_risk] *= 2.5
        late_hazard[high_risk] *= 2.5
        if test:
            censor_mean[high_risk] = 30.0

    cumulative_hazard = rng.exponential(size=size)
    hazard_at_switch = early_hazard * SWITCH_DAY
    loss = rng.exponential(scale=censor_mean)
    failure = np.where(
        cumulative_hazard <= hazard_at_switch,
        cumulative_hazard / early_hazard,
        SWITCH_DAY + (cumulative_hazard - hazard_at_switch) / late_hazard,
    )
    return np.minimum(failure, loss), failure <= loss


def run_simulation(replicates: int, per_arm: int, seed: int) -> pd.DataFrame:
    if replicates < 2 or per_arm < 2:
        raise ValueError("Use at least two replicates and two patients per arm.")
    rows = []
    for scenario, sequence in zip(SCENARIOS, np.random.SeedSequence(seed).spawn(len(SCENARIOS))):
        rng = np.random.default_rng(sequence)
        truth = true_difference(scenario)
        estimates = []
        covered = []
        for _ in range(replicates):
            standard = kaplan_meier(*draw_arm(rng, scenario, per_arm, test=False))
            test = kaplan_meier(*draw_arm(rng, scenario, per_arm, test=True))
            try:
                estimate, _, low, high = rmst_difference(standard, test, HORIZON)
            except UnsupportedHorizon:
                continue
            estimates.append(estimate)
            covered.append(low <= truth <= high)
        coverage = np.mean(covered) if covered else np.nan
        estimate_array = np.asarray(estimates, dtype=float)
        rows.append({
            "scenario": scenario.name,
            "replicates": replicates,
            "estimable": len(estimates),
            "true_difference_days": truth,
            "mean_difference_days": estimate_array.mean() if estimates else np.nan,
            "bias_days": (estimate_array - truth).mean() if estimates else np.nan,
            "rmse_days": np.sqrt(np.mean((estimate_array - truth)**2)) if estimates else np.nan,
            "coverage": coverage,
            "coverage_mcse": np.sqrt(coverage * (1 - coverage) / len(estimates)) if covered else np.nan,
            "mean_mcse_days": (
                estimate_array.std(ddof=1) / np.sqrt(len(estimates)) if len(estimates) > 1 else np.nan
            ),
        })
    return pd.DataFrame.from_records(rows)
