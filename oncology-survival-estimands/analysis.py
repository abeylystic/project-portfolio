"""Reproduce the VA lung cancer trial analysis and simulation results."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import restricted_mean_survival_time

from simulation import run_simulation
from survival import kaplan_meier, rmst_difference
from trial_data import DEFAULT_DATA, load_trial

HORIZONS = (90, 180, 270)
RISK_DAYS = (0, 45, 90, 135, 180, 225, 270)
ARMS = {1: ("Standard", "#163B60"), 2: ("Test", "#B84A21")}


def analyse_trial(trial: pd.DataFrame, output: Path) -> None:
    arms = {arm: trial.loc[trial["trt"] == arm] for arm in ARMS}
    km = {arm: kaplan_meier(data["time"], data["status"]) for arm, data in arms.items()}

    fig, (axis, risk_axis) = plt.subplots(
        2, 1, figsize=(8.4, 5.7), sharex=True, layout="constrained",
        gridspec_kw={"height_ratios": [6, 1]},
    )
    risk_axis.set_title("Number at risk immediately before each day", loc="left", fontsize=10)
    arm_rows = []
    for arm, (label, color) in ARMS.items():
        data = arms[arm]
        reference = KaplanMeierFitter(label=label).fit(data["time"], data["status"])
        reference.plot_survival_function(ax=axis, ci_show=True, ci_alpha=0.10, color=color)
        for horizon in HORIZONS:
            np.testing.assert_allclose(
                km[arm].rmst(horizon).days,
                restricted_mean_survival_time(reference, t=horizon),
                rtol=1e-12,
            )
        row = {
            "arm": label,
            "randomized": len(data),
            "deaths": int(data["status"].sum()),
            "censored": int(data["status"].eq(0).sum()),
            "median_survival_days": reference.median_survival_time_,
        }
        for day in RISK_DAYS:
            risk = int(data["time"].ge(day).sum())
            row[f"at_risk_day_{day}"] = risk
            risk_axis.text(day, 2 - arm, str(risk), ha="center", va="center", color=color)
        arm_rows.append(row)

    axis.axvline(180, linestyle="--", color="#616975", linewidth=1.2)
    axis.set(
        title="Overall survival in the VA lung cancer trial",
        ylabel="Kaplan–Meier survival probability",
        ylim=(0, 1.04),
        xlim=(-8, 282),
    )
    axis.grid(axis="y", color="#DEE4E9", linewidth=0.65)
    axis.legend(loc="upper right", frameon=False)
    axis.spines[["right", "top"]].set_visible(False)
    risk_axis.set(yticks=[1, 0], yticklabels=["Standard", "Test"], ylim=(-0.55, 1.55),
                  xticks=RISK_DAYS, xlabel="Days since randomization")
    risk_axis.tick_params(axis="y", length=0)
    for spine in risk_axis.spines.values():
        spine.set_visible(False)
    fig.savefig(output / "survival.png", dpi=150)
    plt.close(fig)

    pd.DataFrame(arm_rows).to_csv(output / "arm_summary.csv", index=False)
    rows = []
    for horizon in HORIZONS:
        difference, se, low, high = rmst_difference(km[1], km[2], horizon)
        rows.append({
            "horizon_days": horizon,
            "standard_rmst_days": km[1].rmst(horizon).days,
            "test_rmst_days": km[2].rmst(horizon).days,
            "test_minus_standard_days": difference,
            "standard_error_days": se,
            "ci95_low_days": low,
            "ci95_high_days": high,
        })
    pd.DataFrame(rows).to_csv(output / "trial_effects.csv", index=False, float_format="%.8f")

    cox_data = trial[["time", "status"]].assign(test_arm=trial["trt"].eq(2).astype(int))
    cox = CoxPHFitter()
    cox.fit(cox_data, "time", "status")
    summary = cox.summary.loc["test_arm"]
    ph = proportional_hazard_test(cox, cox_data, time_transform="rank").summary.loc["test_arm"]
    pd.DataFrame([{
        "hazard_ratio": summary["exp(coef)"],
        "ci95_low": summary["exp(coef) lower 95%"],
        "ci95_high": summary["exp(coef) upper 95%"],
        "ph_rank_test_statistic": ph["test_statistic"],
        "ph_rank_test_p_value": ph["p"],
    }]).to_csv(output / "cox_model.csv", index=False, float_format="%.8f")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--replicates", type=int, default=3000)
    parser.add_argument("--per-arm", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20261120)
    args = parser.parse_args()
    trial = load_trial(args.data)
    simulation = run_simulation(args.replicates, args.per_arm, args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    analyse_trial(trial, args.output)
    simulation.to_csv(args.output / "simulation.csv", index=False, float_format="%.8f")


if __name__ == "__main__":
    main()
