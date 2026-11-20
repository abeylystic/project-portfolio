# Analysis design

## Clinical question

How many additional days alive within 180 days are associated with assignment to the test regimen compared with assignment to the standard regimen in the historical VA lung cancer trial?

- **Population** includes all 137 archived randomized participants.
- **Treatment comparison** uses assigned test minus assigned standard. The archive supports an assignment-based comparison. It has limited documentation of treatment switching and other intercurrent events.
- **Endpoint** is all-cause mortality. Observed follow-up is the time from randomization used in this analysis.
- **Summary measure** is the between-arm difference in 180-day restricted mean survival time (RMST). Positive values favor test. The demonstrated analysis treats 180 days as its primary horizon. This is a retrospective portfolio analysis of the archived data, with no claim of historical trial prespecification.
- **Sensitivity horizons** are 90 and 270 days. Their intervals are descriptive, with no multiplicity adjustment. The number at risk is reported because declining follow-up support affects the stability of later estimates.

Arm-level RMST measures the area under the overall survival curve from randomization through the chosen horizon. Under independent right censoring within each arm, Kaplan–Meier can estimate that area. That assumption means future event times and censoring times are independent within arm, with sufficient follow-up support to estimate survival through the horizon.

## Inference

Let Ŝ denote the arm's Kaplan–Meier survival function and τ the horizon. Integrate the right-continuous step function exactly to obtain RMST. Failures precede tied censorings in the risk set calculation. A failure at τ changes Ŝ(τ) and contributes zero area.

The Greenwood variance is a sum over event days before τ. Each term multiplies the squared area from that event day to τ by d / [Y(Y − d)], where d is the number of deaths and Y is the risk set immediately before the event. A terminal event with Y = d leaves zero area and contributes zero. Add arm variances and take the square root to get the standard error of the independent-arm contrast. The interval uses a normal approximation.

A treatment-only Cox model with Efron handling of event-time ties provides a complementary hazard ratio over the available follow-up. A rank-time Schoenfeld-residual test assesses evidence of time dependence in the treatment coefficient. The test has limited power here. Under nonproportional hazards, the fitted coefficient can summarize a time-varying pattern with weights driven in part by the event and censoring processes. RMST remains a horizon-specific average in days and has its own censoring assumptions.

## Simulation

Each scenario generates 3,000 independent trials of 120 patients per assigned arm. Exponential event hazards can switch on day 60, and event times are sampled by transforming a unit exponential cumulative hazard. The standard event hazard is 0.008/day. The test arm uses the following hazard multipliers.

| Scenario | Days 0–60 | After day 60 | Censoring |
| --- | ---: | ---: | --- |
| Null | 1.00 | 1.00 | Independent exponential, mean 360 days |
| Proportional benefit | 0.70 | 0.70 | Independent exponential, mean 360 days |
| Early harm, late benefit | 1.65 | 0.42 | Independent exponential, mean 360 days |
| Informative censoring stress test | 1.00 | 1.00 | High-risk test patients have mean 30 days of censoring time |

In the stress test, each arm has the same population mixture, with 45% high-risk patients whose event hazard is 2.5 times baseline. Censoring in the standard arm and among lower-risk test patients has a mean of 360 days. This makes censoring dependent on prognosis within the test arm even though the event-time distributions under both assignments are identical.

The event-time truth comes from analytical integration of the piecewise exponential survival function. The informative-censoring scenario has a true mean difference of zero because the event-time mixtures are identical. The simulation reports bias, root mean squared error, interval coverage, and Monte Carlo standard error. Coverage uses estimable replicates as its denominator, and the output reports the count. A replicate is unestimable when the horizon exceeds observed follow-up while the Kaplan–Meier curve remains above zero.

Simulation scenarios demonstrate estimator behavior under controlled mechanisms. The stress test does not diagnose informative censoring in the VA trial or quantify its actual bias.

## Verification

The trial pipeline compares each arm's exact integration with the `lifelines` Kaplan–Meier RMST at all three horizons. The unit tests check uncensored sample means and variance, tied times, horizon boundaries, unsupported follow-up, contrast direction, the simulation's analytical target, and seed reproducibility.

## References

- [Royston and Parmar, restricted mean survival time in randomized trials](https://doi.org/10.1186/1471-2288-13-152)
- [Lifelines proportional hazards diagnostic documentation](https://lifelines.readthedocs.io/en/stable/lifelines.statistics.html#lifelines.statistics.proportional_hazard_test)
