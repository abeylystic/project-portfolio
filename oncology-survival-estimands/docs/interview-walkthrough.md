# Technical interview walkthrough

## Questions to prepare

1. **What problem did you pose?** Describe the treatment assignment, mortality endpoint, study population, time origin, outcome window, and sign of the comparison before naming the estimator.
2. **Why RMST?** Start with the clinical interpretation in days alive through day 180. Explain why a single hazard ratio answers a different question and how a changing treatment effect complicates that summary.
3. **Why use three horizons?** Discuss the choice of one focal window and the tradeoff between a longer window and smaller risk sets. Explain why the sensitivity estimates are descriptive and what would be required to claim a prospectively specified endpoint.
4. **How did you check your implementation?** Explain the uncensored check against the mean of truncated times, the tied-time risk set test, the comparison against an independent implementation, and coverage in simulations with known truth.
5. **What broke, and why?** Explain how prognostic patients selectively disappearing from the test-arm risk set distorted the later Kaplan–Meier curve despite equal underlying event distributions. Discuss how measured prognostic factors might support inverse probability of censoring weights and why unmeasured causes call for a sensitivity analysis.
6. **What would you change for a real trial?** Agree on the estimand and analysis plan before seeing outcomes, specify how intercurrent events and missing follow-up will be handled, evaluate expected support at the clinical horizon, and simulate design-specific power and coverage.

## Exercises to reproduce

- **Derive the result on paper.** For death days 1 and 2, a tied censoring at day 2, and a final death at day 4 among four patients, draw the Kaplan–Meier curve and calculate the area through day 5. Check it against `test_tied_failure_and_censoring_use_full_risk_set`.
- **Explain the horizon dependence.** Reproduce the arm-specific and contrast estimates and use the risk table to explain what the dataset can support at 90, 180, and 270 days.
- **Run a new seed.** Compare interval coverage with the checked-in output, using Monte Carlo standard error to explain ordinary simulation variability. Calculate how many replicates are needed for a coverage Monte Carlo standard error near 0.5 percentage points.
- **Change one data-generating assumption.** In the informative-censoring simulation, change the mean censoring time for high-risk test patients from 30 to 360 days. State what should happen to bias and coverage before running it. Explain the difference between censoring frequency and dependence between censoring and prognosis.
