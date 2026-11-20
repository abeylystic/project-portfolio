"""Kaplan–Meier integration and Greenwood variance for unadjusted RMST."""

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm


class UnsupportedHorizon(ValueError):
    """RMST cannot be estimated at the requested horizon."""


@dataclass(frozen=True)
class RestrictedMean:
    days: float
    variance: float


@dataclass(frozen=True)
class KaplanMeier:
    event_times: np.ndarray
    survival: np.ndarray
    failures: np.ndarray
    risk: np.ndarray
    last_follow_up: float

    def rmst(self, horizon: float) -> RestrictedMean:
        """Integrate the right-continuous step function and its covariance."""
        if not np.isfinite(horizon) or horizon <= 0:
            raise ValueError("The horizon must be finite and positive.")
        if horizon > self.last_follow_up and (
            not len(self.survival) or self.survival[-1] > 0
        ):
            raise UnsupportedHorizon("The horizon extends beyond observed survival support.")

        # An event exactly at the boundary changes S(horizon), with zero area.
        count = int(np.searchsorted(self.event_times, horizon, side="left"))
        knots = np.r_[0.0, self.event_times[:count], horizon]
        areas = np.diff(knots) * np.r_[1.0, self.survival[:count]]

        # Each Greenwood increment contributes the squared area after its event.
        tail_area = np.cumsum(areas[:0:-1])[::-1]
        at_risk = self.risk[:count]
        failures = self.failures[:count]
        denominator = at_risk * (at_risk - failures)
        increments = np.divide(
            failures,
            denominator,
            out=np.zeros(count, dtype=float),
            where=denominator > 0,
        )
        # If all remaining patients die together, their remaining area is zero.
        variance = float(np.dot(increments, tail_area**2))
        return RestrictedMean(float(areas.sum()), variance)


def kaplan_meier(time: np.ndarray, event: np.ndarray) -> KaplanMeier:
    """Fit a right-censored KM curve, with failures preceding tied censoring."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event)
    if time.ndim != 1 or event.shape != time.shape or not time.size:
        raise ValueError("Supply aligned, nonempty vectors of times and events.")
    if not np.isfinite(time).all() or (time <= 0).any():
        raise ValueError("Follow-up times must be finite and positive.")
    if not np.isin(event, [0, 1]).all():
        raise ValueError("Event indicators must be 0 or 1.")

    times, groups, counts = np.unique(time, return_inverse=True, return_counts=True)
    risk = np.cumsum(counts[::-1])[::-1]
    failures = np.bincount(groups, weights=event.astype(float), minlength=len(times))
    observed = failures > 0
    return KaplanMeier(
        event_times=times[observed],
        survival=np.cumprod(1 - failures[observed] / risk[observed]),
        failures=failures[observed],
        risk=risk[observed],
        last_follow_up=float(times[-1]),
    )


def rmst_difference(
    standard: KaplanMeier, test: KaplanMeier, horizon: float
) -> tuple[float, float, float, float]:
    """Return test minus standard, its standard error, and a 95% Wald CI."""
    control_mean = standard.rmst(horizon)
    test_mean = test.rmst(horizon)
    difference = test_mean.days - control_mean.days
    se = np.sqrt(test_mean.variance + control_mean.variance)
    half_width = norm.ppf(0.975) * se
    return (
        float(difference), float(se),
        float(difference - half_width), float(difference + half_width),
    )
