"""Walk-forward 99% one-day value-at-risk for fixed Treasury exposures."""

import argparse
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import xlogy
from scipy.stats import chi2, norm

URL = "https://home.treasury.gov/system/files/276/yield-curve-rates-1990-2024.csv"
SHA256 = "258a594eb64e9d21d9e8574b0b8580b0a5230dce56bb284ec6aa49a9a260897f"
DV01 = {"2 Yr": 200, "5 Yr": 300, "10 Yr": 500}
MODELS = ("historical", "normal_ewma", "filtered_historical")
QUANTILE = 0.99
WINDOW = 500
DECAY = 0.94


def load_loss(path: Path) -> pd.Series:
    if path.is_file():
        payload = path.read_bytes()
    else:
        with urlopen(URL, timeout=30) as response:
            payload = response.read()
    if hashlib.sha256(payload).hexdigest() != SHA256:
        raise ValueError("Treasury source checksum does not match the reference archive.")
    curve = pd.read_csv(io.BytesIO(payload), usecols=["Date", *DV01])
    curve["Date"] = pd.to_datetime(curve["Date"], format="%m/%d/%y")
    curve = curve.set_index("Date").sort_index()
    if not curve.index.is_unique:
        raise ValueError("Unexpected duplicate Treasury dates.")
    # The archive lists one date with no quotes in any selected maturity.
    curve = curve.dropna(how="all")
    if curve.isna().any().any():
        raise ValueError("Unexpected partial curve in the selected maturities.")
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    changes_bp = 100 * curve.diff().dropna()
    return changes_bp.mul(pd.Series(DV01)).sum(axis=1).rename("loss")


def forecasts(loss: pd.Series, window: int = WINDOW) -> pd.DataFrame:
    past_loss = loss.shift(1)
    variance = past_loss.pow(2).ewm(alpha=1 - DECAY, adjust=False, min_periods=25).mean()
    sigma = np.sqrt(variance).replace(0, np.nan)
    standardized_shocks = loss / sigma
    return pd.DataFrame({
        "loss": loss,
        "historical": past_loss.rolling(window).quantile(QUANTILE),
        "normal_ewma": norm.ppf(QUANTILE) * sigma,
        "filtered_historical": (
            standardized_shocks.shift(1).rolling(window).quantile(QUANTILE) * sigma
        ),
    })


def bernoulli_loglik(hits: int, count: int, p: float) -> float:
    return float(xlogy(hits, p) + xlogy(count - hits, 1 - p))


def coverage_tests(hits: np.ndarray, expected: float = 1 - QUANTILE) -> tuple[float, float]:
    n = len(hits)
    observed = int(hits.sum())
    uc = 2 * (bernoulli_loglik(observed, n, observed / n) -
              bernoulli_loglik(observed, n, expected))
    counts = np.bincount(2 * hits[:-1].astype(int) + hits[1:], minlength=4).reshape(2, 2)
    totals = counts.sum(axis=1)
    if (totals == 0).any():
        return float(chi2.sf(max(uc, 0), 1)), float("nan")
    p01, p11 = counts[:, 1] / totals
    alternative = sum(
        bernoulli_loglik(int(counts[i, 1]), int(totals[i]), p)
        for i, p in enumerate((p01, p11))
    )
    transitions = int(counts[:, 1].sum())
    null = bernoulli_loglik(transitions, n - 1, transitions / (n - 1))
    independence = 2 * (alternative - null)
    return float(chi2.sf(max(uc, 0), 1)), float(chi2.sf(max(independence, 0), 1))


def evaluate(frame: pd.DataFrame, period: str) -> pd.DataFrame:
    rows = []
    for model in MODELS:
        error = frame["loss"] - frame[model]
        hits = error.gt(0).to_numpy()
        coverage_p, independence_p = coverage_tests(hits)
        rows.append({
            "period": period, "model": model, "days": len(frame),
            "exceptions": int(hits.sum()), "exception_rate": hits.mean(),
            "adjacent_exception_pairs": int((hits[:-1] & hits[1:]).sum()),
            "coverage_p": coverage_p, "independence_p": independence_p,
            "mean_var_usd": frame[model].mean(),
            "quantile_loss_usd": np.maximum(QUANTILE * error, (QUANTILE - 1) * error).mean(),
        })
    return pd.DataFrame(rows)


def plot_stress(frame: pd.DataFrame, path: Path) -> None:
    fig, axis = plt.subplots(figsize=(10, 4.5), layout="constrained")
    axis.plot(frame.index, frame["loss"], color="#B4BEC9", linewidth=0.9, label="Realized loss")
    for key, label, color in (
        ("historical", "Historical", "#476884"),
        ("normal_ewma", "Normal EWMA", "#289078"),
        ("filtered_historical", "Filtered historical", "#BD5420"),
    ):
        axis.plot(frame.index, frame[key], label=label, color=color, linewidth=1.2)
    axis.axhline(0, color="#717C87", linewidth=0.65)
    axis.set(title="One-day 99% VaR during the 2022 rate shock", ylabel="Dollars")
    axis.legend(ncol=2, frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=145)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/treasury.csv"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    frame = forecasts(load_loss(args.data)).loc["2015":"2024"].dropna()
    stress = frame.loc["2022"]
    report = pd.concat([
        evaluate(frame, "2015-2024"), evaluate(frame.loc["2020"], "2020"),
        evaluate(stress, "2022"),
    ], ignore_index=True)
    args.output.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.output / "backtest.csv", index=False, float_format="%.8f")
    plot_stress(stress, args.output / "rate-shock.png")


if __name__ == "__main__":
    main()
