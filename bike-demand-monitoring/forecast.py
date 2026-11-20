"""Train, evaluate, and score an hourly bikeshare demand forecast."""

import argparse
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error

URL = "https://archive.ics.uci.edu/static/public/275/bike+sharing+dataset.zip"
SHA256 = "b70182d0d0508e9abbb79306ce5c0cec34869000f8220175ac83d11dbe845401"
LAGS = (1, 2, 3, 24, 48, 168)
TRAIN_END = "2011-12-31"
VALIDATION_START = "2012-01-01"
VALIDATION_END = "2012-06-30"
TEST_START = "2012-07-01"
GRID = ((7, 140), (15, 180), (31, 220))


def validate_counts(history: pd.Series) -> pd.Series:
    if history.index.hasnans or history.index.has_duplicates:
        raise ValueError("History timestamps must be nonmissing and unique.")
    if not history.index.equals(history.index.floor("h")):
        raise ValueError("History timestamps must be aligned to the hour.")
    if not np.isfinite(history).all() or (history < 0).any() or (history % 1 != 0).any():
        raise ValueError("Rental counts must be finite nonnegative integers.")
    return history.sort_index()


def load_history(path: Path) -> pd.Series:
    if path.is_file():
        payload = path.read_bytes()
    else:
        with urlopen(URL, timeout=30) as response:
            payload = response.read()
    if hashlib.sha256(payload).hexdigest() != SHA256:
        raise ValueError("UCI archive checksum does not match the reference release.")
    with ZipFile(io.BytesIO(payload)) as source:
        raw = pd.read_csv(source.open("hour.csv"), usecols=["dteday", "hr", "cnt"])
    timestamp = pd.DatetimeIndex(pd.to_datetime(raw["dteday"]) + pd.to_timedelta(raw["hr"], unit="h"))
    history = validate_counts(pd.Series(raw["cnt"].to_numpy(), index=timestamp, name="count"))
    if len(history) != 17379:
        raise ValueError("Unexpected number of hourly records.")
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    return history


def make_features(history: pd.Series, forecast_times: pd.DatetimeIndex) -> pd.DataFrame:
    """Construct features from the calendar and strictly earlier hourly counts."""
    frame = pd.DataFrame({
        "hour": forecast_times.hour,
        "weekday": forecast_times.dayofweek,
        "month": forecast_times.month,
        "weekend": forecast_times.dayofweek >= 5,
    }, index=forecast_times)
    for hours in LAGS:
        earlier = forecast_times - pd.Timedelta(hours=hours)
        frame[f"lag_{hours}h"] = history.reindex(earlier).to_numpy()
    frame["recent_mean"] = frame[["lag_1h", "lag_2h", "lag_3h"]].mean(axis=1)
    return frame


def baseline_predictions(features: pd.DataFrame, training: pd.Series) -> pd.DataFrame:
    medians = training.groupby([training.index.dayofweek, training.index.hour]).median()
    index = pd.MultiIndex.from_arrays([features.index.dayofweek, features.index.hour])
    fallback = pd.Series(medians.reindex(index).to_numpy(), index=features.index).fillna(training.median())
    return pd.DataFrame({
        "persistence": features["lag_1h"].fillna(fallback),
        "seasonal_naive": features["lag_168h"].fillna(fallback),
    })


def train(history: pd.Series, features: pd.DataFrame):
    training = history.loc[:TRAIN_END]
    validation = history.loc[VALIDATION_START:VALIDATION_END]
    candidates = []
    models = []
    for leaves, iterations in GRID:
        model = HistGradientBoostingRegressor(
            max_leaf_nodes=leaves, max_iter=iterations, learning_rate=0.075,
            early_stopping=False, random_state=2026,
        )
        model.fit(features.loc[training.index], training)
        prediction = model.predict(features.loc[validation.index]).clip(min=0)
        candidates.append({
            "max_leaf_nodes": leaves,
            "max_iter": iterations,
            "validation_mae": mean_absolute_error(validation, prediction),
        })
        models.append(model)
    leaderboard = pd.DataFrame(candidates)
    best = leaderboard["validation_mae"].idxmin()
    leaderboard["selected"] = leaderboard.index == best
    return models[best], leaderboard


def evaluate(actual: pd.Series, predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model in predictions:
        error = predictions[model] - actual
        rows.append({
            "model": model, "hours": len(actual), "mae": error.abs().mean(),
            "rmse": np.sqrt((error**2).mean()), "mean_error": error.mean(),
            "wape": error.abs().sum() / actual.sum(),
        })
    return pd.DataFrame(rows)


def monitor(
    actual: pd.Series, predictions: pd.DataFrame, features: pd.DataFrame, validation_mae: float
) -> pd.DataFrame:
    rows = []
    for month, observed in actual.groupby(actual.index.to_period("M")):
        baseline = (predictions.loc[observed.index, "seasonal_naive"] - observed).abs().mean()
        residual = predictions.loc[observed.index, "boosted_trees"] - observed
        level = 1.25 * validation_mae
        mae = residual.abs().mean()
        rows.append({
            "month": str(month), "expected_hours": month.days_in_month * 24,
            "observed_hours": len(observed), "mean_demand": observed.mean(),
            "missing_previous_hour_rate": features.loc[observed.index, "lag_1h"].isna().mean(),
            "seasonal_naive_mae": baseline, "mae": mae, "mean_error": residual.mean(),
            "review_threshold": level, "review": bool(mae > level),
        })
    return pd.DataFrame(rows)


def monitoring_plot(monitoring: pd.DataFrame, output: Path) -> None:
    fig, axis = plt.subplots(figsize=(8, 4.1), layout="constrained")
    axis.plot(monitoring["month"], monitoring["seasonal_naive_mae"], "o-", color="#9D462B",
              label="Previous-week baseline")
    axis.plot(monitoring["month"], monitoring["mae"], "o-", color="#235B85", label="Boosted trees")
    axis.axhline(monitoring["review_threshold"].iloc[0], color="#636C78", linestyle="--",
                 label="Review threshold")
    axis.set(title="Monthly forecast monitoring on the held-out period",
             ylabel="Mean absolute error (rentals per hour)", xlabel="Test month")
    axis.grid(axis="y", color="#DEE4E9")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "monitoring.png", dpi=145)
    plt.close(fig)


def run(data_path: Path, output: Path, model_out: Path | None) -> None:
    history = load_history(data_path)
    features = make_features(history, history.index)
    model, leaderboard = train(history, features)
    calibration = float(leaderboard.loc[leaderboard["selected"], "validation_mae"].iloc[0])
    actual = history.loc[TEST_START:]
    predictions = baseline_predictions(features.loc[actual.index], history.loc[:TRAIN_END])
    predictions["boosted_trees"] = model.predict(features.loc[actual.index]).clip(min=0)
    monitoring = monitor(actual, predictions, features, calibration)
    output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "validation.csv": leaderboard,
        "test_metrics.csv": evaluate(actual, predictions),
        "monitoring.csv": monitoring,
        "data_audit.csv": pd.DataFrame([{
            "archive_hours": len(history),
            "calendar_hours": len(pd.date_range("2011-01-01", "2012-12-31 23:00", freq="h")),
            "training_hours": len(history.loc[:TRAIN_END]),
            "validation_hours": len(history.loc[VALIDATION_START:VALIDATION_END]),
            "test_hours": len(actual),
        }]),
    }
    for name, table in outputs.items():
        table.to_csv(output / name, index=False, float_format="%.8f")
    monitoring_plot(monitoring, output)
    if model_out:
        model_out.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(model, model_out)


def score(model_path: Path, history_path: Path, target: str) -> None:
    raw = pd.read_csv(history_path, usecols=["timestamp", "count"])
    timestamps = pd.DatetimeIndex(pd.to_datetime(raw["timestamp"], format="ISO8601"))
    horizon = pd.Timestamp(target)
    if timestamps.tz is not None or horizon.tzinfo is not None or horizon != horizon.floor("h"):
        raise ValueError("Use timezone-naive local calendar timestamps aligned to the hour.")
    history = validate_counts(pd.Series(raw["count"].to_numpy(), index=timestamps))
    if history.index.max() >= horizon:
        raise ValueError("Provide only observed rental counts strictly before the forecast hour.")
    features = make_features(history, pd.DatetimeIndex([horizon]))
    # Load only the locally fitted artifact produced by this project.
    model = joblib.load(model_path)
    prediction = max(0.0, float(model.predict(features)[0]))
    print(f"{horizon.isoformat()},{prediction:.2f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="action", required=True)
    training = subcommands.add_parser("train")
    training.add_argument("--data", type=Path, default=Path("data/bike-sharing.zip"))
    training.add_argument("--output", type=Path, default=Path("results"))
    training.add_argument("--model-out", type=Path)
    scoring = subcommands.add_parser("score")
    scoring.add_argument("--model", type=Path, required=True)
    scoring.add_argument("--history", type=Path, required=True)
    scoring.add_argument("--target", required=True)
    args = parser.parse_args()
    if args.action == "train":
        run(args.data, args.output, args.model_out)
    else:
        score(args.model, args.history, args.target)


if __name__ == "__main__":
    main()
