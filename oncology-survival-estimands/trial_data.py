"""Download and validate a pinned copy of the historical VA trial data."""

import hashlib
import io
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

DATA_URL = (
    "https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/"
    "cc03c29690889dbe83089f5206e2422db8c3f71f/csv/survival/veteran.csv"
)
DATA_SHA256 = "3fba5cb9b15a10ab94d54e28b54d2c39d95c27b0fcff1c57aee0eb7edfd93950"
DEFAULT_DATA = Path(__file__).resolve().parent / "data" / "veteran.csv"
EXPECTED_COLUMNS = [
    "rownames", "trt", "celltype", "time", "status", "karno", "diagtime", "age", "prior",
]


def load_trial(path: Path = DEFAULT_DATA) -> pd.DataFrame:
    payload = path.read_bytes() if path.exists() else download()
    if hashlib.sha256(payload).hexdigest() != DATA_SHA256:
        raise ValueError(f"Dataset integrity check failed for {path}.")
    data = pd.read_csv(io.BytesIO(payload))
    if data.columns.tolist() != EXPECTED_COLUMNS or len(data) != 137:
        raise ValueError("Unexpected trial dataset layout.")
    if data.isna().any().any() or not data["rownames"].is_unique:
        raise ValueError("Unexpected missing values or duplicate row identifiers.")
    if data["trt"].value_counts().to_dict() != {1: 69, 2: 68}:
        raise ValueError("Unexpected randomized arm counts.")
    if not data["status"].isin([0, 1]).all() or not data["time"].gt(0).all():
        raise ValueError("Invalid event or follow-up coding.")
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    return data


def download() -> bytes:
    with urlopen(DATA_URL, timeout=30) as response:
        return response.read()
