# Statistical research and data science portfolio

Four reproducible projects connecting biostatistics, applied economics, quantitative risk, and data science. Each project includes its analysis code, assumptions, results, and tests.

| Role | Project | Main finding |
| --- | --- | --- |
| Biostatistician | [Clinical trial survival analysis](oncology-survival-estimands/) | The test regimen had an estimated **13.8 fewer days alive through day 180**. A censoring stress test shows how an apparent benefit can arise under a true null. |
| Quant researcher / risk analyst | [Treasury value-at-risk backtest](treasury-var-backtest/) | Filtered historical VaR had **28 exceptions in 2,501 days**, closest to nominal coverage. Gaussian EWMA earned the lowest quantile loss. |
| Statistician | [Matched-pair cluster experiment design](cluster-experiment-design/) | With a site ICC of 0.15, naive individual-level inference yielded **34.0% false positives** in a nominal 5% test. |
| Data scientist | [Hourly demand forecasting and monitoring](bike-demand-monitoring/) | Boosted trees improved overall MAE by **16.1%** against a weekly baseline but lost to it in two of six test months. |

## Run the projects

Each project can be run from its own directory using the commands in its README. To use a shared environment, install the root dependencies from the repository root with Python 3.12.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

| Directory | Command |
| --- | --- |
| `oncology-survival-estimands` | `python analysis.py --output results/local` |
| `treasury-var-backtest` | `python backtest.py --output results/local` |
| `cluster-experiment-design` | `python design.py --output results/local` |
| `bike-demand-monitoring` | `OMP_NUM_THREADS=2 python forecast.py train --output results/local` |

From a project directory run `python -m unittest discover -s tests -v`. Example results are tracked in each project's `results/` directory. The three empirical projects download checksum-verified public source files on first run. The cluster study generates data from a documented simulation model. Raw data, local outputs, Python environments, and fitted binary models are excluded from version control.
