# Data provenance

The analysis uses the historically published Veterans' Administration lung cancer trial, distributed as `veteran` in the R `survival` package. Its documentation describes a randomized comparison of two lung cancer treatment regimens. The historical dataset is downloaded from the pinned Rdatasets CSV mirror on first run and cached locally. The downloader verifies its SHA-256 hash. The patient-level CSV is excluded from this repository.

- [R survival package documentation](https://search.r-project.org/CRAN/refmans/survival/help/veteran.html)
- [Pinned Rdatasets CSV snapshot](https://github.com/vincentarelbundock/Rdatasets/blob/cc03c29690889dbe83089f5206e2422db8c3f71f/csv/survival/veteran.csv)
- Kalbfleisch JD, Prentice RL. *The Statistical Analysis of Failure Time Data*. Wiley, 1980.

| CSV field | Meaning |
| --- | --- |
| `rownames` | Row identifier added by the mirror |
| `trt` | Randomized arm, 1 standard and 2 test |
| `time` | Observed event or censoring time in days |
| `status` | 1 death, 0 right-censored |
| `celltype` | Tumor cell type |
| `karno` | Karnofsky performance score |
| `diagtime` | Months from diagnosis to randomization |
| `age` | Years |
| `prior` | Prior therapy, 0 no and 10 yes |

The main analysis uses assigned arm, survival time, and event indicator. The data review verifies the expected schema and sample size, complete fields, unique row identifiers, arm sizes, valid status coding, and positive times. Tied event days are expected for integer-day follow-up. The data dictionary provides the treatment labels shown above. The analysis uses these labels throughout.
