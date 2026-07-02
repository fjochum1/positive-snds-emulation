<<<<<<< HEAD
# Safety of endocrine therapy interruption for pregnancy in HR-positive early breast cancer: a nationwide target trial emulation

Analytic code for the target trial emulation of temporary endocrine therapy (ET)
interruption for pregnancy in young women with hormone receptor-positive (HR+)
early breast cancer, using the French National Health Data System (SNDS).

Two target trials are emulated with a clone–censor–weight (CCW) design and
weighted pooled-logistic marginal structural models (MSM):

- **Trial 1** — benchmark emulation of the POSITIVE trial (interruption at
  18–30 months of ET vs uninterrupted ET), overall and by baseline recurrence risk.
- **Trial 2** — extension to nine 6-month interruption windows (6–60 months of ET)
  vs uninterrupted ET for ≥24 or ≥60 months.

## Data availability

The individual-level SNDS data used in this study **cannot be shared publicly**:
French regulation and data-protection law prohibit their redistribution. Any
person or institution may, however, request access from the French Data Protection
Authority (CNIL) to conduct a study of public interest. 

## Repository structure

```
positive-snds-emulation/
├── README.md
└── R/
    ├── 00_setup.R                     # paths, packages, plotting theme (edit paths here)
    ├── 01_trial1_clone_censor_weight.R  # Trial 1: clone → censor → IP weights → MSM
    ├── 02_trial2_clone_censor_weight.R  # Trial 2: clone → censor → IP weights → MSM
    ├── 03_trial1_bootstrap_results.R    # Trial 1: bootstrap CIs, 5-year risk differences, figures
    ├── 04_trial2_bootstrap_results.R    # Trial 2: bootstrap CIs, 5-year risk differences
    ├── functions_msm_bootstrap.R        # shared MSM, bootstrap and CI-extraction functions
    ├── 05_figures.R                     # Figures 2, 3, 4 from the bootstrap outputs
    └── 06_sensitivity_analyses.R        # the single change made for each sensitivity analysis
```

## Input datasets (not shared)

Each trial uses a **person-month** dataset (one row per person per month of
follow-up from time zero) containing, at minimum:

- `person_id`, `month`
- `et`, `et_lag1`, `et_lag2`, `et_lag3` — ET coverage in the month and lags
- `pregnancy`, `had_pregnancy`, `planning_pregnancy_past_6mo`
- `relapse`, `relapse_month` — breast-cancer event (composite) and its month
- baseline covariates: `age_cl`, `previous_pregnancy`, `breast_surgery_2cl`,
  `pnuicc_2cl_modified` (nodal status), `ct_setting_5cl`, `rt`, `ht_type`,
  `her2_status`, `nb_comor_ql`, `Cardiovascular`, `Endocrine_and_metabolism`,
  `Psychiatric_disorders`, `Other`, `depriv_index_quintile`, `period_diag`
- `cumulative_side_effects`, `followup_past6m`

`pop_positive_first.RData` (Trial 1, time zero = 18 months of ET) and
`pop_positive_second.RData` (Trial 2, time zero = ET initiation).

## Software and packages

Analyses were run in **R version 4.5**. Required packages:
`data.table`, `dplyr`, `tidyr`, `purrr`, `splines`, `speedglm`, `broom`,
`ggplot2`, `patchwork`.

## Method summary

At time zero each eligible woman is cloned to every strategy she is eligible to
follow; clones are artificially censored at the first deviation from their
assigned strategy; and inverse-probability-of-censoring weights (from two
pooled-logistic models — monthly probability of being on ET, and monthly
probability of pregnancy) restore the baseline-comparable population. The
per-protocol effect (5-year risk difference under sustained adherence) is
estimated from a weighted pooled-logistic MSM, with 95% confidence intervals
from a nonparametric bootstrap. See the manuscript and Supplementary Methods
for full detail.

## Citation / correspondence

Jochum F, et al. *Safety of endocrine therapy interruption for pregnancy in
hormone receptor-positive early breast cancer: a nationwide target trial
emulation. BMJ 2026* Correspondence: jochum.floriane@gmail.com
=======
# positive-snds-emulation
>>>>>>> 1510ebaf20116867a0285acf22485ccf0ee4d4c4
