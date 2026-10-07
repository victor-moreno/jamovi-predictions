# clinpred — clinical prediction models for jamovi

A [jamovi](https://www.jamovi.org) module to **develop and validate** clinical prediction models,
following the clinical-epidemiology workflow (Steyerberg; TRIPOD+AI 2024; Riley et al.). A
clinician can fit a model, see how optimistic its apparent performance is, validate it in new
patients, and see which predictors drive the predictions. Version 0.1 covers **binary outcomes
(logistic regression)**; linear and Cox models will follow, with the same options and results.

Menu: **Prediction ▸ Binary Outcome (Logistic)**. Needs **jamovi 28.4 or later**.

## Features

- **Predictors**: linear covariates, **restricted cubic splines** (3–7 knots, Harrell's
  placement) for non-linear associations, and factors.
- **Performance**, with confidence intervals where they exist:
  - *Discrimination*: C-statistic (AUC, DeLong CI) and discrimination slope.
  - *Calibration*: calibration-in-the-large, calibration slope, O/E ratio, and ICI, E50, E90
    and Emax from a loess curve.
  - *Overall*: Brier score, scaled Brier score, Nagelkerke and Cox–Snell R².
- **Internal validation**:
  - Harrell's **bootstrap optimism correction**, giving the uniform shrinkage factor too.
  - Or **repeated k-fold cross-validation**, stratified by outcome.
- **Prediction stability** (Riley & Collins 2023): MAPE, plus instability plots of
  individual risks and of calibration curves.
- **Validation in new data**: the development model is applied unchanged to
  - the rows marked by a **validation split** variable (another period or centre), and/or
  - an **external data file** (csv, tsv, txt, omv, xlsx or sav).
  
  For the external file, a variable check reports missing columns, names or levels that
  differ only in case, unknown levels, and values outside the development range. A
  **case-mix** table gives means, percentages, the spread of the linear predictor and the
  membership c-statistic (Debray et al. 2015).
- **Clinical utility**: classification at a fixed, prevalence-based or Youden threshold
  (sensitivity, specificity, PPV, NPV, accuracy with Wilson CIs), and **decision curves**
  with a net-benefit table and interventions avoided.
- **Sample size adequacy**: Riley et al. (2020) criteria, and events per parameter.
- **Variable importance**:
  - **Shapley decomposition of R²** (LMG / general dominance; exact up to 10 predictors).
  - **Permutation importance** (AUC loss, Brier increase), computed in the validation data
    when there is any.
  - Exact **SHAP** values on the log-odds scale: summary (beeswarm) plot, mean |SHAP|, and
    a waterfall explaining one patient's prediction.
- Coefficients (odds ratios) and the **model equation**, with the spline formulae, so the
  model can be reported and transported.
- ROC curves, calibration plots (smoothed curve with CI, grouped observed proportions, risk
  histogram) and risk distributions, one panel per data set.
- Saves the predicted risk and the linear predictor to the spreadsheet, for every row with
  complete predictors (validation rows included).

Example data in jamovi's Data Library: Pima diabetes (`MASS::Pima.tr`), and Pima.tr and Pima.te
stacked with a `sample` variable for split validation. For external validation, load
[`examples/pima_external.csv`](examples/pima_external.csv) (`MASS::Pima.te`).

Planned: linear and Cox models (C-index, time-dependent AUC, calibration at a horizon), validation
of published models, model updating (recalibration), penalised regression and machine-learning
learners through the same validation pipeline, and Spanish and Catalan translations.

## Install

Download the `.jmo` for your CPU from the [Releases](../../releases) page and, in jamovi, choose
**Modules (⊕) ▸ Sideload**:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (28.4 or later, R 4.6.0) | `clinpred_<version>_current_R4.6.0_arm64.jmo` | `clinpred_<version>_current_R4.6.0_x64.jmo` |

The same file works on macOS, Windows and Linux. jamovi *solid* is not supported: the external
validation file uses jamovi's File option, which appeared in 28.4.

## Acknowledgements

The methods follow E. W. Steyerberg, *Clinical Prediction Models* (2019), F. E. Harrell's
bootstrap validation (as in the R package [rms](https://cran.r-project.org/package=rms), used
with [pROC](https://cran.r-project.org/package=pROC) and
[Hmisc](https://cran.r-project.org/package=Hmisc) to check the results), and R. D. Riley and
colleagues' work on sample size and stability. Data from the R package
[MASS](https://cran.r-project.org/package=MASS) (Venables & Ripley). Designed and built together
with Claude Code.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII),
“Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number
FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health
(CIBERESP), action Genrisk.

## Licence

GPL (≥ 2).
