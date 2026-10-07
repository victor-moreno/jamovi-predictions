# clinpred — clinical prediction models for jamovi

A [jamovi](https://www.jamovi.org) module to **develop and validate** clinical prediction models,
following the clinical-epidemiology workflow (Steyerberg; TRIPOD+AI 2024; Riley et al.). A
clinician can fit a model, see how optimistic its apparent performance is, validate it in new
patients, update it, and see which predictors drive the predictions. Three analyses share the
same options, results and validation engine:

| menu: **Prediction ▸** | outcome | model | effects |
| --- | --- | --- | --- |
| **Binary Outcome (Logistic)** | event / no event | logistic regression | odds ratios |
| **Continuous Outcome (Linear)** | a measurement | linear regression | differences (β) |
| **Time to Event (Cox)** | time + event, risk by a **horizon** | Cox regression | hazard ratios |

Needs **jamovi 28.4 or later**.

## Features

- **Predictors**:
  - Linear covariates, with **covariate scaling** (effect per unit, per 1 SD, per 10, 100
    or 1000 units).
  - **Restricted cubic splines** (3–7 knots, Harrell's placement) for non-linear
    associations.
  - Factors, with a chosen **reference level**.
- **Interactions**, two ways: listed in an interaction builder, or all pairs of
  predictors as candidates.
  - A screening table gives each candidate added alone (LR χ², ΔAIC, ΔBIC).
  - The candidates are all kept, or selected forward by AIC, BIC or LR test.
  - The selection is **repeated in every bootstrap or cross-validation fit**, so the
    optimism includes it, and the table reports how often each term was selected.
- **Explaining the model**:
  - A **descriptive table** (by outcome group, or with each predictor's association with a
    continuous outcome), with missing values.
  - **Crude and adjusted effects** (OR, β or HR) with reference rows; splines as Q3 vs Q1.
  - **Likelihood ratio tests** (F tests for the linear model) per predictor, for the
    non-linear part of each spline and for each interaction.
  - **Effects within subgroups** for each interaction kept.
- **Performance**, with confidence intervals where they exist:
  - *Logistic*: C-statistic (AUC, DeLong CI), discrimination slope; calibration-in-the-large,
    calibration slope, O/E, ICI, E50, E90, Emax; Brier, scaled Brier, Nagelkerke and Cox–Snell R².
  - *Linear*: R², RMSE, MAE, calibration-in-the-large and calibration slope.
  - *Cox, at the horizon*: Harrell's and Uno's C, time-dependent AUC (censoring weighted),
    Royston's D and R²_D; O/E, calibration slope, ICI, E50, E90 (Austin et al. 2020);
    censoring-weighted Brier score and the index of prediction accuracy; Nagelkerke and
    Cox–Snell R².
- **Internal validation**:
  - Harrell's **bootstrap optimism correction**, giving the uniform shrinkage factor too.
  - Or **repeated k-fold cross-validation**, stratified by outcome.
- **Prediction stability** (Riley & Collins 2023): MAPE, plus instability plots of
  individual predictions and of calibration curves.
- **Validation in new data**: the development model is applied unchanged to
  - the rows marked by a **validation split** variable (another period or centre), and/or
  - an **external data file** (csv, tsv, txt, omv, xlsx or sav).

  For the external file, a variable check reports missing columns, names or levels that
  differ only in case, unknown levels, and values outside the development range. A
  **case-mix** table gives means, percentages, the spread of the linear predictor and the
  membership c-statistic (Debray et al. 2015).
- **Model updating** (logistic and linear) in each validation set: recalibration in the large,
  recalibration (intercept and slope) and model revision, with the **closed testing
  procedure** of Vergouwe et al. (2017) choosing among them.
- **Clinical utility**: classification at a fixed, prevalence-based or Youden threshold
  (logistic: sensitivity, specificity, PPV, NPV, accuracy with Wilson CIs), and **decision
  curves** with a net-benefit table and interventions avoided (Cox: at the horizon, with
  Kaplan–Meier risks).
- **Sample size adequacy**: Riley et al. criteria for each outcome type (2019, 2020).
- **Variable importance**:
  - **Shapley decomposition of R²** (LMG / general dominance; exact up to 10 predictors).
  - **Permutation importance** (loss in AUC, R² or C; increase in Brier score or RMSE),
    computed in the validation data when there is any.
  - Exact **SHAP** values on the model's linear scale: summary (beeswarm) plot, mean |SHAP|,
    and a waterfall explaining one patient's prediction.
- Coefficients and the **model equation**, with the spline formulae (and, for Cox, the baseline
  survival at the horizon), so the model can be reported and transported.
- Plots, one panel per data set: calibration (smoothed curve with CI, grouped observed values,
  histogram of predictions), ROC curves and risk distributions (logistic), Kaplan–Meier curves
  by risk group (Cox), decision curves.
- Saves the predictions to the spreadsheet (risk and linear predictor; predicted value; risk by
  the horizon and prognostic index), for every row with complete predictors.

Example data in jamovi's Data Library:

- Pima diabetes (`MASS::Pima.tr`), and Pima.tr and Pima.te stacked with a `sample` variable for
  split validation; for external validation load
  [`examples/pima_external.csv`](examples/pima_external.csv) (`MASS::Pima.te`). `glu` serves as
  a continuous outcome.
- Node-positive breast cancer from the Rotterdam tumour bank (`survival::rotterdam`, recurrence-free
  survival in years), and the same stacked with the German Breast Cancer Study Group trial
  (`survival::gbsg`) as in Royston & Altman (2013); for external validation load
  [`examples/gbsg_external.csv`](examples/gbsg_external.csv).

Planned: validation of published models and model export, model updating for Cox, penalised
regression and machine-learning learners through the same validation pipeline, and Spanish and
Catalan translations.

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
bootstrap validation (as in the R package [rms](https://cran.r-project.org/package=rms)), and
R. D. Riley and colleagues' work on sample size and stability. The results are checked against
[rms](https://cran.r-project.org/package=rms), [pROC](https://cran.r-project.org/package=pROC),
[Hmisc](https://cran.r-project.org/package=Hmisc),
[survival](https://cran.r-project.org/package=survival) (T. Therneau),
[timeROC](https://cran.r-project.org/package=timeROC),
[relaimpo](https://cran.r-project.org/package=relaimpo) and
[pmsampsize](https://cran.r-project.org/package=pmsampsize). Data from the R packages
[MASS](https://cran.r-project.org/package=MASS) (Venables & Ripley) and survival. Designed and
built together with Claude Code.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII),
“Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number
FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health
(CIBERESP), action Genrisk.

## Licence

GPL (≥ 2).
