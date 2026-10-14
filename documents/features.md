# clinpred — features and methods in detail

The [README](../README.md) gives the overview. This document lists what each analysis computes,
the options behind it and the methods used. The worked tour of every feature on example data is in
[`examples/README.md`](../examples/README.md), and the
[tutorial](https://victor-moreno.github.io/jamovi-predictions/) explains the methods with worked
examples.

Three development analyses share the same options, results and validation engine; a fourth
validates published models:

| menu: **Prediction ▸** | outcome | model | effects |
| --- | --- | --- | --- |
| **Binary Outcome (Logistic)** | event / no event | logistic regression | odds ratios |
| **Continuous Outcome (Linear)** | a measurement | linear regression | differences (β) |
| **Time to Event (Cox)** | time + event, risk by a **horizon** | Cox regression | hazard ratios |
| **Validate a Published Model** | any of the three | read from a file | as published |

## Predictors

- Linear covariates, with **covariate scaling** (effect per unit, per 1 SD, per 10, 100 or 1000
  units).
- **Restricted cubic splines** (3–7 knots, Harrell's placement) for non-linear associations.
- Factors, with a chosen **reference level**. Ordinal variables are coded as unordered factors
  (level against reference).

## Baseline covariates and selection

- Predictors marked as **baseline** (e.g. institution, age, sex) are always in the model.
- The others are kept, selected **backward** or **forward** (AIC, BIC or likelihood ratio test),
  or penalised by an **elastic net** (λ by cross-validation; checked against glmnet), the selection
  or tuning repeated in every bootstrap or cross-validation fit. Optionally a **relaxed lasso**
  (the unpenalised refit of the predictors the elastic net keeps).
- After any selection, the coefficient tables warn that their CIs ignore the selection, and a
  **bootstrap coefficient table** (Heinze, Wallisch & Dunkler 2018) gives each coefficient's
  inclusion frequency, bootstrap median and percentile interval, RMSD ratio and relative
  conditional bias, with the most frequent selected models.
- Every predictor table and importance plot lists the baseline covariates last, under their own
  heading.
- A penalised (not relaxed) elastic net hides what needs an unpenalised likelihood: standard
  errors, CIs and p values, likelihood ratio tests, the Shapley decomposition and the heuristic
  shrinkage factor.

## Incremental value

Of the other predictors over the baseline model: Δ of every performance measure (apparent,
optimism-corrected on the same resamples, and in validation data), ΔAUC (paired DeLong) or ΔC with
confidence intervals, LR or F test of the added predictors, continuous NRI, NRI over the risk
categories cut at the net-benefit thresholds, IDI, Δ net benefit, and overlaid ROC curves
(logistic) and decision curves of both models.

## Missing values

- An exploration table (missing by variable, by outcome group, with tests; and the share missing in
  each validation set) and the most frequent patterns.
- Methods: complete cases; a **missing category / indicator** (a `(missing)` level for factors; the
  development median plus an indicator predictor for numeric variables); **single imputation**; or
  **multiple imputation** by chained equations (predictive mean matching, written in pure R; the
  outcome in the imputation models, for Cox the event and its Nelson–Aalen cumulative hazard).
- Under multiple imputation: coefficients pooled by Rubin's rules, tests by pooled Wald (D1),
  selection by majority over imputations, performance averaged and optimism estimated within each
  imputation (Wahl et al. 2016); validation data imputed separately.

## Interactions

- Two ways: listed in an interaction builder, or all pairs of predictors as candidates.
- A screening table gives each candidate added alone (LR χ², ΔAIC, ΔBIC).
- The candidates are all kept, or selected forward by AIC, BIC or LR test.
- The selection is **repeated in every bootstrap or cross-validation fit**, so the optimism
  includes it, and the table reports how often each term was selected.

## Explaining the model

- A **descriptive table** (by outcome group, or with each predictor's association with a
  continuous outcome), with missing values.
- **Crude and adjusted effects** (OR, β or HR) with reference rows; splines as Q3 vs Q1; after a
  selection, the model with all predictors next to the selected model; baseline covariates last;
  CI and p columns optional; a **forest plot** of any of the three.
- **Likelihood ratio tests** (F tests for the linear model) per predictor, for the non-linear part
  of each spline and for each interaction.
- **Effects within subgroups** for each interaction kept.
- Coefficients and the **model equation**, with the spline formulae (and, for Cox, the baseline
  survival at the horizon), so the model can be reported and transported.

## Performance

With confidence intervals where they exist:

- *Logistic*: C-statistic (AUC, DeLong CI), discrimination slope; calibration-in-the-large,
  calibration slope, O/E, ICI, E50, E90, Emax; Brier, scaled Brier, Nagelkerke and Cox–Snell R².
- *Linear*: R², RMSE, MAE, calibration-in-the-large and calibration slope.
- *Cox, at the horizon*: Harrell's and Uno's C, time-dependent AUC (censoring weighted), Royston's
  D and R²_D; O/E, calibration slope, ICI, E50, E90 (Austin et al. 2020); censoring-weighted Brier
  score and the index of prediction accuracy; Nagelkerke and Cox–Snell R².

## Internal validation

- Harrell's **bootstrap optimism correction**, giving the uniform shrinkage factor too.
- Or **repeated k-fold cross-validation**, stratified by outcome.
- **Prediction stability** (Riley & Collins 2023): MAPE, plus instability plots of individual
  predictions and of calibration curves.

## Validation in new data

The development model is applied unchanged to

- the rows marked by a **validation split** variable (another period or centre), and/or
- an **external data file** (csv, tsv, txt, omv, xlsx or sav).

For the external file, a variable check reports missing columns, names or levels that differ only
in case, unknown levels, and values outside the development range. A **case-mix** table gives
means, percentages, the spread of the linear predictor and the membership c-statistic (Debray et
al. 2015).

**Model updating** in each validation set: recalibration in the large, recalibration (intercept
and slope) and model revision. For logistic and linear models the **closed testing procedure** of
Vergouwe et al. (2017) chooses among them. For Cox models the baseline risk at the horizon is
re-estimated with the linear predictor as offset; as the partial likelihood does not involve it,
the choice is made by sequential tests from the most complex method down: revision vs
recalibration and slope = 1 by partial likelihood ratio, O/E at the horizon = 1 by a Wald test.

## Clinical utility

Classification at a fixed, prevalence-based or Youden threshold (logistic: sensitivity,
specificity, PPV, NPV, accuracy with Wilson CIs), and **decision curves** with a net-benefit table
and interventions avoided (Cox: at the horizon, with Kaplan–Meier risks).

## Sample size adequacy

Riley et al. criteria for each outcome type (2019, 2020), from the anticipated Cox–Snell (or
adjusted) R²: entered, or this model's (optimism-corrected when validated).

## Variable importance

- **Shapley decomposition of R²** (LMG / general dominance; exact up to 10 predictors, otherwise
  200 random orders); with baseline covariates, also the gain over the baseline model shared among
  the other predictors.
- **Permutation importance** (loss in AUC, R² or C; increase in Brier score or RMSE), computed in
  the validation data when there is any.
- Exact **SHAP** values on the model's linear scale: summary (beeswarm) plot, mean |SHAP|, and a
  waterfall explaining one patient's prediction.

## Plots

One panel per data set: calibration (smoothed curve with CI, grouped observed values, histogram of
predictions), ROC curves and risk distributions (logistic), Kaplan–Meier curves by risk group
(Cox), decision curves, forest plot, elastic-net cross-validation curve, importance bars, SHAP
summary and waterfall, instability plots.

## Output to the spreadsheet and export

- Saves the predictions (risk and linear predictor; predicted value; risk by the horizon and
  prognostic index) for every row with complete predictors (with the missing-indicator method,
  every row; with multiple imputation, imputed development rows get the mean over imputations).
- **Export model (JSON)**: the fitted model (predictors, coding, spline knots, interactions,
  coefficients with SE, and for Cox the baseline survival at the horizon) as a file that *Validate
  a Published Model* reads, or that can accompany a publication.

## Languages

English, Spanish and Catalan (`jamovi/i18n/*.po`): options, tables, notes, messages and plots.

## Multicore

macOS and Linux, opt-in *Cores*: bootstrap and cross-validation fits, the elastic net's λ tuning,
the Shapley decomposition, permutation importance and the imputation chains run on several cores;
random draws are made first, so results do not depend on the number of cores. Long computations
stop when the options change.

## Validating a published model

*Validate a Published Model* applies a model, unchanged, to the current data and reports the same
performance measures, model updating (logistic and linear), classification, net benefit and plots
as the development analyses. The model comes from a clinpred JSON export, or from a table (csv, tsv,
txt or xlsx) with the columns **term**, **estimate** and, optionally, **SE** (also `EE`), e.g.

| term | estimate | SE |
| --- | --- | --- |
| Intercept | −9.10 | 1.37 |
| glu | 0.037 | 0.006 |
| log(bmi) | 2.31 | 0.98 |
| sex = male | 0.42 | 0.20 |
| age * sex = male | 0.011 | 0.005 |

Terms name variables of the data: `age`, `sex = male` (or `sex: male`, or R's `sexmale`),
`log(x)`, `sqrt(x)`, `exp(x)`, `x^2`, `x/10`, and products with `*`. They are parsed, never run as
code. A Cox model adds the rows `baseline survival` (at the horizon), `horizon` and, when the
published prognostic index is centred, `centre`; without a baseline survival, it is estimated from
the data (discrimination and calibration slope are still valid; calibration-in-the-large is not).
Put the variables the model uses in **Predictor variables** (jamovi gives an analysis only the
variables in its boxes).

## Example data

[`examples/README.md`](../examples/README.md) gives a tour of every feature on a simulated heart
cohort (in jamovi's Data Library), with an external file and published model files to open from the
analyses. In the Data Library too:

- Pima diabetes (`MASS::Pima.tr`), and Pima.tr and Pima.te stacked with a `sample` variable for
  split validation; for external validation load
  [`examples/pima_external.csv`](../examples/pima_external.csv) (`MASS::Pima.te`). `glu` serves as
  a continuous outcome.
- Node-positive breast cancer from the Rotterdam tumour bank (`survival::rotterdam`,
  recurrence-free survival in years), and the same stacked with the German Breast Cancer Study
  Group trial (`survival::gbsg`) as in Royston & Altman (2013); for external validation load
  [`examples/gbsg_external.csv`](../examples/gbsg_external.csv).

## Known limitations

- Model updating is not available with multiple imputation.
- With multiple imputation, performance measures are averaged over the imputations without
  confidence intervals, and the incremental-value comparison is not optimism-corrected.
- Optimism-corrected error measures (ICI, E50, E90, Emax, Brier, RMSE, MAE) are floored at 0: a
  noisy measure such as Emax can have an estimated optimism larger than its apparent value.
- The Cox analyses assess the Brier score at the horizon, not integrated over time.

## Planned

Machine-learning learners (penalised and tree-based models, neural networks) through the same
validation pipeline, with a comparison on identical resamples.
