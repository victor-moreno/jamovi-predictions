# Example data for clinpred

Everything here can be regenerated with `Rscript data-raw/make_examples.R` (fixed seed). The data
sets in jamovi's **Data Library** (☰ ▸ Open ▸ Data Library ▸ Clinical Prediction Models) are the
`.omv` files in [`../data`](../data); the files in this folder are opened from an analysis
(the external validation file, or a published model).

## The heart cohort (simulated): every feature

**`heart_cohort.omv`** (Data Library: *Heart cohort, simulated*) — 2,400 patients, each with
three outcomes:

| outcome | analysis | what it is |
| --- | --- | --- |
| `event10` (no / yes) | Binary Outcome (Logistic) | cardiovascular event within 10 years (≈ 16 %) |
| `sbp1y` | Continuous Outcome (Linear) | systolic blood pressure one year later |
| `time` + `status` (censored / event) | Time to Event (Cox) | follow-up in years, horizon 5 |

Predictors: `age` (non-linear effect: use it as a **spline**), `sex`, `smoking` (never / former /
current), `sbp` (scale it **per 10 units**), `chol`, `bmi` (≈ 3 % missing), `diabetes`; `centre`
(A / B / C) is not a predictor. The true model has an **age × sex interaction** (the excess risk of
men falls with age), so the interaction screening and forward selection find it. `period`
(2012–2016 / 2017–2019) splits the patients for **split validation**; the later period has a
slightly lower baseline risk.

**`heart_external.csv`** and **`heart_external.xlsx`** — 1,000 patients from another population
(older, more smokers, a higher baseline risk and weaker effects), with the problems the variable
check exists for: columns named `Age` and `SEX` (another case), sex levels `Female` / `Male`, 12
patients with smoking `unknown`, systolic pressures above the development range, 8 % missing BMI.

### A tour of the features

1. **Binary Outcome (Logistic)**: Outcome `event10` (event `yes`); covariates `sbp`, `chol`, `bmi`;
   spline `age`; factors `sex`, `smoking`, `diabetes`.
   - *Model*: descriptive table, odds ratios (crude and adjusted), likelihood ratio tests,
     coefficients, equation, sample size. *Covariate Scaling*: `sbp` per 10 units. *Reference
     Levels*: e.g. `smoking` = `current`.
   - *Variable Selection*: add `centre` to the factors (it adjusts for the institution without being
     a predictor of interest), then tick `sex` and `centre` as baseline covariates (always in the model),
     choose backward elimination (AIC, BIC or LR test); the selection table shows what was kept,
     and with bootstrap validation how often each predictor was kept in the resampled fits. Or
     choose *Elastic net (penalised)*: λ by cross-validation (plot), baseline covariates unpenalised,
     the tuning repeated in every bootstrap fit. Tick *Incremental value over the baseline model*
     (baseline `age`, `sex`): ΔAUC with a paired CI, the LR test of the added predictors, NRI,
     IDI and Δ net benefit, apparent, optimism-corrected and in the validation sets.
   - *Interactions*: add `age × sex` (and a false one, e.g. `sbp × diabetes`), keep
     *Forward selection, LR test*: `age × sex` enters at step 1, and the table of effects within
     subgroups shows the sex odds ratio at each age quartile.
   - *Validation*: bootstrap (200) with prediction stability; validation split `period` =
     `2017-2019`; external file `heart_external.csv`; model updating. Expect the original model to
     hold in the later period, and calibration problems in the external population (O/E ≈ 1.26,
     slope ≈ 0.80; the closed testing procedure chooses recalibration). The variable check lists
     the renamed columns, the unknown smoking level and the out-of-range pressures; the case-mix
     table shows the membership c-statistic (≈ 0.54 split, ≈ 0.70 external).
   - *Missing Values*: the exploration shows `bmi` and `chol` missing (≈ 3 % and 2 %); compare
     complete cases with the missing category, single and multiple imputation.
   - *Clinical Utility*, *Variable Importance*, *Plots*, *Save*: all of them;
     *Export model (JSON)* writes a file for step 4.
2. **Continuous Outcome (Linear)**: Outcome `sbp1y`; covariates `sbp`, `bmi`; spline `age`; factors
   `smoking`, `diabetes`; cross-validation; the same split and external file; updating.
3. **Time to Event (Cox)**: Time `time`, Event `status` (`event`), horizon `5`; the predictors of
   step 1; bootstrap; split and external file; Kaplan–Meier by risk group, decision curve,
   calibration at 5 years.
4. **Validate a Published Model**: open `heart_external.csv` as the data set (or the cohort), put
   `Age`, `SEX`, `smoking`, `sbp`, `chol`, `bmi`, `diabetes` in *Predictor variables*, and choose a
   model file below.

### Published models (for *Validate a Published Model*)

Fitted on the 2012–2016 patients, as a paper would report them.

| file | model type | shows |
| --- | --- | --- |
| `published_logistic.csv` | binary | term grammar: `age^2`, `sex = male`, `sbp/10`, `log(chol)`, `age * sex = male`; SE column |
| `published_logistic_semicolon.csv` | binary | European export: `;`, decimal commas, R's term names (`sexmale`), an `EE` column |
| `published_logistic.xlsx` | binary | the same table as an Excel file |
| `published_linear.csv` | continuous (outcome `sbp1y`) | a linear model |
| `published_cox.csv` | time to event | `baseline survival` at the `horizon` (5) for a `centre`d prognostic index |
| `published_cox_nobaseline.csv` | time to event | no baseline: set the horizon to 5; the baseline is estimated from the data |
| `heart_logistic_model.json` | binary | a clinpred export: spline for age, sbp per 10, age × sex |
| `heart_cox_model.json` | time to event | a clinpred export with its baseline survival at 5 years |

## Real data

| file | where | use |
| --- | --- | --- |
| `pima_development.omv`, `pima_split.omv` | Data Library | Pima diabetes (`MASS::Pima.tr`, `Pima.te`): logistic (`type`), linear (`glu`), split by `sample` |
| `pima_external.csv` | here | external validation of a Pima model |
| `rotterdam_development.omv`, `breast_split.omv` | Data Library | breast cancer recurrence-free survival (`survival::rotterdam`, `gbsg`), split by `cohort` (Royston & Altman 2013) |
| `gbsg_external.csv` | here | external validation of a Rotterdam Cox model |
