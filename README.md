# clinpred — clinical prediction models for jamovi

A [jamovi](https://www.jamovi.org) module to **develop and validate** clinical prediction models,
following the clinical-epidemiology workflow (Steyerberg; TRIPOD+AI 2024; Riley et al.). A
clinician can fit a model, see how optimistic its apparent performance is, validate it in new
patients, update it, and see which predictors drive the predictions.

| menu: **Prediction ▸** | outcome | model |
| --- | --- | --- |
| **Binary Outcome (Logistic)** | event / no event | logistic regression |
| **Continuous Outcome (Linear)** | a measurement | linear regression |
| **Time to Event (Cox)** | time + event, risk by a horizon | Cox regression |
| **Validate a Published Model** | any of the three | read from a file (clinpred JSON or a coefficient table) |

Needs **jamovi 28.4 or later**.

## What it does

- **Predictors**: linear covariates (with scaling), restricted cubic splines, factors with a
  chosen reference level; interactions, listed or screened among all pairs.
- **Selection**: baseline covariates always in the model; backward or forward selection, or an
  elastic net (optionally relaxed), repeated in every resampled fit; bootstrap coefficient tables
  for inference after selection; incremental value over the baseline model.
- **Missing data**: exploration, complete cases, missing indicator, single or multiple imputation
  by chained equations.
- **Performance**: discrimination, calibration (in the large, slope, O/E, ICI), overall measures
  and clinical utility (classification, decision curves), with confidence intervals.
- **Validation**: bootstrap optimism correction or repeated cross-validation; validation in a
  split of the data or an external file (with variable and case-mix checks); model updating with
  the closed testing procedure (sequential tests for Cox models); stability of individual
  predictions; sample-size adequacy.
- **Explanation**: crude and adjusted effects with a forest plot, likelihood ratio tests, effects
  within subgroups, the model equation, Shapley decomposition of R², permutation importance and
  SHAP.
- **Output**: predictions saved to the spreadsheet, and the fitted model exported as JSON.

Every feature, option and method is described in [`documents/features.md`](documents/features.md).

## Tutorial

A [tutorial website](https://victor-moreno.github.io/jamovi-predictions/) (in `docs/`) introduces
prediction models — overfitting, performance measures, internal and external validation,
updating — and works through every analysis on the example data: a logistic model, selection and
incremental value, validation in new patients, missing data, a Cox model, a linear model, published
models and variable importance. Every table and plot in it is computed by the module itself.

## Languages

The options, tables, notes, messages and plots are available in English, Spanish and Catalan;
jamovi uses the language chosen in its settings. The translations are in `jamovi/i18n/*.po`.

## Example data

The Data Library (☰ ▸ Open ▸ Data Library ▸ Clinical Prediction Models) holds a simulated heart
cohort with a binary, a continuous and a time-to-event outcome, Pima diabetes and Rotterdam/GBSG
breast cancer data. [`examples/README.md`](examples/README.md) is a tour of every feature on them,
with an external validation file and published model files to open from the analyses.

## Install

Download the `.jmo` for your CPU from the [Releases](../../releases) page and, in jamovi, choose
**Modules (⊕) ▸ Sideload**:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (28.4 or later, R 4.6.0) | `clinpred_<version>_current_R4.6.0_arm64.jmo` | `clinpred_<version>_current_R4.6.0_x64.jmo` |

The same file works on macOS, Windows and Linux. jamovi *solid* is not supported: the external
validation file uses jamovi's File option, which appeared in 28.4.

## Build and test (developers)

```sh
bash tools/install.sh desktop        # install into jamovi.app (tools/smoke.R checks the installed module)
bash tools/build-jmo.sh current      # release builds: dist/clinpred_<version>_current_R4.6.0_{arm64,x64}.jmo
bash tools/release.sh                # build and publish one GitHub release
sh tools/tutorial_redo.sh            # render the tutorial (Quarto source in tutorial/, local) into docs/
R CMD INSTALL -l .tmp/lib . && Rscript -e '.libPaths(c(".tmp/lib", .libPaths())); testthat::test_dir("tests/testthat", package = "clinpred", load_package = "installed")'
```

The tests check the results against rms, pROC, Hmisc, survival, timeROC, glmnet, mice, relaimpo
and pmsampsize when those packages are installed (they are skipped otherwise, and kept out of
`Suggests` so that the build does not bundle them).

## Acknowledgements

The methods follow E. W. Steyerberg, *Clinical Prediction Models* (2019), F. E. Harrell's
bootstrap validation (as in the R package [rms](https://cran.r-project.org/package=rms)), and
R. D. Riley and colleagues' work on sample size and stability. The results are checked against
[rms](https://cran.r-project.org/package=rms), [pROC](https://cran.r-project.org/package=pROC),
[Hmisc](https://cran.r-project.org/package=Hmisc),
[survival](https://cran.r-project.org/package=survival) (T. Therneau),
[timeROC](https://cran.r-project.org/package=timeROC), [glmnet](https://cran.r-project.org/package=glmnet),
[mice](https://cran.r-project.org/package=mice),
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
