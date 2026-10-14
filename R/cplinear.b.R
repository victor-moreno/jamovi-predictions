# Continuous outcome prediction model (linear regression): the shared
# orchestration (.cp_private in cp_analysis.R) plus the hooks specific to a
# continuous outcome.

cpLinearClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "cpLinearClass",
    inherit = cpLinearBase,
    private = utils::modifyList(.cp_private, list(
        .famName = "gaussian",
        .outcomeVars = function() self$options$dep,
        .updKeys = c("r2", "rmse"),

        .intro = function() paste0(
            "<p>", .("Develop and validate a prediction model for a continuous outcome (linear regression)."), "</p><ul>",
            "<li>", .("Put the outcome in <b>Outcome</b>."), "</li>",
            "<li>", .("Add predictors: numeric ones as linear covariates or as splines (non-linear), categorical ones as factors."), "</li>",
            "<li>", .("Under <b>Validation</b>, choose bootstrap or cross-validation; optionally a split variable marking validation rows, or an external data file with the same variable names."), "</li></ul>"),

        .setOutcome = function(prep) cp_gaussian_outcome(prep),

        .yOf = function(df, spec, key) {
            y <- as.numeric(df[[spec$outcome$name]])
            if (length(y) < 3 || stats::sd(y) == 0)
                .fail(sprintf(.("%s: the validation data need at least 3 complete rows with a varying outcome (there are %s)."), private$.sets()[[key]], length(y)))
            y
        },

        .sampleValues = function(y, nMissing) list(n = length(y), mean = mean(y), sd = stats::sd(y), excluded = nMissing),

        .outcomeSummary = function(y) sprintf("%.2f (%.2f)", mean(y), stats::sd(y)),

        .text = function(key) switch(key,
            separation = .("Some predicted risks are numerically 0 or 1 (quasi-complete separation): coefficients and their intervals may be unreliable."),
            outcomeRow = self$options$dep,
            outcomeLevel = .("mean (SD)"),
            memberNote = .("Membership c-statistic: how well the predictors and the outcome distinguish validation from development rows (Debray et al. 2015). Near 0.5: similar case-mix (reproducibility); high: different populations (transportability). A larger SD of the predictions means a more heterogeneous case-mix, which by itself raises the R²."),
            samplesEx = .("Rows with a missing outcome or predictor are excluded (complete-case analysis)."),
            coefSpline = .("Spline terms (', '') are the non-linear parts of a restricted cubic spline; their coefficients have no direct interpretation (see the effects, Q3 vs Q1, or SHAP)."),
            orSpline = .("Spline predictors: difference in the mean outcome between the third and the first quartile; the p is the F test of the whole (non-linear) association."),
            orRef = .("The p on a factor's reference row is the F test of the factor."),
            orInt = .("With interactions, an adjusted effect is the effect at the reference level (factor) or median (covariate) of the variables it interacts with; see Effects Within Subgroups."),
            subEx = .("Difference in the mean outcome for each variable of an interaction at the levels (factor) or quartiles (covariate) of the other; other predictors at their reference level or median."),
            cvNote = .("Cross-validation: %d folds × %d repeats, stratified by outcome quintiles; performance of the out-of-fold predictions, averaged over repeats."),
            negNote = .("A negative corrected R² means that, in new patients, the model would predict worse than the mean outcome: it is overfitted."),
            perfEx = .("R²: 1 − residual sum of squares / total sum of squares around the mean of the data in which it is measured (can be negative outside the development data). RMSE: root mean squared error; MAE: mean absolute error, in the units of the outcome (lower is better). Calibration-in-the-large: mean of observed − predicted, 0 is ideal. Calibration slope: of the regression of the observed on the predicted outcome, 1 is ideal; below 1 the predictions are too extreme (overfitting)."),
            valEx = .("The development model is applied unchanged to the validation data. A calibration-in-the-large far from 0 means the mean outcome differs between the populations; a slope below 1, that the predictions are too extreme there."),
            updNote = .("Recalibration re-estimates the intercept and the slope of the original prediction; intercept and slope: of the original prediction in the updated model. Each test compares a method with model revision (F test); the closed testing procedure (Vergouwe et al. 2017) keeps the first method, from the original model up, not significantly worse than revision at α = %s (✓). Performance of an updated model is measured in the data used to update it, so it is optimistic, most of all for revision."),
            stabEx = .("A shrinkage factor well below 1 (e.g. < 0.9) signals overfitting: the coefficients would need to be shrunk towards 0. MAPE: average absolute difference between an individual's prediction from the developed model and from models developed in bootstrap samples (Riley & Collins 2023), in the units of the outcome."),
            shapEx = .("SHAP: the contribution of each predictor to an individual's predicted outcome, relative to the average individual of the development data (exact for this additive model), in the units of the outcome."),
            shapDown = .("Lowers the prediction"),
            shapUp = .("Raises the prediction"),
            forestAxis = .("Difference in the mean outcome"),
            permAxis = .("Loss in R² when permuted (mean ± SD)"),
            shapAxis = .("SHAP value (contribution to the predicted outcome)"),
            instabX = .("Prediction from the developed model"),
            instabY = .("Prediction from bootstrap models"),
            calX = .("Predicted value"),
            calY = .("Observed value (smoothed)")),

        .shapAxis = function(base, lp) c(sprintf(.("Predicted %s"), self$options$dep),
                                         sprintf(.("predicted %.4g"), lp)),

        .equation = function(beta, spec) {
            out <- c(sprintf(.("predicted %s ="), spec$outcome$name), private$.equationTerms(beta))
            paste(c(out, private$.equationVars(spec)), collapse = "\n")
        },

        # descriptive table: each predictor's distribution and its crude
        # association with the outcome (correlation; outcome mean by level)
        .initDesc = function() {
            if (!self$options$descTable) return()
            t <- self$results$descTable
            t$addRow(rowKey = "n", values = list(var = .("Patients"), level = "n"))
            stat <- if (self$options$descStat == "meanSD") .("mean (SD)") else paste(.("median"), "[Q1, Q3]")
            for (v in c(self$options$covs, self$options$splineVars))
                t$addRow(rowKey = v, values = list(var = v, level = stat))
            for (v in self$options$factors)
                for (l in private$.levelsOf(v))
                    t$addRow(rowKey = paste(v, l), values = list(var = v, level = l))
        },

        .fillDesc = function(dd) {
            o <- self$options
            y <- jmvcore::toNumeric(dd[[o$dep]])
            dd$.grp <- "a"
            rows <- cp_describe(dd, ".grp", "a", c(o$covs, o$splineVars), o$factors, o$descStat)
            t <- self$results$descTable
            for (rw in rows) {
                v <- rw$var
                assoc <- if (rw$key == "n") sprintf("%.2f (%.2f)", mean(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE))
                    else if (v %in% o$factors) {
                        yy <- y[!is.na(dd[[v]]) & as.character(dd[[v]]) == rw$level & !is.na(y)]
                        if (length(yy)) sprintf("%.2f (%.2f)", mean(yy), if (length(yy) > 1) stats::sd(yy) else NA) else ""
                    } else sprintf("r = %.3f", stats::cor(dd[[v]], y, use = "pairwise.complete.obs"))
                private$.putRow(t, rw$key, list(var = rw$var, level = rw$level, total = rw$values$total,
                                                outcome = assoc, missing = rw$missing))
            }
            t$setNote("dev", .("Development rows, before excluding missing values. Outcome: its mean (SD) overall and in each factor level, and its Pearson correlation with each numeric predictor."))
        },

        .initExtras = function(sets) NULL,
        .fillExtras = function(all, sets, conf) NULL,
        .extraPlotStates = function(all, sets, conf) NULL,

        .ssRows = function() c(shrinkage = .("1. Shrinkage factor ≥ 0.9"),
                               optimism = .("2. Optimism in R² ≤ 0.05 (apparent − adjusted)"),
                               residual = .("3. Residual SD within × 1.1"),
                               intercept = .("4. Mean outcome within × 1.1"),
                               required = .("Minimum required"),
                               actual = .("This development data")),

        .fillSampleSize = function(y, lp, spec, iv) {
            o <- self$options
            p <- cp_nparams(spec); n <- length(y)
            # entered, or this model's out-of-sample R² (optimism-corrected, else adjusted)
            r2 <- if (o$ssSource == "value") o$ssR2 else
                if (!is.null(iv)) iv$corrected[["r2"]] else
                    1 - (1 - private$.fam$metrics(y, lp)[["r2"]]) * (n - 1) / (n - p - 1)
            if (!is.finite(r2) || r2 <= 0) .fail(.("The anticipated adjusted R² must be positive; enter a value under Sample size adequacy."))
            ss <- cp_gaussian_samplesize(p, r2, mean(y), stats::sd(y))
            t <- self$results$ssTable
            for (k in names(ss$n)) t$setRow(rowKey = k, values = list(n = ss$n[[k]], spp = ss$n[[k]] / p))
            req <- max(ss$n)
            t$setRow(rowKey = "required", values = list(n = req, spp = req / p))
            t$setRow(rowKey = "actual", values = list(n = n, spp = n / p))
            t$setNote("in", sprintf(.("Candidate parameters (excluding the intercept): %d; anticipated adjusted R²: %.3f (%s); outcome mean %.4g, SD %.4g."),
                                    p, r2, if (o$ssSource == "value") .("entered") else if (!is.null(iv)) .("optimism-corrected, this model") else .("adjusted, this model"),
                                    mean(y), stats::sd(y)))
            t$setNote("verdict", if (n >= req) .("The development data meet the minimum sample size.")
                      else sprintf(.("The development data are %d short of the minimum sample size: expect overfitting; consider fewer parameters or penalisation."), req - n))
        },

        .calState = function(all, sets, conf) lapply(names(all), function(k)
            c(cp_gaussian_calcurve(all[[k]]$y, all[[k]]$p, conf), list(label = sets[[k]]))),

        .outputDefs = function() list(predValue = function(lp, beta) lp),

        .plotCal = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_cal_cont(image$state, ggtheme)); TRUE
        }
    ))
)
