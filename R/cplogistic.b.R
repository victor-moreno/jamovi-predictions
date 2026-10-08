# Binary outcome prediction model: the shared orchestration (.cp_private in
# cp_analysis.R) plus the hooks specific to a binary outcome: event coding,
# classification, decision curves, ROC and risk plots, Riley's sample size.

cpLogisticClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "cpLogisticClass",
    inherit = cpLogisticBase,
    private = utils::modifyList(.cp_private, list(
        .famName = "binary",
        .outcomeVars = function() self$options$dep,
        .updKeys = c("auc", "brier", "ici"),

        .intro = function() paste0(
            "<p>", .("Develop and validate a prediction model for a binary outcome."), "</p><ul>",
            "<li>", .("Put the outcome in <b>Outcome</b> and choose the event level."), "</li>",
            "<li>", .("Add predictors: numeric ones as linear covariates or as splines (non-linear), categorical ones as factors."), "</li>",
            "<li>", .("Under <b>Validation</b>, choose bootstrap or cross-validation; optionally a split variable marking validation rows, or an external data file with the same variable names."), "</li></ul>"),

        .setOutcome = function(prep) cp_binary_outcome(prep, self$options$eventLevel),

        .yOf = function(df, spec, key) {
            y <- as.integer(df[[spec$outcome$name]] == spec$outcome$event)
            if (sum(y) < 1 || sum(1 - y) < 1)
                .fail(private$.sets()[[key]], ": the validation data need both events and non-events (there are ",
                      sum(y), " events among ", length(y), " complete rows).")
            y
        },

        .sampleValues = function(y, nMissing) list(n = length(y), events = sum(y), pct = mean(y), excluded = nMissing),

        .outcomeSummary = function(y) sprintf("%d (%.1f%%)", sum(y), 100 * mean(y)),

        .eventLabel = function() private$.level(self$options$dep, self$options$eventLevel),

        .text = function(key) switch(key,
            separation = .("Some predicted risks are numerically 0 or 1 (quasi-complete separation): coefficients and their intervals may be unreliable."),
            outcomeRow = .("Events"),
            outcomeLevel = "n (%)",
            memberNote = .("Membership c-statistic: how well the predictors and the outcome distinguish validation from development rows (Debray et al. 2015). Near 0.5: similar case-mix (reproducibility); high: different populations (transportability). A larger SD of the linear predictor means a more heterogeneous case-mix, which by itself raises the AUC."),
            samplesEx = paste0(.("Event"), ": ", self$options$dep, " = ", private$.eventLabel(), ". ",
                               .("Rows with a missing outcome or predictor are excluded (complete-case analysis).")),
            coefSpline = .("Spline terms (', '') are the non-linear parts of a restricted cubic spline; their coefficients have no direct interpretation (see the odds ratios, Q3 vs Q1, or SHAP)."),
            orSpline = .("Spline predictors: odds ratio between the third and the first quartile; the p is the likelihood ratio test of the whole (non-linear) association."),
            orRef = .("The p on a factor's reference row is the likelihood ratio test of the factor."),
            orInt = .("With interactions, an adjusted odds ratio is the effect at the reference level (factor) or median (covariate) of the variables it interacts with; see Effects Within Subgroups."),
            subEx = .("Odds ratio of each variable of an interaction at the levels (factor) or quartiles (covariate) of the other; other predictors at their reference level or median."),
            cvNote = .("Cross-validation: %d folds × %d repeats, stratified by outcome; performance of the out-of-fold predictions, averaged over repeats."),
            negNote = .("A negative corrected R² or scaled Brier score means that, in new patients, the model would predict worse than giving everyone the overall risk: it is overfitted."),
            perfEx = .("AUC: probability that a random event has a higher predicted risk than a random non-event (0.5 = no discrimination). Calibration-in-the-large: 0 is ideal (log-odds). Calibration slope: 1 is ideal; below 1 the predictions are too extreme (overfitting). O/E: observed / expected events. ICI, E50, E90, Emax: mean, median, 90th percentile and maximum absolute difference between the predicted and the (smoothed) observed risk. Brier: mean squared error of the predicted risks (lower is better); scaled Brier: 1 − Brier / Brier of a model without predictors."),
            valEx = .("The development model is applied unchanged to the validation data. A calibration-in-the-large far from 0 means the overall risk differs between the populations; a slope below 1, that the predictions are too extreme there."),
            updNote = .("Recalibration re-estimates the intercept and the slope of the original prediction; intercept and slope: of the original linear predictor (log-odds) in the updated model. Each test compares a method with model revision (likelihood ratio); the closed testing procedure (Vergouwe et al. 2017) keeps the first method, from the original model up, not significantly worse than revision at α = %s (✓). Performance of an updated model is measured in the data used to update it, so it is optimistic, most of all for revision."),
            stabEx = .("A shrinkage factor well below 1 (e.g. < 0.9) signals overfitting: the coefficients would need to be shrunk towards 0. MAPE: average absolute difference between an individual's risk from the developed model and from models developed in bootstrap samples (Riley & Collins 2023)."),
            shapEx = .("SHAP: the contribution of each predictor to an individual's log-odds, relative to the average individual of the development data (exact for this additive model)."),
            shapDown = .("Lowers risk"),
            shapUp = .("Raises risk"),
            permAxis = .("Loss in AUC when permuted (mean ± SD)"),
            shapAxis = .("SHAP value (contribution to the log-odds)"),
            instabX = .("Risk from the developed model"),
            instabY = .("Risk from bootstrap models"),
            calX = .("Predicted risk"),
            calY = .("Observed proportion (smoothed)")),

        .shapAxis = function(base, lp) c("Log-odds", sprintf("risk %.1f%%", 100 * stats::plogis(lp))),

        .equation = function(beta, spec) {
            out <- c(sprintf("log-odds(%s = %s) =", spec$outcome$name, spec$outcome$event),
                     private$.equationTerms(beta), "", "risk = 1 / (1 + exp(−log-odds))")
            paste(c(out, private$.equationVars(spec)), collapse = "\n")
        },

        # descriptive table: one column per outcome group (event last)
        .initDesc = function() {
            if (!self$options$descTable) return()
            t <- self$results$descTable
            dep <- self$options$dep
            lev <- private$.levelsOf(dep)
            ev <- private$.level(dep, self$options$eventLevel)
            other <- setdiff(lev, ev)
            t$getColumn("g1")$setTitle(if (length(other) == 1) other else .("Other levels"))
            t$getColumn("g2")$setTitle(ev)
            t$addRow(rowKey = "n", values = list(var = .("Patients"), level = "n"))
            stat <- if (self$options$descStat == "meanSD") .("mean (SD)") else .("median [Q1, Q3]")
            for (v in c(self$options$covs, self$options$splineVars))
                t$addRow(rowKey = v, values = list(var = v, level = stat))
            for (v in self$options$factors)
                for (l in private$.levelsOf(v))
                    t$addRow(rowKey = paste(v, l), values = list(var = v, level = l))
        },

        .fillDesc = function(dd) {
            o <- self$options
            ev <- private$.level(o$dep, o$eventLevel)
            yv <- as.character(dd[[o$dep]])
            dd$.grp <- ifelse(is.na(yv), NA, ifelse(yv == ev, "e", "n"))
            rows <- cp_describe(dd, ".grp", c("n", "e"), c(o$covs, o$splineVars), o$factors, o$descStat)
            t <- self$results$descTable
            for (rw in rows)
                private$.putRow(t, rw$key, list(var = rw$var, level = rw$level, g1 = rw$values$g1,
                                                g2 = rw$values$g2, total = rw$values$total, missing = rw$missing))
            t$setNote("dev", .("Development rows, before excluding missing values; percentages of the non-missing values."))
        },

        .initExtras = function(sets) {
            r <- self$results
            if (self$options$classif) {
                meas <- c(sens = .("Sensitivity"), spec = .("Specificity"), ppv = .("PPV"),
                          npv = .("NPV"), acc = .("Accuracy"))
                for (k in names(sets)) {
                    r$confTable$addRow(rowKey = k, values = list(set = sets[[k]]))
                    for (m in names(meas))
                        r$classTable$addRow(rowKey = paste(k, m), values = list(set = sets[[k]], measure = meas[[m]]))
                }
            }
            if (self$options$dca) {
                th <- tryCatch(.parse_numbers(self$options$nbThresholds), error = function(e) numeric(0))
                th <- th[th > 0 & th < 1]
                for (k in names(sets)) for (t in th)
                    r$nbTable$addRow(rowKey = paste(k, t), values = list(set = sets[[k]], threshold = t))
            }
        },

        .fillExtras = function(all, sets, conf) {
            if (self$options$classif) private$.fillClassif(all, sets, conf)
            if (self$options$dca) private$.fillDca(all, sets)
        },

        .fillClassif = function(all, sets, conf) {
            o <- self$options
            thr <- cp_binary_threshold(all$dev$y, all$dev$p, o$thrMethod, o$threshold)
            for (k in names(all)) {
                cl <- cp_binary_classify(all[[k]]$y, all[[k]]$p, thr, conf)
                self$results$confTable$setRow(rowKey = k, values = list(threshold = thr, tp = cl$tp, fp = cl$fp,
                                                                         fn = cl$fn, tn = cl$tn))
                for (m in c("sens", "spec", "ppv", "npv", "acc"))
                    self$results$classTable$setRow(rowKey = paste(k, m), values = list(
                        est = cl[[m]][1], lower = cl[[m]][2], upper = cl[[m]][3]))
            }
            if (o$thrMethod != "value")
                self$results$confTable$setNote("thr", .("The threshold is derived from the development data and applied unchanged to the validation data."))
        },

        .fillDca = function(all, sets) {
            th <- .parse_numbers(self$options$nbThresholds)
            th <- th[th > 0 & th < 1]
            for (k in names(all)) {
                nb <- cp_net_benefit(all[[k]]$y, all[[k]]$p, th)
                for (i in seq_along(th)) {
                    best <- max(nb$all[i], 0)
                    self$results$nbTable$setRow(rowKey = paste(k, th[i]), values = list(
                        model = nb$model[i], all = nb$all[i], gain = nb$model[i] - best,
                        avoided = 100 * (nb$model[i] - nb$all[i]) / (th[i] / (1 - th[i]))))
                }
            }
            if (self$options$explain)
                self$results$nbTable$setNote("ex", .("Net benefit = true positives/N − false positives/N × t/(1 − t), at threshold t. Gain: over the better of treating all and treating none. Interventions avoided per 100: fewer false positives than treating all, at the same benefit."))
        },

        .ssRows = function() c(shrinkage = .("1. Shrinkage factor ≥ 0.9"),
                               optimism = .("2. Optimism in Nagelkerke R² ≤ 0.05"),
                               intercept = .("3. Overall risk within ± 0.05"),
                               required = .("Minimum required"),
                               actual = .("This development data")),

        .fillSampleSize = function(y, lp, spec, iv) {
            o <- self$options
            p <- cp_nparams(spec); phi <- mean(y)
            # entered, or this model's (optimism-corrected when validated)
            r2 <- if (o$ssSource == "value") o$ssR2 else
                if (!is.null(iv)) iv$corrected[["r2cs"]] else private$.fam$metrics(y, lp)[["r2cs"]]
            if (!is.finite(r2) || r2 <= 0) .fail("The anticipated Cox–Snell R² must be positive; enter a value under Sample size adequacy.")
            ss <- cp_binary_samplesize(p, phi, r2)
            t <- self$results$ssTable
            for (k in names(ss$n))
                t$setRow(rowKey = k, values = list(n = ss$n[[k]], events = ss$n[[k]] * phi, epp = ss$n[[k]] * phi / p))
            req <- max(ss$n)
            t$setRow(rowKey = "required", values = list(n = req, events = req * phi, epp = req * phi / p))
            t$setRow(rowKey = "actual", values = list(n = length(y), events = sum(y), epp = sum(y) / p))
            t$setNote("in", sprintf(.("Parameters (excluding the intercept): %d; outcome proportion: %.3f; anticipated Cox–Snell R²: %.3f (%s)."),
                                    p, phi, r2, if (o$ssSource == "value") .("entered") else if (!is.null(iv)) .("optimism-corrected, this model") else .("apparent, this model")))
            t$setNote("verdict", if (length(y) >= req) .("The development data meet the minimum sample size.")
                      else sprintf(.("The development data are %d short of the minimum sample size: expect overfitting; consider fewer parameters or penalisation."), req - length(y)))
        },

        .calState = function(all, sets, conf) lapply(names(all), function(k)
            c(cp_binary_calcurve(all[[k]]$y, all[[k]]$p, self$options$calGroups, conf), list(label = sets[[k]]))),

        .extraPlotStates = function(all, sets, conf) {
            pl <- self$results$plots
            o <- self$options
            if (o$rocPlot && is.null(pl$roc$state))
                pl$roc$setState(lapply(names(all), function(k) {
                    rp <- cp_roc_points(all[[k]]$y, all[[k]]$p)
                    c(rp, list(key = k, label = sets[[k]], auc = .auc_delong(all[[k]]$y, all[[k]]$p)$auc))
                }))
            if (o$riskPlot && is.null(pl$riskDist$state)) {
                lev <- c(paste0("Not ", private$.eventLabel()), private$.eventLabel())
                pl$riskDist$setState(lapply(names(all), function(k)
                    list(label = sets[[k]], p = round(all[[k]]$p, 4),
                         outcome = factor(lev[all[[k]]$y + 1], lev))))
            }
            if (o$dca && is.null(pl$dcaCurves$state)) {
                th <- seq(0.01, o$dcaMax, by = 0.01)
                pl$dcaCurves$setState(lapply(names(all), function(k)
                    list(label = sets[[k]], nb = cp_net_benefit(all[[k]]$y, all[[k]]$p, th))))
            }
        },

        .outputDefs = function() list(predRisk = function(lp, beta) stats::plogis(lp), predLP = function(lp, beta) lp),

        .plotRoc = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_roc(image$state, ggtheme)); TRUE
        },
        .plotCal = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_cal(image$state, ggtheme)); TRUE
        },
        .plotRisk = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_risk(image$state, ggtheme)); TRUE
        },
        .plotDca = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_dca(image$state, ggtheme)); TRUE
        }
    ))
)
