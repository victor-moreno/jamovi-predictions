# Time-to-event prediction model (Cox): the shared orchestration (.cp_private
# in cp_analysis.R) plus the hooks specific to censored outcomes: time and
# event, the prediction horizon, Kaplan–Meier based calibration, net benefit
# and risk groups.

cpCoxClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "cpCoxClass",
    inherit = cpCoxBase,
    private = utils::modifyList(.cp_private, list(
        .famName = "cox",
        .outcomeNames = c("time", "status"),
        .outcomeVars = function() c(self$options$time, self$options$status),
        # NA, not NULL: modifyList drops NULL entries, and R6 cannot add the field later
        .horizon = NA_real_,
        .updKeys = c("c", "oe", "brier"),

        .intro = function() paste0(
            "<p>", .("Develop and validate a prediction model for the risk of an event by a time horizon (Cox regression)."), "</p><ul>",
            "<li>", .("Put the follow-up time in <b>Time</b> and the event indicator in <b>Event</b>, and choose the level that marks an event."), "</li>",
            "<li>", .("Set the <b>prediction horizon</b> in the units of the time variable (0 uses the median follow-up)."), "</li>",
            "<li>", .("Add predictors: numeric ones as linear covariates or as splines (non-linear), categorical ones as factors."), "</li>",
            "<li>", .("Under <b>Validation</b>, choose bootstrap or cross-validation; optionally a split variable marking validation rows, or an external data file with the same variable names."), "</li></ul>"),

        .setOutcome = function(prep) {
            o <- self$options
            prep <- cp_cox_outcome(prep, o$time, o$status, o$eventLevel, o$horizon)
            private$.horizon <- prep$spec$outcome$horizon
            prep
        },

        .makeFamily = function(prep) cp_family("cox", horizon = prep$spec$outcome$horizon),

        .yOf = function(df, spec, key) {
            tm <- as.numeric(df[[spec$outcome$time]])
            st <- as.integer(df[[spec$outcome$name]] == spec$outcome$event)
            h <- spec$outcome$horizon
            if (!any(st == 1 & tm <= h))
                .fail(private$.sets()[[key]], ": no event happens before the horizon (", .num(h), ") in the validation data.")
            if (!any(tm > h))
                .fail(private$.sets()[[key]], ": nobody is followed beyond the horizon (", .num(h), ") in the validation data.")
            cbind(time = tm, status = st)
        },

        .sampleValues = function(y, nMissing) list(n = nrow(y), events = sum(y[, 2]), pct = mean(y[, 2]),
                                                   km = 1 - .km(y[, 1], y[, 2], private$.horizon)$surv,
                                                   excluded = nMissing),

        .outcomeSummary = function(y) sprintf("%d (%.1f%%)", sum(y[, 2]), 100 * mean(y[, 2])),

        .eventLabel = function() private$.level(self$options$status, self$options$eventLevel),

        .hText = function() .num(private$.horizon, 4),

        .text = function(key) switch(key,
            separation = .("The partial likelihood has no maximum for some coefficient (monotone likelihood): a predictor, or a combination of them, separates patients with and without events. Coefficients and their intervals are unreliable."),
            outcomeRow = .("Events"),
            outcomeLevel = "n (%)",
            memberNote = .("Membership c-statistic: how well the predictors and the outcome distinguish validation from development rows (Debray et al. 2015). Near 0.5: similar case-mix (reproducibility); high: different populations (transportability). A larger SD of the linear predictor means a more heterogeneous case-mix, which by itself raises the concordance."),
            samplesEx = paste0(.("Event"), ": ", self$options$status, " = ", private$.eventLabel(), ". ",
                               sprintf(.("Risks are predicted by the horizon t = %s%s."), private$.hText(),
                                       if (isTRUE(self$options$horizon <= 0)) .(" (median follow-up; set a clinically relevant horizon)") else ""), " ",
                               .("Rows with a missing time, event or predictor are excluded (complete-case analysis).")),
            coefSpline = .("Spline terms (', '') are the non-linear parts of a restricted cubic spline; their coefficients have no direct interpretation (see the hazard ratios, Q3 vs Q1, or SHAP)."),
            orSpline = .("Spline predictors: hazard ratio between the third and the first quartile; the p is the likelihood ratio test of the whole (non-linear) association."),
            orRef = .("The p on a factor's reference row is the likelihood ratio test of the factor."),
            orInt = .("With interactions, an adjusted hazard ratio is the effect at the reference level (factor) or median (covariate) of the variables it interacts with; see Effects Within Subgroups."),
            subEx = .("Hazard ratio of each variable of an interaction at the levels (factor) or quartiles (covariate) of the other; other predictors at their reference level or median."),
            cvNote = .("Cross-validation: %d folds × %d repeats, stratified by event; performance of the out-of-fold predictions, averaged over repeats."),
            negNote = .("A negative corrected R² or scaled Brier score means that, in new patients, the model would predict worse than giving everyone the overall risk: it is overfitted."),
            perfEx = sprintf(.("Harrell's C: probability that, of two comparable patients, the one with the event first has the higher predicted risk. Uno's C: the same, weighted for censoring and restricted to the horizon. Time-dependent AUC: discrimination between events by t = %s and patients event-free at t (censoring weighted). Royston's D: separation in log hazard between the halves of the prognostic index. O/E: observed (Kaplan–Meier) / mean predicted risk at t. Calibration slope: 1 is ideal; below 1 the predictions are too extreme. ICI, E50, E90: mean, median and 90th percentile absolute difference between predicted and observed risk at t (smoothed by a Cox model on a spline of the predicted risk). Brier: censoring-weighted mean squared error at t; scaled Brier (index of prediction accuracy): 1 − Brier / Brier of the Kaplan–Meier risk."), private$.hText()),
            valEx = .("The development model, with its baseline hazard, is applied unchanged to the validation data. An O/E far from 1 means the overall risk at the horizon differs between the populations; a slope below 1, that the predictions are too extreme there."),
            updNote = .("Recalibration in the large re-estimates the baseline risk at the horizon (the original linear predictor as offset); recalibration also its slope; revision refits every coefficient. Shift and slope: of the original linear predictor in the updated model. The partial likelihood does not involve the baseline risk, so the closed testing procedure is adapted into sequential tests (each p in the row of the more complex method): revision vs recalibration and slope = 1 by partial likelihood ratio, O/E at the horizon = 1 by a Wald test; the most complex method whose test is significant at α = %s is chosen (✓). Performance of an updated model is measured in the data used to update it, so it is optimistic, most of all for revision."),
            stabEx = .("A shrinkage factor well below 1 (e.g. < 0.9) signals overfitting: the coefficients would need to be shrunk towards 0. MAPE: average absolute difference between an individual's risk from the developed model and from models developed in bootstrap samples (Riley & Collins 2023)."),
            shapEx = .("SHAP: the contribution of each predictor to an individual's log cumulative hazard at the horizon (equivalently, log relative hazard), relative to the average individual of the development data (exact for this additive model)."),
            shapDown = .("Lowers risk"),
            shapUp = .("Raises risk"),
            forestAxis = .("Hazard ratio (log scale)"),
            permAxis = .("Loss in Harrell's C when permuted (mean ± SD)"),
            shapAxis = .("SHAP value (contribution to the log hazard)"),
            instabX = .("Risk from the developed model"),
            instabY = .("Risk from bootstrap models"),
            calX = .("Predicted risk"),
            calY = .("Observed risk (smoothed)")),

        .shapAxis = function(base, lp) c(sprintf("Log cumulative hazard at t = %s", private$.hText()),
                                         sprintf("risk %.1f%% by t = %s", 100 * .cox_risk(lp), private$.hText())),

        .equation = function(beta, spec) {
            out <- c("Prognostic index PI =", private$.equationTerms(beta)[-1], "",
                     sprintf("Risk of %s = %s by t = %s:  1 − S0^exp(PI),  S0 = %.6g",
                             spec$outcome$name, spec$outcome$event, .num(spec$outcome$horizon, 5),
                             exp(-exp(beta[[1]]))),
                     "(S0: event-free probability at t for PI = 0, i.e. every covariate 0 and factors at their reference level)")
            paste(c(out, private$.equationVars(spec)), collapse = "\n")
        },

        # descriptive table: censored and event columns
        .initDesc = function() {
            if (!self$options$descTable) return()
            t <- self$results$descTable
            st <- self$options$status
            lev <- private$.levelsOf(st)
            ev <- private$.level(st, self$options$eventLevel)
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
            ev <- private$.level(o$status, o$eventLevel)
            sv <- as.character(dd[[o$status]])
            dd$.grp <- ifelse(is.na(sv), NA, ifelse(sv == ev, "e", "n"))
            rows <- cp_describe(dd, ".grp", c("n", "e"), c(o$covs, o$splineVars), o$factors, o$descStat)
            t <- self$results$descTable
            for (rw in rows)
                private$.putRow(t, rw$key, list(var = rw$var, level = rw$level, g1 = rw$values$g1,
                                                g2 = rw$values$g2, total = rw$values$total, missing = rw$missing))
            t$setNote("dev", .("Development rows, before excluding missing values; percentages of the non-missing values. Event: during the whole follow-up."))
        },

        .initExtras = function(sets) {
            if (!self$options$dca) return()
            th <- tryCatch(.parse_numbers(self$options$nbThresholds), error = function(e) numeric(0))
            th <- th[th > 0 & th < 1]
            for (k in names(sets)) for (t in th)
                self$results$nbTable$addRow(rowKey = paste(k, t), values = list(set = sets[[k]], threshold = t))
        },

        .fillExtras = function(all, sets, conf) {
            if (!self$options$dca) return()
            th <- .parse_numbers(self$options$nbThresholds)
            th <- th[th > 0 & th < 1]
            for (k in names(all)) {
                nb <- cp_net_benefit_surv(all[[k]]$y, all[[k]]$p, th, private$.horizon)
                for (i in seq_along(th)) {
                    best <- max(nb$all[i], 0)
                    self$results$nbTable$setRow(rowKey = paste(k, th[i]), values = list(
                        model = nb$model[i], all = nb$all[i], gain = nb$model[i] - best,
                        avoided = 100 * (nb$model[i] - nb$all[i]) / (th[i] / (1 - th[i]))))
                }
            }
            if (self$options$explain)
                self$results$nbTable$setNote("ex", sprintf(.("Net benefit at t = %s, with the event risks above each threshold estimated by Kaplan–Meier (censoring accounted for): true positives/N − false positives/N × p/(1 − p) at threshold p. Gain: over the better of treating all and treating none."), private$.hText()))
        },

        .ssRows = function() c(shrinkage = .("1. Shrinkage factor ≥ 0.9"),
                               optimism = .("2. Optimism in Nagelkerke R² ≤ 0.05"),
                               required = .("Minimum required"),
                               actual = .("This development data")),

        .fillSampleSize = function(y, lp, spec, iv) {
            o <- self$options
            p <- cp_nparams(spec); E <- sum(y[, 2])
            r2 <- if (o$ssSource == "value") o$ssR2 else
                if (!is.null(iv)) iv$corrected[["r2cs"]] else private$.fam$metrics(y, lp)[["r2cs"]]
            if (!is.finite(r2) || r2 <= 0) .fail("The anticipated Cox–Snell R² must be positive; enter a value under Sample size adequacy.")
            rate <- E / sum(y[, 1]); fu <- mean(y[, 1])
            ss <- cp_cox_samplesize(p, r2, rate, fu, private$.horizon)
            evRate <- rate * fu
            t <- self$results$ssTable
            for (k in names(ss$n))
                t$setRow(rowKey = k, values = list(n = ss$n[[k]], events = ss$n[[k]] * evRate, epp = ss$n[[k]] * evRate / p))
            req <- max(ss$n)
            t$setRow(rowKey = "required", values = list(n = req, events = req * evRate, epp = req * evRate / p))
            t$setRow(rowKey = "actual", values = list(n = nrow(y), events = E, epp = E / p))
            t$setNote("in", sprintf(.("Candidate parameters: %d; event rate %.4g per unit of time, mean follow-up %.4g; anticipated Cox–Snell R²: %.3f (%s)."),
                                    p, rate, fu, r2, if (o$ssSource == "value") .("entered") else if (!is.null(iv)) .("optimism-corrected, this model") else .("apparent, this model")))
            t$setNote("risk", sprintf(.("3. Overall risk by t = %s: %.3f (exponential), estimated within [%.3f, %.3f] with the minimum required size."),
                                      private$.hText(), ss$risk, ss$riskCI[1], ss$riskCI[2]))
            t$setNote("verdict", if (nrow(y) >= req) .("The development data meet the minimum sample size.")
                      else sprintf(.("The development data are %d short of the minimum sample size: expect overfitting; consider fewer parameters or penalisation."), req - nrow(y)))
        },

        .calState = function(all, sets, conf) lapply(names(all), function(k)
            c(cp_cox_calcurve(all[[k]]$y, all[[k]]$p, private$.horizon, self$options$calGroups, conf),
              list(label = sets[[k]]))),

        .extraPlotStates = function(all, sets, conf) {
            pl <- self$results$plots
            o <- self$options
            h <- private$.horizon
            if (o$kmPlot && is.null(pl$kmCurves$state)) {
                g <- o$kmGroups
                cuts <- unique(stats::quantile(all$dev$p, seq_len(g - 1) / g))
                pl$kmCurves$setState(list(h = h, sets = lapply(names(all), function(k)
                    c(cp_km_groups(all[[k]]$y, all[[k]]$p, cuts, h), list(label = sets[[k]])))))
            }
            if (o$dca && is.null(pl$dcaCurves$state)) {
                th <- seq(0.01, o$dcaMax, by = 0.01)
                pl$dcaCurves$setState(lapply(names(all), function(k)
                    list(label = sets[[k]], nb = cp_net_benefit_surv(all[[k]]$y, all[[k]]$p, th, h))))
            }
        },

        .outputDefs = function() list(predRisk = function(lp, beta) .cox_risk(lp),
                                      predLP = function(lp, beta) lp - beta[[1]]),

        .plotCal = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_cal(image$state, ggtheme, private$.text("calX"), .("Observed risk (Kaplan–Meier)"))); TRUE
        },
        .plotKm = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_km(image$state, ggtheme)); TRUE
        },
        .plotDca = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_dca(image$state, ggtheme)); TRUE
        }
    ))
)
