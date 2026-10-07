# Binary outcome prediction model: orchestration only. Statistics live in the
# family-agnostic cp_*.R files and cp_family_binary.R, so the linear and Cox
# analyses reuse the same structure (sets, performance tables, importance, plots).

cpLogisticClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "cpLogisticClass",
    inherit = cpLogisticBase,
    private = list(
        .fam = NULL,
        .devData = NULL,

        # ---- structure ---------------------------------------------------------
        .preds = function() c(self$options$covs, self$options$splineVars, self$options$factors),

        .ready = function() !is.null(self$options$dep) && length(private$.preds()) > 0,

        # the level chosen in a LevelSelector, or the second level of the variable
        .level = function(var, chosen) {
            if (!is.null(chosen) && nzchar(chosen)) return(chosen)
            lev <- levels(factor(self$data[[var]]))
            if (length(lev) >= 2) lev[2] else if (length(lev)) lev[1] else ""
        },

        # data sets: development, validation rows (split variable), external file
        .sets = function() {
            sets <- c(dev = "Development")
            sv <- self$options$splitVar
            if (!is.null(sv))
                sets["split"] <- paste0("Validation (", sv, " = ", private$.level(sv, self$options$splitLevel), ")")
            f <- self$options$extFile
            if (!is.null(f)) sets["ext"] <- paste0("External (", basename(f$filename %||% f$path), ")")
            sets
        },

        .init = function() {
            r <- self$results
            if (!private$.ready()) {
                r$instructions$setContent(paste0(
                    "<p>", .("Develop and validate a prediction model for a binary outcome."), "</p><ul>",
                    "<li>", .("Put the outcome in <b>Outcome</b> and choose the event level."), "</li>",
                    "<li>", .("Add predictors: numeric ones as linear covariates or as splines (non-linear), categorical ones as factors."), "</li>",
                    "<li>", .("Under <b>Validation</b>, choose bootstrap or cross-validation; optionally a split variable marking validation rows, or an external data file with the same variable names."), "</li></ul>"))
                for (item in r$items) if (item$name != "instructions") item$setVisible(FALSE)
                return()
            }
            r$instructions$setVisible(FALSE)
            sets <- private$.sets()
            defs <- cp_metric_defs("binary")
            val <- self$options$validation

            for (k in names(sets)) r$samples$addRow(rowKey = k, values = list(set = sets[[k]]))

            for (i in seq_len(nrow(defs)))
                r$perfDev$addRow(rowKey = defs$key[i], values = list(group = defs$group[i], metric = defs$label[i]))
            r$perfDev$getColumn("optimism")$setVisible(val != "none")
            r$perfDev$getColumn("corrected")$setVisible(val != "none")
            if (val == "cv") r$perfDev$getColumn("corrected")$setTitle(.("Cross-validated"))
            ciTitle <- paste0(self$options$ciWidth, "% CI")
            for (t in c("perfDev", "perfVal", "coefTable", "classTable"))
                for (cn in c("lower", "upper")) r[[t]]$getColumn(cn)$setSuperTitle(ciTitle)

            vsets <- setdiff(names(sets), "dev")
            if (length(vsets)) {
                r$perfVal$setVisible(TRUE)
                for (k in vsets) for (i in seq_len(nrow(defs)))
                    r$perfVal$addRow(rowKey = paste(k, defs$key[i]),
                                     values = list(set = sets[[k]], group = defs$group[i], metric = defs$label[i]))
                private$.initCaseMix(sets)
            }
            if ("ext" %in% names(sets)) {
                r$extCheck$setVisible(TRUE)
                for (v in c(self$options$dep, private$.preds())) r$extCheck$addRow(rowKey = v, values = list(var = v))
            }

            private$.initCoefs()

            if (val == "boot") {
                r$stabTable$setVisible(TRUE)
                keys <- c("shrink", "heur", if (self$options$stability) c("mape", "maxape"))
                labs <- c(shrink = .("Uniform shrinkage factor (mean bootstrap calibration slope)"),
                          heur = .("Heuristic shrinkage factor (van Houwelingen)"),
                          mape = .("Mean absolute prediction error (MAPE)"),
                          maxape = .("Largest individual absolute prediction error"))
                for (k in keys) r$stabTable$addRow(rowKey = k, values = list(measure = labs[[k]]))
                r$plots$instability$setVisible(self$options$stability)
                r$plots$calInstability$setVisible(self$options$stability)
            }

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
            if (self$options$sampleSize) {
                crit <- c(shrinkage = .("1. Shrinkage factor ≥ 0.9"),
                          optimism = .("2. Optimism in Nagelkerke R² ≤ 0.05"),
                          intercept = .("3. Overall risk within ± 0.05"),
                          required = .("Minimum required"),
                          actual = .("This development data"))
                for (k in names(crit)) r$ssTable$addRow(rowKey = k, values = list(criterion = crit[[k]]))
            }
            impTables <- c(shapley = "shapleyTable", permImp = "permTable", shap = "shapTable")
            for (opt in names(impTables)) if (self$options[[opt]])
                for (v in private$.preds()) r[[impTables[[opt]]]]$addRow(rowKey = v, values = list(var = v))
            r$plots$shapPatient$setVisible(self$options$shap && self$options$shapRow > 0)

            # one facet per data set
            w <- 420 * length(sets) + 60
            for (p in c("calibration", "riskDist", "dcaCurves"))
                if (length(sets) > 1) r$plots[[p]]$setSize(w, 450)
            np <- length(private$.preds())
            for (p in c("shapleyBar", "permBar", "shapSummary", "shapPatient"))
                r$plots[[p]]$setSize(600, max(300, 120 + 32 * np))
        },

        .initCaseMix = function(sets) {
            t <- self$results$caseMix
            t$setVisible(TRUE)
            for (k in setdiff(names(sets), "dev")) t$getColumn(k)$setVisible(TRUE)
            t$addRow(rowKey = "n", values = list(var = .("N")))
            t$addRow(rowKey = "events", values = list(var = .("Events")))
            for (v in c(self$options$covs, self$options$splineVars))
                t$addRow(rowKey = v, values = list(var = v, level = .("mean (SD)")))
            for (v in self$options$factors)
                for (l in levels(factor(self$data[[v]])))
                    t$addRow(rowKey = paste(v, l), values = list(var = v, level = l))
            t$addRow(rowKey = "lp", values = list(var = .("Linear predictor"), level = .("mean (SD)")))
            t$addRow(rowKey = "member", values = list(var = .("Membership c-statistic")))
        },

        .initCoefs = function() {
            if (!self$options$coefs) return()
            t <- self$results$coefTable
            t$addRow(rowKey = "(Intercept)", values = list(term = .("Intercept")))
            for (v in self$options$covs) t$addRow(rowKey = v, values = list(term = v))
            for (v in self$options$splineVars) {
                nm <- c(v, paste0(v, strrep("'", seq_len(self$options$knots - 2))))
                for (n in nm) t$addRow(rowKey = n, values = list(term = n))
            }
            for (v in self$options$factors) {
                lev <- levels(factor(self$data[[v]]))
                if (length(lev)) t$addRow(rowKey = paste0(v, " (ref)"),
                                          values = list(term = paste0(v, ": ", lev[1], " (reference)")))
                for (l in lev[-1]) t$addRow(rowKey = paste0(v, ": ", l), values = list(term = paste0(v, ": ", l)))
            }
        },

        # ---- computation -------------------------------------------------------
        .run = function() .guard(self, private$.analyse()),

        .analyse = function() {
            if (!private$.ready()) return()
            private$.fam <- cp_family("binary")
            fam <- private$.fam
            r <- self$results
            o <- self$options
            conf <- o$ciWidth / 100
            sets <- private$.sets()

            data <- self$data
            for (v in c(o$covs, o$splineVars)) data[[v]] <- jmvcore::toNumeric(data[[v]])

            # development rows: all but the validation rows of the split variable
            # (rows with a missing split value belong to neither and count as excluded)
            isVal <- rep(FALSE, nrow(data)); isDev <- rep(TRUE, nrow(data))
            if (!is.null(o$splitVar)) {
                sl <- private$.level(o$splitVar, o$splitLevel)
                sv <- as.character(data[[o$splitVar]])
                isVal <- !is.na(sv) & sv == sl
                isDev <- !is.na(sv) & !isVal
                if (!any(isDev)) .fail("Every row is a validation row; no data are left to develop the model.")
                if (!any(isVal)) .fail("No row has ", o$splitVar, " = ", sl, ".")
            }
            prep <- cp_prepare(data[isDev, , drop = FALSE], o$dep, o$covs, o$splineVars, o$factors, o$knots)
            prep$nMissing <- prep$nMissing + sum(!isDev & !isVal)
            prep <- cp_binary_outcome(prep, o$eventLevel)
            private$.devData <- prep$data
            spec <- prep$spec
            X <- cp_design(spec, prep$data); y <- prep$y
            beta <- fam$fit(X, y)
            if (anyNA(beta))
                .fail("The predictor column '", names(beta)[is.na(beta)][1], "' is redundant (a combination of ",
                      "other predictors); remove one of the predictors involved.")
            if (!isTRUE(attr(beta, "converged")))
                .fail("The model did not converge. This usually means complete separation: some predictor, or a ",
                      "combination of them, predicts the outcome perfectly. Remove or merge the predictor(s) involved.")
            if (isTRUE(attr(beta, "separation")))
                .show_warning(self, "separation", .("Some predicted risks are numerically 0 or 1 (quasi-complete separation): coefficients and their intervals may be unreliable."))
            lp <- fam$lp(beta, X); p <- fam$risk(lp)

            # validation data sets, each with X, y, p
            vsets <- list()
            if (any(isVal)) {
                sd <- cp_split_data(data[isVal, , drop = FALSE], spec, o$dep)
                vsets$split <- private$.valSet(sd$data, spec, beta, "split")
                vsets$split$nMissing <- sd$nMissing
            }
            if (!is.null(o$extFile)) {
                ext <- cp_read_external(o$extFile$path, o$extFile$filename %||% o$extFile$path)
                h <- cp_harmonise(ext, spec, o$dep, sets[["ext"]])
                for (i in seq_len(nrow(h$check)))
                    r$extCheck$setRow(rowKey = h$check$var[i],
                                      values = list(role = h$check$role[i], status = h$check$status[i],
                                                    details = h$check$details[i]))
                vsets$ext <- private$.valSet(h$data, spec, beta, "ext")
                vsets$ext$nMissing <- h$nMissing
            }

            private$.fillSamples(y, prep$nMissing, vsets)
            if (length(spec$outcome$levels) > 2)
                r$samples$setNote("lev", sprintf(.("The outcome has %d levels: '%s' is modelled against all the others combined."),
                                                 length(spec$outcome$levels), spec$outcome$event))
            private$.fillCoefs(beta, spec, conf)
            if (o$equation) r$equationText$setContent(private$.equation(beta, spec))

            iv <- private$.internal(X, y)
            private$.fillPerfDev(y, lp, iv, conf)
            for (k in names(vsets)) private$.fillPerfVal(k, vsets[[k]], conf)
            if (length(vsets)) private$.fillCaseMix(prep, X, y, lp, vsets)
            if (!is.null(iv) && iv$method == "boot") private$.fillStab(iv, y, lp, X)

            all <- c(list(dev = list(X = X, y = y, lp = lp, p = p, df = prep$data)), vsets)
            if (o$classif) private$.fillClassif(all, sets, conf)
            if (o$dca) private$.fillDca(all, sets)
            if (o$sampleSize) private$.fillSampleSize(y, lp, spec, iv)
            if (o$shapley) private$.fillShapley(X, y, spec)
            if (o$permImp) private$.fillPerm(beta, spec, all, sets)
            if (o$shap) private$.fillShap(beta, spec, X, data)
            private$.plotStates(all, sets, iv)
            private$.saveOutputs(beta, spec, data)
        },

        .valSet = function(df, spec, beta, key) {
            fam <- private$.fam
            y <- as.integer(df[[spec$outcome$name]] == spec$outcome$event)
            if (sum(y) < 1 || sum(1 - y) < 1)
                .fail(private$.sets()[[key]], ": the validation data need both events and non-events (there are ",
                      sum(y), " events among ", length(y), " complete rows).")
            X <- cp_design(spec, df)
            lp <- fam$lp(beta, X)
            list(X = X, y = y, lp = lp, p = fam$risk(lp), df = df)
        },

        .fillSamples = function(y, nMissing, vsets) {
            t <- self$results$samples
            put <- function(k, yy, miss)
                t$setRow(rowKey = k, values = list(n = length(yy), events = sum(yy), pct = mean(yy), excluded = miss))
            put("dev", y, nMissing)
            for (k in names(vsets)) put(k, vsets[[k]]$y, vsets[[k]]$nMissing)
            if (self$options$explain)
                t$setNote("ev", paste0(.("Event"), ": ", self$options$dep, " = ",
                                       private$.level(self$options$dep, self$options$eventLevel), ". ",
                                       .("Rows with a missing outcome or predictor are excluded (complete-case analysis).")))
        },

        .fillCoefs = function(beta, spec, conf) {
            if (!self$options$coefs) return()
            t <- self$results$coefTable
            V <- .glm_vcov(attr(beta, "fit"))
            se <- sqrt(diag(V))
            z <- stats::qnorm(1 - (1 - conf) / 2)
            for (i in seq_along(beta)) {
                key <- names(beta)[i]
                if (!(key %in% t$rowKeys)) next
                b <- unname(beta[i])
                t$setRow(rowKey = key, values = list(
                    est = b, se = se[i],
                    or = if (i > 1) exp(b) else "", lower = if (i > 1) exp(b - z * se[i]) else "",
                    upper = if (i > 1) exp(b + z * se[i]) else "",
                    p = 2 * stats::pnorm(-abs(b / se[i]))))
            }
            if (length(self$options$splineVars))
                t$setNote("spl", .("Spline terms (', '') are the non-linear parts of a restricted cubic spline; their odds ratios have no direct interpretation. Use the model equation, or SHAP, to see the shape of the association."))
        },

        .equation = function(beta, spec) {
            b <- beta
            terms <- vapply(seq_along(b)[-1], function(i)
                sprintf("%s %.5g × [%s]", if (b[i] < 0) "−" else "+", abs(b[i]), names(b)[i]), "")
            out <- c(sprintf("log-odds(%s = %s) =", spec$outcome$name, spec$outcome$event),
                     sprintf("    %.5g", b[1]), paste0("    ", terms), "",
                     "risk = 1 / (1 + exp(−log-odds))")
            for (info in spec$vars) {
                if (info$role == "factor")
                    out <- c(out, "", sprintf("[%s: level] = 1 if %s is that level, else 0 (reference: %s)",
                                              info$name, info$name, info$levels[1]))
                if (info$role == "spline") {
                    kn <- info$knots; k <- length(kn)
                    out <- c(out, "", sprintf("Spline of %s, knots t1..t%d = %s", info$name, k,
                                              paste(signif(kn, 5), collapse = ", ")),
                             sprintf("  [%s'] ... : ((x − tj)₊³ − (x − t%d)₊³ (t%d − tj)/(t%d − t%d) + (x − t%d)₊³ (t%d − tj)/(t%d − t%d)) / (t%d − t1)²,  j = 1..%d",
                                     info$name, k - 1, k, k, k - 1, k, k - 1, k, k - 1, k, k - 2))
                }
            }
            paste(out, collapse = "\n")
        },

        # internal validation, cached in the performance table's state
        .internal = function(X, y) {
            o <- self$options
            if (o$validation == "none") return(NULL)
            iv <- self$results$perfDev$state
            if (is.null(iv)) {
                iv <- cp_internal(X, y, private$.fam, o$validation, B = o$nBoot, k = o$nFolds,
                                  R = o$nRepeats, seed = o$seed, stability = o$validation == "boot")
                self$results$perfDev$setState(iv)
            }
            iv
        },

        .fillPerfDev = function(y, lp, iv, conf) {
            t <- self$results$perfDev
            app <- private$.fam$metrics(y, lp)
            ci <- cp_binary_ci(y, lp, conf)
            for (k in t$rowKeys) {
                # apparent CITL is 0 by construction; print it as 0, not as rounding noise
                a <- if (abs(app[[k]]) < 1e-8) 0 else app[[k]]
                v <- list(apparent = a, lower = ci[[k]][1] %||% "", upper = ci[[k]][2] %||% "")
                if (!is.null(iv)) { v$optimism <- iv$optimism[[k]]; v$corrected <- iv$corrected[[k]] }
                t$setRow(rowKey = k, values = v)
            }
            if (!is.null(iv)) {
                t$setNote("iv", if (iv$method == "boot")
                    sprintf(.("Bootstrap optimism correction (Harrell): %d samples (%d usable); corrected = apparent − optimism."), iv$B, iv$nOk)
                else sprintf(.("Cross-validation: %d folds × %d repeats, stratified by outcome; performance of the out-of-fold predictions, averaged over repeats."), iv$k, iv$R))
            } else {
                t$setNote("iv", .("Apparent performance is optimistic: it is measured in the data used to fit the model. Choose bootstrap or cross-validation under Validation."))
            }
            if (self$options$explain)
                t$setNote("ex", .("AUC: probability that a random event has a higher predicted risk than a random non-event (0.5 = no discrimination). Calibration-in-the-large: 0 is ideal (log-odds). Calibration slope: 1 is ideal; below 1 the predictions are too extreme (overfitting). O/E: observed / expected events. ICI, E50, E90, Emax: mean, median, 90th percentile and maximum absolute difference between the predicted and the (smoothed) observed risk. Brier: mean squared error of the predicted risks (lower is better); scaled Brier: 1 − Brier / Brier of a model without predictors."))
        },

        .fillPerfVal = function(key, vs, conf) {
            t <- self$results$perfVal
            m <- private$.fam$metrics(vs$y, vs$lp)
            ci <- cp_binary_ci(vs$y, vs$lp, conf)
            for (k in names(m)) {
                rk <- paste(key, k)
                if (!(rk %in% t$rowKeys)) next
                t$setRow(rowKey = rk, values = list(est = m[[k]], lower = ci[[k]][1] %||% "",
                                                    upper = ci[[k]][2] %||% ""))
            }
            if (self$options$explain)
                t$setNote("ex", .("The development model is applied unchanged to the validation data. A calibration-in-the-large far from 0 means the overall risk differs between the populations; a slope below 1, that the predictions are too extreme there."))
        },

        .fillCaseMix = function(prep, X, y, lp, vsets) {
            t <- self$results$caseMix
            fmt <- function(m, s) sprintf("%.2f (%.2f)", m, s)
            col <- function(k, df, yy, lpk, Xk) {
                vals <- list()
                vals$n <- as.character(length(yy))
                vals$events <- sprintf("%d (%.1f%%)", sum(yy), 100 * mean(yy))
                for (info in prep$spec$vars) {
                    x <- df[[info$name]]
                    if (info$role == "factor") {
                        for (l in info$levels) vals[[paste(info$name, l)]] <- sprintf("%.1f%%", 100 * mean(x == l))
                    } else vals[[info$name]] <- fmt(mean(x), stats::sd(x))
                }
                vals$lp <- fmt(mean(lpk), stats::sd(lpk))
                vals$member <- if (k == "dev") "" else sprintf("%.3f", cp_membership_c(X, y, Xk, yy))
                vals
            }
            cols <- list(dev = col("dev", prep$data, y, lp, X))
            for (k in names(vsets)) cols[[k]] <- col(k, vsets[[k]]$df, vsets[[k]]$y, vsets[[k]]$lp, vsets[[k]]$X)
            for (rk in t$rowKeys) {
                v <- lapply(cols, function(cl) cl[[rk]] %||% "")
                t$setRow(rowKey = rk, values = v)
            }
            t$setNote("member", .("Membership c-statistic: how well the predictors and the outcome distinguish validation from development rows (Debray et al. 2015). Near 0.5: similar case-mix (reproducibility); high: different populations (transportability). A larger SD of the linear predictor means a more heterogeneous case-mix, which by itself raises the AUC."))
        },

        .fillStab = function(iv, y, lp, X) {
            t <- self$results$stabTable
            n <- length(y); ybar <- mean(y)
            ll <- .binary_ll(y, stats::plogis(lp)); ll0 <- n * (ybar * log(ybar) + (1 - ybar) * log(1 - ybar))
            lr <- 2 * (ll - ll0)
            t$setRow(rowKey = "shrink", values = list(value = iv$shrinkage))
            t$setRow(rowKey = "heur", values = list(value = (lr - ncol(X)) / lr))
            if (!is.null(iv$stability) && "mape" %in% t$rowKeys) {
                t$setRow(rowKey = "mape", values = list(value = iv$stability$mape))
                t$setRow(rowKey = "maxape", values = list(value = iv$stability$maxape))
            }
            if (self$options$explain)
                t$setNote("ex", .("A shrinkage factor well below 1 (e.g. < 0.9) signals overfitting: the coefficients would need to be shrunk towards 0. MAPE: average absolute difference between an individual's risk from the developed model and from models developed in bootstrap samples (Riley & Collins 2023)."))
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

        .fillShapley = function(X, y, spec) {
            t <- self$results$shapleyTable
            sh <- t$state
            if (is.null(sh)) {
                sh <- cp_shapley_r2(X, y, spec, private$.fam, self$options$shapleyMeasure, seed = self$options$seed)
                t$setState(sh)
            }
            for (i in seq_len(nrow(sh)))
                t$setRow(rowKey = sh$var[i], values = list(shapley = sh$shapley[i], share = sh$share[i], rank = sh$rank[i]))
            lab <- if (self$options$shapleyMeasure == "r2mf") "McFadden R²" else "Nagelkerke R²"
            t$setNote("tot", sprintf(.("%s of the full model: %.4f, shared among the predictors by their average gain over all orders of entry (Shapley / LMG)%s."),
                                     lab, attr(sh, "total"), if (isTRUE(attr(sh, "exact"))) "" else .("; estimated from 200 random orders")))
            img <- self$results$plots$shapleyBar
            if (is.null(img$state)) img$setState(list(df = sh[, c("var", "shapley")], label = lab))
        },

        .fillPerm = function(beta, spec, all, sets) {
            t <- self$results$permTable
            key <- if ("ext" %in% names(all)) "ext" else if ("split" %in% names(all)) "split" else "dev"
            pi <- t$state
            if (is.null(pi)) {
                pi <- cp_perm_importance(beta, spec, all[[key]]$df, all[[key]]$y, private$.fam,
                                         nPerm = self$options$nPerm, seed = self$options$seed)
                t$setState(pi)
            }
            for (i in seq_len(nrow(pi)))
                t$setRow(rowKey = pi$var[i], values = as.list(pi[i, -1]))
            t$setNote("data", sprintf(.("Computed in: %s. Each predictor is shuffled %d times; larger losses mean the model relies more on it."),
                                      sets[[key]], self$options$nPerm))
            img <- self$results$plots$permBar
            if (is.null(img$state)) img$setState(pi)
        },

        .fillShap = function(beta, spec, X, data) {
            phi <- cp_shap(beta, spec, X)
            ma <- colMeans(abs(phi))
            t <- self$results$shapTable
            for (v in names(ma)) t$setRow(rowKey = v, values = list(meanAbs = ma[[v]], share = ma[[v]] / sum(ma)))
            t$setNote("ex", .("SHAP: the contribution of each predictor to an individual's log-odds, relative to the average individual of the development data (exact for this additive model)."))
            img <- self$results$plots$shapSummary
            if (is.null(img$state)) {
                n <- nrow(phi)
                idx <- if (n > 1000) .with_seed(self$options$seed, sort(sample.int(n, 1000))) else seq_len(n)
                df <- private$.devData
                long <- do.call(rbind, lapply(colnames(phi), function(v) {
                    x <- df[[v]][idx]
                    val <- if (is.factor(x)) (as.integer(x) - 1) / max(1, nlevels(x) - 1) else
                        (rank(x) - 1) / max(1, length(x) - 1)
                    data.frame(var = v, phi = round(phi[idx, v], 4), value = round(val, 3))
                }))
                img$setState(long)
            }
            row <- self$options$shapRow
            img <- self$results$plots$shapPatient
            if (row > 0 && is.null(img$state)) {
                nd <- cp_newdata(data, spec)
                rk <- as.character(row)
                if (!(rk %in% rownames(nd)))
                    .fail("Row ", row, " has missing predictor values (or does not exist), so its prediction cannot be explained.")
                Xi <- cp_design(spec, nd[rk, , drop = FALSE])
                ph <- cp_shap(beta, spec, Xi, Xbg = X)
                labs <- vapply(names(spec$vars), function(v) {
                    x <- nd[rk, v]
                    paste0(v, " = ", if (is.factor(x)) as.character(x) else .num(x, 4))
                }, "")
                lpi <- private$.fam$lp(beta, Xi)
                img$setState(list(label = unname(labs), phi = unname(ph[1, ]), base = attr(ph, "base"),
                                  lp = lpi, risk = private$.fam$risk(lpi),
                                  title = sprintf(.("Row %d"), row)))
            }
        },

        .plotStates = function(all, sets, iv) {
            pl <- self$results$plots
            o <- self$options
            conf <- o$ciWidth / 100
            if (o$rocPlot && is.null(pl$roc$state))
                pl$roc$setState(lapply(names(all), function(k) {
                    rp <- cp_roc_points(all[[k]]$y, all[[k]]$p)
                    c(rp, list(key = k, label = sets[[k]], auc = .auc_delong(all[[k]]$y, all[[k]]$p)$auc))
                }))
            if (o$calPlot && is.null(pl$calibration$state))
                pl$calibration$setState(lapply(names(all), function(k)
                    c(cp_binary_calcurve(all[[k]]$y, all[[k]]$p, o$calGroups, conf), list(label = sets[[k]]))))
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
            if (!is.null(iv$stability)) {
                if (is.null(pl$instability$state)) pl$instability$setState(iv$stability)
                if (is.null(pl$calInstability$state)) pl$calInstability$setState(iv$stability)
            }
        },

        .eventLabel = function() private$.level(self$options$dep, self$options$eventLevel),

        .saveOutputs = function(beta, spec, data) {
            o <- self$options
            if (!(o$predRisk || o$predLP)) return()
            nd <- cp_newdata(data, spec)
            lp <- rep(NA_real_, nrow(data)); names(lp) <- rownames(data)
            lp[rownames(nd)] <- private$.fam$lp(beta, cp_design(spec, nd))
            if (o$predRisk && self$results$predRisk$isNotFilled()) {
                self$results$predRisk$setRowNums(rownames(data))
                self$results$predRisk$setValues(unname(private$.fam$risk(lp)))
            }
            if (o$predLP && self$results$predLP$isNotFilled()) {
                self$results$predLP$setRowNums(rownames(data))
                self$results$predLP$setValues(unname(lp))
            }
        },

        # ---- plots -------------------------------------------------------------
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
        },
        .plotShapley = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_shapley(image$state$df, image$state$label, ggtheme)); TRUE
        },
        .plotPerm = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_perm(image$state, ggtheme)); TRUE
        },
        .plotShapSummary = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_shap_summary(image$state, ggtheme)); TRUE
        },
        .plotShapPatient = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_shap_patient(image$state, ggtheme)); TRUE
        },
        .plotInstability = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_instability(image$state, ggtheme)); TRUE
        },
        .plotCalInstability = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_cal_instability(image$state, ggtheme)); TRUE
        }
    )
)
