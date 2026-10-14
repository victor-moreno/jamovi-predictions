# Validate a published model: a model read from a file (clinpred JSON, or a
# table of coefficients) is applied unchanged to the current data. Reuses the
# family functions, the updating and classification engine, and the plots of
# the development analyses.

# performance columns of the updating table, per model type
.pub_upd_keys <- list(binary = c("auc", "brier", "ici"), gaussian = c("r2", "rmse"), cox = c("c", "oe", "brier"))

cpPublishedClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "cpPublishedClass",
    inherit = cpPublishedBase,
    private = list(
        .model = NULL,

        .ready = function() {
            o <- self$options
            !is.null(o$modelFile) && !is.null(o$dep) && length(o$covs) > 0 && (o$family != "cox" || !is.null(o$time))
        },

        .levelsOf = function(var) {
            x <- self$data[[var]]
            if (is.factor(x)) levels(x) else levels(factor(x))
        },
        .level = function(var, chosen) {
            if (!is.null(chosen) && nzchar(chosen)) return(chosen)
            lev <- private$.levelsOf(var)
            if (length(lev) >= 2) lev[2] else if (length(lev)) lev[1] else ""
        },

        # the model file, read once per analysis object (it is small)
        .readModel = function() {
            if (is.null(private$.model)) {
                f <- self$options$modelFile
                m <- cp_model_read(f$path, f$filename %||% f$path)
                if (m$kind == "terms") m$family <- self$options$family
                private$.model <- m
            }
            private$.model
        },

        .init = function() {
            r <- self$results
            o <- self$options
            if (!private$.ready()) {
                r$instructions$setContent(paste0(
                    "<p>", .("Validate a published prediction model in the current data."), "</p><ul>",
                    "<li>", .("<b>Model file</b>: a model exported by clinpred (Export model, JSON), or a table (csv, tsv, txt, xlsx) with the columns <b>term</b>, <b>estimate</b> and, optionally, <b>SE</b>."), "</li>",
                    "<li>", .("Terms name variables of the data: <code>age</code>, <code>sex = male</code>, <code>log(crp)</code>, <code>age^2</code>, <code>age/10</code>, <code>age * sex = male</code>; a row <code>Intercept</code>. For Cox models, rows <code>baseline survival</code> (at the horizon), <code>horizon</code> and, when the published prognostic index is centred, <code>centre</code>."), "</li>",
                    "<li>", .("Choose the model type, put the outcome (for Cox, the event and the time) in their boxes, and the variables the model uses in <b>Predictor variables</b>."), "</li></ul>"))
                for (item in r$items) if (item$name != "instructions") item$setVisible(FALSE)
                return()
            }
            r$instructions$setVisible(FALSE)
            fam <- o$family
            ciTitle <- paste0(o$ciWidth, "% CI")
            for (t in c("modelTable", "perf", "classTable"))
                for (cn in c("lower", "upper")) r[[t]]$getColumn(cn)$setSuperTitle(ciTitle)

            # model rows, so that they exist after .init (a read error is reported by .run)
            m <- tryCatch(private$.readModel(), error = function(e) NULL)
            t <- r$modelTable
            if (fam != "gaussian") {
                t$getColumn("ratio")$setVisible(TRUE)
                t$getColumn("ratio")$setTitle(if (fam == "cox") .("HR") else .("OR"))
            }
            if (!is.null(m)) {
                if (fam != "cox") t$addRow(rowKey = "int", values = list(term = .("Intercept")))
                for (i in seq_along(m$terms)) t$addRow(rowKey = paste0("t", i), values = list(term = m$terms[i]))
                if (fam == "cox") t$addRow(rowKey = "s0", values = list(term = .("Baseline survival")))
            }

            defs <- cp_metric_defs(fam)
            for (i in seq_len(nrow(defs)))
                r$perf$addRow(rowKey = defs$key[i], values = list(group = defs$group[i], metric = defs$label[i]))
            vis <- switch(fam, binary = c("events", "pct"), gaussian = c("mean", "sd"), cox = c("events", "pct", "km"))
            for (cn in vis) r$samples$getColumn(cn)$setVisible(TRUE)

            if (o$updating) {
                u <- r$updTable
                u$setVisible(TRUE)
                for (cn in .pub_upd_keys[[fam]]) u$getColumn(cn)$setVisible(TRUE)
                if (fam == "gaussian") u$getColumn("stat")$setTitle("F")
                if (fam == "cox") for (cn in c("stat", "df", "p")) u$getColumn(cn)$setSuperTitle(.("Sequential test"))
                meth <- c(original = .("Original model"), citl = .("Recalibration in the large"),
                          recal = .("Recalibration"), revision = .("Model revision"))
                for (k in names(meth)) u$addRow(rowKey = k, values = list(method = meth[[k]]))
            }
            if (o$classif && fam == "binary") {
                r$confTable$setVisible(TRUE); r$classTable$setVisible(TRUE)
                meas <- c(sens = .("Sensitivity"), spec = .("Specificity"), ppv = .("PPV"), npv = .("NPV"), acc = .("Accuracy"))
                for (k in names(meas)) r$classTable$addRow(rowKey = k, values = list(measure = meas[[k]]))
            }
            if (o$dca && fam != "gaussian") {
                r$nbTable$setVisible(TRUE)
                th <- tryCatch(.parse_numbers(o$nbThresholds), error = function(e) numeric(0))
                for (x in th[th > 0 & th < 1]) r$nbTable$addRow(rowKey = as.character(x), values = list(threshold = x))
                r$plots$dcaCurves$setVisible(TRUE)
            }
            r$plots$roc$setVisible(o$rocPlot && fam == "binary")
            r$plots$kmCurves$setVisible(o$kmPlot && fam == "cox")
        },

        .run = function() .guard(self, private$.analyse()),

        .analyse = function() {
            if (!private$.ready()) return()
            o <- self$options
            r <- self$results
            famName <- o$family
            m <- private$.readModel()
            if (m$kind == "spec" && m$family != famName)
                .fail(sprintf(.("The model file holds a model of type %s: choose that model type."), .cp_family_title()[[m$family]]))
            data <- self$data
            des <- cp_model_design(m, data)
            okX <- stats::complete.cases(des$X)

            # outcome
            ev <- private$.level(o$dep, o$eventLevel)
            if (famName == "gaussian") {
                yv <- jmvcore::toNumeric(data[[o$dep]])
                if (is.factor(yv)) yv <- suppressWarnings(as.numeric(as.character(yv)))
                okY <- !is.na(yv)
            } else {
                sv <- as.character(data[[o$dep]])
                if (!(ev %in% sv)) .fail(sprintf(.("The event level '%s' does not occur in '%s'."), ev, o$dep))
                st <- ifelse(is.na(sv), NA_integer_, as.integer(sv == ev))
                okY <- !is.na(st)
                if (famName == "cox") {
                    tm <- jmvcore::toNumeric(data[[o$time]])
                    okY <- okY & !is.na(tm)
                }
            }
            ok <- okX & okY
            if (sum(ok) < 3) .fail(.("Fewer than 3 rows have the outcome and every predictor of the model."))
            X <- des$X[ok, , drop = FALSE]
            y <- switch(famName, gaussian = yv[ok], binary = st[ok], cox = cbind(time = tm[ok], status = st[ok]))
            if (famName == "binary" && (sum(y) < 1 || sum(1 - y) < 1))
                .fail(sprintf(.("The data need both events and non-events (there are %s events among %s rows)."), sum(y), length(y)))
            if (famName == "gaussian" && stats::sd(y) == 0) .fail(.("The outcome does not vary."))

            # Cox: the horizon, and the baseline (from the model, else from these data)
            h <- NULL; baseNote <- NULL
            beta <- cp_model_beta(m)
            if (famName == "cox") {
                if (!any(y[, 2] == 1)) .fail(.("There are no events in the data."))
                if (o$horizon > 0 && !is.null(m$baseline) && !is.null(m$horizon) && !isTRUE(all.equal(o$horizon, m$horizon)))
                    .fail(sprintf(.("The model's baseline survival is for t = %s; set the horizon to 0 (from the model file) or to %s."), .num(m$horizon), .num(m$horizon)))
                h <- if (o$horizon > 0) o$horizon else m$horizon %||% .median_followup(y[, 1], y[, 2])
                if (!is.null(m$baseline) && is.null(m$horizon) && o$horizon <= 0)
                    .fail(.("The model gives a baseline survival but not its time: add a row 'horizon' to the model table, or set the prediction horizon."))
                if (!any(y[, 2] == 1 & y[, 1] <= h)) .fail(sprintf(.("No event happens before the horizon (%s)."), .num(h)))
                if (!any(y[, 1] > h)) .fail(sprintf(.("Nobody is followed beyond the horizon (%s)."), .num(h)))
                if (is.na(beta[1])) {
                    eta <- drop(X %*% beta[-1])
                    beta[1] <- log(.cox_H0(y[, 1], y[, 2], eta, h))
                    baseNote <- .("The model gives no baseline survival: it was estimated from these data, so the O/E ratio is close to 1 by construction and calibration-in-the-large is not assessed.")
                }
            }
            fam <- cp_family(famName, horizon = h)
            lp <- fam$lp(beta, X); p <- fam$risk(lp)
            conf <- o$ciWidth / 100

            private$.fillModel(m, beta, famName, h, conf, baseNote)
            private$.fillCheck(des$check)
            sm <- switch(famName,
                binary = list(n = length(y), events = sum(y), pct = mean(y)),
                gaussian = list(n = length(y), mean = mean(y), sd = stats::sd(y)),
                cox = list(n = nrow(y), events = sum(y[, 2]), pct = mean(y[, 2]), km = 1 - .km(y[, 1], y[, 2], h)$surv))
            r$samples$setRow(rowNo = 1, values = c(sm, list(excluded = nrow(data) - sum(ok))))
            if (famName == "cox")
                r$samples$setNote("h", sprintf(.("Risks by the horizon t = %s%s."), .num(h, 4),
                                               if (o$horizon > 0) "" else if (!is.null(m$horizon)) .(" (from the model file)") else .(" (median follow-up)")))

            mt <- fam$metrics(y, lp); ci <- fam$ci(y, lp, conf)
            for (k in r$perf$rowKeys)
                # rounding noise (e.g. validating in the development data) printed as 0
                r$perf$setRow(rowKey = k, values = list(est = if (isTRUE(abs(mt[[k]]) < 1e-8)) 0 else mt[[k]], lower = ci[[k]][1] %||% "", upper = ci[[k]][2] %||% ""))
            r$perf$setNote("ex", .("The published model is applied unchanged. Calibration-in-the-large (or O/E) far from its ideal means that the overall risk differs from the population where the model was developed; a calibration slope below 1, that its predictions are too extreme here."))

            if (o$updating) private$.fillUpd(fam, X, y, lp)
            if (o$classif && famName == "binary") private$.fillClassif(y, p, conf)
            if (o$dca && famName != "gaussian") private$.fillDca(y, p, h)
            private$.plotStates(famName, y, p, h, conf)
            private$.saveOutputs(fam, beta, des$X, okX)
        },

        .fillModel = function(m, beta, famName, h, conf, baseNote) {
            t <- self$results$modelTable
            z <- stats::qnorm(1 - (1 - conf) / 2)
            row <- function(key, est, se, ratio) {
                v <- list(est = est, se = if (is.na(se)) "" else se,
                          p = if (is.na(se) || se <= 0) "" else 2 * stats::pnorm(-abs(est / se)))
                if (ratio) v <- c(v, list(ratio = exp(est), lower = if (is.na(se)) "" else exp(est - z * se),
                                          upper = if (is.na(se)) "" else exp(est + z * se)))
                else v <- c(v, list(ratio = "", lower = if (is.na(se)) "" else est - z * se, upper = if (is.na(se)) "" else est + z * se))
                if (key %in% t$rowKeys) t$setRow(rowKey = key, values = v)
            }
            if (famName != "cox") {
                row("int", beta[[1]], m$interceptSE %||% NA_real_, FALSE)
                if (is.null(m$intercept) && m$kind == "terms")
                    t$setNote("int", .("The model table has no intercept: 0 is used, so calibration-in-the-large will be poor. Add a row 'Intercept'."))
            }
            for (i in seq_along(m$terms)) row(paste0("t", i), m$estimate[i], m$se[i], famName != "gaussian")
            if (famName == "cox" && "s0" %in% t$rowKeys) {
                s0 <- exp(-exp(beta[[1]]))
                t$setRow(rowKey = "s0", values = list(est = s0, se = "", lower = "", upper = "", p = ""))
                t$setNote("s0", if (!is.null(baseNote)) baseNote else
                    sprintf(.("Baseline survival at t = %s, for a prognostic index of 0%s."), .num(h, 4),
                            if (isTRUE(m$centre != 0)) sprintf(.(" after subtracting the centre %s"), .num(m$centre, 5)) else ""))
            }
            t$setNote("src", sprintf(.("Read from %s (%s)."), m$source,
                                     if (m$kind == "spec") .("a model exported by clinpred") else .("a table of coefficients")))
        },

        .fillCheck = function(check) {
            t <- self$results$check
            o <- self$options
            add <- function(key, v) if (key %in% t$rowKeys) t$setRow(rowKey = key, values = v) else t$addRow(rowKey = key, values = v)
            add(".dep", list(var = o$dep, role = if (o$family == "cox") .("Event") else .("Outcome"), status = "OK",
                             details = if (o$family != "gaussian") paste0(.("event"), ": ", private$.level(o$dep, o$eventLevel)) else ""))
            if (o$family == "cox") add(".time", list(var = o$time, role = .("Time"), status = "OK", details = ""))
            for (i in seq_len(nrow(check)))
                add(check$var[i], list(var = check$var[i], role = .check_label(check$role[i]), status = .check_label(check$status[i]), details = check$details[i]))
        },

        .fillUpd = function(fam, X, y, lp) {
            t <- self$results$updTable
            u <- cp_update(fam, X, y, lp, self$options$updAlpha)
            keys <- .pub_upd_keys[[fam$name]]
            for (k in t$rowKeys) {
                mm <- u$methods[[k]]
                if (is.null(mm)) next
                perf <- fam$metrics(y, mm$lp)[keys]
                t$setRow(rowKey = k, values = c(list(
                    a = if (is.na(mm$a)) "" else mm$a, b = if (is.na(mm$b)) "" else mm$b, k = mm$k,
                    stat = mm$test[["stat"]] %||% "", df = mm$test[["df"]] %||% "", p = mm$test[["p"]] %||% "",
                    chosen = if (identical(u$chosen, k)) "✓" else ""), as.list(perf)))
            }
            if (is.na(u$chosen))
                t$setNote("fail", .("The revised model could not be fitted in these data, so the closed testing procedure was not run."))
            t$setNote("ex", if (fam$name == "cox") sprintf(.("Recalibration in the large re-estimates the baseline risk at the horizon; recalibration also the slope of the published linear predictor; revision re-estimates every coefficient. Sequential tests (each p in the row of the more complex method): revision vs recalibration and slope = 1 by partial likelihood ratio, O/E at the horizon = 1 by a Wald test; the most complex method whose test is significant at α = %s is chosen (✓). Performance of an updated model is measured in the data used to update it, so it is optimistic."), .num(self$options$updAlpha)) else sprintf(.("Recalibration re-estimates the intercept and the slope of the published linear predictor; revision re-estimates every coefficient. Each test compares a method with revision; the closed testing procedure (Vergouwe et al. 2017) keeps the first method, from the original model up, not significantly worse than revision at α = %s (✓). Performance of an updated model is measured in the data used to update it, so it is optimistic."),
                                    .num(self$options$updAlpha)))
        },

        .fillClassif = function(y, p, conf) {
            o <- self$options
            thr <- if (o$thrMethod == "prevalence") mean(y) else o$threshold
            cl <- cp_binary_classify(y, p, thr, conf)
            self$results$confTable$setRow(rowNo = 1, values = list(threshold = thr, tp = cl$tp, fp = cl$fp, fn = cl$fn, tn = cl$tn))
            for (k in c("sens", "spec", "ppv", "npv", "acc"))
                self$results$classTable$setRow(rowKey = k, values = list(est = cl[[k]][1], lower = cl[[k]][2], upper = cl[[k]][3]))
        },

        .fillDca = function(y, p, h) {
            t <- self$results$nbTable
            keys <- unlist(t$rowKeys)
            th <- as.numeric(keys)
            nb <- if (is.null(h)) cp_net_benefit(y, p, th) else cp_net_benefit_surv(y, p, th, h)
            for (i in seq_along(th))
                t$setRow(rowKey = keys[i], values = list(model = nb$model[i], all = nb$all[i],
                    gain = nb$model[i] - max(nb$all[i], 0),
                    avoided = 100 * (nb$model[i] - nb$all[i]) / (th[i] / (1 - th[i]))))
        },

        .plotStates = function(famName, y, p, h, conf) {
            pl <- self$results$plots
            o <- self$options
            lab <- .("These data")
            if (o$calPlot && is.null(pl$calibration$state)) {
                cc <- switch(famName,
                    binary = cp_binary_calcurve(y, p, o$calGroups, conf),
                    gaussian = cp_gaussian_calcurve(y, p, conf),
                    cox = cp_cox_calcurve(y, p, h, o$calGroups, conf))
                pl$calibration$setState(list(family = famName, sets = list(c(cc, list(label = lab)))))
            }
            if (o$rocPlot && famName == "binary" && is.null(pl$roc$state))
                pl$roc$setState(list(c(cp_roc_points(y, p), list(key = "ext", label = lab, auc = .auc_delong(y, p)$auc))))
            if (o$kmPlot && famName == "cox" && is.null(pl$kmCurves$state)) {
                g <- o$kmGroups
                cuts <- unique(stats::quantile(p, seq_len(g - 1) / g))
                pl$kmCurves$setState(list(h = h, sets = list(c(cp_km_groups(y, p, cuts, h), list(label = lab)))))
            }
            if (o$dca && famName != "gaussian" && is.null(pl$dcaCurves$state)) {
                th <- seq(0.01, o$dcaMax, by = 0.01)
                nb <- if (famName == "cox") cp_net_benefit_surv(y, p, th, h) else cp_net_benefit(y, p, th)
                pl$dcaCurves$setState(list(list(label = lab, nb = nb)))
            }
        },

        # predictions for every row with complete predictors
        .saveOutputs = function(fam, beta, Xall, okX) {
            o <- self$options
            if (!(isTRUE(o$pred) || isTRUE(o$predLP))) return()
            lp <- rep(NA_real_, nrow(Xall))
            lp[okX] <- fam$lp(beta, Xall[okX, , drop = FALSE])
            rows <- rownames(self$data)
            if (isTRUE(o$pred) && self$results$pred$isNotFilled()) {
                self$results$pred$setRowNums(rows)
                self$results$pred$setValues(fam$risk(lp))
            }
            if (isTRUE(o$predLP) && self$results$predLP$isNotFilled()) {
                self$results$predLP$setRowNums(rows)
                self$results$predLP$setValues(lp)
            }
        },

        .plotCal = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(switch(st$family,
                gaussian = cp_plot_cal_cont(st$sets, ggtheme),
                cox = cp_plot_cal(st$sets, ggtheme, .("Predicted risk"), .("Observed risk (Kaplan–Meier)")),
                cp_plot_cal(st$sets, ggtheme)))
            TRUE
        },
        .plotRoc = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_roc(image$state, ggtheme)); TRUE
        },
        .plotKm = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_km(image$state, ggtheme)); TRUE
        },
        .plotDca = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(cp_plot_dca(image$state, ggtheme)); TRUE
        }
    )
)
