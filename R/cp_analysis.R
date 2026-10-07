# The orchestration shared by every analysis (cpLogistic, cpLinear, ...): a
# list of R6 private methods. Each analysis builds its class from
# modifyList(.cp_private, <its hooks>), where the hooks are the few methods
# that depend on the outcome type:
#   .famName, .outcomeVars(), .intro(), .setOutcome(prep), .yOf(df, spec, key),
#   .sampleValues(y, nMissing), .outcomeSummary(y), .initDesc(), .fillDesc(dd),
#   .equation(beta, spec), .text(key), .initExtras(sets),
#   .fillExtras(all, sets, conf), .extraPlotStates(all, sets, conf),
#   .calState(all, sets, conf), .ssRows(), .fillSampleSize(y, lp, spec, iv),
#   .outputDefs(), .updKeys
# Statistics live in the cp_*.R files; this code only fills results.

.cp_private <- list(
    .fam = NULL,
    .devData = NULL,

    # ---- structure ---------------------------------------------------------
    .preds = function() c(self$options$covs, self$options$splineVars, self$options$factors),

    .ready = function() length(private$.outcomeVars()) == length(private$.outcomeNames) &&
        length(private$.preds()) > 0,

    # names of the outcome options (their values: .outcomeVars())
    .outcomeNames = "dep",

    # the level chosen in a LevelSelector, or the second level of the variable
    .level = function(var, chosen) {
        if (!is.null(chosen) && nzchar(chosen)) return(chosen)
        lev <- levels(factor(self$data[[var]]))
        if (length(lev) >= 2) lev[2] else if (length(lev)) lev[1] else ""
    },

    # Reference Levels and Covariate Scaling choices, as named lists
    .refs = function() {
        out <- list()
        for (it in self$options$refLevels) if (!is.null(it$ref) && nzchar(it$ref)) out[[it$var]] <- it$ref
        out
    },
    .scales = function() {
        out <- list()
        for (it in self$options$covScales) out[[it$var]] <- it$scale %||% "unit"
        out
    },
    .refOf = function(v) private$.refs()[[v]] %||% levels(factor(self$data[[v]]))[1],

    # candidate interactions: the listed two-way terms, or every pair of predictors
    .candidates = function() {
        preds <- private$.preds()
        if (self$options$intAllPairs)
            return(if (length(preds) < 2) list() else utils::combn(preds, 2, simplify = FALSE))
        ints <- lapply(self$options$interactions, function(t) unique(unlist(t)))
        Filter(function(t) length(t) == 2 && all(t %in% preds), ints)
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

    # the updating methods, labelled for this family
    .updMethods = function() c(original = .("Original model"), citl = .("Recalibration in the large"),
                               recal = private$.text("recalLabel"), revision = .("Model revision")),

    .init = function() {
        r <- self$results
        if (!private$.ready()) {
            r$instructions$setContent(private$.intro())
            for (item in r$items) if (item$name != "instructions") item$setVisible(FALSE)
            return()
        }
        r$instructions$setVisible(FALSE)
        sets <- private$.sets()
        defs <- cp_metric_defs(private$.famName)
        val <- self$options$validation

        for (k in names(sets)) r$samples$addRow(rowKey = k, values = list(set = sets[[k]]))

        for (i in seq_len(nrow(defs)))
            r$perfDev$addRow(rowKey = defs$key[i], values = list(group = defs$group[i], metric = defs$label[i]))
        r$perfDev$getColumn("optimism")$setVisible(val != "none")
        r$perfDev$getColumn("corrected")$setVisible(val != "none")
        if (val == "cv") r$perfDev$getColumn("corrected")$setTitle(.("Cross-validated"))
        ciTitle <- paste0(self$options$ciWidth, "% CI")
        for (t in c("perfDev", "perfVal", "coefTable", "classTable"))
            if (!is.null(r[[t]])) for (cn in c("lower", "upper")) r[[t]]$getColumn(cn)$setSuperTitle(ciTitle)

        vsets <- setdiff(names(sets), "dev")
        if (length(vsets)) {
            r$perfVal$setVisible(TRUE)
            for (k in vsets) for (i in seq_len(nrow(defs)))
                r$perfVal$addRow(rowKey = paste(k, defs$key[i]),
                                 values = list(set = sets[[k]], group = defs$group[i], metric = defs$label[i]))
            private$.initCaseMix(sets)
            if (isTRUE(self$options$updating)) {
                r$updTable$setVisible(TRUE)
                meth <- private$.updMethods()
                for (k in vsets) for (m in names(meth))
                    r$updTable$addRow(rowKey = paste(k, m), values = list(set = sets[[k]], method = meth[[m]]))
            }
        }
        if ("ext" %in% names(sets)) {
            r$extCheck$setVisible(TRUE)
            for (v in c(private$.outcomeVars(), private$.preds())) r$extCheck$addRow(rowKey = v, values = list(var = v))
        }

        private$.initDesc()
        private$.initOR()
        cands <- private$.candidates()
        if (length(cands)) {
            t <- r$intTable
            t$setVisible(TRUE)
            for (ab in cands) t$addRow(rowKey = paste(ab, collapse = ":"),
                                       values = list(term = paste(ab, collapse = " × ")))
            t$getColumn("freq")$setVisible(val != "none" && self$options$intMode != "all")
            if (val == "cv") t$getColumn("freq")$setTitle(.("CV fits"))
            else t$getColumn("freq")$setTitle(.("Bootstrap fits"))
            r$subTable$setVisible(TRUE)
        }
        for (cn in c("lower", "upper")) r$subTable$getColumn(cn)$setSuperTitle(ciTitle)

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

        private$.initExtras(sets)
        if (self$options$sampleSize) {
            crit <- private$.ssRows()
            for (k in names(crit)) r$ssTable$addRow(rowKey = k, values = list(criterion = crit[[k]]))
        }
        impTables <- c(shapley = "shapleyTable", permImp = "permTable", shap = "shapTable")
        for (opt in names(impTables)) if (self$options[[opt]])
            for (v in private$.preds()) r[[impTables[[opt]]]]$addRow(rowKey = v, values = list(var = v))
        r$plots$shapPatient$setVisible(self$options$shap && self$options$shapRow > 0)

        # one facet per data set
        w <- 420 * length(sets) + 60
        for (p in c("calibration", "riskDist", "dcaCurves"))
            if (length(sets) > 1 && !is.null(r$plots[[p]])) r$plots[[p]]$setSize(w, 450)
        np <- length(private$.preds())
        for (p in c("shapleyBar", "permBar", "shapSummary", "shapPatient"))
            r$plots[[p]]$setSize(600, max(300, 120 + 32 * np))
    },

    .initCaseMix = function(sets) {
        t <- self$results$caseMix
        t$setVisible(TRUE)
        for (k in setdiff(names(sets), "dev")) t$getColumn(k)$setVisible(TRUE)
        t$addRow(rowKey = "n", values = list(var = .("N")))
        t$addRow(rowKey = "events", values = list(var = private$.text("outcomeRow")))
        for (v in c(self$options$covs, self$options$splineVars))
            t$addRow(rowKey = v, values = list(var = v, level = .("mean (SD)")))
        for (v in self$options$factors)
            for (l in levels(factor(self$data[[v]])))
                t$addRow(rowKey = paste(v, l), values = list(var = v, level = l))
        t$addRow(rowKey = "lp", values = list(var = .("Linear predictor"), level = .("mean (SD)")))
        t$addRow(rowKey = "member", values = list(var = .("Membership c-statistic")))
    },

    # effects table: one row per covariate, reference + one per level of a factor
    .initOR = function() {
        if (!self$options$oddsRatios) return()
        t <- self$results$orTable
        sc <- private$.scales()
        labs <- c(unit = .("per unit"), sd = .("per 1 SD"), ten = .("per 10 units"),
                  hundred = .("per 100 units"), thousand = .("per 1000 units"))
        for (v in self$options$covs) t$addRow(rowKey = v, values = list(var = v, level = labs[[sc[[v]] %||% "unit"]]))
        for (v in self$options$splineVars) t$addRow(rowKey = v, values = list(var = v, level = .("Q3 vs Q1")))
        for (v in self$options$factors) {
            ref <- private$.refOf(v)
            lev <- levels(factor(self$data[[v]]))
            t$addRow(rowKey = paste0(v, " (ref)"), values = list(var = v, level = ref))
            for (l in setdiff(lev, ref)) t$addRow(rowKey = paste0(v, ": ", l), values = list(var = v, level = l))
        }
    },

    # write a row by key, adding it when .init could not know it
    .putRow = function(tab, key, values) {
        if (key %in% tab$rowKeys) tab$setRow(rowKey = key, values = values)
        else tab$addRow(rowKey = key, values = values)
    },

    # ---- computation -------------------------------------------------------
    .run = function() .guard(self, private$.analyse()),

    .analyse = function() {
        if (!private$.ready()) return()
        private$.fam <- cp_family(private$.famName)
        fam <- private$.fam
        r <- self$results
        o <- self$options
        conf <- o$ciWidth / 100
        sets <- private$.sets()
        outVars <- private$.outcomeVars()

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
        if (o$descTable) private$.fillDesc(data[isDev, , drop = FALSE])
        cands <- private$.candidates()
        prep <- cp_prepare(data[isDev, , drop = FALSE], outVars, o$covs, o$splineVars, o$factors, o$knots,
                           refLevels = private$.refs(), scales = private$.scales(), interactions = cands)
        prep$nMissing <- prep$nMissing + sum(!isDev & !isVal)
        prep <- private$.setOutcome(prep)
        private$.devData <- prep$data
        y <- prep$y
        # candidate interactions (specC, XC); the model keeps the selected ones
        specC <- prep$spec
        XC <- cp_design(specC, prep$data)
        sel <- cp_select_interactions(XC, y, specC, fam, o$intMode, o$intAlpha)
        spec <- cp_with_interactions(specC, attr(sel, "selected"))
        X <- cp_design(spec, prep$data)
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

        # validation data sets, each with X, y, lp, p
        vsets <- list()
        if (any(isVal)) {
            sd <- cp_split_data(data[isVal, , drop = FALSE], spec, outVars)
            vsets$split <- private$.valSet(sd$data, spec, beta, "split")
            vsets$split$nMissing <- sd$nMissing
        }
        if (!is.null(o$extFile)) {
            ext <- cp_read_external(o$extFile$path, o$extFile$filename %||% o$extFile$path)
            h <- cp_harmonise(ext, spec, outVars, sets[["ext"]])
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
        lrt <- if (o$oddsRatios || o$predTests) cp_lr_tests(spec, X, y, fam) else NULL
        if (o$oddsRatios) private$.fillOR(spec, beta, prep, y, lrt, conf)
        if (o$predTests) private$.fillTests(lrt)

        iv <- private$.internal(XC, y, specC)
        if (length(cands)) private$.fillInt(XC, y, specC, sel, iv, conf, spec, beta)
        private$.fillPerfDev(y, lp, iv, conf)
        for (k in names(vsets)) private$.fillPerfVal(k, vsets[[k]], conf)
        if (length(vsets)) private$.fillCaseMix(prep, X, y, lp, vsets)
        if (isTRUE(o$updating)) for (k in names(vsets)) private$.fillUpd(k, vsets[[k]])
        if (!is.null(iv) && iv$method == "boot") private$.fillStab(iv, y, lp, X)

        all <- c(list(dev = list(X = X, y = y, lp = lp, p = p, df = prep$data)), vsets)
        private$.fillExtras(all, sets, conf)
        if (o$sampleSize) private$.fillSampleSize(y, lp, spec, iv)
        if (o$shapley) private$.fillShapley(X, y, spec)
        if (o$permImp) private$.fillPerm(beta, spec, all, sets)
        if (o$shap) private$.fillShap(beta, spec, X, data)
        private$.plotStates(all, sets, iv)
        private$.saveOutputs(beta, spec, data)
    },

    .valSet = function(df, spec, beta, key) {
        fam <- private$.fam
        y <- private$.yOf(df, spec, key)
        X <- cp_design(spec, df)
        lp <- fam$lp(beta, X)
        list(X = X, y = y, lp = lp, p = fam$risk(lp), df = df)
    },

    .fillSamples = function(y, nMissing, vsets) {
        t <- self$results$samples
        t$setRow(rowKey = "dev", values = private$.sampleValues(y, nMissing))
        for (k in names(vsets)) t$setRow(rowKey = k, values = private$.sampleValues(vsets[[k]]$y, vsets[[k]]$nMissing))
        if (self$options$explain) t$setNote("ev", private$.text("samplesEx"))
    },

    .fillCoefs = function(beta, spec, conf) {
        if (!self$options$coefs) return()
        fam <- private$.fam
        t <- self$results$coefTable
        V <- fam$vcov(beta)
        se <- sqrt(diag(V))
        df <- fam$dfResid(beta)
        z <- stats::qt(1 - (1 - conf) / 2, df)
        for (i in seq_along(beta)) {
            b <- unname(beta[i])
            v <- list(term = if (i == 1) .("Intercept") else names(beta)[i], est = b, se = se[i],
                      p = 2 * stats::pt(-abs(b / se[i]), df))
            # ratio column (odds ratio) only on a log scale; the CI then is of the ratio
            if (fam$exp) v <- c(v, list(or = if (i > 1) exp(b) else "", lower = if (i > 1) exp(b - z * se[i]) else "",
                                        upper = if (i > 1) exp(b + z * se[i]) else ""))
            else v <- c(v, list(lower = b - z * se[i], upper = b + z * se[i]))
            private$.putRow(t, names(beta)[i], v)
        }
        refs <- vapply(Filter(function(i) i$role == "factor", spec$vars),
                       function(i) paste0(i$name, " = ", i$levels[1]), "")
        if (length(refs)) t$setNote("ref", paste0(.("Reference levels"), ": ", paste(refs, collapse = "; "), "."))
        if (length(self$options$splineVars)) t$setNote("spl", private$.text("coefSpline"))
    },

    .fillOR = function(spec, beta, prep, y, lrt, conf) {
        t <- self$results$orTable
        fam <- private$.fam
        adj <- cp_effects(spec, beta, fam$vcov(beta), conf, exp = fam$exp, df = fam$dfResid(beta))
        cr <- if (self$options$crude) cp_crude(spec, prep$data, y, fam, conf) else NULL
        pAdj <- stats::setNames(lrt$p, lrt$key)
        ref <- if (fam$exp) 1 else 0
        for (info in spec$vars) {
            v <- info$name
            if (info$role == "factor")
                private$.putRow(t, paste0(v, " (ref)"), list(var = v, level = info$levels[1],
                    or_u = ref, lower_u = "", upper_u = "", p_u = if (!is.null(cr)) attr(cr, "pOverall")[[v]] else "",
                    or = ref, lower = "", upper = "", p = pAdj[[v]]))
            for (i in which(adj$var == v)) {
                a <- adj[i, ]
                vals <- list(var = v, level = a$level, or = a$est, lower = a$lower, upper = a$upper,
                             p = if (info$role == "spline") pAdj[[v]] else a$p)
                if (!is.null(cr)) {
                    cc <- cr[cr$key == a$key, ]
                    vals <- c(vals, list(or_u = cc$est, lower_u = cc$lower, upper_u = cc$upper,
                                         p_u = if (info$role == "spline") attr(cr, "pOverall")[[v]] else cc$p))
                }
                private$.putRow(t, a$key, vals)
            }
        }
        if (self$options$crude)
            t$setNote("crude", .("Crude: each predictor alone, on the same rows. Adjusted: in the prediction model."))
        if (length(self$options$splineVars)) t$setNote("spl", private$.text("orSpline"))
        if (length(self$options$factors)) t$setNote("ref", private$.text("orRef"))
        if (length(spec$interactions)) t$setNote("int", private$.text("orInt"))
    },

    .fillTests = function(lrt) {
        t <- self$results$testTable
        for (i in seq_len(nrow(lrt)))
            private$.putRow(t, lrt$key[i], list(term = lrt$term[i], chisq = lrt$chisq[i],
                                                df = lrt$df[i], df2 = lrt$df2[i], p = lrt$p[i]))
        t$setNote("ex", .("Each test drops the term from the model (a predictor with its non-linear and interaction terms)."))
    },

    # screening of each candidate added alone to the main effects, the
    # selection, and how often resampled fits selected it
    .fillInt = function(XC, y, specC, sel, iv, conf, spec, beta) {
        o <- self$options
        fam <- private$.fam
        t <- self$results$intTable
        screen <- attr(sel, "screen") %||%
            attr(cp_select_interactions(XC, y, specC, fam, "lrt", 0), "screen")
        kept <- attr(sel, "selected")
        for (i in seq_len(NROW(screen))) {
            k <- screen$term[i]
            step <- match(k, kept)
            private$.putRow(t, k, list(chisq = screen$chisq[i], df = screen$df[i], p = screen$p[i],
                daic = screen$daic[i], dbic = screen$dbic[i],
                kept = if (o$intMode == "all") .("yes") else if (is.na(step)) .("no") else sprintf(.("step %d"), step),
                freq = if (!is.null(iv$selection)) iv$selection[[k]] %||% 0 else ""))
        }
        t$setNote("ex", .("χ², ΔAIC and ΔBIC: the term added alone to the main-effects model (negative Δ: the term improves the criterion)."))
        if (o$intMode != "all") {
            t$setNote("sel", switch(o$intMode,
                aic = .("Forward selection: at each step the term that lowers the AIC most enters, while the AIC decreases."),
                bic = .("Forward selection: at each step the term that lowers the BIC most enters, while the BIC decreases."),
                lrt = sprintf(.("Forward selection: at each step the term with the smallest likelihood ratio p enters, while p < %s."), .num(o$intAlpha))))
            if (o$validation != "none")
                t$setNote("res", .("The selection is repeated in every resampled fit, so the internal validation includes its optimism; the last column is the share of fits that selected the term."))
        }
        sub <- if (length(spec$interactions))
            cp_subgroups(spec, beta, fam$vcov(beta), conf, exp = fam$exp, df = fam$dfResid(beta)) else NULL
        st <- self$results$subTable
        if (is.null(sub)) {
            st$setNote("none", .("No interaction is in the model."))
            return()
        }
        for (i in seq_len(nrow(sub)))
            private$.putRow(st, sub$key[i], list(effect = sub$effect[i], within = sub$within[i],
                or = sub$est[i], lower = sub$lower[i], upper = sub$upper[i], p = sub$p[i]))
        st$setNote("ex", private$.text("subEx"))
    },

    # how each predictor enters the equation (scaling, dummies, spline basis)
    .equationVars = function(spec) {
        out <- character(0)
        for (info in spec$vars) {
            if (info$role == "numeric" && info$scale != 1)
                out <- c(out, "", sprintf("[%s] = %s / %s  (%s)", info$name, info$name, .num(info$scale, 5), info$scaleLabel))
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
        out
    },

    # the terms of a linear predictor, one per line
    .equationTerms = function(beta) {
        b <- beta
        c(sprintf("    %.5g", b[1]),
          paste0("    ", vapply(seq_along(b)[-1], function(i)
              sprintf("%s %.5g × [%s]", if (b[i] < 0) "−" else "+", abs(b[i]), names(b)[i]), "")))
    },

    # internal validation, cached in the performance table's state
    .internal = function(X, y, specC) {
        o <- self$options
        if (o$validation == "none") return(NULL)
        iv <- self$results$perfDev$state
        if (is.null(iv)) {
            fam <- private$.fam
            if (length(specC$interactions) && o$intMode != "all") {
                base <- fam
                fam$fit <- function(X, y) cp_select_interactions(X, y, specC, base, o$intMode, o$intAlpha)
            }
            iv <- cp_internal(X, y, fam, o$validation, B = o$nBoot, k = o$nFolds,
                              R = o$nRepeats, seed = o$seed, stability = o$validation == "boot")
            self$results$perfDev$setState(iv)
        }
        iv
    },

    .fillPerfDev = function(y, lp, iv, conf) {
        t <- self$results$perfDev
        app <- private$.fam$metrics(y, lp)
        ci <- private$.fam$ci(y, lp, conf)
        for (k in t$rowKeys) {
            # apparent calibration-in-the-large is 0 by construction; print it as 0, not as rounding noise
            a <- if (abs(app[[k]]) < 1e-8) 0 else app[[k]]
            v <- list(apparent = a, lower = ci[[k]][1] %||% "", upper = ci[[k]][2] %||% "")
            if (!is.null(iv)) { v$optimism <- iv$optimism[[k]]; v$corrected <- iv$corrected[[k]] }
            t$setRow(rowKey = k, values = v)
        }
        if (!is.null(iv)) {
            t$setNote("iv", if (iv$method == "boot")
                sprintf(.("Bootstrap optimism correction (Harrell): %d samples (%d usable); corrected = apparent − optimism."), iv$B, iv$nOk)
            else sprintf(private$.text("cvNote"), iv$k, iv$R))
        } else {
            t$setNote("iv", .("Apparent performance is optimistic: it is measured in the data used to fit the model. Choose bootstrap or cross-validation under Validation."))
        }
        neg <- intersect(c("r2n", "r2cs", "sbrier", "r2"), names(app))
        if (!is.null(iv) && any(iv$corrected[neg] < 0, na.rm = TRUE)) t$setNote("neg", private$.text("negNote"))
        if (self$options$explain) t$setNote("ex", private$.text("perfEx"))
    },

    .fillPerfVal = function(key, vs, conf) {
        t <- self$results$perfVal
        m <- private$.fam$metrics(vs$y, vs$lp)
        ci <- private$.fam$ci(vs$y, vs$lp, conf)
        for (k in names(m)) {
            rk <- paste(key, k)
            if (!(rk %in% t$rowKeys)) next
            t$setRow(rowKey = rk, values = list(est = m[[k]], lower = ci[[k]][1] %||% "",
                                                upper = ci[[k]][2] %||% ""))
        }
        if (self$options$explain) t$setNote("ex", private$.text("valEx"))
    },

    # model updating in one validation set, with the closed testing procedure
    .fillUpd = function(key, vs) {
        t <- self$results$updTable
        u <- cp_update(private$.fam, vs$X, vs$y, vs$lp, self$options$updAlpha)
        for (m in names(private$.updMethods())) {
            mm <- u$methods[[m]]
            if (is.null(mm)) next
            perf <- private$.fam$metrics(vs$y, mm$lp)[private$.updKeys]
            vals <- c(list(a = mm$a, b = mm$b, k = mm$k,
                           stat = mm$test[["stat"]] %||% "", df = mm$test[["df"]] %||% "", p = mm$test[["p"]] %||% "",
                           chosen = if (identical(u$chosen, m)) "✓" else ""),
                      as.list(perf))
            t$setRow(rowKey = paste(key, m), values = vals)
        }
        if (is.na(u$chosen))
            t$setNote("fail", .("The revised model could not be fitted in the validation data (too few events or separation), so the closed testing procedure was not run."))
        t$setNote("ex", sprintf(private$.text("updNote"), .num(self$options$updAlpha)))
    },

    .fillCaseMix = function(prep, X, y, lp, vsets) {
        t <- self$results$caseMix
        fmt <- function(m, s) sprintf("%.2f (%.2f)", m, s)
        col <- function(k, df, yy, lpk, Xk) {
            vals <- list()
            vals$n <- as.character(NROW(yy))
            vals$events <- private$.outcomeSummary(yy)
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
        fam <- private$.fam
        t <- self$results$stabTable
        lr <- 2 * (fam$nullLoss(y) - fam$loss(y, lp))
        t$setRow(rowKey = "shrink", values = list(value = iv$shrinkage))
        t$setRow(rowKey = "heur", values = list(value = (lr - ncol(X)) / lr))
        if (!is.null(iv$stability) && "mape" %in% t$rowKeys) {
            t$setRow(rowKey = "mape", values = list(value = iv$stability$mape))
            t$setRow(rowKey = "maxape", values = list(value = iv$stability$maxape))
        }
        if (self$options$explain) t$setNote("ex", private$.text("stabEx"))
    },

    .fillShapley = function(X, y, spec) {
        t <- self$results$shapleyTable
        sh <- t$state
        measure <- self$options$shapleyMeasure %||% "r2"
        if (is.null(sh)) {
            sh <- cp_shapley_r2(X, y, spec, private$.fam, measure, seed = self$options$seed)
            t$setState(sh)
        }
        for (i in seq_len(nrow(sh)))
            t$setRow(rowKey = sh$var[i], values = list(shapley = sh$shapley[i], share = sh$share[i], rank = sh$rank[i]))
        lab <- switch(measure, r2mf = "McFadden R²", r2n = "Nagelkerke R²", "R²")
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
        fam <- private$.fam
        phi <- cp_shap(beta, spec, X, intercept = fam$intercept)
        ma <- colMeans(abs(phi))
        t <- self$results$shapTable
        for (v in names(ma)) t$setRow(rowKey = v, values = list(meanAbs = ma[[v]], share = ma[[v]] / sum(ma)))
        t$setNote("ex", private$.text("shapEx"))
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
            ph <- cp_shap(beta, spec, Xi, Xbg = X, intercept = fam$intercept)
            labs <- vapply(names(spec$vars), function(v) {
                x <- nd[rk, v]
                paste0(v, " = ", if (is.factor(x)) as.character(x) else .num(x, 4))
            }, "")
            lpi <- fam$lp(beta, Xi)
            base <- attr(ph, "base")
            img$setState(list(label = unname(labs), phi = unname(ph[1, ]), base = base,
                              lp = lpi, risk = fam$risk(lpi), axis = private$.shapAxis(base, lpi),
                              title = sprintf(.("Row %d"), row)))
        }
    },

    .plotStates = function(all, sets, iv) {
        pl <- self$results$plots
        o <- self$options
        conf <- o$ciWidth / 100
        if (o$calPlot && is.null(pl$calibration$state)) pl$calibration$setState(private$.calState(all, sets, conf))
        private$.extraPlotStates(all, sets, conf)
        if (!is.null(iv$stability)) {
            if (is.null(pl$instability$state)) pl$instability$setState(iv$stability)
            if (is.null(pl$calInstability$state)) pl$calInstability$setState(iv$stability)
        }
    },

    # predictions written to the spreadsheet, for every row with complete predictors
    .saveOutputs = function(beta, spec, data) {
        defs <- private$.outputDefs()
        on <- Filter(function(k) isTRUE(self$options[[k]]), names(defs))
        if (!length(on)) return()
        nd <- cp_newdata(data, spec)
        lp <- rep(NA_real_, nrow(data)); names(lp) <- rownames(data)
        lp[rownames(nd)] <- private$.fam$lp(beta, cp_design(spec, nd))
        for (k in on) {
            out <- self$results[[k]]
            if (!out$isNotFilled()) next
            out$setRowNums(rownames(data))
            out$setValues(unname(defs[[k]](lp)))
        }
    },

    # ---- plots -------------------------------------------------------------
    .plotShapley = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_shapley(image$state$df, image$state$label, ggtheme)); TRUE
    },
    .plotPerm = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        fam <- cp_family(private$.famName)
        print(cp_plot_perm(image$state, ggtheme, fam$permKeys[1], private$.text("permAxis"))); TRUE
    },
    .plotShapSummary = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_shap_summary(image$state, ggtheme, private$.text("shapAxis"))); TRUE
    },
    .plotShapPatient = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_shap_patient(image$state, ggtheme)); TRUE
    },
    .plotInstability = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_instability(image$state, ggtheme, private$.text("instabX"), private$.text("instabY"))); TRUE
    },
    .plotCalInstability = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_cal_instability(image$state, ggtheme, private$.text("calX"), private$.text("calY"))); TRUE
    }
)
