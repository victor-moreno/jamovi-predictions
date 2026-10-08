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

    # an optional switch: FALSE when this analysis has no such option
    .on = function(name) self$options$has(name) && isTRUE(self$options$get(name)),

    # the family object; Cox needs the horizon, known once the outcome is read
    .makeFamily = function(prep) cp_family(private$.famName),

    # Levels of a factor in the data. In jamovi, .init() sees only the data
    # header (0 rows), and factor() would drop every level of an empty factor.
    .levelsOf = function(var) {
        x <- self$data[[var]]
        if (is.factor(x)) levels(x) else levels(factor(x))
    },

    # the level chosen in a LevelSelector, or the second level of the variable
    .level = function(var, chosen) {
        if (!is.null(chosen) && nzchar(chosen)) return(chosen)
        lev <- private$.levelsOf(var)
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
    .refOf = function(v) private$.refs()[[v]] %||% private$.levelsOf(v)[1],

    # candidate interactions: the listed two-way terms, or every pair of predictors
    .candidates = function() {
        preds <- private$.preds()
        if (self$options$intAllPairs)
            return(if (length(preds) < 2) list() else utils::combn(preds, 2, simplify = FALSE))
        ints <- lapply(self$options$interactions, function(t) unique(unlist(t)))
        Filter(function(t) length(t) == 2 && all(t %in% preds), ints)
    },

    # the fitting procedure's options: baseline (forced) covariates, selection
    # of the others, interaction selection
    .selOpts = function() {
        o <- self$options
        forced <- character(0)
        for (it in o$baseline) if (isTRUE(it$forced)) forced <- c(forced, it$var)
        list(mode = o$selMode, crit = o$selCrit, alpha = o$selAlpha,
             forced = intersect(forced, private$.preds()), intMode = o$intMode, intAlpha = o$intAlpha,
             enetAlpha = o$enetAlpha, enetRule = o$enetRule, enetFolds = o$enetFolds, seed = o$seed,
             relax = isTRUE(o$enetRelax))
    },
    .selecting = function() {
        s <- private$.selOpts()
        s$mode != "none" || (length(private$.candidates()) > 0 && s$intMode != "all")
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
                               recal = .("Recalibration"), revision = .("Model revision")),

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
            if (private$.on("updating")) {
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

        if (self$options$selMode != "none") {
            st <- r$selTable
            st$setVisible(TRUE)
            for (v in private$.preds()) st$addRow(rowKey = v, values = list(var = v))
            st$getColumn("freq")$setVisible(val != "none")
            st$getColumn("freq")$setTitle(if (val == "cv") .("CV fits") else .("Bootstrap fits"))
            r$plots$enetCv$setVisible(self$options$selMode == "enet")
            for (cn in c("step", "chisq", "df", "p")) st$getColumn(cn)$setVisible(self$options$selMode != "enet")
        }

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

        if (self$options$incremental) private$.initInc(sets)
        if (self$options$missExplore) private$.initMiss()
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
        # every cell of the two heading columns gets a value: jamovi spans a heading
        # cell over blank ones below it but still draws them, shifting those rows
        t$addRow(rowKey = "n", values = list(var = .("Patients"), level = "n"))
        t$addRow(rowKey = "events", values = list(var = private$.text("outcomeRow"), level = private$.text("outcomeLevel")))
        for (v in c(self$options$covs, self$options$splineVars))
            t$addRow(rowKey = v, values = list(var = v, level = .("mean (SD)")))
        for (v in self$options$factors)
            for (l in private$.levelsOf(v))
                t$addRow(rowKey = paste(v, l), values = list(var = v, level = l))
        t$addRow(rowKey = "lp", values = list(var = .("Linear predictor"), level = .("mean (SD)")))
        t$addRow(rowKey = "member", values = list(var = .("Membership c-statistic"), level = "c"))
    },

    # effects table: one row per covariate, reference + one per level of a factor
    # predictors in table order: baseline covariates last
    .orderedPreds = function() {
        f <- private$.selOpts()$forced
        c(setdiff(private$.preds(), f), intersect(private$.preds(), f))
    },

    # effects table: one row per covariate, reference + one per level of a factor;
    # baseline covariates last, below a line; crude, adjusted and (after a
    # selection) selected-model blocks; CI and p columns optional
    .initOR = function() {
        o <- self$options
        if (!o$oddsRatios) return()
        t <- self$results$orTable
        sc <- private$.scales()
        labs <- c(unit = .("per unit"), sd = .("per 1 SD"), ten = .("per 10 units"),
                  hundred = .("per 100 units"), thousand = .("per 1000 units"))
        forced <- private$.selOpts()$forced
        first <- TRUE; nEff <- 0
        for (v in private$.orderedPreds()) {
            keys <- character(0)
            if (v %in% o$covs) { t$addRow(rowKey = v, values = list(var = v, level = labs[[sc[[v]] %||% "unit"]])); keys <- v; nEff <- nEff + 1 }
            else if (v %in% o$splineVars) { t$addRow(rowKey = v, values = list(var = v, level = .("Q3 vs Q1"))); keys <- v; nEff <- nEff + 1 }
            else {
                ref <- private$.refOf(v)
                lev <- private$.levelsOf(v)
                t$addRow(rowKey = paste0(v, " (ref)"), values = list(var = v, level = ref)); keys <- paste0(v, " (ref)")
                for (l in setdiff(lev, ref)) t$addRow(rowKey = paste0(v, ": ", l), values = list(var = v, level = l))
                nEff <- nEff + max(1, length(lev) - 1)
            }
            if (v %in% forced && first) {
                if (length(setdiff(private$.preds(), forced))) t$addFormat(rowKey = keys[1], col = 1, format = jmvcore::Cell.BEGIN_GROUP)
                first <- FALSE
            }
        }
        if (length(forced)) t$setNote("forced", .("Below the line: baseline covariates (always in the model)."))
        sel <- private$.selecting()
        penal <- o$selMode == "enet" && !isTRUE(o$enetRelax)
        vis <- function(cn, on) t$getColumn(cn)$setVisible(on)
        for (b in c("_u", "", "_s")) {
            on <- switch(b, `_u` = o$crude, `_s` = sel, TRUE)
            vis(paste0("or", b), on)
            vis(paste0("lower", b), on && o$orCI && !(b == "_s" && penal))
            vis(paste0("upper", b), on && o$orCI && !(b == "_s" && penal))
            vis(paste0("p", b), on && o$orP && !(b == "_s" && penal))
        }
        if (sel) for (cn in c("or", "lower", "upper", "p")) t$getColumn(cn)$setSuperTitle(.("Adjusted (all predictors)"))
        if (o$forest) self$results$plots$forestPlot$setSize(620, max(300, 120 + 26 * nEff * (1 + sel)))
    },


    # write a row by key, adding it when .init could not know it
    .putRow = function(tab, key, values) {
        # a missing number is a blank cell (jamovi would show NaN)
        values <- lapply(values, function(v) if (length(v) == 1 && !is.character(v) && is.na(v)) "" else v)
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
        if (o$missExplore) private$.fillMiss(data[isDev, , drop = FALSE])
        if (o$missMethod == "mi") return(private$.analyseMI(data, isDev, isVal, sets, conf))
        cands <- private$.candidates()
        # missing values: complete cases, missing category / indicator, or single imputation
        devData <- data[isDev, , drop = FALSE]
        covs <- o$covs
        private$.mh <- NULL
        if (o$missMethod == "category") {
            private$.mh <- cp_miss_category(devData, c(o$covs, o$splineVars), o$factors)
            devData <- cp_miss_category_apply(private$.mh, devData)
            covs <- c(covs, unname(private$.mh$indicator))
        } else if (o$missMethod == "single") {
            private$.nImputed <- private$.countImputed(devData)
            devData <- private$.impSet(devData, 1, o$seed)[[1]]
        }
        prep <- cp_prepare(devData, outVars, covs, o$splineVars, o$factors, o$knots,
                           refLevels = private$.refs(), scales = private$.scales(), interactions = cands)
        prep$nMissing <- prep$nMissing + sum(isDev) - nrow(devData) + sum(!isDev & !isVal)
        prep <- private$.setOutcome(prep)
        fam <- private$.fam <- private$.makeFamily(prep)
        private$.devData <- prep$data
        y <- prep$y
        # candidate interactions (specC, XC); the model keeps the selected ones
        specC <- prep$spec
        XC <- cp_design(specC, prep$data)
        sel <- cp_fit_procedure(XC, y, specC, fam, private$.selOpts())
        if (!length(attr(sel, "vars")))
            .fail("The selection removed every predictor. Mark baseline covariates under Variable Selection, ",
                  "or use a less strict criterion.")
        spec <- cp_with_terms(specC, attr(sel, "vars"), attr(sel, "ints"))
        penalized <- isTRUE(attr(sel, "penalized"))
        X <- cp_design(spec, prep$data)
        # a penalised model is its penalised coefficients; any other is refitted on the kept terms
        beta <- if (penalized) {
            b <- sel[c("(Intercept)", spec$columns)]
            attributes(b) <- c(attributes(b)["names"], attributes(sel)[c("converged", "penalized", "lambda", "cv")])
            b
        } else fam$fit(X, y)
        if (anyNA(beta))
            .fail("The predictor column '", names(beta)[is.na(beta)][1], "' is redundant (a combination of ",
                  "other predictors); remove one of the predictors involved.")
        if (!isTRUE(attr(beta, "converged")))
            .fail("The model did not converge. This usually means complete separation: some predictor, or a ",
                  "combination of them, predicts the outcome perfectly. Remove or merge the predictor(s) involved.")
        if (isTRUE(attr(beta, "separation")))
            .show_warning(self, "separation", private$.text("separation"))
        lp <- fam$lp(beta, X); p <- fam$risk(lp)

        # validation data sets, each with X, y, lp, p
        vsets <- list()
        cc <- o$missMethod == "cc"
        specRaw <- private$.specRaw(spec)
        if (any(isVal)) {
            sd <- cp_split_data(data[isVal, , drop = FALSE], specRaw, outVars, complete = cc)
            sd2 <- cp_split_data(private$.missVal(sd$data), spec, outVars)
            vsets$split <- private$.valSet(sd2$data, spec, beta, "split")
            vsets$split$nMissing <- sum(isVal) - nrow(sd2$data)
        }
        if (!is.null(o$extFile)) {
            ext <- cp_read_external(o$extFile$path, o$extFile$filename %||% o$extFile$path)
            h <- cp_harmonise(ext, specRaw, outVars, sets[["ext"]], complete = cc)
            for (i in seq_len(nrow(h$check)))
                r$extCheck$setRow(rowKey = h$check$var[i],
                                  values = list(role = h$check$role[i], status = h$check$status[i],
                                                details = h$check$details[i]))
            e2 <- cp_split_data(private$.missVal(h$data), spec, outVars)
            vsets$ext <- private$.valSet(e2$data, spec, beta, "ext")
            vsets$ext$nMissing <- h$nTotal - nrow(e2$data)
        }

        private$.fillSamples(y, prep$nMissing, vsets)
        private$.missNote()
        if (length(spec$outcome$levels) > 2)
            r$samples$setNote("lev", sprintf(.("The outcome has %d levels: '%s' is modelled against all the others combined."),
                                             length(spec$outcome$levels), spec$outcome$event))
        private$.fillCoefs(beta, spec, conf)
        if (o$equation) r$equationText$setContent(private$.equation(beta, spec))
        lrt <- if ((o$oddsRatios || o$predTests) && !penalized) cp_lr_tests(spec, X, y, fam) else NULL
        if (o$oddsRatios) private$.fillOR(spec, beta, prep, y, lrt, conf,
                                          full = if (private$.selecting()) private$.fullModel(XC, y, specC) else NULL)
        if (o$predTests) private$.fillTests(lrt)

        iv <- private$.internal(XC, y, specC, attr(sel, "lambda"))
        if (length(cands)) private$.fillInt(XC, y, specC, sel, iv, conf, spec, beta)
        if (o$selMode != "none") private$.fillSel(sel, iv)
        private$.postSelection(penalized, iv, beta, spec)
        private$.fillPerfDev(y, lp, iv, conf)
        for (k in names(vsets)) private$.fillPerfVal(k, vsets[[k]], conf)
        if (length(vsets)) private$.fillCaseMix(prep, X, y, lp, vsets)
        if (private$.on("updating")) for (k in names(vsets)) private$.fillUpd(k, vsets[[k]])
        # shrinkage and sample size count the candidate parameters (selection included)
        if (o$incremental) private$.fillInc(prep, specC, spec, y, lp, iv, vsets, conf, beta)
        if (!is.null(iv) && iv$method == "boot") private$.fillStab(iv, y, lp, XC)

        all <- c(list(dev = list(X = X, y = y, lp = lp, p = p, df = prep$data)), vsets)
        private$.fillExtras(all, sets, conf)
        if (o$sampleSize) private$.fillSampleSize(y, lp, specC, iv)
        if (o$shapley) private$.fillShapley(X, y, spec, penalized)
        if (o$permImp) private$.fillPerm(beta, spec, all, sets)
        if (o$shap) private$.fillShap(beta, spec, X, data)
        dropped <- setdiff(private$.preds(), names(spec$vars))
        if (length(dropped)) private$.zeroImportance(dropped)
        private$.plotStates(all, sets, iv)
        private$.saveOutputs(beta, spec, data)
        # a failed export must not hide the results
        if (private$.on("exportModel"))
            tryCatch(private$.exportModel(beta, spec, y, lp, iv),
                     error = function(e) .show_warning(self, "export", .("The model could not be exported.")))
    },

    # the fitted model as JSON, handed to the host (desktop: opened; cloud: downloaded)
    .exportModel = function(beta, spec, y, lp, iv) {
        fam <- private$.fam
        app <- fam$metrics(y, lp)
        dev <- list(n = NROW(y),
                    events = if (is.matrix(y)) sum(y[, 2]) else if (fam$name == "binary") sum(y) else NULL,
                    validation = if (is.null(iv)) "none" else iv$method,
                    apparent = as.list(app), corrected = if (!is.null(iv)) as.list(iv$corrected))
        json <- cp_model_json(spec, beta, fam$name, se = sqrt(pmax(0, diag(fam$vcov(beta)))),
                              version = as.character(utils::packageVersion("clinpred")), development = dev)
        file <- paste0("clinpred_", gsub("[^A-Za-z0-9_.-]", "_", spec$outcome$name), "_", fam$name, ".json")
        self$options$option("exportModel")$perform(function(action) {
            writeLines(json, action$params$fullPath, useBytes = TRUE)
            list(filename = file)
        })
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
        V <- private$.vcov(beta)
        se <- if (is.null(V)) rep(NA_real_, length(beta)) else sqrt(diag(V))
        df <- fam$dfResid(beta)
        z <- stats::qt(1 - (1 - conf) / 2, df)
        for (i in seq_along(beta)) {
            # the Cox "intercept" (log baseline cumulative hazard) is not a model coefficient
            if (i == 1 && identical(fam$showIntercept, FALSE)) next
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
        if (is.null(V)) t$setNote("pen", .("Penalised (elastic net) estimates: they have no standard errors, confidence intervals or p values; their uncertainty is in the internal validation."))
    },

    # full: list(spec, beta, V, lrt) of the model with every candidate predictor,
    # when a selection made the final model (spec, beta, lrt) a reduced one
    .fillOR = function(spec, beta, prep, y, lrt, conf, crude = NULL, full = NULL) {
        o <- self$options
        t <- self$results$orTable
        fam <- private$.fam
        ref <- if (fam$exp) 1 else 0
        cr <- if (o$crude || (o$forest && o$forestCrude)) crude %||% cp_crude(prep$spec, prep$data, y, fam, conf) else NULL
        blockEff <- function(sp, b, V, lr) list(eff = cp_effects(sp, b, V, conf, exp = fam$exp, df = fam$dfResid(b)),
                                                p = if (is.null(lr)) list() else as.list(stats::setNames(lr$p, lr$key)), vars = names(sp$vars))
        sel <- private$.selecting()
        fin <- blockEff(spec, beta, private$.vcov(beta), lrt)
        adj <- if (sel && !is.null(full)) blockEff(full$spec, full$beta, full$V, full$lrt) else if (sel) NULL else fin
        forced <- private$.selOpts()$forced
        rows <- list()
        put <- function(key, vals) private$.putRow(t, key, vals)
        # table order: baseline covariates last (each panel of the forest plot then holds its own rows)
        ord <- c(setdiff(names(prep$spec$vars), forced), intersect(names(prep$spec$vars), forced))
        for (info in prep$spec$vars[ord]) {
            v <- info$name
            if (info$role == "factor")
                put(paste0(v, " (ref)"), list(var = v, level = info$levels[1],
                    or_u = ref, lower_u = "", upper_u = "", p_u = if (!is.null(cr)) attr(cr, "pOverall")[[v]] else "",
                    or = if (!is.null(adj) && v %in% adj$vars) ref else "", lower = "", upper = "", p = if (!is.null(adj)) adj$p[[v]] %||% "" else "",
                    or_s = if (sel && v %in% fin$vars) ref else "", lower_s = "", upper_s = "", p_s = if (sel && v %in% fin$vars) fin$p[[v]] %||% "" else ""))
            for (d in cp_effect_defs(prep$spec, v)) {
                vals <- list(var = v, level = d$level)
                pick <- function(bl, suffix) {
                    if (is.null(bl) || !(v %in% bl$vars)) return(stats::setNames(list("", "", "", ""), paste0(c("or", "lower", "upper", "p"), suffix)))
                    a <- bl$eff[bl$eff$key == d$key, ]
                    stats::setNames(list(a$est, a$lower, a$upper, if (info$role == "spline") bl$p[[v]] %||% "" else a$p),
                                    paste0(c("or", "lower", "upper", "p"), suffix))
                }
                vals <- c(vals, pick(adj, ""), if (sel) pick(fin, "_s"))
                if (!is.null(cr)) {
                    cc <- cr[cr$key == d$key, ]
                    vals <- c(vals, list(or_u = cc$est, lower_u = cc$lower, upper_u = cc$upper,
                                         p_u = if (info$role == "spline") attr(cr, "pOverall")[[v]] else cc$p))
                }
                put(d$key, vals)
                # rows of the forest plot
                lab <- if (info$role == "factor") sprintf("%s: %s (vs %s)", v, d$level, info$levels[1]) else sprintf("%s (%s)", v, d$level)
                add <- function(block, est, lo, hi) if (length(est) && !identical(est, "") && is.finite(as.numeric(est)))
                    rows[[length(rows) + 1]] <<- data.frame(label = lab, block = block, est = as.numeric(est),
                        lower = suppressWarnings(as.numeric(lo)), upper = suppressWarnings(as.numeric(hi)),
                        forced = v %in% forced, stringsAsFactors = FALSE)
                if (!is.null(cr)) add("Crude", vals$or_u, vals$lower_u, vals$upper_u)
                add(if (sel) "Adjusted (all predictors)" else "Adjusted", vals$or, vals$lower, vals$upper)
                if (sel) add("Selected model", vals$or_s, vals$lower_s, vals$upper_s)
            }
        }
        if (o$crude)
            t$setNote("crude", .("Crude: each predictor alone, on the same rows. Adjusted: in the prediction model."))
        if (sel) {
            t$setNote("crude", .("Crude: each predictor alone. Adjusted (all predictors): the model with every candidate predictor, before selection. Selected model: the prediction model."))
            if (is.null(full)) t$setNote("full", .("The model with all predictors could not be fitted (e.g. separation)."))
        }
        if (length(self$options$splineVars)) t$setNote("spl", private$.text("orSpline"))
        if (length(self$options$factors)) t$setNote("ref", private$.text("orRef"))
        if (length(spec$interactions)) t$setNote("int", private$.text("orInt"))
        img <- self$results$plots$forestPlot
        if (o$forest && length(rows) && is.null(img$state)) {
            df <- do.call(rbind, rows)
            img$setState(list(rows = df, exp = fam$exp, order = unique(df$label),
                              xlab = private$.text("forestAxis")))
        }
    },


    .fillTests = function(lrt) {
        t <- self$results$testTable
        if (is.null(lrt)) {
            t$setNote("ex", .("Likelihood ratio tests do not apply to a penalised (elastic net) model."))
            return()
        }
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
        # candidates are screened among the predictors the selection kept
        specK <- cp_with_terms(specC, attr(sel, "vars"), names(specC$interactions))
        screen <- attr(sel, "screen") %||%
            attr(cp_select_interactions(XC[, specK$columns, drop = FALSE], y, specK, fam, "lrt", 0), "screen")
        kept <- attr(sel, "ints")
        for (k in setdiff(names(specC$interactions), names(specK$interactions)))
            private$.putRow(t, k, list(chisq = "", df = "", p = "", daic = "", dbic = "", kept = .("main effect removed"),
                                       freq = if (!is.null(iv$selection)) .sel_freq(iv$selection, k) else ""))
        for (i in seq_len(NROW(screen))) {
            k <- screen$term[i]
            step <- match(k, kept)
            private$.putRow(t, k, list(chisq = screen$chisq[i], df = screen$df[i], p = screen$p[i],
                daic = screen$daic[i], dbic = screen$dbic[i],
                kept = if (o$intMode == "all") .("yes") else if (is.na(step)) .("no") else sprintf(.("step %d"), step),
                freq = if (!is.null(iv$selection)) .sel_freq(iv$selection, k) else ""))
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
            cp_subgroups(spec, beta, private$.vcov(beta), conf, exp = fam$exp, df = fam$dfResid(beta)) else NULL
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

    # selection of the non-baseline predictors: decision, step, test against the
    # final main-effects model, share of resampled fits that kept it
    .fillSel = function(sel, iv) {
        o <- self$options
        t <- self$results$selTable
        tt <- attr(sel, "mainTests")
        forced <- private$.selOpts()$forced
        for (v in forced) private$.putRow(t, v, list(var = v, status = .("baseline"), step = "", chisq = "", df = "", p = "",
                                                     freq = if (!is.null(iv$selection)) 1 else ""))
        for (i in seq_len(NROW(tt))) {
            v <- tt$var[i]
            st <- if (o$selMode == "enet") { if (tt$kept[i]) .("kept") else .("shrunk to 0") } else
                if (tt$kept[i]) { if (is.na(tt$step[i])) .("kept") else .("added") } else
                { if (is.na(tt$step[i])) .("not added") else .("removed") }
            private$.putRow(t, v, list(var = v, status = st, step = if (is.na(tt$step[i])) "" else tt$step[i],
                chisq = tt$chisq[i], df = tt$df[i], p = tt$p[i],
                freq = if (!is.null(iv$selection)) .sel_freq(iv$selection, v) else ""))
        }
        if (o$selMode == "enet") {
            lam <- attr(sel, "lambda"); cv <- attr(sel, "cv")
            if (isTRUE(attr(sel, "relaxed")))
                t$setNote("relax", .("Relaxed lasso: the model is the unpenalised refit of the terms the elastic net kept (also in every resampled fit)."))
            t$setNote("how", sprintf(.("Elastic net, mixing α = %s, λ = %s (%s, %d-fold cross-validation); baseline covariates are not penalised. A predictor is kept when any of its coefficients is not shrunk to 0; all candidate interactions enter the penalised model."),
                                     .num(o$enetAlpha), .num(lam, 4),
                                     if (o$enetRule == "1se") .("largest within 1 SE of the minimum deviance") else .("minimum deviance"),
                                     o$enetFolds))
            if (identical(private$.famName, "cox"))
                t$setNote("ties", .("The penalised Cox model uses Breslow's handling of tied times (the unpenalised one, Efron's)."))
            if (!is.null(iv)) t$setNote("res", if (o$enetTune)
                .("λ is tuned again by cross-validation in every resampled fit, so the internal validation includes the tuning (slow: untick it for a quick look); the last column is the share of fits that kept the predictor.")
                else .("Resampled fits use the λ of the final model (the tuning itself is not validated); the last column is the share of fits that kept the predictor."))
            img <- self$results$plots$enetCv
            if (!is.null(cv) && is.null(img$state))
                img$setState(list(lambda = cv$lambda, cvm = cv$cvm, cvsd = cv$cvsd, imin = cv$imin, i1se = cv$i1se,
                                  nz = colSums(cv$path$beta != 0), chosen = lam,
                                  ylab = if (identical(private$.famName, "gaussian")) "Mean squared error" else
                                      if (identical(private$.famName, "cox")) "Partial likelihood deviance" else "Binomial deviance"))
            return()
        }
        crit <- switch(o$selCrit, aic = "AIC", bic = "BIC", lrt = sprintf(.("LR test, α = %s"), .num(o$selAlpha)))
        t$setNote("how", sprintf(if (o$selMode == "backward")
            .("Backward elimination (%s): at each step the predictor whose removal improves the criterion most leaves. Baseline covariates are never removed.")
            else .("Forward selection (%s): from the baseline covariates, at each step the predictor that improves the criterion most enters."), crit))
        t$setNote("test", .("χ², df, p: each predictor dropped from (if kept) or added to (if not) the final main-effects model; a factor or a spline is one term."))
        if (!is.null(iv))
            t$setNote("res", .("The selection is repeated in every resampled fit, so the internal validation includes its optimism; the last column is the share of fits that kept the predictor."))
    },

    # inference after selection: cautions, and the bootstrap coefficient tables
    # (Heinze, Wallisch & Dunkler 2018) from the internal-validation bootstrap
    .postSelection = function(penalized, iv, beta, spec) {
        o <- self$options; r <- self$results
        if (private$.selecting() && !penalized) {
            txt <- .("The model was selected in these data: its confidence intervals and p values do not account for the selection, so they are too narrow (see the bootstrap coefficients).")
            if (o$coefs) r$coefTable$setNote("post", txt)
            if (o$oddsRatios) r$orTable$setNote("post", txt)
            if (o$predTests) r$testTable$setNote("post", txt)
            if (length(spec$interactions)) r$subTable$setNote("post", txt)
        }
        if (!o$bootCoefs) return()
        t <- r$bootCoefTable; mt <- r$modelFreqTable
        t$setVisible(TRUE)
        if (is.null(iv) || iv$method != "boot" || is.null(iv$coefs)) {
            t$setNote("need", .("Choose Internal validation: Bootstrap; the bootstrap fits repeat the whole procedure (selection, tuning)."))
            return()
        }
        full <- stats::setNames(numeric(ncol(iv$coefs)), colnames(iv$coefs))
        full[intersect(names(beta), names(full))] <- beta[intersect(names(beta), names(full))]
        se <- NULL
        V <- private$.vcov(beta)
        if (!is.null(V)) {
            se <- stats::setNames(rep(NA_real_, length(full)), names(full))
            se[names(beta)] <- sqrt(pmax(0, diag(V)))
        }
        bc <- cp_boot_coefs(iv$coefs, full, se)
        if (identical(private$.famName, "cox")) bc <- bc[bc$term != "(Intercept)", ]
        for (i in seq_len(nrow(bc)))
            private$.putRow(t, bc$term[i], list(term = if (bc$term[i] == "(Intercept)") .("Intercept") else bc$term[i],
                est = bc$est[i], median = bc$median[i], lower = bc$lower[i], upper = bc$upper[i], incl = bc$incl[i],
                rmsd = bc$rmsd[i], relbias = bc$relbias[i]))
        t$setNote("ex", sprintf(.("%d bootstrap fits, each repeating the whole procedure. Median and 2.5th-97.5th percentiles count a coefficient left out of a fit as 0; inclusion: share of fits with the coefficient not 0. RMSD ratio: root mean squared difference of the bootstrap estimates from the final one, over its model-based SE (above 1: the selection adds variability the SE ignores). Relative conditional bias: mean estimate where selected over the final estimate, minus 1 (Heinze, Wallisch & Dunkler 2018)."),
                                nrow(iv$coefs)))
        t$setNote("cave", .("Percentile intervals describe the variability of the whole procedure; they are not exact confidence intervals, least of all for coefficients often shrunk or selected to 0."))
        mf <- cp_model_freq(iv$sets)
        if (!is.null(mf)) {
            mt$setVisible(TRUE)
            for (i in seq_len(nrow(mf)))
                private$.putRow(mt, paste0("m", i), list(model = mf$model[i], n = mf$n[i], pct = mf$pct[i], cum = mf$cum[i]))
            mt$setNote("ex", sprintf(.("The %d most frequent sets of selected terms among the bootstrap fits (%d distinct sets)."),
                                     nrow(mf), length(unique(stats::na.omit(iv$sets)))))
        }
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
    .internal = function(X, y, specC, lambda = NULL) {
        o <- self$options
        if (o$validation == "none") return(NULL)
        iv <- self$results$perfDev$state
        if (is.null(iv)) {
            fam <- private$.fam
            if (private$.selecting()) {
                base <- fam; so <- private$.selOpts()
                # elastic net: tune lambda again in each fit, or keep the final model's
                if (so$mode == "enet" && !o$enetTune) so$lambda <- lambda
                fam$fit <- function(X, y) cp_fit_procedure(X, y, specC, base, so)
            }
            iv <- cp_internal(X, y, fam, o$validation, B = o$nBoot, k = o$nFolds,
                              R = o$nRepeats, seed = o$seed, stability = o$validation == "boot")
            self$results$perfDev$setState(iv)
        }
        iv
    },

    # app: apparent measures already averaged over imputations (no CIs then)
    .fillPerfDev = function(y, lp, iv, conf, app = NULL) {
        t <- self$results$perfDev
        ci <- if (is.null(app)) private$.fam$ci(y, lp, conf) else list()
        app <- app %||% private$.fam$metrics(y, lp)
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
        # multiple imputation: averaged over the imputations, no CIs
        if (!is.null(vs$lps)) {
            m <- rowMeans(vapply(vs$lps, function(l) private$.fam$metrics(vs$y, l), private$.fam$metrics(vs$y, vs$lps[[1]])))
            ci <- list()
        } else {
            m <- private$.fam$metrics(vs$y, vs$lp)
            ci <- private$.fam$ci(vs$y, vs$lp, conf)
        }
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
            # revision has no single intercept and slope: blank, not NA (jamovi shows NaN)
            vals <- c(list(a = if (is.na(mm$a)) "" else mm$a, b = if (is.na(mm$b)) "" else mm$b, k = mm$k,
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
                # validation data hold the final model's predictors only (selection)
                if (is.null(x)) next
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
        t$setNote("member", private$.text("memberNote"))
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

    # a predictor the selection removed has no importance in the model: exactly 0
    .zeroImportance = function(dropped) {
        o <- self$options; r <- self$results
        keys <- private$.fam$permKeys
        note <- .("Predictors removed by the selection are not in the model: their importance is 0.")
        for (v in dropped) {
            if (o$shapley && v %in% r$shapleyTable$rowKeys) r$shapleyTable$setRow(rowKey = v, values = list(shapley = 0, share = 0, rank = ""))
            if (o$permImp && v %in% r$permTable$rowKeys)
                r$permTable$setRow(rowKey = v, values = as.list(stats::setNames(rep(0, 2 * length(keys)), c(paste0("d_", keys), paste0("sd_", keys)))))
            if (o$shap && v %in% r$shapTable$rowKeys) r$shapTable$setRow(rowKey = v, values = list(meanAbs = 0, share = 0))
        }
        for (t in c("shapleyTable", "permTable", "shapTable")) if (o[[c(shapleyTable = "shapley", permTable = "permImp", shapTable = "shap")[[t]]]])
            r[[t]]$setNote("drop", note)
    },

    # ---- incremental value of the other predictors over the baseline model ----
    .incKeys = function() switch(private$.famName,
        binary = c("auc", "dslope", "slope", "brier", "sbrier", "r2n"),
        gaussian = c("r2", "rmse", "mae", "slope"),
        cox = c("c", "uno", "auc", "d", "brier", "sbrier", "r2d", "r2n")),

    .incThresholds = function() {
        if (identical(private$.famName, "gaussian")) return(numeric(0))
        th <- tryCatch(.parse_numbers(self$options$nbThresholds), error = function(e) numeric(0))
        th[th > 0 & th < 1]
    },

    .initInc = function(sets) {
        t <- self$results$incTable
        t$setVisible(TRUE)
        defs <- cp_metric_defs(private$.famName)
        lab <- stats::setNames(defs$label, defs$key)
        for (k in names(sets)) {
            for (m in private$.incKeys()) t$addRow(rowKey = paste(k, m), values = list(set = sets[[k]], measure = lab[[m]]))
            if (identical(private$.famName, "binary")) {
                t$addRow(rowKey = paste(k, "nri"), values = list(set = sets[[k]], measure = .("NRI (continuous)")))
                t$addRow(rowKey = paste(k, "idi"), values = list(set = sets[[k]], measure = .("IDI")))
            }
            for (x in private$.incThresholds())
                t$addRow(rowKey = paste(k, "nb", x), values = list(set = sets[[k]], measure = sprintf(.("Net benefit at %s"), .num(x))))
            if (k == "dev") t$addRow(rowKey = "dev lr", values = list(set = sets[[k]], measure = .("Test of the added predictors")))
        }
        if (self$options$validation != "none")
            for (cn in c("baseCorr", "fullCorr", "deltaCorr")) t$getColumn(cn)$setVisible(TRUE)
        for (cn in c("lower", "upper")) t$getColumn(cn)$setSuperTitle(paste0("Δ ", self$options$ciWidth, "% CI"))
        # the decision curves of both models only when decision curves are asked for
        if (!is.null(self$results$plots[["incDca"]])) self$results$plots[["incDca"]]$setVisible(isTRUE(self$options$dca))
    },

    .fillInc = function(prep, specC, spec, y, lp, iv, vsets, conf, beta) {
        o <- self$options
        t <- self$results$incTable
        fam <- private$.fam
        forced <- private$.selOpts()$forced
        if (!length(forced)) {
            t$setNote("need", .("Tick the baseline covariates under Variable Selection: they form the baseline model."))
            return()
        }
        if (!length(setdiff(names(spec$vars), forced))) {
            t$setNote("need", .("The full model has no predictor beyond the baseline covariates (the selection removed all of them)."))
            return()
        }
        specB <- cp_with_terms(specC, forced, character(0))
        XB <- cp_design(specB, prep$data)
        bB <- fam$fit(XB, y)
        lpB <- fam$lp(bB, XB)
        keys <- private$.incKeys()
        th <- private$.incThresholds()
        nbf <- function(yy, p) if (identical(fam$name, "cox")) cp_net_benefit_surv(yy, p, th, fam$horizon) else cp_net_benefit(yy, p, th)
        fillSet <- function(k, yy, lb, lf) {
            mb <- fam$metrics(yy, lb); mf <- fam$metrics(yy, lf)
            for (m in keys) {
                dl <- mf[[m]] - mb[[m]]
                # rounding noise (e.g. apparent calibration slopes of 1) printed as 0
                v <- list(base = mb[[m]], full = mf[[m]], delta = if (isTRUE(abs(dl) < 1e-8)) 0 else dl)
                ci <- if (m == "auc" && fam$name == "binary") cp_delta_auc(yy, fam$risk(lb), fam$risk(lf), conf)
                      else if (m == "c" && fam$name == "cox") cp_delta_c(yy, lb, lf, conf) else NULL
                if (!is.null(ci)) v <- c(v, list(lower = ci[["lower"]], upper = ci[["upper"]], p = ci[["p"]]))
                private$.putRow(t, paste(k, m), v)
            }
            if (fam$name == "binary") {
                ni <- cp_nri_idi(yy, fam$risk(lb), fam$risk(lf), conf)
                for (z in c("nri", "idi"))
                    private$.putRow(t, paste(k, z), list(delta = ni[[z]][["est"]], lower = ni[[z]][["lower"]],
                                                         upper = ni[[z]][["upper"]], p = ni[[z]][["p"]]))
            }
            if (length(th)) {
                b <- nbf(yy, fam$risk(lb)); f <- nbf(yy, fam$risk(lf))
                for (i in seq_along(th))
                    private$.putRow(t, paste(k, "nb", th[i]), list(base = b$model[i], full = f$model[i], delta = f$model[i] - b$model[i]))
            }
        }
        fillSet("dev", y, lpB, lp)
        for (k in names(vsets)) fillSet(k, vsets[[k]]$y, fam$lp(bB, cp_design(specB, vsets[[k]]$df)), vsets[[k]]$lp)

        # nested test of the added predictors (not for a penalised model)
        if (isTRUE(attr(beta, "penalized"))) {
            private$.putRow(t, "dev lr", list(measure = .("Test of the added predictors: not for a penalised model")))
        } else {
            kF <- sum(beta != 0 & !is.na(beta)); kB <- sum(bB != 0 & !is.na(bB))
            n <- NROW(y)
            tst <- .nested_test(fam, fam$loss(y, lpB), fam$loss(y, lp), kF - kB, n, n - kF)
            private$.putRow(t, "dev lr", list(
                measure = if (fam$name == "gaussian") sprintf(.("F test of the added predictors (df = %d, %d)"), kF - kB, n - kF)
                          else sprintf(.("LR test of the added predictors (χ², df = %d)"), kF - kB),
                delta = tst[["stat"]], p = tst[["p"]]))
            if (o$selMode %in% c("backward", "forward"))
                t$setNote("sel", .("The test ignores that the added predictors were selected: its p is too small."))
        }

        # optimism-corrected, the baseline model validated on the same resamples
        if (!is.null(iv)) {
            ivB <- t$state
            if (is.null(ivB)) {
                ivB <- cp_internal(XB, y, fam, o$validation, B = o$nBoot, k = o$nFolds, R = o$nRepeats,
                                   seed = o$seed, stability = FALSE)
                t$setState(ivB)
            }
            for (m in keys)
                private$.putRow(t, paste("dev", m), list(baseCorr = ivB$corrected[[m]], fullCorr = iv$corrected[[m]],
                                                         deltaCorr = iv$corrected[[m]] - ivB$corrected[[m]]))
            t$setNote("corr", .("Optimism-corrected: each model validated on the same resamples; the full model's validation repeats its selection."))
        }
        t$setNote("base", sprintf(.("Baseline model: %s. Full model: the prediction model. Δ = full − baseline (for Brier, RMSE and MAE, negative is better)."),
                                  paste(forced, collapse = ", ")))
        if (fam$name == "binary")
            t$setNote("nri", .("NRI and IDI are reported for completeness; they can favour miscalibrated models (Pepe et al. 2015) and say less than ΔAUC, the LR test and Δ net benefit."))

        img <- self$results$plots[["incDca"]]
        if (!is.null(img) && isTRUE(o$dca) && is.null(img$state) && fam$name != "gaussian") {
            grid <- seq(0.01, o$dcaMax, by = 0.01)
            nbg <- function(yy, p) if (fam$name == "cox") cp_net_benefit_surv(yy, p, grid, fam$horizon) else cp_net_benefit(yy, p, grid)
            mk <- function(label, yy, lb, lf) { f <- nbg(yy, fam$risk(lf)); f$base <- nbg(yy, fam$risk(lb))$model; list(label = label, nb = f) }
            sets <- private$.sets()
            st <- list(mk(sets[["dev"]], y, lpB, lp))
            for (k in names(vsets)) st[[length(st) + 1]] <- mk(sets[[k]], vsets[[k]]$y, fam$lp(bB, cp_design(specB, vsets[[k]]$df)), vsets[[k]]$lp)
            img$setState(st)
        }
    },

    .plotIncDca = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_dca(image$state, ggtheme)); TRUE
    },

    # ---- missing values ----------------------------------------------------
    .mh = NULL,
    .nImputed = 0,

    # rows (with the outcome) in which at least one predictor is missing
    .countImputed = function(dd) {
        out <- private$.outcomeVars(); pr <- intersect(private$.preds(), names(dd))
        ok <- stats::complete.cases(dd[, out, drop = FALSE])
        sum(ok & !stats::complete.cases(dd[, pr, drop = FALSE]))
    },

    # the spec without the missing indicators (they are made from the data, not read)
    .specRaw = function(spec) {
        if (is.null(private$.mh) || !length(private$.mh$indicator)) return(spec)
        s <- spec; s$vars <- spec$vars[setdiff(names(spec$vars), private$.mh$indicator)]
        s
    },

    # validation data under the missing-data method (indicators with the
    # development medians, or their own single imputation)
    .missVal = function(df) switch(self$options$missMethod,
        category = cp_miss_category_apply(private$.mh, df),
        single = private$.impSet(df, 1, self$options$seed + 1)[[1]],
        df),

    # m imputed copies of the rows with the outcome; the outcome (for Cox, the
    # event and its Nelson-Aalen cumulative hazard; White & Royston 2009) is in
    # every imputation model
    .impSet = function(dd, m, seed) {
        o <- self$options
        outVars <- private$.outcomeVars()
        preds <- intersect(private$.preds(), names(dd))
        d <- dd[, c(outVars, preds), drop = FALSE]
        for (v in intersect(c(o$covs, o$splineVars), preds)) {
            x <- d[[v]]; d[[v]] <- if (is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else as.numeric(x)
        }
        for (v in intersect(o$factors, preds)) d[[v]] <- factor(d[[v]], ordered = FALSE)
        d <- d[stats::complete.cases(d[, outVars, drop = FALSE]), , drop = FALSE]
        if (identical(private$.famName, "cox")) {
            tm <- as.numeric(d[[o$time]])
            st <- as.integer(as.character(d[[o$status]]) == private$.level(o$status, o$eventLevel))
            ut <- sort(unique(tm[st == 1]))
            na <- cumsum(tabulate(match(tm[st == 1], ut), length(ut)) / (length(tm) - findInterval(ut, sort(tm), left.open = TRUE)))
            d[[".H"]] <- c(0, na)[findInterval(tm, ut) + 1]
            d[[o$time]] <- NULL
            imps <- cp_impute(d, preds, m = m, maxit = o$miIter, seed = seed)
            return(lapply(imps, function(x) { x[[o$time]] <- tm; x[[".H"]] <- NULL; x }))
        }
        if (identical(private$.famName, "gaussian")) d[[outVars]] <- jmvcore::toNumeric(d[[outVars]])
        cp_impute(d, preds, m = m, maxit = o$miIter, seed = seed)
    },

    .missNote = function() {
        o <- self$options
        t <- self$results$samples
        if (o$missMethod == "category" && !is.null(private$.mh) && length(c(private$.mh$num, private$.mh$fac)))
            t$setNote("miss", sprintf(.("Missing values: a '(missing)' level for %s; for %s, the development median and an indicator predictor."),
                                      if (length(private$.mh$fac)) paste(private$.mh$fac, collapse = ", ") else "—",
                                      if (length(private$.mh$num)) paste(private$.mh$num, collapse = ", ") else "—"))
        if (o$missMethod == "single")
            t$setNote("miss", sprintf(.("Missing predictor values imputed once by chained equations (predictive mean matching, %d iterations; the outcome in the imputation models): %d development rows. Validation data are imputed in the same way, separately. Rows with a missing outcome are excluded."),
                                      o$miIter, private$.nImputed))
    },

    .initMiss = function() {
        t <- self$results$missTable
        for (v in private$.outcomeVars()) t$addRow(rowKey = v, values = list(var = v, role = .("Outcome")))
        for (v in private$.preds()) t$addRow(rowKey = v, values = list(var = v, role = .("Predictor")))
        if (identical(private$.famName, "gaussian")) {
            for (cn in c("a_mean", "b_mean")) t$getColumn(cn)$setVisible(TRUE)
        } else {
            for (cn in c("a_pct", "b_pct")) t$getColumn(cn)$setVisible(TRUE)
            cox <- identical(private$.famName, "cox")
            t$getColumn("a_pct")$setTitle(if (cox) .("Censored") else .("Non-events"))
            t$getColumn("b_pct")$setTitle(if (cox) .("Events") else .("Events"))
        }
    },

    .fillMiss = function(dd) {
        o <- self$options
        outVars <- private$.outcomeVars()
        g <- if (identical(private$.famName, "gaussian")) jmvcore::toNumeric(dd[[o$dep]]) else {
            sv <- if (identical(private$.famName, "cox")) o$status else o$dep
            ev <- private$.level(sv, if (identical(private$.famName, "cox")) o$eventLevel else o$eventLevel)
            x <- as.character(dd[[sv]]); ifelse(is.na(x), NA, as.integer(x == ev))
        }
        ex <- cp_miss_explore(dd, c(outVars, private$.preds()), g, identical(private$.famName, "gaussian"))
        t <- self$results$missTable
        for (i in seq_len(nrow(ex))) {
            isOut <- ex$var[i] %in% outVars
            private$.putRow(t, ex$var[i], list(n = ex$n[i], pct = ex$pct[i],
                a_pct = if (isOut) "" else ex$a[i], b_pct = if (isOut) "" else ex$b[i],
                a_mean = if (isOut) "" else ex$a[i], b_mean = if (isOut) "" else ex$b[i], p = if (isOut) "" else ex$p[i]))
        }
        t$setNote("ex", if (identical(private$.famName, "gaussian"))
            .("Outcome mean in the rows where the predictor is observed and where it is missing (Welch's t test): a difference means that missingness depends on the outcome (not missing completely at random).")
            else .("Share of rows with the predictor missing among non-events and events (chi-square test): a difference means that missingness depends on the outcome (not missing completely at random)."))
        pt <- cp_miss_patterns(dd, c(outVars, private$.preds()))
        tp <- self$results$missPatterns
        for (i in seq_len(nrow(pt)))
            private$.putRow(tp, paste0("p", i), list(pattern = if (pt$pattern[i] == "(complete)") .("none (complete)") else pt$pattern[i],
                                                     nvars = pt$nvars[i], n = pt$n[i], pct = pt$pct[i]))
        tp$setNote("ex", .("The most frequent combinations of missing variables in the development rows (at most 10)."))
    },

    # the model with every candidate predictor (unpenalised), for the effects table
    .fullModel = function(XC, y, specC) {
        fam <- private$.fam
        b <- tryCatch(fam$fit(XC, y), error = function(e) NULL)
        if (is.null(b) || anyNA(b) || !isTRUE(attr(b, "converged"))) return(NULL)
        list(spec = specC, beta = b, V = fam$vcov(b), lrt = cp_lr_tests(specC, XC, y, fam))
    },

    .plotForest = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        o <- self$options
        blocks <- c(if (o$forestCrude) "Crude", if (o$forestAdj) c("Adjusted", "Adjusted (all predictors)"),
                    if (o$forestSel) "Selected model")
        p <- cp_plot_forest(image$state, blocks, ggtheme)
        if (is.null(p)) return(FALSE)
        print(p); TRUE
    },

    # the covariance of the coefficients; none for a penalised model
    .vcov = function(beta) if (isTRUE(attr(beta, "penalized"))) NULL else
        if (!is.null(attr(beta, "pool"))) attr(beta, "pool")$vcov else private$.fam$vcov(beta),

    # Xs: design matrices of the imputed data sets (the decomposition is averaged)
    .fillShapley = function(X, y, spec, penalized = FALSE, Xs = NULL) {
        t <- self$results$shapleyTable
        if (penalized) {
            t$setNote("tot", .("The Shapley decomposition refits unpenalised models for every subset of predictors, so it does not describe a penalised (elastic net) model; use SHAP or permutation importance."))
            return()
        }
        sh <- t$state
        # the linear analysis has a single R² (no shapleyMeasure option)
        measure <- if (self$options$has("shapleyMeasure")) self$options$shapleyMeasure else "r2"
        if (is.null(sh)) {
            sh <- if (is.null(Xs)) cp_shapley_r2(X, y, spec, private$.fam, measure, seed = self$options$seed) else {
                shs <- lapply(Xs, function(Xm) cp_shapley_r2(Xm, y, spec, private$.fam, measure, seed = self$options$seed))
                out <- shs[[1]]
                out$shapley <- rowMeans(vapply(shs, `[[`, numeric(nrow(out)), "shapley"))
                tot <- mean(vapply(shs, function(z) attr(z, "total"), 0))
                out$share <- out$shapley / tot; out$rank <- rank(-out$shapley, ties.method = "min")
                attr(out, "total") <- tot
                out
            }
            t$setState(sh)
        }
        for (i in seq_len(nrow(sh)))
            private$.putRow(t, sh$var[i], list(var = sh$var[i], shapley = sh$shapley[i], share = sh$share[i], rank = sh$rank[i]))
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
            private$.putRow(t, pi$var[i], as.list(pi[i, ]))
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
        for (v in names(ma)) private$.putRow(t, v, list(var = v, meanAbs = ma[[v]], share = ma[[v]] / sum(ma)))
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
                              signs = c(private$.text("shapDown"), private$.text("shapUp")),
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
    # imputed: list(rows, lp) of imputed development rows (mean over imputations)
    .saveOutputs = function(beta, spec, data, imputed = NULL) {
        defs <- private$.outputDefs()
        on <- Filter(function(k) isTRUE(self$options[[k]]), names(defs))
        if (!length(on)) return()
        nd <- cp_newdata(data, spec)
        lp <- rep(NA_real_, nrow(data)); names(lp) <- rownames(data)
        lp[rownames(nd)] <- private$.fam$lp(beta, cp_design(spec, nd))
        if (!is.null(imputed)) {
            fill <- imputed$rows %in% names(lp)[is.na(lp)]
            lp[imputed$rows[fill]] <- imputed$lp[fill]
        }
        for (k in on) {
            out <- self$results[[k]]
            if (!out$isNotFilled()) next
            out$setRowNums(rownames(data))
            out$setValues(unname(defs[[k]](lp, beta)))
        }
    },

    # ---- plots -------------------------------------------------------------
    .plotEnet = function(image, ggtheme, theme, ...) {
        if (is.null(image$state)) return(FALSE)
        print(cp_plot_enet_cv(image$state, ggtheme)); TRUE
    },
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
