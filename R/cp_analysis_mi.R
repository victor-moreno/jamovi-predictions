# Multiple imputation: the development (and validation) data are imputed m
# times by chained equations; the model is chosen and fitted in each imputed
# data set and pooled:
#   - predictors and interactions kept by the majority of the imputations
#     (Wood, White & Royston 2008); baseline covariates always;
#   - coefficients and their covariance by Rubin's rules (elastic net: the
#     average of the penalised coefficients); tests by pooled Wald (D1);
#   - performance of the pooled model averaged over the imputations, optimism
#     by bootstrap or cross-validation within each imputation (Wahl et al. 2016);
#   - plots, classification, net benefit, case-mix, permutation importance and
#     SHAP on the imputations stacked; the Shapley decomposition averaged.
# Appended to the shared private methods (this file collates after cp_analysis.R).

.cp_private <- c(.cp_private, list(

    .analyseMI = function(data, isDev, isVal, sets, conf) {
        o <- self$options
        r <- self$results
        outVars <- private$.outcomeVars()
        M <- o$miM
        cands <- private$.candidates()
        devRaw <- data[isDev, , drop = FALSE]
        private$.nImputed <- private$.countImputed(devRaw)
        imps <- private$.impSet(devRaw, M, o$seed)
        n <- nrow(imps[[1]])
        if (n < 3) .fail("Fewer than 3 development rows have the outcome.")

        # one spec for every imputation (knots and quartiles from them all)
        prepS <- cp_prepare(do.call(rbind, imps), outVars, o$covs, o$splineVars, o$factors, o$knots,
                            refLevels = private$.refs(), scales = private$.scales(), interactions = cands)
        prep <- cp_prepare(imps[[1]], outVars, o$covs, o$splineVars, o$factors, o$knots,
                           refLevels = private$.refs(), scales = private$.scales(), interactions = cands)
        prep$spec <- prepS$spec
        prep <- private$.setOutcome(prep)
        prep$nMissing <- sum(isDev) - n + sum(!isDev & !isVal)
        fam <- private$.fam <- private$.makeFamily(prep)
        y <- prep$y
        specC <- prep$spec
        Ds <- lapply(imps, function(im) cp_newdata(im, specC))
        # stacked imputations: SHAP values (computed on them) are coloured by them
        private$.devData <- do.call(rbind, Ds)
        XCs <- lapply(Ds, function(D) cp_design(specC, D))

        # the procedure in each imputation; terms kept by the majority
        so <- private$.selOpts()
        sels <- lapply(XCs, function(XC) cp_fit_procedure(XC, y, specC, fam, so))
        penalized <- isTRUE(attr(sels[[1]], "penalized"))
        kept <- table(unlist(lapply(sels, function(s) c(attr(s, "vars"), attr(s, "ints"))))) / M
        maj <- names(kept)[kept >= 0.5]
        vars <- names(specC$vars)[names(specC$vars) %in% c(so$forced, maj)]
        if (!length(vars))
            .fail("The selection removed every predictor. Mark baseline covariates under Variable Selection, ",
                  "or use a less strict criterion.")
        spec <- cp_with_terms(specC, vars, intersect(names(specC$interactions), maj))
        Xs <- lapply(Ds, function(D) cp_design(spec, D))
        if (penalized) {
            beta <- rowMeans(vapply(sels, function(s) s[c("(Intercept)", spec$columns)], numeric(length(spec$columns) + 1)))
            attr(beta, "converged") <- TRUE
            attr(beta, "penalized") <- TRUE
            attr(beta, "lambda") <- mean(vapply(sels, function(s) attr(s, "lambda"), 0))
            attr(beta, "cv") <- attr(sels[[1]], "cv")
            pool <- NULL
        } else {
            fits <- lapply(Xs, function(X) fam$fit(X, y))
            if (any(vapply(fits, function(b) anyNA(b), TRUE)))
                .fail("A predictor column is redundant (a combination of other predictors) in some imputed data; remove one of the predictors involved.")
            if (!all(vapply(fits, function(b) isTRUE(attr(b, "converged")), TRUE)))
                .fail("The model did not converge in some imputed data (complete separation); remove or merge the predictor(s) involved.")
            pool <- cp_rubin(lapply(fits, function(b) stats::setNames(as.numeric(b), names(b))), lapply(fits, fam$vcov))
            beta <- pool$est
            attr(beta, "converged") <- TRUE
            attr(beta, "pool") <- pool
        }
        lps <- lapply(Xs, function(X) fam$lp(beta, X))
        stackY <- function(yy, k) if (is.matrix(yy)) do.call(rbind, rep(list(yy), k)) else rep(yy, k)
        Xst <- do.call(rbind, Xs); yst <- stackY(y, M); lpst <- unlist(lps)

        # validation data: imputed separately, M times
        vsets <- list()
        specRaw <- spec
        mkSet <- function(df, key, nTotal) {
            vi <- private$.impSet(df, M, o$seed + 1)
            vd <- lapply(vi, function(x) cp_split_data(x, spec, outVars)$data)
            yv <- private$.yOf(vd[[1]], spec, key)
            Xv <- lapply(vd, function(x) cp_design(spec, x))
            lv <- lapply(Xv, function(X) fam$lp(beta, X))
            list(y = yv, lps = lv, Xs = Xv, X = do.call(rbind, Xv), lp = unlist(lv), p = fam$risk(unlist(lv)),
                 ys = stackY(yv, M), df = do.call(rbind, vd), nMissing = nTotal - nrow(vd[[1]]))
        }
        if (any(isVal)) {
            sd <- cp_split_data(data[isVal, , drop = FALSE], specRaw, outVars, complete = FALSE)
            vsets$split <- mkSet(sd$data, "split", sum(isVal))
        }
        if (!is.null(o$extFile)) {
            ext <- cp_read_external(o$extFile$path, o$extFile$filename %||% o$extFile$path)
            h <- cp_harmonise(ext, specRaw, outVars, sets[["ext"]], complete = FALSE)
            for (i in seq_len(nrow(h$check)))
                r$extCheck$setRow(rowKey = h$check$var[i], values = list(role = h$check$role[i], status = h$check$status[i],
                                                                         details = h$check$details[i]))
            vsets$ext <- mkSet(h$data, "ext", h$nTotal)
        }
        # stacked views of the validation sets for the functions that take one data set
        stacked <- lapply(vsets, function(v) list(X = v$X, y = v$ys, lp = v$lp, p = v$p, df = v$df))

        private$.fillSamples(y, prep$nMissing, vsets)
        r$samples$setNote("miss", sprintf(.("Multiple imputation: %d imputations by chained equations (predictive mean matching, %d iterations; the outcome in the imputation models) of %d development rows with a missing predictor; validation data imputed separately. Coefficients pooled by Rubin's rules, performance averaged over the imputations. Rows with a missing outcome are excluded."),
                                          M, o$miIter, private$.nImputed))
        private$.fillCoefs(beta, spec, conf)
        if (o$equation) r$equationText$setContent(private$.equation(beta, spec))
        lrt <- if ((o$oddsRatios || o$predTests) && !penalized) cp_lr_tests_mi(spec, pool) else NULL
        if (o$oddsRatios) {
            prepD <- prep; prepD$data <- Ds[[1]]
            fullMI <- NULL
            if (private$.selecting()) {
                fs <- tryCatch(lapply(XCs, function(X) fam$fit(X, y)), error = function(e) NULL)
                if (!is.null(fs) && !any(vapply(fs, anyNA, TRUE))) {
                    pf <- cp_rubin(lapply(fs, function(b) stats::setNames(as.numeric(b), names(b))), lapply(fs, fam$vcov))
                    bF <- pf$est; attr(bF, "converged") <- TRUE
                    fullMI <- list(spec = specC, beta = bF, V = pf$vcov, lrt = cp_lr_tests_mi(specC, pf))
                }
            }
            private$.fillOR(spec, beta, prepD, y, lrt, conf,
                            crude = if (o$crude || (o$forest && o$forestCrude)) private$.crudeMI(specC, Ds, y, conf) else NULL,
                            full = fullMI)
        }
        if (o$predTests) {
            private$.fillTests(lrt)
            if (!is.null(lrt)) r$testTable$setNote("mi", .("Pooled Wald tests (D1) over the imputations; χ² = D1 × df, p from F(df, df2)."))
        }
        if (o$oddsRatios && o$crude) r$orTable$setNote("mi", .("Crude and adjusted effects pooled over the imputations by Rubin's rules."))

        iv <- private$.internalMI(XCs, y, specC, so)
        if (length(cands)) {
            private$.fillInt(XCs[[1]], y, specC, sels[[1]], iv, conf, spec, beta)
            r$intTable$setNote("mi", .("Screening shown in the first imputed data set; the terms in the model are those kept in most imputations."))
        }
        if (o$selMode != "none") {
            private$.fillSel(sels[[1]], iv)
            for (v in names(specC$vars)) if (v %in% r$selTable$rowKeys && !(v %in% so$forced))
                private$.putRow(r$selTable, v, list(status = if (v %in% vars) .("kept") else .("removed")))
            r$selTable$setNote("mi", sprintf(.("Kept: selected in at least half of the %d imputations (tests from the first). Selected in each imputation: %s."),
                                             M, paste(sprintf("%s %.0f%%", names(kept), 100 * as.numeric(kept)), collapse = ", ")))
        }
        private$.postSelection(penalized, iv, beta, spec)
        app <- rowMeans(vapply(lps, function(l) fam$metrics(y, l), fam$metrics(y, lps[[1]])))
        private$.fillPerfDev(y, lps[[1]], iv, conf, app = app)
        r$perfDev$setNote("mi", .("Averaged over the imputations (no confidence intervals); internal validation within each imputation."))
        for (k in names(vsets)) private$.fillPerfVal(k, vsets[[k]], conf)
        if (length(vsets)) {
            prepS2 <- prep; prepS2$data <- do.call(rbind, Ds)
            private$.fillCaseMix(prepS2, Xst, yst, lpst, stacked)
        }
        if (private$.on("updating"))
            r$updTable$setNote("mi", .("Model updating is not available with multiple imputation (its tests are not defined across imputations); use single imputation or complete cases."))
        if (o$incremental) private$.fillIncMI(prep, specC, spec, y, Ds, lps, vsets, beta, pool)
        if (!is.null(iv) && iv$method == "boot") private$.fillStab(iv, y, lps[[1]], XCs[[1]])

        all <- c(list(dev = list(X = Xst, y = yst, lp = lpst, p = fam$risk(lpst), df = do.call(rbind, Ds))), stacked)
        private$.fillExtras(all, sets, conf)
        if (o$sampleSize) private$.fillSampleSize(y, lps[[1]], specC, iv)
        if (o$shapley) private$.fillShapley(Xs[[1]], y, spec, penalized, Xs = Xs)
        if (o$permImp) private$.fillPerm(beta, spec, all, sets)
        if (o$shap) private$.fillShap(beta, spec, Xst, data)
        dropped <- setdiff(private$.preds(), names(spec$vars))
        if (length(dropped)) private$.zeroImportance(dropped)
        private$.plotStates(all, sets, iv)
        # predictions: rows with complete predictors, and imputed development rows (mean over imputations)
        private$.saveOutputs(beta, spec, data, imputed = list(rows = rownames(imps[[1]]), lp = rowMeans(do.call(cbind, lps))))
        if (private$.on("exportModel"))
            tryCatch(private$.exportModel(beta, spec, y, lps[[1]], iv),
                     error = function(e) .show_warning(self, "export", .("The model could not be exported.")))
    },

    # internal validation within each imputation; optimism and selection
    # frequencies averaged (bootstrap samples shared out over the imputations)
    .internalMI = function(XCs, y, specC, so) {
        o <- self$options
        if (o$validation == "none") return(NULL)
        iv <- self$results$perfDev$state
        if (!is.null(iv)) return(iv)
        fam <- private$.fam
        if (private$.selecting()) {
            base <- fam
            fam$fit <- function(X, y) cp_fit_procedure(X, y, specC, base, so)
        }
        M <- length(XCs)
        Bm <- max(10, ceiling(o$nBoot / M))
        ivs <- lapply(seq_len(M), function(m) cp_internal(XCs[[m]], y, fam, o$validation, B = Bm, k = o$nFolds,
                                                           R = o$nRepeats, seed = o$seed + m, stability = m == 1 && o$validation == "boot"))
        avg <- function(f) rowMeans(vapply(ivs, f, ivs[[1]]$apparent))
        iv <- ivs[[1]]
        iv$apparent <- avg(function(z) z$apparent)
        iv$optimism <- avg(function(z) z$optimism)
        iv$corrected <- iv$apparent - iv$optimism
        iv$B <- Bm * M; iv$nOk <- sum(vapply(ivs, function(z) z$nOk, 0))
        if (!is.null(iv$shrinkage)) iv$shrinkage <- mean(vapply(ivs, function(z) z$shrinkage, 0), na.rm = TRUE)
        iv$coefs <- do.call(rbind, lapply(ivs, function(z) z$coefs))
        iv$sets <- unlist(lapply(ivs, function(z) z$sets))
        sels <- lapply(ivs, function(z) z$selection)
        if (!all(vapply(sels, is.null, TRUE))) {
            terms <- unique(unlist(lapply(sels, names)))
            iv$selection <- stats::setNames(vapply(terms, function(t) mean(vapply(sels, function(s) .sel_freq(s, t), 0)), 0), terms)
        }
        self$results$perfDev$setState(iv)
        iv
    },

    # crude effects pooled over imputations (Rubin's rules on the log scale for ratios)
    .crudeMI = function(specC, Ds, y, conf) {
        fam <- private$.fam
        z <- stats::qnorm(1 - (1 - conf) / 2)
        crs <- lapply(Ds, function(D) cp_crude(specC, D, y, fam, conf))
        tr <- if (fam$exp) log else identity; inv <- if (fam$exp) exp else identity
        est <- vapply(crs, function(cr) tr(cr$est), numeric(nrow(crs[[1]])))
        se <- vapply(crs, function(cr) (tr(cr$upper) - tr(cr$lower)) / (2 * z), numeric(nrow(crs[[1]])))
        est <- matrix(est, nrow = nrow(crs[[1]])); se <- matrix(se, nrow = nrow(crs[[1]]))
        M <- length(Ds)
        q <- rowMeans(est)
        tv <- rowMeans(se^2) + (1 + 1 / M) * apply(est, 1, function(x) if (M > 1) stats::var(x) else 0)
        out <- crs[[1]]
        out$est <- inv(q); out$lower <- inv(q - z * sqrt(tv)); out$upper <- inv(q + z * sqrt(tv))
        out$p <- 2 * stats::pnorm(-abs(q / sqrt(tv)))
        # overall p of each predictor: the median over imputations (a pooled test needs the fits)
        attr(out, "pOverall") <- lapply(names(attr(crs[[1]], "pOverall")), function(v)
            stats::median(vapply(crs, function(cr) attr(cr, "pOverall")[[v]], 0)))
        names(attr(out, "pOverall")) <- names(attr(crs[[1]], "pOverall"))
        out
    },

    # incremental value under multiple imputation: measures averaged over the
    # imputations, Wald (D1) test of the added predictors
    .fillIncMI = function(prep, specC, spec, y, Ds, lps, vsets, beta, pool) {
        o <- self$options
        t <- self$results$incTable
        fam <- private$.fam
        forced <- private$.selOpts()$forced
        if (!length(forced)) { t$setNote("need", .("Tick the baseline covariates under Variable Selection: they form the baseline model.")); return() }
        if (!length(setdiff(names(spec$vars), forced))) {
            t$setNote("need", .("The full model has no predictor beyond the baseline covariates (the selection removed all of them).")); return()
        }
        specB <- cp_with_terms(specC, forced, character(0))
        XBs <- lapply(Ds, function(D) cp_design(specB, D))
        fitsB <- lapply(XBs, function(X) fam$fit(X, y))
        bB <- cp_rubin(lapply(fitsB, function(b) stats::setNames(as.numeric(b), names(b))), lapply(fitsB, fam$vcov))$est
        lpBs <- lapply(XBs, function(X) fam$lp(bB, X))
        keys <- private$.incKeys(); th <- private$.incThresholds()
        nbf <- function(yy, p) if (identical(fam$name, "cox")) cp_net_benefit_surv(yy, p, th, fam$horizon) else cp_net_benefit(yy, p, th)
        fillSet <- function(k, yy, lbs, lfs) {
            mb <- rowMeans(vapply(lbs, function(l) fam$metrics(yy, l), fam$metrics(yy, lbs[[1]])))
            mf <- rowMeans(vapply(lfs, function(l) fam$metrics(yy, l), fam$metrics(yy, lfs[[1]])))
            for (m in keys) private$.putRow(t, paste(k, m), list(base = mb[[m]], full = mf[[m]], delta = mf[[m]] - mb[[m]]))
            if (fam$name == "binary") {
                nis <- lapply(seq_along(lbs), function(i) cp_nri_idi(yy, fam$risk(lbs[[i]]), fam$risk(lfs[[i]])))
                for (z in c("nri", "idi"))
                    private$.putRow(t, paste(k, z), list(delta = mean(vapply(nis, function(x) x[[z]][["est"]], 0))))
            }
            if (length(th)) for (i in seq_along(th)) {
                b <- mean(vapply(lbs, function(l) nbf(yy, fam$risk(l))$model[i], 0))
                f <- mean(vapply(lfs, function(l) nbf(yy, fam$risk(l))$model[i], 0))
                private$.putRow(t, paste(k, "nb", th[i]), list(base = b, full = f, delta = f - b))
            }
        }
        fillSet("dev", y, lpBs, lps)
        for (k in names(vsets)) {
            vs <- vsets[[k]]
            fillSet(k, vs$y, lapply(vs$Xs, function(X) fam$lp(bB, X[, specB$columns, drop = FALSE])), vs$lps)
        }
        if (is.null(pool)) {
            private$.putRow(t, "dev lr", list(measure = .("Test of the added predictors: not for a penalised model")))
        } else {
            added <- which(!(spec$colTerm %in% forced)) + 1
            d <- cp_d1(pool, added)
            private$.putRow(t, "dev lr", list(measure = sprintf(.("Pooled Wald test of the added predictors (D1, df = %d, %.1f)"), d[["df1"]], d[["df2"]]),
                                              delta = d[["stat"]], p = d[["p"]]))
        }
        t$setNote("mi", .("Averaged over the imputations, without confidence intervals; the optimism-corrected comparison is not computed with multiple imputation."))
        t$setNote("base", sprintf(.("Baseline model: %s. Full model: the prediction model. Δ = full − baseline (for Brier, RMSE and MAE, negative is better)."),
                                  paste(forced, collapse = ", ")))
    }
))
