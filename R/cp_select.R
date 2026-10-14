# The model-fitting procedure: selection of the predictors that are not
# forced (baseline covariates), then of the interactions among the kept ones.
# A predictor is one term (all the columns of a factor or a spline). The same
# procedure fits the final model and every bootstrap or cross-validation fit,
# so the optimism includes the selection.

# The spec with some predictors and interaction terms only
cp_with_terms <- function(spec, vars, ints = character(0)) {
    spec$vars <- spec$vars[intersect(names(spec$vars), vars)]
    keep <- vapply(spec$interactions, function(ab) all(ab %in% vars), TRUE)
    spec$interactions <- spec$interactions[keep & names(spec$interactions) %in% ints]
    .spec_columns(spec)
}

# Columns of X (of the full spec) that belong to the main effects of vars
.main_cols <- function(spec, vars) spec$colTerm %in% vars

# Stepwise selection of main effects. Returns the kept predictors and, for
# each non-forced predictor, the likelihood ratio test of the term against the
# final main-effects model (dropping it if kept, adding it if not) and the
# step at which it entered or left.
cp_select_main <- function(X, y, spec, fam, mode = "none", crit = "aic", alpha = 0.05,
                           forced = character(0)) {
    vars <- names(spec$vars)
    forced <- intersect(forced, vars)
    free <- setdiff(vars, forced)
    n <- NROW(X)
    fitv <- function(v) .fit_cols(X, y, fam, .main_cols(spec, v))
    # change of the criterion when a term of chi-square chi on df is added
    # (AIC, BIC: lower is better; LR test: its p). Dropping a term changes
    # AIC and BIC by minus that; terms without columns (df < 1) never move.
    addScore <- function(chi, df) if (df < 1) Inf else switch(crit,
        aic = -chi + 2 * df, bic = -chi + log(n) * df, lrt = stats::pchisq(chi, df, lower.tail = FALSE))
    step <- stats::setNames(rep(NA_integer_, length(free)), free)
    kept <- if (mode == "forward") forced else vars
    if (mode != "none" && length(free)) {
        cur <- fitv(kept); k <- 0
        repeat {
            cand <- if (mode == "backward") intersect(free, kept) else setdiff(free, kept)
            if (!length(cand)) break
            st <- lapply(cand, function(v) {
                f <- fitv(if (mode == "backward") setdiff(kept, v) else c(kept, v))
                big <- if (mode == "backward") cur else f; small <- if (mode == "backward") f else cur
                list(fit = f, s = addScore(max(0, 2 * (small$loss - big$loss)), big$df - small$df))
            })
            s <- vapply(st, `[[`, 0, "s")
            if (mode == "backward") {
                # drop the term whose removal improves the criterion most (largest p)
                sb <- replace(s, !is.finite(s), -Inf)
                best <- which.max(sb)
                ok <- if (crit == "lrt") sb[best] > alpha else sb[best] > 0
            } else {
                best <- which.min(s)
                ok <- if (crit == "lrt") s[best] < alpha else s[best] < 0
            }
            if (!length(best) || !isTRUE(ok)) break
            k <- k + 1
            step[cand[best]] <- k
            kept <- if (mode == "backward") setdiff(kept, cand[best]) else c(kept, cand[best])
            cur <- st[[best]]$fit
        }
    }
    kept <- vars[vars %in% kept]
    # each free predictor's test against the final main-effects model
    final <- fitv(kept)
    tests <- do.call(rbind, lapply(free, function(v) {
        other <- if (v %in% kept) fitv(setdiff(kept, v)) else fitv(c(kept, v))
        big <- if (v %in% kept) final else other; small <- if (v %in% kept) other else final
        chi <- max(0, 2 * (small$loss - big$loss)); df <- big$df - small$df
        data.frame(var = v, kept = v %in% kept, step = step[[v]], chisq = chi, df = df,
                   p = stats::pchisq(chi, df, lower.tail = FALSE), stringsAsFactors = FALSE)
    }))
    list(kept = kept, tests = tests)
}

# The whole procedure. sel: list(mode, crit, alpha, forced, intMode, intAlpha).
# Returns full-length coefficients (0 for the terms left out) with attributes
# vars (kept predictors), ints (kept interactions), selected (both, for the
# resampling frequencies), mainTests and screen (interaction screening).
cp_fit_procedure <- function(X, y, spec, fam, sel) {
    mode <- sel$mode %||% "none"
    if (mode == "enet") return(.fit_enet(X, y, spec, fam, sel))
    ms <- cp_select_main(X, y, spec, fam, mode, sel$crit %||% "aic", sel$alpha %||% 0.05, sel$forced %||% character(0))
    specK <- cp_with_terms(spec, ms$kept, names(spec$interactions))
    XK <- X[, specK$columns, drop = FALSE]
    b <- cp_select_interactions(XK, y, specK, fam, sel$intMode %||% "all", sel$intAlpha %||% 0.05)
    ints <- attr(b, "selected")
    full <- stats::setNames(numeric(ncol(X)), colnames(X))
    bb <- b; bb[is.na(bb)] <- 0
    full[specK$columns] <- if (fam$intercept) bb[-1] else bb
    out <- if (fam$intercept) c(`(Intercept)` = unname(bb[1]), full) else full
    attr(out, "converged") <- attr(b, "converged")
    attr(out, "separation") <- attr(b, "separation")
    attr(out, "vars") <- ms$kept
    attr(out, "ints") <- ints
    attr(out, "selected") <- c(if (mode != "none") ms$kept, ints)
    attr(out, "mainTests") <- ms$tests
    attr(out, "screen") <- attr(b, "screen")
    out
}

# Elastic net over every candidate column (main effects and all candidate
# interactions); the columns of baseline covariates are not penalised. The
# penalised coefficients are the model (no refit). A predictor is kept when
# any of its columns, or of an interaction it is part of, is not zero.
.fit_enet <- function(X, y, spec, fam, sel) {
    forced <- sel$forced %||% character(0)
    pf <- ifelse(spec$colTerm %in% forced, 0, 1)
    # nothing to penalise: baseline covariates are never penalised
    if (all(pf == 0))
        .fail(.("The elastic net has nothing to select: every predictor is a baseline covariate. Choose no selection, or leave some predictors out of the baseline."))
    b <- cp_enet_fit(X, y, fam, alpha = sel$enetAlpha %||% 0.5, pf = pf, rule = sel$enetRule %||% "min",
                     nfolds = sel$enetFolds %||% 10, seed = sel$seed %||% 1, lambda = sel$lambda)
    nz <- spec$colTerm[b[-1] != 0]
    # relaxed lasso (Meinshausen 2007, with the unpenalised refit): the
    # selected terms refitted without penalty; this refit is the model
    if (isTRUE(sel$relax)) {
        keep <- b[-1] != 0 | spec$colTerm %in% forced
        rf <- fam$fit(X[, keep, drop = FALSE], y)
        full <- stats::setNames(numeric(ncol(X)), colnames(X))
        rr <- rf; rr[is.na(rr)] <- 0
        full[keep] <- rr[-1]
        out <- c(`(Intercept)` = unname(rr[1]), full)
        attr(out, "converged") <- attr(rf, "converged")
        attr(out, "separation") <- attr(rf, "separation")
        attr(out, "lambda") <- attr(b, "lambda"); attr(out, "cv") <- attr(b, "cv")
        attr(out, "relaxed") <- TRUE
        b <- out
    }
    ints <- intersect(names(spec$interactions), nz)
    vars <- names(spec$vars)[names(spec$vars) %in% c(forced, nz, unlist(spec$interactions[ints]))]
    free <- setdiff(names(spec$vars), forced)
    tests <- data.frame(var = free, kept = free %in% vars, step = NA_integer_, chisq = NA_real_, df = NA_real_,
                        p = NA_real_, stringsAsFactors = FALSE)
    attr(b, "vars") <- vars
    attr(b, "ints") <- ints
    attr(b, "selected") <- c(vars, ints)
    attr(b, "mainTests") <- tests
    b
}
