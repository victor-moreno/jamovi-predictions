# Internal validation, the same for every model family: Harrell's bootstrap
# optimism correction or repeated k-fold cross-validation, plus the stability
# of individual predictions across bootstrap models (Riley & Collins 2023).
# y may be a vector or a matrix (e.g. time and status), hence .ysub().

.ysub <- function(y, i) if (is.matrix(y)) y[i, , drop = FALSE] else y[i]

# Fold labels, balanced within strata (the outcome for a binary model)
.folds <- function(strata, k) {
    f <- integer(length(strata))
    for (s in unique(strata)) {
        i <- which(strata == s)
        f[i] <- sample(rep_len(seq_len(k), length(i)))
    }
    f
}

.try_fit <- function(fam, X, y) {
    beta <- tryCatch(fam$fit(X, y), error = function(e) NULL)
    if (is.null(beta) || !isTRUE(attr(beta, "converged"))) return(NULL)
    beta
}

# method: "boot" or "cv". Returns apparent, optimism and corrected metrics;
# for the bootstrap also the mean calibration slope of the bootstrap models in
# the original data (a uniform shrinkage factor) and stability summaries.
cp_internal <- function(X, y, fam, method = "boot", B = 200, k = 10, R = 1, seed = 1,
                        stability = TRUE, maxCurves = 100, maxIndiv = 400) {
    beta <- fam$fit(X, y)
    lp0 <- fam$lp(beta, X)
    apparent <- fam$metrics(y, lp0)
    n <- NROW(X)

    .with_seed(seed, {
        if (method == "boot") {
            # the bootstrap samples are drawn first (in this order, as a plain
            # loop would), then fitted, possibly on several cores
            idxs <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
            res <- .cp_lapply(idxs, function(idx) {
                yb <- .ysub(y, idx)
                fb <- .try_fit(fam, X[idx, , drop = FALSE], yb)
                if (is.null(fb)) return(NULL)
                mApp <- tryCatch(fam$metrics(yb, fam$lp(fb, X[idx, , drop = FALSE])), error = function(e) NULL)
                lpTest <- fam$lp(fb, X)
                mTest <- tryCatch(fam$metrics(y, lpTest), error = function(e) NULL)
                if (is.null(mApp) || is.null(mTest)) return(NULL)
                list(opt = mApp - mTest, slope = mTest[["slope"]], selected = attr(fb, "selected"),
                     coef = as.numeric(fb),
                     set = if (!is.null(attr(fb, "vars"))) paste(c(sort(attr(fb, "vars")), sort(attr(fb, "ints"))), collapse = ", ") else NA_character_,
                     risk = if (stability) fam$risk(lpTest) else NULL)
            })
            ok <- !vapply(res, is.null, TRUE)
            opt <- matrix(NA_real_, B, length(apparent), dimnames = list(NULL, names(apparent)))
            coefs <- matrix(NA_real_, B, length(beta), dimnames = list(NULL, names(beta)))
            testSlope <- rep(NA_real_, B); sets <- rep(NA_character_, B); selected <- vector("list", B)
            risk <- if (stability) matrix(NA_real_, n, B) else NULL
            for (b in which(ok)) {
                z <- res[[b]]
                opt[b, ] <- z$opt; testSlope[b] <- z$slope; coefs[b, ] <- z$coef; sets[b] <- z$set
                if (!is.null(z$selected)) selected[[b]] <- z$selected
                if (stability) risk[, b] <- z$risk
            }
            # each measure averaged over the samples where it could be computed
            optimism <- colMeans(opt[ok, , drop = FALSE], na.rm = TRUE)
            out <- list(method = "boot", B = B, nOk = sum(ok),
                        shrinkage = mean(testSlope[ok], na.rm = TRUE),
                        selection = .selection_freq(selected[ok]),
                        coefs = coefs[ok, , drop = FALSE], sets = sets[ok])
            if (stability && sum(ok) > 1)
                out$stability <- .stability(y, fam$risk(lp0), risk[, ok, drop = FALSE],
                                            maxCurves, maxIndiv, fam$smooth %||% .cal_smooth,
                                            bounded = isTRUE(fam$bounded))
        } else {
            strata <- fam$strata(y)
            # all folds of all repeats drawn first, then every fold fitted
            folds <- lapply(seq_len(R), function(r) .folds(strata, k))
            tasks <- expand.grid(f = seq_len(k), r = seq_len(R))
            fits <- .cp_lapply(seq_len(nrow(tasks)), function(i) {
                test <- folds[[tasks$r[i]]] == tasks$f[i]
                fb <- .try_fit(fam, X[!test, , drop = FALSE], .ysub(y, !test))
                if (is.null(fb)) return(NULL)
                list(lp = fam$lp(fb, X[test, , drop = FALSE]), selected = attr(fb, "selected"))
            })
            cvm <- matrix(NA_real_, R, length(apparent), dimnames = list(NULL, names(apparent)))
            selected <- list()
            for (r in seq_len(R)) {
                oof <- rep(NA_real_, n)
                for (i in which(tasks$r == r)) {
                    z <- fits[[i]]
                    if (is.null(z)) next
                    oof[folds[[r]] == tasks$f[i]] <- z$lp
                    selected[[length(selected) + 1]] <- z$selected
                }
                if (!anyNA(oof))
                    cvm[r, ] <- tryCatch(fam$metrics(y, oof), error = function(e) NA_real_)
            }
            ok <- rowSums(!is.na(cvm)) > 0
            if (!any(ok)) .fail("Cross-validation failed in every repeat: the folds are too small for this model. Use fewer folds or the bootstrap.")
            optimism <- apparent - colMeans(cvm[ok, , drop = FALSE], na.rm = TRUE)
            out <- list(method = "cv", k = k, R = R, nOk = sum(ok),
                        selection = .selection_freq(selected))
        }
    })
    out$apparent <- apparent
    out$optimism <- optimism
    out$corrected <- apparent - optimism
    out
}

# How often each selected term (e.g. an interaction) was chosen across the
# resampled fits; NULL when the fitting procedure selects nothing
.selection_freq <- function(sel) {
    if (!length(sel) || all(vapply(sel, is.null, TRUE))) return(NULL)
    terms <- unique(unlist(sel))
    stats::setNames(vapply(terms, function(t) mean(vapply(sel, function(s) t %in% s, TRUE)), 0), terms)
}

# The share of resampled fits that selected term k (0 when none did)
.sel_freq <- function(freq, k) if (k %in% names(freq)) freq[[k]] else 0

# Stability of individual risks across bootstrap models applied to the original
# data: MAPE (mean absolute prediction error versus the original model's risk),
# a sample of points for the instability plot, and calibration curves.
# Predictions are stored as integers on a 1/1000 scale of [lo, hi]: [0, 1] for
# risks, the range of the predictions for an unbounded outcome.
.stability <- function(y, p0, risk, maxCurves, maxIndiv, smooth = .cal_smooth, bounded = TRUE) {
    n <- length(p0); B <- ncol(risk)
    ape <- rowMeans(abs(risk - p0))
    ind <- if (n > maxIndiv) sort(sample.int(n, maxIndiv)) else seq_len(n)
    bs <- if (B > 50) sort(sample.int(B, 50)) else seq_len(B)
    grid <- seq(stats::quantile(p0, 0.01), stats::quantile(p0, 0.99), length.out = 50)
    cb <- if (B > maxCurves) sort(sample.int(B, maxCurves)) else seq_len(B)
    curves <- vapply(cb, function(b) smooth(y, risk[, b], grid), numeric(length(grid)))
    q <- apply(risk, 1, stats::quantile, probs = c(0.025, 0.975))
    # stored compactly as integers (results state is size-limited)
    lim <- if (bounded) c(0, 1) else range(p0, risk[ind, bs])
    enc <- function(v) round((v - lim[1]) / (lim[2] - lim[1]) * 1000)
    list(mape = mean(ape), maxape = max(ape), lim = lim,
         p0 = enc(p0[ind]), pts = matrix(as.integer(enc(risk[ind, bs])), length(ind)),
         lo = enc(q[1, ind]), hi = enc(q[2, ind]),
         grid = grid, curves = if (bounded) round(curves, 4) else signif(curves, 5))
}

# Bootstrap summaries of a model built with selection (Heinze, Wallisch &
# Dunkler 2018): for each coefficient the share of fits in which it is not 0
# (inclusion frequency), the bootstrap median and 2.5th-97.5th percentiles
# (zeros of the fits that left it out included), the RMSD ratio (root mean
# squared difference of the bootstrap estimates from the final one, over its
# model-based SE) and the relative conditional bias (mean of the bootstrap
# estimates where selected, over the final estimate, minus 1).
# final: full-length coefficients (0 for terms left out); se: their SEs or NULL.
cp_boot_coefs <- function(coefs, final, se = NULL) {
    cols <- intersect(colnames(coefs), names(final))
    do.call(rbind, lapply(cols, function(cn) {
        b <- coefs[, cn]; b[is.na(b)] <- 0
        est <- final[[cn]]
        sel <- b != 0
        q <- stats::quantile(b, c(0.5, 0.025, 0.975), names = FALSE)
        data.frame(term = cn, est = est, median = q[1], lower = q[2], upper = q[3], incl = mean(sel),
                   rmsd = if (!is.null(se) && is.finite(se[[cn]]) && se[[cn]] > 0) sqrt(mean((b - est)^2)) / se[[cn]] else NA_real_,
                   relbias = if (any(sel) && est != 0) mean(b[sel]) / est - 1 else NA_real_,
                   stringsAsFactors = FALSE)
    }))
}

# The most frequent sets of selected terms among the bootstrap fits
cp_model_freq <- function(sets, top = 10) {
    sets <- sets[!is.na(sets)]
    if (!length(sets)) return(NULL)
    tab <- sort(table(sets), decreasing = TRUE)
    k <- utils::head(names(tab), top)
    data.frame(model = k, n = as.integer(tab[k]), pct = as.numeric(tab[k]) / length(sets),
               cum = cumsum(as.numeric(tab[k])) / length(sets), stringsAsFactors = FALSE)
}
