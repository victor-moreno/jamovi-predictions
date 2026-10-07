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
            opt <- matrix(NA_real_, B, length(apparent), dimnames = list(NULL, names(apparent)))
            testSlope <- rep(NA_real_, B)
            risk <- if (stability) matrix(NA_real_, n, B) else NULL
            for (b in seq_len(B)) {
                idx <- sample.int(n, n, replace = TRUE)
                yb <- .ysub(y, idx)
                fb <- .try_fit(fam, X[idx, , drop = FALSE], yb)
                if (is.null(fb)) next
                mApp <- tryCatch(fam$metrics(yb, fam$lp(fb, X[idx, , drop = FALSE])),
                                 error = function(e) NULL)
                lpTest <- fam$lp(fb, X)
                mTest <- tryCatch(fam$metrics(y, lpTest), error = function(e) NULL)
                if (is.null(mApp) || is.null(mTest)) next
                opt[b, ] <- mApp - mTest
                testSlope[b] <- mTest[["slope"]]
                if (stability) risk[, b] <- fam$risk(lpTest)
            }
            ok <- stats::complete.cases(opt)
            optimism <- colMeans(opt[ok, , drop = FALSE])
            out <- list(method = "boot", B = B, nOk = sum(ok),
                        shrinkage = mean(testSlope[ok]))
            if (stability && sum(ok) > 1)
                out$stability <- .stability(y, fam$risk(lp0), risk[, ok, drop = FALSE],
                                            maxCurves, maxIndiv)
        } else {
            strata <- fam$strata(y)
            cvm <- matrix(NA_real_, R, length(apparent), dimnames = list(NULL, names(apparent)))
            for (r in seq_len(R)) {
                fold <- .folds(strata, k)
                oof <- rep(NA_real_, n)
                for (f in seq_len(k)) {
                    test <- fold == f
                    fb <- .try_fit(fam, X[!test, , drop = FALSE], .ysub(y, !test))
                    if (!is.null(fb)) oof[test] <- fam$lp(fb, X[test, , drop = FALSE])
                }
                if (!anyNA(oof))
                    cvm[r, ] <- tryCatch(fam$metrics(y, oof), error = function(e) NA_real_)
            }
            ok <- stats::complete.cases(cvm)
            if (!any(ok)) .fail("Cross-validation failed in every repeat: the folds are too small for this model. Use fewer folds or the bootstrap.")
            optimism <- apparent - colMeans(cvm[ok, , drop = FALSE])
            out <- list(method = "cv", k = k, R = R, nOk = sum(ok))
        }
    })
    out$apparent <- apparent
    out$optimism <- optimism
    out$corrected <- apparent - optimism
    out
}

# Stability of individual risks across bootstrap models applied to the original
# data: MAPE (mean absolute prediction error versus the original model's risk),
# a sample of points for the instability plot, and calibration curves.
.stability <- function(y, p0, risk, maxCurves, maxIndiv) {
    n <- length(p0); B <- ncol(risk)
    ape <- rowMeans(abs(risk - p0))
    ind <- if (n > maxIndiv) sort(sample.int(n, maxIndiv)) else seq_len(n)
    bs <- if (B > 50) sort(sample.int(B, 50)) else seq_len(B)
    grid <- seq(stats::quantile(p0, 0.01), stats::quantile(p0, 0.99), length.out = 50)
    cb <- if (B > maxCurves) sort(sample.int(B, maxCurves)) else seq_len(B)
    curves <- vapply(cb, function(b) .cal_smooth(y, risk[, b], grid), numeric(length(grid)))
    q <- apply(risk, 1, stats::quantile, probs = c(0.025, 0.975))
    # stored compactly: risks in thousandths as integers (results state is size-limited)
    list(mape = mean(ape), maxape = max(ape),
         p0 = round(p0[ind] * 1000), pts = matrix(as.integer(round(risk[ind, bs] * 1000)), length(ind)),
         lo = round(q[1, ind] * 1000), hi = round(q[2, ind] * 1000),
         grid = grid, curves = round(curves, 4))
}
