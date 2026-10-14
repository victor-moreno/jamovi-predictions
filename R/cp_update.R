# Model updating in validation data (Steyerberg 2019, ch. 20; Vergouwe et al.
# 2017), family-agnostic: the original model, recalibration in the large (new
# intercept), recalibration (new intercept and slope of the linear predictor)
# and model revision (every coefficient re-estimated). The closed testing
# procedure tests the revised model against the original, then against
# recalibration in the large, then against recalibration, and keeps the first
# method that is not significantly worse than revision.

cp_update <- function(fam, X, y, lp, alpha = 0.05) {
    if (identical(fam$name, "cox")) return(.cp_update_cox(fam, X, y, lp, alpha))
    n <- NROW(X)
    # calibration-in-the-large only (the full metric set would fit the loess smoother too)
    a1 <- if (identical(fam$name, "gaussian")) mean(y - lp) else unname(.cal_glm(y, lp)$citl)
    f2 <- fam$fit(matrix(lp, ncol = 1, dimnames = list(NULL, "lp")), y)
    br <- tryCatch(fam$fit(X, y), error = function(e) NULL)
    revOk <- !is.null(br) && isTRUE(attr(br, "converged"))
    methods <- list(
        original = list(a = 0, b = 1, lp = lp, k = 0),
        citl = list(a = a1, b = 1, lp = lp + a1, k = 1),
        recal = list(a = unname(f2[1]), b = unname(f2[2]), lp = fam$lp(f2, matrix(lp, ncol = 1)), k = 2),
        revision = if (revOk) list(a = NA_real_, b = NA_real_, lp = fam$lp(br, X),
                                   k = sum(!is.na(br)), beta = br) else NULL)
    methods <- Filter(Negate(is.null), methods)
    for (m in names(methods)) methods[[m]]$loss <- fam$loss(y, methods[[m]]$lp)

    chosen <- NA_character_
    if (revOk) {
        rev <- methods$revision
        dfResid <- n - rev$k
        for (m in c("original", "citl", "recal")) {
            t <- .nested_test(fam, methods[[m]]$loss, rev$loss, rev$k - methods[[m]]$k, n, dfResid)
            methods[[m]]$test <- t
            if (is.na(chosen) && t[["p"]] >= alpha) chosen <- m
        }
        if (is.na(chosen)) chosen <- "revision"
    }
    list(methods = methods, chosen = chosen)
}

# Cox models (van Houwelingen 2000; Steyerberg 2019, ch. 20): the baseline
# hazard at the horizon is re-estimated with the linear predictor as offset
# (recalibration in the large), with a re-estimated slope (recalibration), or
# every coefficient is refitted (revision). The partial likelihood does not
# involve the baseline hazard, so the closed testing procedure is adapted into
# sequential tests from the most complex method down: revision vs recalibration
# (partial likelihood ratio, df = parameters - 1), slope = 1 (partial likelihood
# ratio, 1 df), observed/expected risk at the horizon = 1 (Wald test of log O/E,
# Kaplan-Meier variance); the first significant test picks its method.
.cp_update_cox <- function(fam, X, y, lp, alpha = 0.05) {
    h <- fam$horizon; tm <- y[, 1]; st <- y[, 2]
    a1 <- log(.cox_H0(tm, st, lp, h))
    b <- .cox_slope(y, lp)[1]
    br <- tryCatch(fam$fit(X, y), error = function(e) NULL)
    revOk <- !is.null(br) && isTRUE(attr(br, "converged")) && ncol(X) > 1
    methods <- list(original = list(a = 0, b = 1, lp = lp, k = 0),
                    citl = list(a = a1, b = 1, lp = lp + a1, k = 1))
    if (is.finite(b)) {
        a2 <- log(.cox_H0(tm, st, b * lp, h))
        methods$recal <- list(a = a2, b = b, lp = a2 + b * lp, k = 2)
    }
    if (revOk) methods$revision <- list(a = NA_real_, b = NA_real_, lp = fam$lp(br, X), k = sum(!is.na(br)), beta = br)
    for (m in names(methods)) methods[[m]]$loss <- fam$loss(y, methods[[m]]$lp)
    lrt <- function(big, small, df) {
        chi <- max(0, 2 * (methods[[small]]$loss - methods[[big]]$loss))
        c(stat = chi, df = df, df2 = NA_real_, p = stats::pchisq(chi, df, lower.tail = FALSE))
    }
    if (revOk && !is.null(methods$recal)) methods$revision$test <- lrt("revision", "recal", ncol(X) - 1)
    if (!is.null(methods$recal)) methods$recal$test <- lrt("recal", "original", 1)
    km <- .km(tm, st, h); O <- 1 - km$surv; E <- mean(.cox_risk(lp))
    z <- if (O > 0 && km$var > 0) log(O / E) / (sqrt(km$var) / O) else NA_real_
    methods$citl$test <- c(stat = z^2, df = 1, df2 = NA_real_, p = stats::pchisq(z^2, 1, lower.tail = FALSE))
    sig <- function(m) !is.null(methods[[m]]$test) && isTRUE(methods[[m]]$test[["p"]] < alpha)
    chosen <- if (sig("revision")) "revision" else if (sig("recal")) "recal" else
        if (sig("citl")) "citl" else if (is.na(z)) NA_character_ else "original"
    list(methods = methods, chosen = chosen)
}
