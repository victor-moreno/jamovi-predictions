# Explaining the model, family-agnostic: interaction selection, effects as
# contrasts between two covariate profiles (crude and adjusted), likelihood
# ratio tests per predictor, effects within subgroups, descriptive table.

# Full-length coefficients from a fit on a subset of the columns (the others 0)
.expand_beta <- function(b, keep, cols, intercept) {
    full <- stats::setNames(numeric(length(cols)), cols)
    full[keep] <- if (intercept) b[-1] else b
    out <- if (intercept) c(`(Intercept)` = unname(b[1]), full) else full
    attr(out, "converged") <- attr(b, "converged")
    attr(out, "separation") <- attr(b, "separation")
    out
}

.fit_cols <- function(X, y, fam, keep) {
    b <- fam$fit(X[, keep, drop = FALSE], y)
    b[is.na(b)] <- 0
    list(beta = b, loss = fam$loss(y, fam$lp(b, X[, keep, drop = FALSE])), df = sum(b != 0))
}

# Interaction terms: screening of every candidate added alone to the
# main-effects model, then forward selection (AIC, BIC or likelihood ratio test
# at alpha) or all candidates kept. Returns full-length coefficients (0 for the
# terms left out), with the selected terms and the screening table as
# attributes. Used as the fitting procedure in the bootstrap, so that the
# selection is repeated in every bootstrap sample (honest optimism).
cp_select_interactions <- function(X, y, spec, fam, method = "all", alpha = 0.05) {
    terms <- names(spec$interactions)
    cols <- colnames(X)
    if (!length(terms) || method == "all") {
        b <- fam$fit(X, y)
        attr(b, "selected") <- terms
        return(b)
    }
    n <- NROW(X)
    inTerm <- function(t) spec$colTerm == t
    current <- !(spec$colTerm %in% terms)
    base <- .fit_cols(X, y, fam, current)
    screen <- NULL
    selected <- character(0)
    remaining <- terms
    repeat {
        stats1 <- lapply(remaining, function(t) {
            f <- tryCatch(.fit_cols(X, y, fam, current | inTerm(t)), error = function(e) NULL)
            if (is.null(f)) return(NULL)
            chi <- max(0, 2 * (base$loss - f$loss)); df <- f$df - base$df
            if (df < 1) return(NULL)
            data.frame(term = t, chisq = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE),
                       daic = -chi + 2 * df, dbic = -chi + log(n) * df, stringsAsFactors = FALSE)
        })
        st <- do.call(rbind, stats1)
        if (is.null(screen)) screen <- st
        if (is.null(st) || !nrow(st)) break
        score <- switch(method, aic = st$daic, bic = st$dbic, lrt = st$p)
        best <- which.min(score)
        ok <- switch(method, aic = score[best] < 0, bic = score[best] < 0, lrt = score[best] < alpha)
        if (!ok) break
        t <- st$term[best]
        selected <- c(selected, t)
        current <- current | inTerm(t)
        base <- .fit_cols(X, y, fam, current)
        remaining <- setdiff(remaining, t)
        if (!length(remaining)) break
    }
    out <- .expand_beta(base$beta, current, cols, fam$intercept)
    attr(out, "selected") <- selected
    attr(out, "screen") <- screen
    out
}

# One-row covariate profile: numeric predictors at their median, factors at
# their reference level, except the values given
.profile <- function(spec, set = list()) {
    df <- data.frame(row.names = 1)
    for (info in spec$vars)
        df[[info$name]] <- if (info$role == "factor")
            factor(info$levels[1], levels = info$levels) else info$quartiles[2]
    for (v in names(set))
        df[[v]] <- if (spec$vars[[v]]$role == "factor")
            factor(set[[v]], levels = spec$vars[[v]]$levels) else set[[v]]
    df
}

# Log effect (log-odds ratio, log hazard ratio) of profile p1 versus p0, with SE
cp_contrast <- function(spec, beta, V, p1, p0, intercept = TRUE) {
    L <- cp_design(spec, p1) - cp_design(spec, p0)
    if (intercept) L <- cbind(0, L)
    b <- beta; b[is.na(b)] <- 0
    est <- sum(L * b)
    se <- sqrt(max(0, drop(L %*% V %*% t(L))))
    c(est = est, se = se)
}

# The effect rows of a variable: factor levels against the reference, a
# covariate per its scale unit (from its median), a spline Q3 versus Q1
cp_effect_defs <- function(spec, v) {
    info <- spec$vars[[v]]
    switch(info$role,
        factor = lapply(info$levels[-1], function(l)
            list(key = paste0(v, ": ", l), level = l, set1 = stats::setNames(list(l), v),
                 set0 = stats::setNames(list(info$levels[1]), v))),
        numeric = list(list(key = v, level = info$scaleLabel,
                            set1 = stats::setNames(list(info$quartiles[2] + info$scale), v),
                            set0 = stats::setNames(list(info$quartiles[2]), v))),
        spline = list(list(key = v,
                           level = sprintf("Q3 vs Q1 (%s vs %s)", .num(info$quartiles[3], 4), .num(info$quartiles[1], 4)),
                           set1 = stats::setNames(list(info$quartiles[3]), v),
                           set0 = stats::setNames(list(info$quartiles[1]), v))))
}

# Effects of every predictor in a fitted model, at the reference profile
cp_effects <- function(spec, beta, V, conf = 0.95, intercept = TRUE, vars = names(spec$vars)) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    rows <- list()
    for (v in vars) for (d in cp_effect_defs(spec, v)) {
        cc <- cp_contrast(spec, beta, V, .profile(spec, d$set1), .profile(spec, d$set0), intercept)
        rows[[length(rows) + 1]] <- data.frame(key = d$key, var = v, level = d$level,
            est = exp(cc[["est"]]), lower = exp(cc[["est"]] - z * cc[["se"]]),
            upper = exp(cc[["est"]] + z * cc[["se"]]),
            p = 2 * stats::pnorm(-abs(cc[["est"]] / cc[["se"]])), stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

# A spec with a single predictor (for crude effects)
.spec_single <- function(spec, v) {
    s <- spec; s$vars <- spec$vars[v]; s$interactions <- list()
    .spec_columns(s)
}

# Crude effects: each predictor alone, on the same rows, with the same coding;
# the overall p of each predictor is a likelihood ratio test against the null model
cp_crude <- function(spec, df, y, fam, conf = 0.95) {
    rows <- list(); pAll <- c()
    for (v in names(spec$vars)) {
        sv <- .spec_single(spec, v)
        Xv <- cp_design(sv, df)
        b <- fam$fit(Xv, y)
        V <- fam$vcov(b)
        rows[[v]] <- cp_effects(sv, b, V, conf, fam$intercept)
        chi <- 2 * (fam$nullLoss(y) - fam$loss(y, fam$lp(b, Xv)))
        pAll[[v]] <- stats::pchisq(chi, ncol(Xv), lower.tail = FALSE)
    }
    out <- do.call(rbind, rows)
    attr(out, "pOverall") <- pAll
    out
}

# Likelihood ratio tests in the final model: each predictor overall (with its
# interactions), the non-linear part of each spline, each interaction term
cp_lr_tests <- function(spec, X, y, fam) {
    full <- .fit_cols(X, y, fam, rep(TRUE, ncol(X)))
    test <- function(drop) {
        red <- .fit_cols(X, y, fam, !drop)
        chi <- max(0, 2 * (red$loss - full$loss)); df <- full$df - red$df
        c(chisq = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE))
    }
    rows <- list()
    add <- function(key, label, drop)
        rows[[length(rows) + 1]] <<- c(list(key = key, term = label), as.list(test(drop)))
    for (info in spec$vars) {
        v <- info$name
        involves <- vapply(spec$colVars, function(vs) v %in% vs, TRUE)
        add(v, v, involves)
        if (info$role == "spline") {
            nl <- paste0(v, strrep("'", seq_len(length(info$knots) - 2)))
            dropNl <- spec$columns %in% nl |
                vapply(spec$colParts, function(ab) any(ab %in% nl), TRUE)
            add(paste(v, "nonlinear"), paste0("   ", v, ": non-linear"), dropNl)
        }
    }
    for (key in names(spec$interactions))
        add(key, paste(spec$interactions[[key]], collapse = " × "), spec$colTerm == key)
    do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# Effects of each variable of an interaction within subgroups of the other:
# factor levels, or the quartiles of a numeric variable
cp_subgroups <- function(spec, beta, V, conf = 0.95, intercept = TRUE) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    rows <- list()
    for (key in names(spec$interactions)) {
        ab <- spec$interactions[[key]]
        for (k in 1:2) {
            a <- ab[k]; b <- ab[3 - k]
            ib <- spec$vars[[b]]
            within <- if (ib$role == "factor") stats::setNames(as.list(ib$levels), paste(b, "=", ib$levels))
                else stats::setNames(as.list(ib$quartiles),
                                     sprintf("%s = %s (%s)", b, .num(ib$quartiles, 4), c("Q1", "median", "Q3")))
            for (d in cp_effect_defs(spec, a)) for (w in names(within)) {
                s1 <- c(d$set1, stats::setNames(list(within[[w]]), b))
                s0 <- c(d$set0, stats::setNames(list(within[[w]]), b))
                cc <- cp_contrast(spec, beta, V, .profile(spec, s1), .profile(spec, s0), intercept)
                rows[[length(rows) + 1]] <- data.frame(
                    key = paste(key, d$key, w), term = key, effect = paste0(a, ": ", d$level), within = w,
                    est = exp(cc[["est"]]), lower = exp(cc[["est"]] - z * cc[["se"]]),
                    upper = exp(cc[["est"]] + z * cc[["se"]]),
                    p = 2 * stats::pnorm(-abs(cc[["est"]] / cc[["se"]])), stringsAsFactors = FALSE)
            }
        }
    }
    if (length(rows)) do.call(rbind, rows) else NULL
}

# Descriptive table by outcome group: mean (SD) or median [Q1, Q3] for numeric
# predictors, n (%) per level for factors, and the missing values of each
cp_describe <- function(data, group, groupLevels, numVars, facVars, stat = "meanSD") {
    g <- factor(as.character(data[[group]]), levels = groupLevels)
    sets <- c(lapply(groupLevels, function(l) which(!is.na(g) & g == l)), list(seq_len(nrow(data))))
    names(sets) <- c(paste0("g", seq_along(groupLevels)), "total")
    fmtNum <- function(x) {
        x <- x[!is.na(x)]
        if (!length(x)) return("")
        if (stat == "meanSD") sprintf("%.2f (%.2f)", mean(x), if (length(x) > 1) stats::sd(x) else NA)
        else { q <- stats::quantile(x, c(0.5, 0.25, 0.75)); sprintf("%.2f [%.2f, %.2f]", q[1], q[2], q[3]) }
    }
    rows <- list(list(key = "n", var = "N", level = "",
                      values = lapply(sets, function(i) as.character(length(i))), missing = ""))
    for (v in numVars)
        rows[[length(rows) + 1]] <- list(key = v, var = v, level = if (stat == "meanSD") "mean (SD)" else "median [Q1, Q3]",
                                         values = lapply(sets, function(i) fmtNum(data[[v]][i])),
                                         missing = as.character(sum(is.na(data[[v]]))))
    for (v in facVars) {
        x <- factor(data[[v]])
        for (j in seq_along(levels(x))) {
            l <- levels(x)[j]
            rows[[length(rows) + 1]] <- list(key = paste(v, l), var = v, level = l,
                values = lapply(sets, function(i) {
                    xi <- x[i]; nn <- sum(!is.na(xi))
                    if (!nn) "" else sprintf("%d (%.1f%%)", sum(xi == l, na.rm = TRUE), 100 * mean(xi == l, na.rm = TRUE))
                }),
                missing = if (j == 1) as.character(sum(is.na(x))) else "")
        }
    }
    rows
}
