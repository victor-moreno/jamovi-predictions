# Development data and the model "spec": everything needed to rebuild the design
# matrix on new data (bootstrap samples, a validation split, an external file).
# Family-agnostic: the outcome part of the spec is filled by the family.

# Harrell's default knot placement (Hmisc::rcspline.eval, without its
# tied-extremes rule): quantiles, outer knots at the 5th smallest/largest
# value when n < 100
.rcs_knots <- function(x, nk, name) {
    x <- x[!is.na(x)]
    nxu <- length(unique(x))
    if (nxu - 2 <= nk)
        .fail(sprintf(.("'%s' has only %s distinct values, too few for a spline with %s knots. Use fewer knots or enter it as a linear covariate."), name, nxu, nk))
    outer <- if (nk > 6) 0.025 else if (nk > 3) 0.05 else 0.1
    knots <- unname(stats::quantile(x, seq(outer, 1 - outer, length.out = nk)))
    if (length(x) < 100) {
        xs <- sort(x)
        knots[1] <- xs[5]
        knots[nk] <- xs[length(xs) - 4]
    }
    knots <- sort(unique(knots))
    if (length(knots) < 3)
        .fail(sprintf(.("Could not place 3 distinct knots for '%s': most of its values are tied."), name))
    knots
}

# Restricted cubic spline basis (non-linear part only; the linear term is x
# itself), normalised as Hmisc's norm = 2
.rcs_basis <- function(x, knots) {
    k <- length(knots)
    tk <- knots[k]; tk1 <- knots[k - 1]
    pp <- function(u) pmax(u, 0)^3
    out <- vapply(seq_len(k - 2), function(j) {
        tj <- knots[j]
        (pp(x - tj) - pp(x - tk1) * (tk - tj) / (tk - tk1) +
             pp(x - tk) * (tk1 - tj) / (tk - tk1)) / (tk - knots[1])^2
    }, numeric(length(x)))
    matrix(out, nrow = length(x))
}

# Covariate scaling, as in the prognosis module: the covariate is divided by
# the scale before fitting, so its odds ratio is per that many units
.scale_value <- function(key, x) switch(key %||% "unit",
    unit = 1, sd = stats::sd(x), ten = 10, hundred = 100, thousand = 1000, 1)
.scale_label <- function(key, value) switch(key %||% "unit",
    unit = .("per unit"), sd = sprintf(.("per 1 SD (%s)"), .num(value, 3)), ten = .("per 10 units"),
    hundred = .("per 100 units"), thousand = .("per 1000 units"), .("per unit"))

# The roles a variable can take and their checks.
# covs: linear numeric; splineVars: numeric with RCS; factors: treatment dummies
# against a reference level (refLevels[[v]], else the first level of the data).
# scales: named list of scale keys for covs. interactions: list of pairs of
# predictor names (candidates; cp_with_interactions() keeps a subset).
cp_prepare <- function(data, outcome, covs = character(0), splineVars = character(0),
                       factors = character(0), nknots = 4, refLevels = list(),
                       scales = list(), interactions = list()) {
    preds <- c(covs, splineVars, factors)
    if (length(preds) == 0) .fail(.("Add at least one predictor."))
    dup <- unique(preds[duplicated(preds)])
    if (length(dup))
        .fail(sprintf(.("'%s' is entered more than once; each predictor goes in one box only."), dup[1]))
    if (any(outcome %in% preds))
        .fail(.("The outcome cannot also be a predictor."))

    keep <- c(outcome, preds)
    df <- data[, keep, drop = FALSE]
    for (v in c(covs, splineVars)) {
        x <- df[[v]]
        if (is.factor(x)) x <- suppressWarnings(as.numeric(as.character(x)))
        df[[v]] <- as.numeric(x)
    }
    for (v in factors) df[[v]] <- factor(df[[v]], ordered = FALSE)

    complete <- stats::complete.cases(df)
    nMissing <- sum(!complete)
    df <- df[complete, , drop = FALSE]
    for (v in factors) df[[v]] <- droplevels(df[[v]])

    vars <- list()
    for (v in preds) {
        x <- df[[v]]
        role <- if (v %in% factors) "factor" else if (v %in% splineVars) "spline" else "numeric"
        info <- list(name = v, role = role)
        if (role == "factor") {
            if (nlevels(x) < 2)
                .fail(sprintf(.("Factor '%s' has a single level in the complete cases; remove it."), v))
            lev <- levels(x)
            ref <- refLevels[[v]]
            if (!is.null(ref) && ref %in% lev) lev <- c(ref, setdiff(lev, ref))
            info$levels <- lev
        } else {
            if (stats::sd(x) == 0) .fail(sprintf(.("'%s' does not vary in the complete cases; remove it."), v))
            info$range <- range(x)
            info$mean <- mean(x)
            info$sd <- stats::sd(x)
            info$quartiles <- unname(stats::quantile(x, c(0.25, 0.5, 0.75)))
            if (role == "spline") info$knots <- .rcs_knots(x, nknots, v)
            if (role == "numeric") {
                key <- scales[[v]] %||% "unit"
                info$scale <- .scale_value(key, x)
                info$scaleLabel <- .scale_label(key, info$scale)
            }
        }
        vars[[v]] <- info
    }
    ints <- list()
    for (t in interactions) {
        t <- unique(unlist(t))
        if (length(t) != 2 || !all(t %in% preds)) next
        key <- paste(t, collapse = ":")
        if (!(key %in% names(ints)) && !(paste(rev(t), collapse = ":") %in% names(ints)))
            ints[[key]] <- t
    }
    spec <- list(outcome = list(name = outcome), vars = vars, interactions = ints)
    spec <- .spec_columns(spec)
    list(data = df, spec = spec, nMissing = nMissing, rows = rownames(data)[complete])
}

# The same spec with only some of its interaction terms (by key "a:b")
cp_with_interactions <- function(spec, keys) {
    spec$interactions <- spec$interactions[intersect(names(spec$interactions), keys)]
    .spec_columns(spec)
}

# Columns of the design matrix. For each column: its name, its term ("glu" or
# "glu:sex"), the variables it involves, and for interaction columns the two
# main-effect columns whose product it is. Products of two non-linear spline
# columns are left out (as rms's %ia%).
.spec_columns <- function(spec) {
    cols <- character(0); colTerm <- character(0); colVars <- list(); colParts <- list()
    nonlin <- character(0)
    for (info in spec$vars) {
        nm <- switch(info$role,
            numeric = info$name,
            spline = c(info$name, paste0(info$name, strrep("'", seq_len(length(info$knots) - 2)))),
            factor = paste0(info$name, ": ", info$levels[-1]))
        if (info$role == "spline") nonlin <- c(nonlin, nm[-1])
        cols <- c(cols, nm)
        colTerm <- c(colTerm, rep(info$name, length(nm)))
        colVars <- c(colVars, rep(list(info$name), length(nm)))
        colParts <- c(colParts, rep(list(NULL), length(nm)))
    }
    for (key in names(spec$interactions)) {
        ab <- spec$interactions[[key]]
        ca <- cols[colTerm == ab[1]]; cb <- cols[colTerm == ab[2]]
        for (a in ca) for (b in cb) {
            if (a %in% nonlin && b %in% nonlin) next
            cols <- c(cols, paste(a, "×", b))
            colTerm <- c(colTerm, key)
            colVars <- c(colVars, list(ab))
            colParts <- c(colParts, list(c(a, b)))
        }
    }
    spec$columns <- cols
    spec$colTerm <- colTerm
    spec$colVar <- colTerm
    spec$colVars <- colVars
    spec$colParts <- colParts
    spec
}

# Design matrix (no intercept) for a data frame whose predictor columns are
# already of the spec's types (numeric; factor with the spec's levels or a subset)
cp_design <- function(spec, df) {
    blocks <- lapply(spec$vars, function(info) {
        x <- df[[info$name]]
        switch(info$role,
            numeric = matrix(as.numeric(x) / (info$scale %||% 1), ncol = 1),
            spline = cbind(as.numeric(x), .rcs_basis(as.numeric(x), info$knots)),
            factor = {
                xc <- as.character(x)
                vapply(info$levels[-1], function(l) as.numeric(xc == l), numeric(length(xc)))
            })
    })
    X <- do.call(cbind, lapply(blocks, function(b) matrix(b, nrow = nrow(df))))
    nMain <- ncol(X)
    colnames(X) <- spec$columns[seq_len(nMain)]
    if (length(spec$columns) > nMain) {
        P <- vapply(spec$colParts[-seq_len(nMain)], function(ab) X[, ab[1]] * X[, ab[2]],
                    numeric(nrow(X)))
        X <- cbind(X, matrix(P, nrow = nrow(X)))
        colnames(X) <- spec$columns
    }
    X
}

# Number of model parameters, excluding the intercept
cp_nparams <- function(spec) length(spec$columns)

# Predictor data of any rows (e.g. the whole spreadsheet, for saved
# predictions): the spec's types and levels, complete predictors only; unknown
# factor levels become missing
cp_newdata <- function(data, spec) {
    out <- data.frame(row.names = rownames(data))
    for (info in spec$vars) {
        x <- data[[info$name]]
        out[[info$name]] <- if (info$role == "factor")
            factor(as.character(x), levels = info$levels)
        else if (is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else as.numeric(x)
    }
    out[stats::complete.cases(out), , drop = FALSE]
}
