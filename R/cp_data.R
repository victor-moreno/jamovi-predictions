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
        .fail("'", name, "' has only ", nxu, " distinct values, too few for a spline with ",
              nk, " knots. Use fewer knots or enter it as a linear covariate.")
    outer <- if (nk > 6) 0.025 else if (nk > 3) 0.05 else 0.1
    knots <- unname(stats::quantile(x, seq(outer, 1 - outer, length.out = nk)))
    if (length(x) < 100) {
        xs <- sort(x)
        knots[1] <- xs[5]
        knots[nk] <- xs[length(xs) - 4]
    }
    knots <- sort(unique(knots))
    if (length(knots) < 3)
        .fail("Could not place 3 distinct knots for '", name, "': most of its values are tied.")
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

# The roles a variable can take and their checks.
# covs: linear numeric; splineVars: numeric with RCS; factors: treatment dummies,
# first level as reference (the order of levels in the data).
cp_prepare <- function(data, outcome, covs = character(0), splineVars = character(0),
                       factors = character(0), nknots = 4) {
    preds <- c(covs, splineVars, factors)
    if (length(preds) == 0) .fail("Add at least one predictor.")
    dup <- unique(preds[duplicated(preds)])
    if (length(dup))
        .fail("'", dup[1], "' is entered more than once; each predictor goes in one box only.")
    if (any(outcome %in% preds))
        .fail("The outcome cannot also be a predictor.")

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
                .fail("Factor '", v, "' has a single level in the complete cases; remove it.")
            info$levels <- levels(x)
        } else {
            if (stats::sd(x) == 0) .fail("'", v, "' does not vary in the complete cases; remove it.")
            info$range <- range(x)
            info$mean <- mean(x)
            info$sd <- stats::sd(x)
            if (role == "spline") info$knots <- .rcs_knots(x, nknots, v)
        }
        vars[[v]] <- info
    }
    spec <- list(outcome = list(name = outcome), vars = vars)
    spec <- .spec_columns(spec)
    list(data = df, spec = spec, nMissing = nMissing, rows = rownames(data)[complete])
}

# Column names of the design matrix and the variable each column belongs to
.spec_columns <- function(spec) {
    cols <- character(0); colVar <- character(0)
    for (info in spec$vars) {
        nm <- switch(info$role,
            numeric = info$name,
            spline = c(info$name, paste0(info$name, strrep("'", seq_len(length(info$knots) - 2)))),
            factor = paste0(info$name, ": ", info$levels[-1]))
        cols <- c(cols, nm)
        colVar <- c(colVar, rep(info$name, length(nm)))
    }
    spec$columns <- cols
    spec$colVar <- colVar
    spec
}

# Design matrix (no intercept) for a data frame whose predictor columns are
# already of the spec's types (numeric; factor with the spec's levels or a subset)
cp_design <- function(spec, df) {
    blocks <- lapply(spec$vars, function(info) {
        x <- df[[info$name]]
        switch(info$role,
            numeric = matrix(as.numeric(x), ncol = 1),
            spline = cbind(as.numeric(x), .rcs_basis(as.numeric(x), info$knots)),
            factor = {
                xc <- as.character(x)
                vapply(info$levels[-1], function(l) as.numeric(xc == l), numeric(length(xc)))
            })
    })
    X <- do.call(cbind, lapply(blocks, function(b) matrix(b, nrow = nrow(df))))
    colnames(X) <- spec$columns
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
