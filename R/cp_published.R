# Published models: a model exported by clinpred (JSON: the full spec, with
# splines, scaling and levels) or a table of coefficients (CSV and friends:
# term, estimate, optional SE) is read, matched to the data, and gives the
# design matrix and coefficients the family functions need. Terms in a table
# follow a small grammar that is parsed, never evaluated as R code:
#   age                  numeric variable
#   sex = male           indicator of one level (also "sex: male", or R's "sexmale")
#   log(x), ln(x), log2(x), log10(x), sqrt(x), exp(x), x^2, x/10
#   age * sex = male     product of the above ("×" also works)
# and special rows: Intercept (constant); for Cox, "baseline survival" (S0 at
# the horizon), "horizon", and "centre" (subtracted from the prognostic index).

# model types as named in the panel (a function: translated when called)
.cp_family_title <- function() c(binary = .("Binary (logistic)"), gaussian = .("Continuous (linear)"), cox = .("Time to event (Cox)"))

# ---- export ----------------------------------------------------------------

# JSON text of a fitted clinpred model. Arrays are wrapped in I() so that
# length-1 vectors stay arrays.
# missing: the missing-category coding (cp_miss_category), so that the indicator
# predictors can be rebuilt from other data
cp_model_json <- function(spec, beta, family, se = NULL, version = "", development = NULL, missing = NULL) {
    vars <- lapply(spec$vars, function(info) {
        out <- list(name = info$name, role = info$role)
        if (info$role == "factor") out$levels <- I(info$levels)
        else {
            out$range <- I(info$range); out$mean <- info$mean; out$sd <- info$sd
            out$quartiles <- I(info$quartiles)
            if (info$role == "spline") out$knots <- I(info$knots)
            if (info$role == "numeric") { out$scale <- info$scale; out$scaleLabel <- info$scaleLabel }
        }
        out
    })
    b <- beta; b[is.na(b)] <- 0
    coefs <- lapply(seq_along(spec$columns), function(i) list(term = spec$columns[i], estimate = unname(b[i + 1])))
    if (!is.null(se)) for (i in seq_along(coefs)) coefs[[i]]$se <- unname(se[i + 1])
    oc <- spec$outcome
    out <- list(format = "clinpred-model", formatVersion = 1, clinpred = version, family = family,
                created = format(Sys.Date()),
                outcome = c(list(name = oc$name), if (!is.null(oc$event)) list(event = oc$event, levels = I(oc$levels)),
                            if (!is.null(oc$time)) list(time = oc$time)),
                predictors = unname(vars),
                interactions = unname(lapply(spec$interactions, I)))
    if (family == "cox") {
        out$baseline <- list(horizon = oc$horizon, survival = exp(-exp(unname(beta[1]))))
    } else {
        out$intercept <- list(estimate = unname(beta[1]), se = if (!is.null(se)) unname(se[1]))
    }
    out$coefficients <- coefs
    if (!is.null(missing) && length(c(missing$num, missing$fac)))
        out$missing <- list(method = "indicator", medians = as.list(missing$median),
                            factors = I(missing$fac))
    if (!is.null(development)) out$development <- development
    jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE, null = "null")
}

# ---- import ----------------------------------------------------------------

# A model from a file: list(kind = "spec"|"terms", family, ...)
cp_model_read <- function(path, filename = path) {
    ext <- tolower(tools::file_ext(filename))
    if (ext == "json") return(.model_from_json(path, filename))
    tab <- cp_read_external(path, filename)
    .model_from_table(tab, basename(filename))
}

.model_from_json <- function(path, filename) {
    js <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
                   error = function(e) .fail(sprintf(.("The file '%s' is not valid JSON."), basename(filename))))
    if (!identical(js$format, "clinpred-model"))
        .fail(sprintf(.("The file '%s' is not a model exported by clinpred. A published model can also be entered as a table with the columns term, estimate and (optionally) SE."), basename(filename)))
    vars <- list()
    for (p in js$predictors) {
        info <- list(name = p$name, role = p$role)
        if (p$role == "factor") info$levels <- unlist(p$levels)
        else {
            info$range <- unlist(p$range); info$mean <- p$mean; info$sd <- p$sd
            info$quartiles <- unlist(p$quartiles)
            if (p$role == "spline") info$knots <- unlist(p$knots)
            if (p$role == "numeric") { info$scale <- p$scale %||% 1; info$scaleLabel <- p$scaleLabel %||% "per unit" }
        }
        vars[[p$name]] <- info
    }
    ints <- list()
    for (ab in js$interactions) { ab <- unlist(ab); ints[[paste(ab, collapse = ":")]] <- ab }
    spec <- .spec_columns(list(outcome = js$outcome, vars = vars, interactions = ints))
    terms <- vapply(js$coefficients, function(cf) cf$term, "")
    if (!identical(terms, spec$columns))
        .fail(sprintf(.("The model file '%s' is damaged: its coefficients do not match its predictors."), basename(filename)))
    est <- vapply(js$coefficients, function(cf) as.numeric(cf$estimate), 0)
    se <- vapply(js$coefficients, function(cf) as.numeric(cf$se %||% NA), 0)
    out <- list(kind = "spec", family = js$family, spec = spec, terms = terms, estimate = est, se = se,
                source = basename(filename), outcome = js$outcome)
    # missing-indicator coding: development medians and the factors with a "(missing)" level
    if (!is.null(js$missing)) {
        med <- unlist(js$missing$medians)
        out$missing <- list(num = names(med), fac = unlist(js$missing$factors), median = med,
                            indicator = stats::setNames(paste(names(med), "(missing)"), names(med)))
    }
    if (js$family == "cox") {
        out$horizon <- js$baseline$horizon
        out$baseline <- js$baseline$survival
        out$centre <- 0
    } else {
        out$intercept <- js$intercept$estimate
        out$interceptSE <- as.numeric(js$intercept$se %||% NA)
    }
    out
}

.find_col <- function(nms, candidates) {
    k <- match(candidates, tolower(trimws(nms)))
    k <- k[!is.na(k)]
    if (length(k)) nms[k[1]] else NA_character_
}

.model_from_table <- function(tab, source) {
    nms <- names(tab)
    tcol <- .find_col(nms, c("term", "terms", "variable", "predictor", "parameter", "name"))
    if (is.na(tcol)) tcol <- nms[1]
    ecol <- .find_col(nms, c("estimate", "coef", "coefficient", "beta", "b", "β", "est", "log(or)", "logor", "log(hr)", "loghr"))
    if (is.na(ecol)) {
        if (length(nms) < 2) .fail(sprintf(.("The model table '%s' needs the columns term and estimate."), source))
        ecol <- setdiff(nms, tcol)[1]
    }
    scol <- .find_col(nms, c("se", "ee", "std.error", "std. error", "stderr", "standard error", "s.e.", "e.e.", "sd"))
    term <- trimws(as.character(tab[[tcol]]))
    num <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(x)), fixed = TRUE)))
    est <- num(tab[[ecol]])
    se <- if (!is.na(scol)) num(tab[[scol]]) else rep(NA_real_, length(term))
    keep <- nzchar(term) & !is.na(term)
    term <- term[keep]; est <- est[keep]; se <- se[keep]
    if (anyNA(est))
        .fail(sprintf(.("The model table '%s': the estimate of '%s' is not a number."), source, term[is.na(est)][1]))
    low <- tolower(term)
    special <- function(pattern) which(grepl(pattern, low))
    iInt <- special("^\\(?(intercept|constant|const|_cons)\\)?$")
    iS0 <- special("^(baseline survival|baseline|s0)\\b")
    iH <- special("^(horizon|time horizon)$")
    iC <- special("^(centre|center|mean pi|pi mean|mean lp|lp mean)$")
    rest <- setdiff(seq_along(term), c(iInt, iS0, iH, iC))
    out <- list(kind = "terms", family = NA_character_, source = source,
                terms = term[rest], estimate = est[rest], se = se[rest],
                intercept = if (length(iInt)) est[iInt[1]] else NULL,
                interceptSE = if (length(iInt)) se[iInt[1]] else NA_real_,
                baseline = if (length(iS0)) est[iS0[1]] else NULL,
                horizon = if (length(iH)) est[iH[1]] else NULL,
                centre = if (length(iC)) est[iC[1]] else 0)
    if (!length(out$terms)) .fail(sprintf(.("The model table '%s' has no predictor terms."), source))
    if (!is.null(out$baseline) && !(out$baseline > 0 && out$baseline < 1))
        .fail(sprintf(.("The baseline survival must be a probability between 0 and 1 (it is %s)."), .num(out$baseline)))
    out
}

# ---- terms against the data ----------------------------------------------

# One component of a term, resolved against the data columns
.parse_component <- function(s, data) {
    s <- trimws(s)
    cols <- names(data)
    col <- function(v) .match_column(trimws(v), cols)
    m <- regmatches(s, regexec("^(log|ln|log2|log10|sqrt|exp)\\s*\\((.+)\\)$", s, ignore.case = TRUE))[[1]]
    if (length(m)) return(list(type = "fn", fn = tolower(m[2]), var = col(m[3]), name = trimws(m[3])))
    m <- regmatches(s, regexec("^(.+?)\\s*\\^\\s*([0-9.]+)$", s))[[1]]
    if (length(m)) return(list(type = "pow", power = as.numeric(m[3]), var = col(m[2]), name = trimws(m[2])))
    m <- regmatches(s, regexec("^(.+?)\\s*/\\s*([0-9.eE+-]+)$", s))[[1]]
    if (length(m)) return(list(type = "div", by = as.numeric(m[3]), var = col(m[2]), name = trimws(m[2])))
    m <- regmatches(s, regexec("^(.+?)\\s*[=:]\\s*(.+)$", s))[[1]]
    if (length(m) && !is.na(col(m[2]))) return(list(type = "level", var = col(m[2]), level = trimws(m[3]), name = trimws(m[2])))
    if (!is.na(col(s))) return(list(type = "num", var = col(s), name = s))
    # R's naming of factor dummies: the variable name followed by the level
    for (v in cols) {
        x <- data[[v]]
        if (!(is.factor(x) || is.character(x)) || !startsWith(tolower(s), tolower(v))) next
        rest <- substring(s, nchar(v) + 1)
        lev <- if (is.factor(x)) levels(x) else unique(as.character(x))
        hit <- lev[tolower(lev) == tolower(rest)]
        if (length(hit)) return(list(type = "level", var = v, level = hit[1], name = v))
    }
    list(type = "missing", var = NA_character_, name = s)
}

.component_values <- function(cp, data) {
    x <- data[[cp$var]]
    if (cp$type == "level") {
        xs <- as.character(x)
        lv <- .match_levels(cp$level, unique(xs[!is.na(xs)]))
        return(ifelse(is.na(xs), NA_real_, as.numeric(!is.na(lv) & xs == lv)))
    }
    xn <- if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
    if (any(!is.na(x)) && all(is.na(xn[!is.na(x)])))
        .fail(sprintf(.("'%s' is categorical in the data: name the level in the model, e.g. '%s = %s'."), cp$name, cp$name, if (is.factor(x)) levels(x)[min(2, nlevels(x))] else "level"))
    switch(cp$type,
        num = xn,
        div = xn / cp$by,
        pow = xn^cp$power,
        fn = switch(cp$fn, log = , ln = log(xn), log2 = log2(xn), log10 = log10(xn), sqrt = sqrt(xn), exp = exp(xn)))
}

# Design matrix of a coefficient table, with a check of every variable used
cp_terms_design <- function(model, data) {
    comps <- lapply(model$terms, function(t) lapply(strsplit(t, "\\s*[*×]\\s*")[[1]], .parse_component, data = data))
    miss <- unique(unlist(lapply(comps, function(cs) vapply(Filter(function(c) c$type == "missing", cs), `[[`, "", "name"))))
    if (length(miss))
        .fail(sprintf(.("No variable matches the model terms %s. Add the variables the model uses to Predictor variables, and name each term after one of them (e.g. 'age', 'sex = male', 'log(crp)', 'age * sex = male')."),
                      paste0("'", miss, "'", collapse = ", ")))
    X <- vapply(comps, function(cs) Reduce(`*`, lapply(cs, .component_values, data = data)), numeric(nrow(data)))
    X <- matrix(X, nrow = nrow(data), dimnames = list(NULL, model$terms))
    # one check row per variable
    vars <- unique(unlist(lapply(comps, function(cs) vapply(cs, `[[`, "", "var"))))
    check <- do.call(rbind, lapply(vars, function(v) {
        notes <- character(0); status <- "OK"
        cs <- Filter(function(c) identical(c$var, v), unlist(comps, recursive = FALSE))
        if (any(vapply(cs, function(c) c$name != v, TRUE))) notes <- c(notes, sprintf(.("column '%s'"), v))
        lv <- Filter(function(c) c$type == "level", cs)
        if (length(lv)) {
            xs <- as.character(data[[v]])
            inData <- unique(xs[!is.na(xs)])
            used <- character(0)
            for (c in lv) {
                m <- .match_levels(c$level, inData)
                if (is.na(m)) {
                    status <- "Warning"
                    notes <- c(notes, sprintf(.("level '%s' not in the data (its indicator is always 0)"), c$level))
                } else used <- c(used, m)
            }
            # data levels the model does not name get every indicator 0: the reference
            ref <- setdiff(inData, used)
            if (length(ref)) notes <- c(notes, sprintf(.("as reference: %s"), paste(sort(ref), collapse = ", ")))
            # several levels the model does not name all fall into the reference: one of them
            # is probably a category the model never saw (e.g. 'unknown')
            if (length(ref) > 1) status <- "Warning"
        }
        mis <- mean(is.na(data[[v]]))
        if (mis > 0) notes <- c(notes, sprintf(.("%s missing"), .pct(mis)))
        data.frame(var = v, role = "Predictor", status = status, details = paste(unique(notes), collapse = "; "),
                   stringsAsFactors = FALSE)
    }))
    bad <- colSums(!is.finite(X) & !is.na(X)) > 0
    if (any(bad))
        .fail(sprintf(.("The term '%s' gives values that are not finite (e.g. the log of 0 or of a negative value)."), model$terms[bad][1]))
    list(X = X, check = check)
}

# Design matrix and checks for either kind of model
cp_model_design <- function(model, data) {
    if (model$kind == "terms") return(cp_terms_design(model, data))
    # indicator predictors are made from the data, as in development
    if (!is.null(model$missing)) data <- cp_miss_category_apply(model$missing, data)
    miss <- names(model$spec$vars)[is.na(vapply(names(model$spec$vars), .match_column, "", available = names(data)))]
    if (length(miss))
        .fail(sprintf(.("The model uses %s: add the variables to Predictor variables."), paste0("'", miss, "'", collapse = ", ")))
    h <- cp_harmonise(data, model$spec, character(0), "Data")
    X <- matrix(NA_real_, nrow(data), length(model$spec$columns), dimnames = list(NULL, model$spec$columns))
    # cp_harmonise keeps the complete rows only, named by their row number
    if (nrow(h$data)) X[as.integer(rownames(h$data)), ] <- cp_design(model$spec, h$data)
    list(X = X, check = h$check)
}

# Coefficients with the intercept first (log H0(h) for Cox, NA when the
# baseline must come from the data)
cp_model_beta <- function(model) {
    b <- model$estimate
    names(b) <- model$terms
    # Cox: lp = log H0(h) + PI − centre, with H0(h) = −log S0(h)
    c(`(Intercept)` = if (model$family == "cox") {
        if (is.null(model$baseline)) NA_real_ else log(-log(model$baseline)) - (model$centre %||% 0)
    } else model$intercept %||% 0, b)
}
