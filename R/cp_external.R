# Validation data: an external file or a subset of the current data. The data
# are matched to the development spec (names, types, levels), every mismatch is
# reported, and the case-mix is compared with the development data.

# Read csv / tsv / txt / omv / xlsx / sav. `filename` is the user's name for the
# file (the session copy is named by a hash, so the extension comes from it).
cp_read_external <- function(path, filename = path) {
    ext <- tolower(tools::file_ext(filename))
    df <- tryCatch(switch(ext,
        csv = , tsv = , txt = .read_delim(path, ext),
        omv = jmvReadWrite::read_omv(path),
        xlsx = openxlsx::read.xlsx(path, sheet = 1, check.names = FALSE),
        sav = suppressWarnings(foreign::read.spss(path, to.data.frame = TRUE,
                                                  use.value.labels = TRUE)),
        .fail("Files of type '.", ext, "' cannot be read. Use csv, tsv, txt, omv, xlsx or sav.")),
        cp_user_error = function(e) stop(e),
        error = function(e) .fail("The file '", basename(filename), "' could not be read as ", ext, "."))
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(df)) .fail("The file '", basename(filename), "' has no rows.")
    df
}

# Delimited text: separator sniffed from the header line; a semicolon implies
# decimal commas (European spreadsheet export)
.read_delim <- function(path, ext) {
    first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
    counts <- vapply(c(`\t` = "\t", `;` = ";", `,` = ","),
                     function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))), 0)
    sep <- if (ext == "tsv") "\t" else names(counts)[which.max(counts)]
    utils::read.table(path, header = TRUE, sep = sep, dec = if (sep == ";") "," else ".",
                      quote = "\"'", na.strings = c("", "NA", "."), stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = "", strip.white = TRUE,
                      fileEncoding = "UTF-8")
}

# Find a column: exact name, otherwise ignoring case and surrounding spaces
.match_column <- function(name, available) {
    if (name %in% available) return(name)
    hit <- available[tolower(trimws(available)) == tolower(trimws(name))]
    if (length(hit)) hit[1] else NA_character_
}

# Match factor values to the development levels: exact, then ignoring case and
# spaces, then numerically ("1" and "1.0")
.match_levels <- function(x, levels) {
    xs <- trimws(as.character(x))
    out <- ifelse(xs %in% levels, xs, NA_character_)
    miss <- is.na(out) & !is.na(xs)
    if (any(miss)) {
        lo <- match(tolower(xs[miss]), tolower(trimws(levels)))
        out[miss] <- levels[lo]
    }
    miss <- is.na(out) & !is.na(xs)
    if (any(miss)) {
        nx <- suppressWarnings(as.numeric(xs[miss])); nl <- suppressWarnings(as.numeric(levels))
        if (!all(is.na(nl))) out[miss] <- levels[match(nx, nl)]
    }
    out
}

# Harmonise validation data to the spec. Returns the data (outcome + predictors,
# complete cases) and a check table with one row per variable.
# outcomeCols: names of the outcome columns in the spec (family-specific)
cp_harmonise <- function(ext, spec, outcomeCols, setName, complete = TRUE) {
    avail <- names(ext)
    rows <- list()
    out <- data.frame(row.names = seq_len(nrow(ext)))
    add <- function(var, role, status, details)
        rows[[length(rows) + 1]] <<- data.frame(var = var, role = role, status = status,
                                                details = details, stringsAsFactors = FALSE)
    missingCols <- character(0)

    for (v in outcomeCols) {
        col <- .match_column(v, avail)
        if (is.na(col)) { missingCols <- c(missingCols, v); add(v, "Outcome", "Missing", "not in the file"); next }
        x <- ext[[col]]
        note <- if (col != v) paste0("matched to column '", col, "'") else ""
        lev <- spec$outcome$levels
        if (!is.null(lev) && v == spec$outcome$name) {
            m <- .match_levels(x, lev)
            bad <- sum(is.na(m) & !is.na(x))
            unk <- unique(trimws(as.character(x[is.na(m) & !is.na(x)])))
            if (!(spec$outcome$event %in% m))
                add(v, "Outcome", "Problem", paste0("event level '", spec$outcome$event, "' never occurs",
                                                    if (bad) paste0("; unknown values: ", paste(utils::head(unk, 5), collapse = ", ")) else ""))
            else add(v, "Outcome", if (bad) "Warning" else "OK",
                     paste(c(note, if (bad) paste0(bad, " rows with unknown values (", paste(utils::head(unk, 5), collapse = ", "), ") set to missing")), collapse = "; "))
            out[[v]] <- factor(m, levels = lev)
        } else {
            # decimal commas, as for the predictors
            xn <- if (is.numeric(x)) x else suppressWarnings(as.numeric(sub(",", ".", as.character(x), fixed = TRUE)))
            bad <- sum(is.na(xn) & !is.na(x))
            add(v, "Outcome", if (bad) "Warning" else "OK",
                paste(c(note, if (bad) paste0(bad, " non-numeric values set to missing")), collapse = "; "))
            out[[v]] <- xn
        }
    }

    for (info in spec$vars) {
        v <- info$name
        col <- .match_column(v, avail)
        if (is.na(col)) { missingCols <- c(missingCols, v); add(v, "Predictor", "Missing", "not in the file"); next }
        x <- ext[[col]]
        notes <- if (col != v) paste0("matched to column '", col, "'") else character(0)
        status <- "OK"
        if (info$role == "factor") {
            m <- .match_levels(x, info$levels)
            bad <- is.na(m) & !is.na(x)
            if (any(bad)) {
                unk <- unique(trimws(as.character(x[bad])))
                status <- "Warning"
                notes <- c(notes, paste0(sum(bad), " rows with levels not in development (",
                                         paste(utils::head(unk, 5), collapse = ", "), ") set to missing"))
            }
            absent <- setdiff(info$levels, m)
            if (length(absent)) notes <- c(notes, paste0("levels not observed: ", paste(absent, collapse = ", ")))
            out[[v]] <- factor(m, levels = info$levels)
        } else {
            xn <- if (is.numeric(x)) x else suppressWarnings(as.numeric(sub(",", ".", as.character(x), fixed = TRUE)))
            bad <- is.na(xn) & !is.na(x)
            if (any(bad)) {
                status <- "Warning"
                notes <- c(notes, paste0(sum(bad), " non-numeric values set to missing"))
            }
            outside <- mean(xn < info$range[1] | xn > info$range[2], na.rm = TRUE)
            if (is.finite(outside) && outside > 0) {
                if (outside > 0.05) status <- "Warning"
                notes <- c(notes, paste0(.pct(outside), " outside the development range [",
                                         .num(info$range[1]), ", ", .num(info$range[2]), "]"))
            }
            out[[v]] <- xn
        }
        mis <- mean(is.na(out[[v]]))
        if (mis > 0) notes <- c(notes, paste0(.pct(mis), " missing"))
        if (mis > 0.2 && status == "OK") status <- "Warning"
        add(v, "Predictor", status, paste(notes, collapse = "; "))
    }

    check <- do.call(rbind, rows)
    if (length(missingCols))
        .fail(setName, ": ", if (length(missingCols) == 1) "the variable " else "the variables ",
              paste0("'", missingCols, "'", collapse = ", "),
              if (length(missingCols) == 1) " is" else " are",
              " not in the validation data. Rename the columns to the names used in the development data.")
    # complete = FALSE keeps incomplete rows (for a missing-data method to fill)
    keep <- if (complete) stats::complete.cases(out) else rep(TRUE, nrow(out))
    list(data = out[keep, , drop = FALSE], check = check, nMissing = sum(!keep),
         nTotal = nrow(ext), all = out)
}

# Validation rows of the current data defined by a split variable. The data are
# already of the right types (jamovi), but levels must match the spec.
cp_split_data <- function(data, spec, outcomeCols, complete = TRUE) {
    out <- data.frame(row.names = rownames(data))
    for (v in outcomeCols) {
        x <- data[[v]]
        out[[v]] <- if (!is.null(spec$outcome$levels) && v == spec$outcome$name)
            factor(as.character(x), levels = spec$outcome$levels) else as.numeric(x)
    }
    for (info in spec$vars) {
        x <- data[[info$name]]
        out[[info$name]] <- if (info$role == "factor")
            factor(as.character(x), levels = info$levels) else as.numeric(x)
    }
    keep <- if (complete) stats::complete.cases(out) else rep(TRUE, nrow(out))
    list(data = out[keep, , drop = FALSE], nMissing = sum(!keep), nTotal = nrow(data),
         rows = rownames(data)[keep])
}

# Case-mix: per predictor, mean (SD) or level percentages in each data set
cp_casemix <- function(spec, sets) {
    rows <- list()
    for (info in spec$vars) {
        if (info$role == "factor") {
            for (l in info$levels) {
                vals <- vapply(sets, function(d) mean(d[[info$name]] == l), 0)
                rows[[length(rows) + 1]] <- list(var = info$name, level = l, type = "pct", values = vals)
            }
        } else {
            m <- vapply(sets, function(d) mean(d[[info$name]]), 0)
            s <- vapply(sets, function(d) stats::sd(d[[info$name]]), 0)
            rows[[length(rows) + 1]] <- list(var = info$name, level = "", type = "mean", values = m, sd = s)
        }
    }
    rows
}

# Relatedness of the validation data to the development data (Debray et al.
# 2015): c-statistic of a logistic "membership" model of validation vs
# development on the predictors and the outcome; 0.5 = same case-mix.
cp_membership_c <- function(Xdev, ydev, Xval, yval) {
    X <- rbind(cbind(Xdev, .y = ydev), cbind(Xval, .y = yval))
    g <- c(rep(0L, nrow(Xdev)), rep(1L, nrow(Xval)))
    keep <- apply(X, 2, function(col) stats::sd(col) > 0)
    fit <- suppressWarnings(stats::glm.fit(cbind(1, X[, keep, drop = FALSE]), g, family = stats::binomial()))
    .auc_delong(g, fit$fitted.values)$auc
}
