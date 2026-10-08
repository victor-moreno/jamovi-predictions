# Shared helpers: the module's own error messages, small formatting utilities.

# A problem the user can fix: raised with a message written for the user
.fail <- function(...) {
    stop(structure(class = c("cp_user_error", "error", "condition"),
                   list(message = paste0(...), call = NULL)))
}

# Runs an analysis: its own messages reach jamovi as written; any other R error
# is replaced by a plain message, so R's wording never reaches the user. The
# message is a Notice and the analysis ends normally (an error from .run would
# carry the R call stack into jamovi's developer view and copied results).
# The message is picked inside tryCatch and shown outside it, because tryCatch
# nests its handlers.
.guard <- function(analysis, expr) {
    # jamovi's request to restart (new options while computing) must reach jamovi
    restart <- NULL
    msg <- tryCatch({ expr; NULL },
        cp_user_error = function(e) conditionMessage(e),
        error = function(e) if (identical(e$code, "restart")) { restart <<- e; NULL } else paste(
            "The model could not be computed with these data and options.",
            "Check that the outcome has enough events, that each factor level",
            "has observations in both outcome groups, and that the predictors vary."))
    if (!is.null(restart)) stop(restart)
    if (!is.null(msg)) .show_problem(analysis, msg)
    invisible(msg)
}

# The problem as an error Notice above the results, which are hidden
.show_problem <- function(analysis, msg) {
    for (item in analysis$results$items) item$setVisible(FALSE)
    notice <- jmvcore::Notice$new(options = analysis$options, name = "problem",
                                  type = jmvcore::NoticeType$ERROR)
    notice$setContent(msg)
    analysis$results$insert(1, notice)
}

# A non-blocking remark (warning) above the results
.show_warning <- function(analysis, name, msg, index = 1) {
    notice <- jmvcore::Notice$new(options = analysis$options, name = name,
                                  type = jmvcore::NoticeType$WARNING)
    notice$setContent(msg)
    analysis$results$insert(index, notice)
}

# Numbers for messages (jmvcore's format() shadows base::format() in the module)
.num <- function(x, digits = 3) base::format(signif(x, digits), scientific = FALSE, trim = TRUE)

.pct <- function(x, digits = 1) paste0(base::format(round(100 * x, digits), nsmall = digits, trim = TRUE), "%")

# Parse a list of numbers typed by the user ("0.05, 0.1; 0.2")
.parse_numbers <- function(s) {
    if (is.null(s) || !nzchar(trimws(s))) return(numeric(0))
    parts <- trimws(strsplit(s, "[,; ]+")[[1]])
    parts <- parts[nzchar(parts)]
    x <- suppressWarnings(as.numeric(sub(",", ".", parts, fixed = TRUE)))
    if (anyNA(x)) .fail("Could not read '", s, "' as a list of numbers separated by commas.")
    x
}

# Run expr with a fixed seed and restore the caller's random stream afterwards
.with_seed <- function(seed, expr) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
        if (is.null(old)) rm(".Random.seed", envir = globalenv())
        else assign(".Random.seed", old, envir = globalenv())
    })
    set.seed(seed)
    expr
}
