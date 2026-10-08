# jamovi draws the leading text columns of a table (title '' or combineBelow)
# as row headings, and spans a heading cell over blank cells below it; in a
# column without combineBelow it still draws those blank cells, which shifts
# the row one column to the right (client/resultsview/table.ts). So such
# columns must have a value in every row.
heading_problems <- function(results) {
    out <- character(0)
    walk <- function(g) for (it in g$items) {
        if (inherits(it, "Group")) { walk(it); next }
        if (!inherits(it, "Table") || !it$visible || it$rowCount == 0) next
        cols <- Filter(function(cl) cl$visible, it$columns)
        df <- it$asDF
        for (cl in cols) {
            heading <- (cl$title == "" && cl$type == "text") || isTRUE(cl$combineBelow)
            if (!heading) break
            if (isTRUE(cl$combineBelow)) next
            v <- df[[cl$name]]
            bad <- is.na(v) | v == ""
            if (any(bad)) out <- c(out, paste0(it$name, "$", cl$name, " rows ", paste(which(bad), collapse = ",")))
        }
    }
    walk(results)
    out
}

test_that("no table has a blank row-heading cell (logistic, linear, Cox)", {
    d <- pima_stacked()
    r <- cpLogistic(data = d, dep = "type", eventLevel = "Yes", covs = c("bp", "bmi"), splineVars = "glu",
                    factors = "agegrp", splitVar = "sample", splitLevel = "test", descTable = TRUE,
                    coefs = TRUE, predTests = TRUE, classif = TRUE, dca = TRUE, updating = TRUE,
                    interactions = list(list("bmi", "agegrp")), validation = "boot", nBoot = 10,
                    stability = TRUE, sampleSize = TRUE, shapley = TRUE, permImp = TRUE, nPerm = 2, shap = TRUE)
    expect_equal(heading_problems(r), character(0))
    r <- cpLinear(data = d, dep = "glu", covs = c("bp", "bmi"), factors = c("type", "agegrp"),
                  splitVar = "sample", splitLevel = "test", descTable = TRUE, coefs = TRUE, predTests = TRUE,
                  updating = TRUE, sampleSize = TRUE, shapley = TRUE, permImp = TRUE, nPerm = 2, shap = TRUE)
    expect_equal(heading_problems(r), character(0))
    g <- survival::gbsg
    g$status <- factor(g$status, 0:1, c("censored", "event"))
    g$grade <- factor(g$grade)
    g$half <- factor(rep(c("a", "b"), length.out = nrow(g)))
    r <- cpCox(data = g, time = "rfstime", status = "status", eventLevel = "event", horizon = 1826,
               covs = c("age", "nodes"), factors = "grade", splitVar = "half", splitLevel = "b",
               descTable = TRUE, coefs = TRUE, predTests = TRUE, dca = TRUE, sampleSize = TRUE,
               shapley = TRUE, permImp = TRUE, nPerm = 2, shap = TRUE)
    expect_equal(heading_problems(r), character(0))
})

test_that("factor rows exist after .init with jamovi's 0-row data header", {
    d <- pima_stacked()
    opts <- cpLogisticOptions$new(dep = "type", covs = "glu", factors = "agegrp", splitVar = "sample",
                                  splitLevel = "test", descTable = TRUE)
    an <- cpLogisticClass$new(options = opts, data = d[0, ], datasetId = "d", analysisId = "1", revision = 1)
    an$init()
    expect_true(all(c("agegrp <30", "agegrp 30-45", "agegrp 45+") %in% an$results$caseMix$rowKeys))
    expect_true(all(c("agegrp: 30-45", "agegrp: 45+") %in% an$results$orTable$rowKeys))
    expect_true("agegrp 45+" %in% an$results$descTable$rowKeys)
})
