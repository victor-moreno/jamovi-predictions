test_that("effects table: baseline covariates last below a line, adjusted = full model, selected = reduced model", {
    f <- test_path("..", "..", "data", "heart_cohort.omv"); skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f); set.seed(3); d$noise <- stats::rnorm(nrow(d))
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    opts <- cpLogisticOptions$new(dep = "event10", eventLevel = "yes", covs = c("sbp", "chol", "noise"), factors = c("sex", "smoking", "centre"),
                                  baseline = list(list(var = "sex", forced = TRUE), list(var = "centre", forced = TRUE)),
                                  selMode = "backward", selCrit = "bic", orCI = FALSE, orP = TRUE, forest = TRUE, forestCrude = TRUE)
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1); an$init(); an$run()
    t <- an$results$orTable
    keys <- unlist(t$rowKeys)
    expect_equal(keys[1:3], c("sbp", "chol", "noise"))
    expect_equal(tail(keys, 6), c(".baseline", "sex (ref)", "sex: male", "centre (ref)", "centre: B", "centre: C"))
    expect_equal(t$getCell(col = 1, rowKey = ".baseline")$value, "Baseline covariates")
    expect_true(match("sex (ref)", keys) > match("smoking: current", keys))
    vis <- vapply(t$columns, function(cl) cl$visible, TRUE); names(vis) <- vapply(t$columns, function(cl) cl$name, "")
    expect_false(any(vis[c("lower", "upper", "lower_u", "upper_s")]))
    expect_true(all(vis[c("or", "p", "or_s", "p_s", "or_u")]))
    dd <- d[stats::complete.cases(d[, c("event10", "sbp", "chol", "noise", "sex", "smoking", "centre")]), ]
    full <- glm(event10 ~ sbp + chol + noise + sex + smoking + centre, binomial, dd)
    df <- t$asDF
    expect_equal(df$or[keys == "noise"], unname(exp(coef(full)[["noise"]])), tolerance = 1e-8)
    expect_false(any(c("sex", "centre") %in% an$results$selTable$asDF$var))
    kept <- c(an$results$selTable$asDF$var[an$results$selTable$asDF$status == "kept"], "sex", "centre")
    red <- glm(stats::reformulate(kept, "event10"), binomial, dd)
    expect_equal(df$or_s[keys == "smoking: current"], unname(exp(coef(red)[["smokingcurrent"]])), tolerance = 1e-8)
    expect_true(is.na(df$or_s[keys == "noise"]))
    expect_equal(heading_problems(an$results), character(0))
    st <- an$results$plots$forestPlot$state
    expect_setequal(unique(st$rows$block), c("Crude", "Adjusted (all predictors)", "Selected model"))
    expect_true(all(st$rows$forced[st$rows$label == "centre: B (vs A)"]))
    grDevices::png(file.path(tempdir(), "forest.png"))
    expect_true(an$.__enclos_env__$private$.plotForest(an$results$plots$forestPlot, ggtheme = ggplot2::theme_grey(), theme = NULL))
    grDevices::dev.off()
    rl <- cpLinear(data = d, dep = "sbp1y", covs = c("sbp", "bmi"), factors = "sex", forest = TRUE)
    expect_null(rl$get("problem"))
    expect_false(rl$plots$forestPlot$state$exp)
    expect_false(any(vapply(rl$orTable$columns, function(cl) cl$name %in% c("or_s", "p_s") && cl$visible, TRUE)))
    rc <- cpCox(data = d, time = "time", status = "status", eventLevel = "event", horizon = 5, covs = c("sbp", "noise"),
                factors = "sex", baseline = list(list(var = "sex", forced = TRUE)), selMode = "enet", enetFolds = 5, forest = TRUE)
    expect_null(rc$get("problem"))
    vis <- vapply(rc$orTable$columns, function(cl) cl$visible, TRUE); names(vis) <- vapply(rc$orTable$columns, function(cl) cl$name, "")
    expect_true(vis[["or_s"]] && !vis[["lower_s"]] && !vis[["p_s"]])
})
