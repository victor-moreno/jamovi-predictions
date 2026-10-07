full_options <- function(...) list(
    dep = "type", eventLevel = "Yes", covs = c("npreg", "bp", "bmi", "ped"), splineVars = "glu",
    factors = "agegrp", splitVar = "sample", splitLevel = "test", validation = "boot", nBoot = 20,
    stability = TRUE, equation = TRUE, sampleSize = TRUE, classif = TRUE, dca = TRUE, shapley = TRUE,
    permImp = TRUE, nPerm = 5, shap = TRUE, shapRow = 3, rocPlot = TRUE, calPlot = TRUE,
    riskPlot = TRUE, explain = TRUE, refLevels = list(list(var = "agegrp", ref = "30-45")),
    covScales = list(list(var = "bmi", scale = "ten")), interactions = list(list("bmi", "agegrp")),
    intMode = "aic", descTable = TRUE, predTests = TRUE, coefs = TRUE, ...)

test_that("the analysis runs with every option and matches the engine", {
    d <- pima_stacked()
    r <- do.call(cpLogistic, c(list(data = d), full_options()))
    perf <- r$perfDev$asDF
    prep <- pima_prep(factors = "agegrp", data = d[d$sample == "train", ])
    X <- cp_design(prep$spec, prep$data)
    m <- cp_binary_metrics(prep$y, cp_lp(cp_binary_fit(X, prep$y), X))
    expect_equal(perf$apparent[1], m[["auc"]])
    expect_true(all(is.finite(perf$corrected)))
    expect_equal(r$samples$asDF$n, c(200, 332))
    expect_equal(nrow(r$perfVal$asDF), 13)
    expect_equal(r$shapleyTable$asDF$var, c("npreg", "bp", "bmi", "ped", "glu", "agegrp"))
})

test_that("the analysis serialises as jamovi does, with state under the size limit", {
    d <- pima_stacked()
    opts <- do.call(cpLogisticOptions$new, full_options())
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    jmvcore:::initProtoBuf()
    pb <- RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)
    expect_gt(length(pb), 1000)
    st <- an$results$perfDev$state
    expect_lt(length(memCompress(serialize(st, NULL), "gzip")), 500000)
})

test_that("problems reach the user as the module's own message", {
    d <- pima_stacked()
    r <- cpLogistic(data = d, dep = "type", covs = c("glu", "glu"))
    expect_match(r$get("problem")$content, "more than once")
    d$const <- 1
    r <- cpLogistic(data = d, dep = "type", covs = c("glu", "const"))
    expect_match(r$get("problem")$content, "does not vary")
    d$sep <- as.numeric(d$type == "Yes")
    r <- cpLogistic(data = d, dep = "type", covs = c("sep"))
    expect_match(r$get("problem")$content, "separation")
})

test_that("saved predictions cover the validation rows too", {
    d <- pima_stacked()
    opts <- cpLogisticOptions$new(dep = "type", covs = c("glu", "bmi"), splitVar = "sample", splitLevel = "test")
    opts <- enable_output(opts, "predRisk")
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    v <- unlist(an$results$predRisk$.__enclos_env__$private$.values[[1]])
    expect_equal(length(v), 532)
    fit <- glm(type ~ glu + bmi, binomial, d[d$sample == "train", ])
    expect_equal(unname(v[400]), unname(predict(fit, d[400, ], type = "response")), tolerance = 1e-8)
})
