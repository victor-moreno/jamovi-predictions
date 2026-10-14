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

test_that("missing category: saved predictions and the SHAP waterfall use the indicator coding", {
    d <- pima_stacked()
    d$bmi[c(3, 10, 20)] <- NA
    opts <- cpLogisticOptions$new(dep = "type", covs = c("glu", "bmi"), missMethod = "category", shap = TRUE, shapRow = 3)
    opts <- enable_output(opts, "predRisk")
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    expect_null(an$results$get("problem"))
    v <- unlist(an$results$predRisk$.__enclos_env__$private$.values[[1]])
    expect_false(anyNA(v))
    d$bmiMiss <- as.numeric(is.na(d$bmi))
    d$bmi[is.na(d$bmi)] <- median(d$bmi, na.rm = TRUE)
    fit <- glm(type ~ glu + bmi + bmiMiss, binomial, d)
    expect_equal(unname(v[3]), unname(predict(fit, d[3, ], type = "response")), tolerance = 1e-8)
    expect_false(is.null(an$results$plots$shapPatient$state))
})

test_that("a SHAP row with missing predictors is a warning, not a failed analysis", {
    d <- pima_stacked()
    d$bmi[3] <- NA
    r <- cpLogistic(data = d, dep = "type", covs = c("glu", "bmi"), shap = TRUE, shapRow = 3)
    expect_null(r$get("problem"))
    expect_match(r$get("shapRow")$content, "Row 3")
    expect_true(all(is.finite(r$shapTable$asDF$meanAbs)))
})

test_that("elastic net and multiply imputed models can be exported", {
    d <- MASS::Pima.tr
    d$bmi[c(3, 10)] <- NA
    tmp <- tempfile("export"); dir.create(tmp)
    for (sm in c("enet", "none")) for (mm in c("cc", "mi")) {
        opts <- cpLogisticOptions$new(dep = "type", eventLevel = "Yes", covs = c("npreg", "glu", "bmi"),
                                      selMode = sm, missMethod = mm, miM = 2)
        ea <- opts$option("exportModel"); ea$value <- TRUE
        an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
        unlockBinding(".getSessionTemp", an)
        an$.getSessionTemp <- function() tmp
        an$init(); an$run()
        expect_null(an$results$get("export"))
        expect_length(an$results$get("exportModel"), 1)
    }
    files <- list.files(tmp, "[.]json$", recursive = TRUE, full.names = TRUE)
    expect_gt(length(files), 0)
    for (f in files) expect_equal(jsonlite::fromJSON(f)$format, "clinpred-model")
})
