test_that("exploration: shares missing by outcome group with a chi-square test, patterns", {
    d <- data.frame(y = rep(0:1, each = 50), a = c(rep(NA, 10), rnorm(90)), b = c(rnorm(40), rep(NA, 30), rnorm(30)))
    ex <- cp_miss_explore(d, c("a", "b"), d$y)
    expect_equal(ex$n, c(10, 30))
    expect_equal(ex$a[2], mean(is.na(d$b[d$y == 0])))
    expect_equal(ex$p[2], suppressWarnings(chisq.test(table(d$y, is.na(d$b)), correct = FALSE)$p.value))
    pt <- cp_miss_patterns(d, c("a", "b"))
    expect_equal(pt$n[pt$pattern == "(complete)"], 60)
    expect_equal(sum(pt$n), 100)
    ec <- cp_miss_explore(d, "a", d$y + rnorm(100), numericOutcome = TRUE)
    expect_true(is.finite(ec$p))
})

test_that("missing category: a '(missing)' level, an indicator and the development median", {
    d <- data.frame(x = c(1, NA, 3, 5), f = factor(c("a", NA, "b", "a")))
    h <- cp_miss_category(d, "x", "f")
    out <- cp_miss_category_apply(h, d)
    expect_equal(out$x, c(1, 3, 3, 5))
    expect_equal(out[["x (missing)"]], c(0, 1, 0, 0))
    expect_equal(levels(out$f), c("a", "b", "(missing)"))
    v <- cp_miss_category_apply(h, data.frame(x = c(NA, 2), f = factor(c("b", NA))))
    expect_equal(v$x, c(3, 2)); expect_equal(as.character(v$f), c("b", "(missing)"))
})

test_that("chained equations keep observed values, impute from them, are reproducible and recover a MAR association", {
    set.seed(10)
    n <- 1500
    x1 <- rnorm(n); f <- factor(ifelse(runif(n) < plogis(x1), "b", "a")); g <- factor(sample(c("u", "v", "w"), n, TRUE))
    y <- 1 + 0.8 * x1 + 0.5 * (f == "b") + rnorm(n)
    d <- data.frame(y, x1, f, g)
    d$x1[runif(n) < plogis(-1 + y)] <- NA          # MAR on y
    d$f[runif(n) < 0.15] <- NA; d$g[runif(n) < 0.1] <- NA
    imps <- cp_impute(d, c("x1", "f", "g"), m = 5, maxit = 5, seed = 3)
    expect_length(imps, 5)
    for (im in imps) {
        expect_false(anyNA(im))
        expect_equal(im$x1[!is.na(d$x1)], d$x1[!is.na(d$x1)])
        expect_true(all(im$x1[is.na(d$x1)] %in% d$x1[!is.na(d$x1)]))
        expect_setequal(levels(im$g), c("u", "v", "w"))
    }
    expect_identical(cp_impute(d, c("x1", "f", "g"), m = 2, maxit = 2, seed = 3)[[1]],
                     cp_impute(d, c("x1", "f", "g"), m = 2, maxit = 2, seed = 3)[[1]])
    fits <- lapply(imps, function(im) lm(y ~ x1 + f, im))
    pool <- cp_rubin(lapply(fits, coef), lapply(fits, vcov))
    skip_if_not_installed("mice")
    mi <- mice::mice(d, m = 5, maxit = 5, printFlag = FALSE, seed = 3)
    est <- summary(mice::pool(with(mi, lm(y ~ x1 + f))))
    # PMM by mice and here agree (with 56 % of x1 missing both depart from 0.8 alike)
    expect_equal(unname(pool$est), est$estimate, tolerance = 0.1)
    expect_equal(unname(sqrt(diag(pool$vcov))), est$std.error, tolerance = 0.3)
})

test_that("Rubin's rules and the D1 test equal mice on the same fits", {
    skip_if_not_installed("mice")
    set.seed(11)
    d <- mice::nhanes2
    imp <- mice::mice(d, m = 5, printFlag = FALSE, seed = 1)
    f1 <- with(imp, lm(chl ~ age + bmi + hyp)); f0 <- with(imp, lm(chl ~ age))
    fits1 <- mice::getfit(f1); fits0 <- mice::getfit(f0)
    pool <- cp_rubin(lapply(fits1, coef), lapply(fits1, vcov))
    mp <- mice::pool(f1)$pooled
    expect_equal(unname(pool$est), mp$estimate, tolerance = 1e-10)
    expect_equal(unname(diag(pool$vcov)), mp$t, tolerance = 1e-10)
    d1 <- cp_d1(pool, grep("^(bmi|hyp)", names(pool$est)))
    # mitml's D1 without a complete-data df (Li, Raghunathan & Rubin 1991)
    md <- mitml::testModels(mice::as.mitml.result(f1), mice::as.mitml.result(f0), method = "D1")$test
    expect_equal(d1[["stat"]], md[1, "F.value"], tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(d1[["df1"]], md[1, "df1"], ignore_attr = TRUE)
    expect_equal(d1[["df2"]], md[1, "df2"], tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(d1[["p"]], md[1, "P(>F)"], tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("the analyses handle missing values: category, single and multiple imputation", {
    f <- test_path("..", "..", "data", "heart_cohort.omv"); skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f)
    ext <- list(path = test_path("..", "..", "examples", "heart_external.csv"), filename = "ext.csv")
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    args <- list(data = d, dep = "event10", eventLevel = "yes", covs = c("sbp", "chol", "bmi"), factors = c("sex", "smoking"),
                 splitVar = "period", splitLevel = "2017-2019", extFile = ext, missExplore = TRUE, coefs = TRUE, predTests = TRUE)
    r <- do.call(cpLogistic, c(args, list(missMethod = "category")))
    expect_null(r$get("problem"))
    expect_true("bmi (missing)" %in% r$coefTable$asDF$term)
    expect_equal(r$samples$asDF$n, c(1600, 800, 988))
    expect_equal(heading_problems(r), character(0))
    r <- do.call(cpLogistic, c(args, list(missMethod = "single", updating = TRUE)))
    expect_null(r$get("problem"))
    expect_equal(r$samples$asDF$n, c(1600, 800, 1000))
    expect_equal(nrow(r$updTable$asDF), 8)
    # multiple imputation: pooled coefficients are Rubin's rules over glm fits of the same imputations
    opts <- do.call(cpLogisticOptions$new, c(args[-1], list(missMethod = "mi", miM = 3, miIter = 3)))
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    rm <- an$results
    expect_null(rm$get("problem"))
    imps <- an$.__enclos_env__$private$.impSet(d[d$period == "2012-2016", ], 3, 1234)
    fits <- lapply(imps, function(im) glm(event10 ~ sbp + chol + bmi + sex + smoking, binomial, im))
    pool <- cp_rubin(lapply(fits, coef), lapply(fits, vcov))
    expect_equal(rm$coefTable$asDF$est, unname(pool$est), tolerance = 1e-6)
    expect_equal(rm$coefTable$asDF$se, unname(sqrt(diag(pool$vcov))), tolerance = 1e-5)
    expect_equal(heading_problems(rm), character(0))
    jmvcore:::initProtoBuf()
    expect_gt(length(RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)), 1000)
    rl <- cpLinear(data = d, dep = "sbp1y", covs = c("sbp", "bmi"), factors = "smoking", missMethod = "mi", miM = 2, miIter = 2,
                   validation = "cv", nFolds = 3, nRepeats = 1, shapley = TRUE)
    expect_null(rl$get("problem"))
    rc <- cpCox(data = d, time = "time", status = "status", eventLevel = "event", horizon = 5, covs = c("chol", "bmi"),
                factors = "sex", missMethod = "mi", miM = 2, miIter = 2, kmPlot = TRUE, dca = TRUE)
    expect_null(rc$get("problem"))
    expect_equal(heading_problems(rc), character(0))
})
