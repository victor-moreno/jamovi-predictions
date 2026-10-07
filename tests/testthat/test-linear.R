lin_prep <- function(data = pima_age(), ...) {
    prep <- cp_prepare(data, "glu", covs = c("npreg", "bp", "ped"), splineVars = "bmi",
                       factors = c("type", "agegrp"), refLevels = list(agegrp = "30-45"),
                       scales = list(bp = "ten"), ...)
    cp_gaussian_outcome(prep)
}
pima_age <- function() {
    d <- MASS::Pima.tr
    d$agegrp <- cut(d$age, c(0, 30, 45, 100), labels = c("<30", "30-45", "45+"))
    d
}
lin_lm <- function(prep) {
    d <- pima_age(); d$ag2 <- relevel(d$agegrp, "30-45"); d$bp10 <- d$bp / 10
    k <- prep$spec$vars$bmi$knots
    lm(glu ~ npreg + bp10 + ped + rms::rcs(bmi, k) + type + ag2, d)
}

test_that("coefficients, SEs and F tests reproduce lm and anova", {
    skip_if_not_installed("rms")
    prep <- lin_prep(); fam <- cp_family("gaussian")
    X <- cp_design(prep$spec, prep$data)
    b <- fam$fit(X, prep$y)
    m <- lin_lm(prep)
    expect_equal(as.numeric(b), unname(coef(m)), tolerance = 1e-8)
    expect_equal(sqrt(diag(fam$vcov(b))), unname(sqrt(diag(vcov(m)))), tolerance = 1e-8)
    lrt <- cp_lr_tests(prep$spec, X, prep$y, fam)
    a <- drop1(m, test = "F")
    expect_equal(lrt$chisq[lrt$key == "agegrp"], a["ag2", "F value"], tolerance = 1e-8)
    expect_equal(lrt$p[lrt$key == "type"], a["type", "Pr(>F)"], tolerance = 1e-8)
    expect_equal(lrt$df2[1], m$df.residual)
})

test_that("adjusted and crude effects are lm differences with t intervals", {
    prep <- lin_prep(); fam <- cp_family("gaussian")
    X <- cp_design(prep$spec, prep$data)
    b <- fam$fit(X, prep$y)
    eff <- cp_effects(prep$spec, b, fam$vcov(b), 0.95, exp = FALSE, df = fam$dfResid(b))
    m <- lin_lm(prep)
    ci <- confint(m)
    expect_equal(eff$est[eff$key == "type: Yes"], unname(coef(m)["typeYes"]), tolerance = 1e-8)
    expect_equal(unlist(eff[eff$key == "type: Yes", c("lower", "upper")]), ci["typeYes", ],
                 tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(eff$est[eff$key == "bp"], unname(coef(m)["bp10"]), tolerance = 1e-8)
    cr <- cp_crude(prep$spec, prep$data, prep$y, fam)
    m1 <- lm(glu ~ npreg, pima_age())
    expect_equal(cr$est[cr$key == "npreg"], unname(coef(m1)[2]), tolerance = 1e-8)
    expect_equal(cr$p[cr$key == "npreg"], summary(m1)$coefficients[2, 4], tolerance = 1e-8)
    expect_equal(attr(cr, "pOverall")[["agegrp"]], anova(lm(glu ~ agegrp, pima_age()))[1, 5], tolerance = 1e-8)
})

test_that("metrics: apparent R² is lm's, calibration slope interval is lm's", {
    prep <- lin_prep(); fam <- cp_family("gaussian")
    X <- cp_design(prep$spec, prep$data)
    lp <- fam$lp(fam$fit(X, prep$y), X)
    m <- fam$metrics(prep$y, lp)
    expect_equal(m[["r2"]], summary(lin_lm(prep))$r.squared, tolerance = 1e-10)
    expect_equal(m[["slope"]], 1, tolerance = 1e-10)
    yv <- MASS::Pima.te$glu; lv <- yv * 0.8 + rnorm(length(yv), 20, 10)
    ci <- cp_gaussian_ci(yv, lv)
    expect_equal(ci$slope, unname(confint(lm(yv ~ lv))[2, ]), tolerance = 1e-8)
    expect_equal(cp_gaussian_metrics(yv, lv)[["slope"]], unname(coef(lm(yv ~ lv))[2]), tolerance = 1e-10)
})

test_that("bootstrap optimism matches rms::validate(ols) with the same seed", {
    skip_if_not_installed("rms")
    prep <- cp_gaussian_outcome(cp_prepare(MASS::Pima.tr, "glu", covs = c("npreg", "bp", "bmi", "ped"),
                                           splineVars = "age"))
    X <- cp_design(prep$spec, prep$data)
    g <- rms::ols(prep$y ~ X, x = TRUE, y = TRUE)
    set.seed(11); vv <- rms::validate(g, B = 30)
    iv <- cp_internal(X, prep$y, cp_family("gaussian"), "boot", B = 30, seed = 11, stability = TRUE)
    expect_equal(iv$optimism[["r2"]], vv["R-square", "optimism"], tolerance = 1e-6)
    expect_equal(iv$optimism[["slope"]], vv["Slope", "optimism"], tolerance = 1e-6)
    expect_gt(iv$stability$mape, 0)
    expect_true(all(iv$stability$pts >= 0 & iv$stability$pts <= 1000))
})

test_that("Shapley R² equals relaimpo's LMG", {
    skip_if_not_installed("relaimpo")
    d <- MASS::Pima.tr
    prep <- cp_gaussian_outcome(cp_prepare(d, "glu", covs = c("npreg", "bp", "bmi", "ped", "age")))
    X <- cp_design(prep$spec, prep$data)
    sh <- cp_shapley_r2(X, prep$y, prep$spec, cp_family("gaussian"), "r2")
    ri <- relaimpo::calc.relimp(lm(glu ~ npreg + bp + bmi + ped + age, d), type = "lmg", rela = FALSE)
    expect_equal(sh$shapley, unname(ri@lmg), tolerance = 1e-10)
})

test_that("sample size criteria 1-3 match pmsampsize", {
    skip_if_not_installed("pmsampsize")
    ss <- cp_gaussian_samplesize(10, 0.2, intercept = 120, sd = 30)
    pm <- pmsampsize::pmsampsize(type = "c", rsquared = 0.2, parameters = 10, intercept = 120, sd = 30)
    expect_equal(unname(ss$n[1:3]), unname(pm$results_table[1:3, 1]))
    # criterion 4 solved directly: the upper limit of the mean is within × 1.1
    n4 <- ss$n[["intercept"]]
    expect_lte(120 + qt(0.975, n4 - 11) * sqrt(30^2 * 0.8 / n4), 1.1 * 120)
})

test_that("updating F tests reproduce anova against the revised model", {
    prep <- lin_prep(); fam <- cp_family("gaussian")
    X <- cp_design(prep$spec, prep$data)
    b <- fam$fit(X, prep$y)
    te <- MASS::Pima.te; te$agegrp <- cut(te$age, c(0, 30, 45, 100), labels = c("<30", "30-45", "45+"))
    v <- cp_split_data(te, prep$spec, "glu")$data
    Xv <- cp_design(prep$spec, v); yv <- v$glu
    lp <- fam$lp(b, Xv)
    u <- cp_update(fam, Xv, yv, lp)
    rev <- lm(yv ~ Xv)
    expect_equal(u$methods$original$test[["stat"]], anova(lm(yv ~ 0 + offset(lp)), rev)$F[2], tolerance = 1e-8)
    expect_equal(u$methods$citl$test[["stat"]], anova(lm(yv ~ 1 + offset(lp)), rev)$F[2], tolerance = 1e-8)
    expect_equal(u$methods$recal$test[["stat"]], anova(lm(yv ~ lp), rev)$F[2], tolerance = 1e-8)
})

test_that("the linear analysis runs with every option, saves predictions and serialises", {
    d <- pima_stacked()
    opts <- cpLinearOptions$new(dep = "glu", covs = c("npreg", "bp", "ped"), splineVars = "bmi",
        factors = c("type", "agegrp"), splitVar = "sample", splitLevel = "test", validation = "boot", nBoot = 20,
        stability = TRUE, equation = TRUE, sampleSize = TRUE, shapley = TRUE, permImp = TRUE, nPerm = 3,
        shap = TRUE, shapRow = 3, calPlot = TRUE, explain = TRUE, interactions = list(list("bmi", "type")),
        descTable = TRUE, predTests = TRUE, coefs = TRUE, updating = TRUE)
    opts <- enable_output(opts, "predValue")
    an <- cpLinearClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    r <- an$results
    expect_null(r$get("problem"))
    expect_equal(r$samples$asDF$n, c(200, 332))
    expect_true(all(is.finite(r$perfDev$asDF$corrected)))
    v <- unlist(an$results$predValue$.__enclos_env__$private$.values[[1]])
    tr <- d[d$sample == "train", ]
    k <- cp_prepare(tr, "glu", splineVars = "bmi")$spec$vars$bmi$knots
    m <- lm(glu ~ npreg + bp + ped + rms::rcs(bmi, k) * type + agegrp, tr)
    expect_equal(unname(v[400]), unname(predict(m, d[400, ])), tolerance = 1e-8)
    jmvcore:::initProtoBuf()
    pb <- RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)
    expect_gt(length(pb), 1000)
    render <- c(calibration = ".plotCal", instability = ".plotInstability", calInstability = ".plotCalInstability",
                shapSummary = ".plotShapSummary", shapPatient = ".plotShapPatient", permBar = ".plotPerm",
                shapleyBar = ".plotShapley")
    priv <- an$.__enclos_env__$private
    for (p in names(render)) {
        f <- file.path(tempdir(), paste0(p, ".png"))
        grDevices::png(f)
        ok <- priv[[render[[p]]]](r$plots[[p]], ggtheme = ggplot2::theme_grey(), theme = NULL)
        grDevices::dev.off()
        expect_true(ok, label = p)
    }
})

test_that("an external file validates the linear model", {
    f <- file.path(tempdir(), "pima_te_lin.csv")
    utils::write.csv(MASS::Pima.te, f, row.names = FALSE)
    r <- cpLinear(data = MASS::Pima.tr, dep = "glu", covs = c("bp", "bmi"), factors = "type",
                  extFile = list(path = f, filename = "pima_te.csv"), updating = TRUE, permImp = TRUE, nPerm = 2)
    expect_null(r$get("problem"))
    pv <- r$perfVal$asDF
    m <- lm(glu ~ bp + bmi + type, MASS::Pima.tr)
    pr <- predict(m, MASS::Pima.te); y <- MASS::Pima.te$glu
    expect_equal(pv$est[1], 1 - sum((y - pr)^2) / sum((y - mean(y))^2), tolerance = 1e-8)
    expect_equal(r$extCheck$asDF$status, rep("OK", 4))
    expect_equal(nrow(r$updTable$asDF), 4)
})
