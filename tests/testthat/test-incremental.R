inc_data <- function() {
    f <- test_path("..", "..", "data", "heart_cohort.omv")
    skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f); d[stats::complete.cases(d), ]
}

test_that("Delta AUC equals pROC's paired DeLong test; NRI and IDI equal Hmisc::improveProb", {
    skip_if_not_installed("pROC")
    d <- inc_data(); y <- as.integer(d$event10 == "yes")
    p1 <- fitted(glm(y ~ age + sex, binomial, d)); p2 <- fitted(glm(y ~ age + sex + smoking + sbp + chol + diabetes, binomial, d))
    da <- cp_delta_auc(y, p1, p2)
    rt <- pROC::roc.test(pROC::roc(y, p2, quiet = TRUE, direction = "<"), pROC::roc(y, p1, quiet = TRUE, direction = "<"),
                         method = "delong", paired = TRUE)
    expect_equal(da[["p"]], rt$p.value, tolerance = 1e-8)
    expect_equal(da[["est"]], unname(diff(rev(rt$estimate))), tolerance = 1e-10)
    ip <- Hmisc::improveProb(p1, p2, y)
    ni <- cp_nri_idi(y, p1, p2)
    expect_equal(ni$nri[["est"]], unname(ip$nri), tolerance = 1e-10)
    expect_equal(ni$idi[["est"]], unname(ip$idi), tolerance = 1e-10)
    expect_equal((ni$nri[["upper"]] - ni$nri[["est"]]) / stats::qnorm(0.975), unname(ip$se.nri), tolerance = 1e-6)
    expect_equal((ni$idi[["upper"]] - ni$idi[["est"]]) / stats::qnorm(0.975), unname(ip$se.idi), tolerance = 1e-6)
})

test_that("Delta C of two Cox predictors uses survival's joint concordance variance", {
    d <- inc_data(); y <- cbind(d$time, as.integer(d$status == "event"))
    S <- survival::Surv(y[, 1], y[, 2])
    l1 <- predict(survival::coxph(S ~ age + sex, d)); l2 <- predict(survival::coxph(S ~ age + sex + smoking + sbp, d))
    dc <- cp_delta_c(y, l1, l2)
    c1 <- survival::concordance(S ~ l1, reverse = TRUE)$concordance; c2 <- survival::concordance(S ~ l2, reverse = TRUE)$concordance
    expect_equal(dc[["est"]], unname(c2 - c1), tolerance = 1e-10)
    expect_true(dc[["lower"]] < dc[["est"]] && dc[["est"]] < dc[["upper"]])
})

test_that("the analyses report incremental value: LR test as anova, validation of both models, rendering", {
    f <- test_path("..", "..", "data", "heart_cohort.omv"); skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f)
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    base <- list(list(var = "age", forced = TRUE), list(var = "sex", forced = TRUE))
    opts <- cpLogisticOptions$new(dep = "event10", eventLevel = "yes", covs = c("age", "sbp", "chol"), factors = c("sex", "smoking"),
                                  baseline = base, incremental = TRUE, validation = "boot", nBoot = 10, dca = TRUE,
                                  splitVar = "period", splitLevel = "2017-2019")
    an <- cpLogisticClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1); an$init(); an$run()
    r <- an$results
    expect_null(r$get("problem"))
    dev <- d[d$period == "2012-2016", ]; dev <- dev[stats::complete.cases(dev[, c("event10", "age", "sbp", "chol", "sex", "smoking")]), ]
    g0 <- glm(event10 ~ age + sex, binomial, dev); g1 <- glm(event10 ~ age + sbp + chol + sex + smoking, binomial, dev)
    it <- r$incTable$asDF; keys <- unlist(r$incTable$rowKeys)
    expect_equal(it$delta[keys == "dev lr"], anova(g0, g1)$Deviance[2], tolerance = 1e-6)
    expect_equal(it$delta[keys == "dev auc"], cp_delta_auc(as.integer(dev$event10 == "yes"), fitted(g0), fitted(g1))[["est"]], tolerance = 1e-8)
    X0 <- cbind(age = dev$age, `sex: male` = as.numeric(dev$sex == "male"))
    iv0 <- cp_internal(X0, as.integer(dev$event10 == "yes"), cp_family("binary"), "boot", B = 10, seed = 1234, stability = FALSE)
    expect_equal(it$baseCorr[keys == "dev auc"], iv0$corrected[["auc"]], tolerance = 1e-8)
    expect_equal(heading_problems(r), character(0))
    grDevices::png(file.path(tempdir(), "incdca.png"))
    expect_true(an$.__enclos_env__$private$.plotIncDca(r$plots$incDca, ggtheme = ggplot2::theme_grey(), theme = NULL))
    grDevices::dev.off()
    rc <- cpCox(data = d, time = "time", status = "status", eventLevel = "event", horizon = 5, covs = c("age", "sbp"),
                factors = c("sex", "smoking"), baseline = base, incremental = TRUE)
    expect_null(rc$get("problem"))
    expect_true(is.finite(rc$incTable$asDF$lower[1]))
    rl <- cpLinear(data = d, dep = "sbp1y", covs = c("age", "sbp", "bmi"), factors = "sex", baseline = base, incremental = TRUE)
    expect_null(rl$get("problem"))
    m0 <- lm(sbp1y ~ age + sex, d[stats::complete.cases(d[, c("sbp1y", "age", "sbp", "bmi", "sex")]), ])
    m1 <- lm(sbp1y ~ age + sbp + bmi + sex, d[stats::complete.cases(d[, c("sbp1y", "age", "sbp", "bmi", "sex")]), ])
    kl <- unlist(rl$incTable$rowKeys)
    expect_equal(rl$incTable$asDF$delta[kl == "dev lr"], anova(m0, m1)$F[2], tolerance = 1e-6)
    rn <- cpLogistic(data = d, dep = "event10", eventLevel = "yes", covs = "sbp", incremental = TRUE)
    expect_match(rn$incTable$notes$need$note, "baseline covariates")
})
