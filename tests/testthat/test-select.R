heart_dev <- function() {
    f <- test_path("..", "..", "data", "heart_cohort.omv")
    skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f)
    d$noise1 <- stats::rnorm(nrow(d)); d$noise2 <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
    d[d$period == "2012-2016" & stats::complete.cases(d), ]
}

test_that("backward and forward selection by AIC equal MASS::stepAIC with a forced lower scope", {
    set.seed(4)
    d <- heart_dev()
    prep <- cp_binary_outcome(cp_prepare(d, "event10", covs = c("sbp", "chol", "bmi", "noise1"),
                                         factors = c("sex", "smoking", "diabetes", "noise2", "centre")), "yes")
    fam <- cp_family("binary"); X <- cp_design(prep$spec, prep$data)
    forced <- c("sex", "centre")
    full <- glm(event10 ~ sbp + chol + bmi + noise1 + sex + smoking + diabetes + noise2 + centre, binomial, d)
    st <- MASS::stepAIC(full, scope = list(lower = ~ sex + centre), direction = "backward", trace = 0)
    sb <- cp_select_main(X, prep$y, prep$spec, fam, "backward", "aic", forced = forced)
    expect_setequal(sb$kept, attr(terms(st), "term.labels"))
    null <- glm(event10 ~ sex + centre, binomial, d)
    sf <- MASS::stepAIC(null, scope = list(lower = ~ sex + centre, upper = formula(full)), direction = "forward", trace = 0)
    ff <- cp_select_main(X, prep$y, prep$spec, fam, "forward", "aic", forced = forced)
    expect_setequal(ff$kept, attr(terms(sf), "term.labels"))
    expect_true(all(forced %in% sb$kept) && all(forced %in% ff$kept))
    sbic <- MASS::stepAIC(full, scope = list(lower = ~ sex + centre), direction = "backward", trace = 0, k = log(nrow(d)))
    expect_setequal(cp_select_main(X, prep$y, prep$spec, fam, "backward", "bic", forced = forced)$kept,
                    attr(terms(sbic), "term.labels"))
})

test_that("the LR test of each term is against the final model, and selection runs inside the bootstrap", {
    set.seed(5)
    d <- heart_dev()
    prep <- cp_binary_outcome(cp_prepare(d, "event10", covs = c("sbp", "chol", "noise1"), splineVars = "age",
                                         factors = c("sex", "noise2"), interactions = list(c("age", "sex"), c("noise1", "sex"))), "yes")
    fam <- cp_family("binary"); X <- cp_design(prep$spec, prep$data)
    sel <- list(mode = "backward", crit = "lrt", alpha = 0.05, forced = "sex", intMode = "lrt", intAlpha = 0.05)
    b <- cp_fit_procedure(X, prep$y, prep$spec, fam, sel)
    expect_true(all(c("sbp", "chol", "age", "sex") %in% attr(b, "vars")))
    tt <- attr(b, "mainTests")
    final <- glm(event10 ~ sbp + chol + rms::rcs(age, prep$spec$vars$age$knots) + sex, binomial, d)
    for (v in tt$var[!tt$kept]) {
        ext <- update(final, as.formula(paste(". ~ . +", v)))
        expect_equal(tt$chisq[tt$var == v], 2 * (as.numeric(logLik(ext)) - as.numeric(logLik(final))), tolerance = 1e-6)
    }
    spec <- cp_with_terms(prep$spec, attr(b, "vars"), attr(b, "ints"))
    expect_equal(unname(b[c("(Intercept)", spec$columns)]), unname(fam$fit(cp_design(spec, prep$data), prep$y)), tolerance = 1e-6, ignore_attr = TRUE)
    fit <- function(X, y) cp_fit_procedure(X, y, prep$spec, fam, sel)
    f2 <- fam; f2$fit <- fit
    iv <- cp_internal(X, prep$y, f2, "boot", B = 5, seed = 1, stability = FALSE)
    expect_true(all(iv$selection[c("sex")] == 1))
})

test_that("the analyses select predictors with baseline covariates and report it", {
    set.seed(6)
    d <- jmvReadWrite::read_omv(test_path("..", "..", "data", "heart_cohort.omv"))
    d$noise <- stats::rnorm(nrow(d))
    base <- list(list(var = "sex", forced = TRUE), list(var = "centre", forced = TRUE))
    r <- cpLogistic(data = d, dep = "event10", eventLevel = "yes", covs = c("sbp", "chol", "noise"), splineVars = "age",
                    factors = c("sex", "centre", "smoking"), baseline = base, selMode = "backward", selCrit = "lrt",
                    interactions = list(list("noise", "sex"), list("age", "sex")), validation = "boot", nBoot = 10,
                    shapley = TRUE, permImp = TRUE, nPerm = 2, shap = TRUE, sampleSize = TRUE)
    expect_null(r$get("problem"))
    st <- r$selTable$asDF
    expect_equal(st$status[st$var == "sex"], "baseline")
    expect_equal(st$status[st$var == "noise"], "removed")
    expect_true(all(st$freq[st$var %in% c("sex", "centre")] == 1))
    orr <- r$orTable$asDF
    expect_true(is.na(orr$or[orr$var == "noise"]) && is.finite(orr$or_u[orr$var == "noise"]))
    expect_match(r$intTable$asDF$kept[r$intTable$asDF$term == "noise × sex"], "main effect removed")
    expect_equal(r$shapleyTable$asDF$shapley[r$shapleyTable$asDF$var == "noise"], 0)
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    expect_equal(heading_problems(r), character(0))
    rc <- cpCox(data = d, time = "time", status = "status", eventLevel = "event", horizon = 5,
                covs = c("sbp", "noise"), factors = c("sex", "smoking"), baseline = list(list(var = "sex", forced = TRUE)),
                selMode = "forward", selCrit = "bic")
    expect_null(rc$get("problem"))
    expect_equal(rc$selTable$asDF$status[rc$selTable$asDF$var == "noise"], "not added")
    rl <- cpLinear(data = d, dep = "sbp1y", covs = c("sbp", "bmi", "noise"), factors = "smoking", selMode = "backward")
    expect_null(rl$get("problem"))
})
