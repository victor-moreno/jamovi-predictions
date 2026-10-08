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
    # removed by the selection: blank in the selected model, estimated in the full and crude models
    expect_true(is.na(orr$or_s[orr$var == "noise"]) && is.finite(orr$or[orr$var == "noise"]) &&
                is.finite(orr$or_u[orr$var == "noise"]))
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

test_that("bootstrap coefficient summaries are computed as Heinze et al. define them", {
    coefs <- cbind(`(Intercept)` = c(1, 1.2, 0.8, 1), a = c(0.5, 0, 0.7, 0.6), b = c(0, 0, 0.2, 0))
    final <- c(`(Intercept)` = 1, a = 0.6, b = 0); se <- c(`(Intercept)` = 0.1, a = 0.2, b = NA)
    bc <- cp_boot_coefs(coefs, final, se)
    expect_equal(bc$incl, c(1, 0.75, 0.25))
    expect_equal(bc$median[2], median(c(0.5, 0, 0.7, 0.6)))
    expect_equal(bc$lower[2], unname(quantile(c(0.5, 0, 0.7, 0.6), 0.025)))
    expect_equal(bc$rmsd[2], sqrt(mean((c(0.5, 0, 0.7, 0.6) - 0.6)^2)) / 0.2)
    expect_equal(bc$relbias[2], mean(c(0.5, 0.7, 0.6)) / 0.6 - 1)
    expect_true(is.na(bc$relbias[3]) && is.na(bc$rmsd[3]))
    mf <- cp_model_freq(c("a, b", "a", "a", NA))
    expect_equal(mf$model[1], "a"); expect_equal(mf$pct[1], 2 / 3)
})

test_that("the relaxed lasso is the unpenalised refit of the kept terms; bootstrap tables and cautions appear", {
    f <- test_path("..", "..", "data", "heart_cohort.omv"); skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f); set.seed(3); d$noise <- stats::rnorm(nrow(d))
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    r <- cpLogistic(data = d, dep = "event10", eventLevel = "yes", covs = c("sbp", "chol", "noise"), factors = c("sex", "smoking"),
                    baseline = list(list(var = "sex", forced = TRUE)), selMode = "enet", enetRelax = TRUE, enetFolds = 5,
                    validation = "boot", nBoot = 10, bootCoefs = TRUE, coefs = TRUE, oddsRatios = TRUE)
    expect_null(r$get("problem"))
    ct <- r$coefTable$asDF
    kept <- r$selTable$asDF$var[r$selTable$asDF$status %in% c("kept", "baseline")]
    dd <- d[stats::complete.cases(d[, c("event10", "sbp", "chol", "noise", "sex", "smoking")]), ]
    g <- glm(stats::reformulate(kept, "event10"), binomial, dd)
    expect_equal(sort(unname(ct$est)), sort(unname(coef(g))), tolerance = 1e-6)
    expect_true(all(is.finite(ct$se)))
    expect_match(r$coefTable$notes$post$note, "too narrow")
    expect_match(r$selTable$notes$relax$note, "Relaxed lasso")
    expect_gt(nrow(r$bootCoefTable$asDF), 3)
    expect_gt(nrow(r$modelFreqTable$asDF), 0)
    expect_equal(heading_problems(r), character(0))
    rn <- cpLogistic(data = d, dep = "event10", eventLevel = "yes", covs = c("sbp", "noise"), selMode = "backward", bootCoefs = TRUE)
    expect_match(rn$bootCoefTable$notes$need$note, "Bootstrap")
    rm <- cpCox(data = d, time = "time", status = "status", eventLevel = "event", horizon = 5, covs = c("sbp", "chol", "noise"),
                factors = "sex", selMode = "backward", missMethod = "mi", miM = 2, miIter = 2, validation = "boot", nBoot = 20,
                bootCoefs = TRUE)
    expect_null(rm$get("problem"))
    expect_false("Intercept" %in% rm$bootCoefTable$asDF$term)
})

test_that("results do not depend on the number of cores; a restart request passes through the guard", {
    skip_on_os("windows")
    f <- test_path("..", "..", "data", "heart_cohort.omv"); skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f)
    run <- function(cores) cpLogistic(data = d, dep = "event10", eventLevel = "yes", covs = c("sbp", "chol"), splineVars = "age",
                                      factors = "smoking", selMode = "enet", enetFolds = 3, validation = "boot", nBoot = 10,
                                      shapley = TRUE, permImp = TRUE, nPerm = 3, bootCoefs = TRUE, cores = cores)
    r1 <- run(1); r3 <- run(3)
    for (t in c("perfDev", "shapleyTable", "permTable", "bootCoefTable", "coefTable"))
        expect_equal(r3[[t]]$asDF, r1[[t]]$asDF, label = t)
    e <- jmvcore::createError("restarting", "restart")
    expect_error(.guard(list(), stop(e)), "restarting")
})
