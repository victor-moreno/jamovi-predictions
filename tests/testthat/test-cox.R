rott <- function() {
    d <- survival::rotterdam
    d$rfs <- pmax(d$recur, d$death)
    d$years <- d$rtime / 365.25
    i <- d$rfs == 1 & d$recur == 0
    d$years[i] <- d$dtime[i] / 365.25
    d <- d[d$nodes > 0, ]
    d$status <- factor(d$rfs, 0:1, c("censored", "event"))
    d$grade <- factor(ifelse(d$grade == 3, "3", "1-2"))
    d$hormon <- factor(d$hormon)
    rownames(d) <- NULL
    d
}
cox_prep <- function(h = 5) {
    prep <- cp_prepare(rott(), c("years", "status"), covs = c("age", "nodes"), splineVars = "pgr",
                       factors = c("size", "grade", "hormon"))
    cp_cox_outcome(prep, "years", "status", "event", h)
}

test_that("fit, partial likelihood, baseline hazard and risks reproduce coxph and survfit", {
    prep <- cox_prep(); X <- cp_design(prep$spec, prep$data); y <- prep$y
    fam <- cp_family("cox", horizon = 5)
    b <- fam$fit(X, y)
    fit <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ X)
    expect_equal(as.numeric(b[-1]), unname(coef(fit)), tolerance = 1e-8)
    expect_equal(sqrt(diag(fam$vcov(b)))[-1], unname(sqrt(diag(vcov(fit)))), tolerance = 1e-8)
    expect_equal(-fam$loss(y, fam$lp(b, X)), fit$loglik[2], tolerance = 1e-8)
    expect_equal(-fam$nullLoss(y), fit$loglik[1], tolerance = 1e-8)
    sf <- survival::survfit(fit, newdata = list(X = X[1:5, , drop = FALSE]))
    expect_equal(fam$risk(fam$lp(b, X[1:5, ])), unname(1 - summary(sf, times = 5)$surv[1, 1:5]), tolerance = 1e-8)
})

test_that("discrimination, D and overall measures reproduce survival and timeROC", {
    prep <- cox_prep(); X <- cp_design(prep$spec, prep$data); y <- prep$y
    fam <- cp_family("cox", horizon = 5)
    lp <- fam$lp(fam$fit(X, y), X)
    m <- fam$metrics(y, lp)
    fit <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ X)
    expect_equal(m[["c"]], unname(survival::concordance(fit)$concordance), tolerance = 1e-10)
    ro <- survival::royston(fit)
    expect_equal(m[["d"]], unname(ro["D"]), tolerance = 1e-6)
    expect_equal(m[["r2d"]], unname(ro["R.D"]), tolerance = 1e-6)
    expect_equal(m[["r2n"]], unname(ro["R.N"]), tolerance = 1e-6)
    # survival::brier divides by the sum of the censoring weights, Graf (and this
    # module, as riskRegression) by n: equal up to that normalisation
    br <- survival::brier(fit, times = 5)
    expect_equal(m[["brier"]], br$brier, tolerance = 1e-4)
    expect_equal(m[["sbrier"]], br$rsquared, tolerance = 1e-4)
    # (left after normalising: survival::brier moves censorings at t just after t)
    w <- .ipcw(y[, 1], y[, 2], 5)
    expect_equal(m[["brier"]] * nrow(y) / sum(w$case + w$ctrl), br$brier, tolerance = 1e-5)
    expect_equal(m[["slope"]], 1, tolerance = 1e-6)
    skip_if_not_installed("timeROC")
    suppressPackageStartupMessages(library(survival))
    tr <- timeROC::timeROC(y[, 1], y[, 2], marker = lp, cause = 1, weighting = "marginal", times = 5, iid = FALSE)
    # timeROC (prodlim) differs in the tie convention of the censoring weights only
    expect_equal(m[["auc"]], unname(tr$AUC[2]), tolerance = 1e-6)
})

test_that("calibration at the horizon is a Cox model on a spline of the cloglog risk", {
    prep <- cox_prep(); X <- cp_design(prep$spec, prep$data); y <- prep$y
    fam <- cp_family("cox", horizon = 5)
    lp <- fam$lp(fam$fit(X, y), X)
    cal <- .cox_calfit(y, lp, 5)
    k <- .rcs_knots(lp, 3, "lp")
    cf <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ rms::rcs(lp, k))
    nd <- data.frame(lp = lp[1:4])
    sf <- survival::survfit(cf, newdata = nd)
    expect_equal(cal(lp[1:4])$fit, unname(1 - summary(sf, times = 5)$surv[1, ]), tolerance = 1e-6)
    r <- fam$risk(lp)
    expect_equal(fam$metrics(y, lp)[["ici"]], mean(abs(cal(lp)$fit - r)), tolerance = 1e-10)
    expect_equal(fam$metrics(y, lp)[["oe"]],
                 (1 - summary(survival::survfit(survival::Surv(y[, 1], y[, 2]) ~ 1), times = 5)$surv) / mean(r),
                 tolerance = 1e-10)
})

test_that("net benefit at the horizon uses Kaplan–Meier above the threshold", {
    prep <- cox_prep(); X <- cp_design(prep$spec, prep$data); y <- prep$y
    fam <- cp_family("cox", horizon = 5)
    p <- fam$risk(fam$lp(fam$fit(X, y), X))
    nb <- cp_net_benefit_surv(y, p, 0.6, 5)
    pos <- p >= 0.6
    s <- summary(survival::survfit(survival::Surv(y[pos, 1], y[pos, 2]) ~ 1), times = 5)$surv
    expect_equal(nb$model, mean(pos) * ((1 - s) - s * 0.6 / 0.4), tolerance = 1e-10)
})

test_that("bootstrap optimism equals an explicit loop with coxph and concordance", {
    prep <- cox_prep(); X <- cp_design(prep$spec, prep$data); y <- prep$y
    fam <- cp_family("cox", horizon = 5)
    iv <- cp_internal(X, y, fam, "boot", B = 5, seed = 3, stability = TRUE)
    set.seed(3)
    oc <- os <- numeric(5)
    S <- survival::Surv(y[, 1], y[, 2])
    for (b in 1:5) {
        i <- sample(nrow(y), nrow(y), replace = TRUE)
        f <- survival::coxph(S[i] ~ X[i, ])
        lpb <- drop(X[i, ] %*% coef(f)); lpt <- drop(X %*% coef(f))
        oc[b] <- survival::concordance(S[i] ~ lpb, reverse = TRUE)$concordance -
            survival::concordance(S ~ lpt, reverse = TRUE)$concordance
        os[b] <- 1 - unname(coef(survival::coxph(S ~ lpt)))
    }
    expect_equal(iv$optimism[["c"]], mean(oc), tolerance = 1e-8)
    expect_equal(iv$optimism[["slope"]], mean(os), tolerance = 1e-6)
    expect_true(all(iv$stability$pts >= 0 & iv$stability$pts <= 1000))
})

test_that("sample size criteria match pmsampsize (survival)", {
    skip_if_not_installed("pmsampsize")
    ss <- cp_cox_samplesize(10, 0.1, rate = 0.065, meanfup = 2.07, horizon = 2)
    pm <- pmsampsize::pmsampsize(type = "s", csrsquared = 0.1, parameters = 10, rate = 0.065,
                                 timepoint = 2, meanfup = 2.07)
    expect_equal(unname(ss$n), unname(pm$results_table[1:2, 1]))
    expect_equal(round(ss$riskCI, 3), c(pm$int_lci, pm$int_uci))
})

test_that("the Cox analysis runs with every option, validates externally, saves risks and serialises", {
    d <- rott()
    f <- file.path(tempdir(), "gbsg_ext.csv")
    g <- survival::gbsg
    g$years <- g$rfstime / 365.25
    g$status <- ifelse(g$status == 1, "event", "censored")
    g$size <- cut(g$size, c(-Inf, 20, 50, Inf), labels = c("<=20", "20-50", ">50"))
    g$grade <- ifelse(g$grade == 3, "3", "1-2")
    utils::write.csv(g, f, row.names = FALSE)
    opts <- cpCoxOptions$new(time = "years", status = "status", eventLevel = "event", horizon = 5,
        covs = c("age", "nodes"), splineVars = "pgr", factors = c("size", "grade", "hormon"),
        extFile = list(path = f, filename = "gbsg.csv"), validation = "boot", nBoot = 10,
        stability = TRUE, equation = TRUE, sampleSize = TRUE, shapley = TRUE, permImp = TRUE, nPerm = 2,
        shap = TRUE, shapRow = 3, calPlot = TRUE, kmPlot = TRUE, dca = TRUE, explain = TRUE,
        interactions = list(list("age", "hormon")), intMode = "aic", descTable = TRUE, predTests = TRUE,
        coefs = TRUE)
    opts <- enable_output(opts, "predRisk")
    an <- cpCoxClass$new(options = opts, data = d, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    r <- an$results
    expect_null(r$get("problem"))
    expect_equal(r$samples$asDF$n, c(1546, 686))
    expect_equal(r$extCheck$asDF$status[r$extCheck$asDF$var == "status"], "OK")
    expect_true(all(is.finite(r$perfDev$asDF$corrected)))
    expect_false("(Intercept)" %in% rownames(r$coefTable$asDF))
    v <- unlist(an$results$predRisk$.__enclos_env__$private$.values[[1]])
    k <- cp_prepare(d, "years", splineVars = "pgr")$spec$vars$pgr$knots
    fit <- survival::coxph(survival::Surv(years, status == "event") ~ age + nodes + rms::rcs(pgr, k) +
                               size + grade + hormon, d)
    expect_equal(unname(v[10]), unname(1 - summary(survival::survfit(fit, newdata = d[10, ]), times = 5)$surv),
                 tolerance = 1e-6)
    jmvcore:::initProtoBuf()
    pb <- RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)
    expect_gt(length(pb), 1000)
    expect_lt(length(memCompress(serialize(r$perfDev$state, NULL), "gzip")), 500000)
    render <- c(calibration = ".plotCal", kmCurves = ".plotKm", dcaCurves = ".plotDca",
                instability = ".plotInstability", calInstability = ".plotCalInstability",
                shapSummary = ".plotShapSummary", shapPatient = ".plotShapPatient")
    priv <- an$.__enclos_env__$private
    for (p in names(render)) {
        grDevices::png(file.path(tempdir(), paste0(p, ".png")))
        ok <- priv[[render[[p]]]](r$plots[[p]], ggtheme = ggplot2::theme_grey(), theme = NULL)
        grDevices::dev.off()
        expect_true(ok, label = p)
    }
})

test_that("the default horizon is the median follow-up and problems reach the user", {
    d <- rott()
    r <- cpCox(data = d, time = "years", status = "status", eventLevel = "event", covs = "nodes", explain = TRUE)
    expect_match(r$samples$notes$ev$note, sprintf("t = %s", .num(median(d$years), 4)), fixed = TRUE)
    r <- cpCox(data = d, time = "years", status = "status", eventLevel = "event", covs = "nodes", horizon = 100)
    expect_match(r$get("problem")$content, "Nobody is followed beyond")
})
