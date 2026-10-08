enet_data <- function() {
    f <- test_path("..", "..", "data", "heart_cohort.omv")
    skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f)
    d <- d[stats::complete.cases(d), ]
    X <- stats::model.matrix(~ age + sex + smoking + sbp + chol + bmi + diabetes, d)[, -1]
    list(d = d, X = X)
}

test_that("the elastic net path equals glmnet (linear, logistic) with an unpenalised column", {
    skip_if_not_installed("glmnet")
    e <- enet_data(); X <- e$X; d <- e$d
    pf <- c(0, rep(1, ncol(X) - 1))
    for (fam in c("gaussian", "binary")) {
        y <- if (fam == "gaussian") d$sbp1y else as.integer(d$event10 == "yes")
        g <- glmnet::glmnet(X, y, family = if (fam == "gaussian") "gaussian" else "binomial", alpha = 0.5,
                            penalty.factor = pf, control = list(thresh = 1e-12))
        ours <- cp_enet_path(X, y, fam, alpha = 0.5, pf = pf, lambda = g$lambda[c(5, 20, 40)])
        expect_equal(rbind(ours$a0, ours$beta), as.matrix(stats::coef(g, s = g$lambda[c(5, 20, 40)])),
                     tolerance = 1e-5, ignore_attr = TRUE)
        expect_equal(cp_enet_path(X, y, fam, alpha = 0.5, pf = pf)$lambda[1], g$lambda[1], tolerance = 1e-8)
    }
})

test_that("the Cox elastic net reaches at least glmnet's penalised partial likelihood", {
    skip_if_not_installed("glmnet")
    e <- enet_data(); X <- e$X; d <- e$d
    y <- cbind(d$time, as.integer(d$status == "event"))
    g <- glmnet::glmnet(X, survival::Surv(y[, 1], y[, 2]), family = "cox", alpha = 0.5, cox.ties = "breslow",
                        control = list(thresh = 1e-10))
    lam <- g$lambda[c(3, length(g$lambda))]
    ours <- cp_enet_path(X, y, "cox", alpha = 0.5, lambda = lam)
    sdv <- sqrt(colMeans(sweep(X, 2, colMeans(X))^2))
    obj <- function(b, l) .enet_loss("cox", y, drop(X %*% b)) / nrow(X) + l * sum(0.25 * (b * sdv)^2 + 0.5 * abs(b * sdv))
    for (j in 1:2) {
        gb <- as.matrix(stats::coef(g, s = lam[j]))[, 1]
        expect_lte(obj(ours$beta[, j], lam[j]), obj(gb, lam[j]) + 1e-10)
        expect_equal(unname(ours$beta[, j]), unname(gb), tolerance = 0.01)
    }
    expect_equal(cp_enet_path(X, y, "cox", alpha = 0.5)$lambda[1], g$lambda[1], tolerance = 1e-8)
})

test_that("cross-validated deviance equals cv.glmnet with the same folds", {
    skip_if_not_installed("glmnet")
    e <- enet_data(); X <- e$X; d <- e$d
    y <- as.integer(d$event10 == "yes")
    cv <- cp_enet_cv(X, y, "binary", alpha = 1, pf = rep(1, ncol(X)), nfolds = 5, seed = 3, strata = y)
    foldid <- .with_seed(3, .folds(y, 5))
    g <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1, foldid = foldid, lambda = cv$lambda,
                           type.measure = "deviance", control = list(thresh = 1e-12))
    expect_equal(cv$cvm, g$cvm, tolerance = 1e-5)
    expect_equal(cv$lambda[cv$imin], g$lambda.min)
    expect_equal(cv$lambda[cv$i1se], g$lambda.1se)
})

test_that("the penalised fit keeps baseline columns, and its Cox baseline is at the horizon", {
    e <- enet_data(); X <- e$X; d <- e$d
    y <- cbind(d$time, as.integer(d$status == "event"))
    fam <- cp_family("cox", horizon = 5)
    pf <- c(0, rep(1, ncol(X) - 1))
    b <- cp_enet_fit(X, y, fam, alpha = 1, pf = pf, rule = "1se", nfolds = 5, seed = 1)
    expect_true(isTRUE(attr(b, "penalized")))
    expect_true(b[["age"]] != 0)
    expect_equal(b[[1]], log(.cox_H0(y[, 1], y[, 2], drop(X %*% b[-1]), 5)))
})

test_that("the three analyses fit an elastic net with baseline covariates, validated and rendered", {
    f <- test_path("..", "..", "data", "heart_cohort.omv")
    skip_if_not(file.exists(f))
    d <- jmvReadWrite::read_omv(f); set.seed(2); d$noise <- stats::rnorm(nrow(d))
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    base <- list(list(var = "sex", forced = TRUE))
    common <- list(data = d, covs = c("sbp", "chol", "noise"), splineVars = "age", factors = c("sex", "smoking"),
                   baseline = base, selMode = "enet", enetFolds = 5, coefs = TRUE, predTests = TRUE, shapley = TRUE,
                   shap = TRUE, permImp = TRUE, nPerm = 2, splitVar = "period", splitLevel = "2017-2019",
                   validation = "boot", nBoot = 10, enetTune = FALSE)
    runs <- list(
        logistic = do.call(cpLogisticClass$new, list(options = do.call(cpLogisticOptions$new, c(common[-1], list(dep = "event10", eventLevel = "yes"))),
                                                     data = d, datasetId = "d", analysisId = "1", revision = 1)),
        linear = do.call(cpLinearClass$new, list(options = do.call(cpLinearOptions$new, c(common[-1], list(dep = "sbp1y"))),
                                                 data = d, datasetId = "d", analysisId = "1", revision = 1)),
        cox = do.call(cpCoxClass$new, list(options = do.call(cpCoxOptions$new, c(common[-1], list(time = "time", status = "status",
                                                                                                 eventLevel = "event", horizon = 5))),
                                           data = d, datasetId = "d", analysisId = "1", revision = 1)))
    for (nm in names(runs)) {
        an <- runs[[nm]]; an$init(); an$run()
        r <- an$results
        expect_null(r$get("problem"), label = nm)
        expect_equal(r$selTable$asDF$status[r$selTable$asDF$var == "sex"], "baseline", label = nm)
        expect_true(all(is.na(r$coefTable$asDF$se)), label = nm)
        expect_true(all(is.finite(r$perfDev$asDF$corrected)), label = nm)
        expect_equal(heading_problems(r), character(0), label = nm)
        grDevices::png(file.path(tempdir(), "enetcv.png"))
        expect_true(an$.__enclos_env__$private$.plotEnet(r$plots$enetCv, ggtheme = ggplot2::theme_grey(), theme = NULL), label = nm)
        grDevices::dev.off()
    }
})
