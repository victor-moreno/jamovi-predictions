test_that("spline knots and basis match Hmisc", {
    x <- MASS::Pima.tr$glu
    k <- .rcs_knots(x, 4, "glu")
    expect_equal(k, unname(Hmisc::rcspline.eval(x, nk = 4, knots.only = TRUE)))
    expect_equal(.rcs_basis(x, k), unname(Hmisc::rcspline.eval(x, knots = k)), ignore_attr = TRUE)
    small <- x[1:60]
    expect_equal(.rcs_knots(small, 5, "glu"), unname(Hmisc::rcspline.eval(small, nk = 5, knots.only = TRUE)))
})

test_that("logistic fit matches rms::lrm", {
    prep <- pima_prep()
    beta <- cp_binary_fit(cp_design(prep$spec, prep$data), prep$y)
    k <- prep$spec$vars$glu$knots
    g <- rms::lrm(type ~ npreg + bp + bmi + ped + rms::rcs(glu, k), data = pima_dev())
    expect_equal(unname(beta), unname(coef(g)), tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("apparent metrics match rms::val.prob and pROC", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data); y <- prep$y
    lp <- cp_lp(cp_binary_fit(X, y), X); p <- plogis(lp)
    m <- cp_binary_metrics(y, lp)
    vp <- rms::val.prob(p, y, pl = FALSE)
    expect_equal(m[["auc"]], unname(vp["C (ROC)"]), tolerance = 1e-7)
    expect_equal(m[["brier"]], unname(vp["Brier"]), tolerance = 1e-7)
    expect_equal(m[["r2n"]], unname(vp["R2"]), tolerance = 1e-7)
    expect_equal(m[["slope"]], 1, tolerance = 1e-7)
    expect_equal(m[["citl"]], 0, tolerance = 1e-7)
    roc <- pROC::roc(y, p, quiet = TRUE, direction = "<")
    expect_equal(.auc_delong(y, p)$se, sqrt(pROC::var(roc)), tolerance = 1e-7)
})

test_that("validation metrics match val.prob for a miscalibrated model", {
    prep <- pima_prep()
    beta <- cp_binary_fit(cp_design(prep$spec, prep$data), prep$y)
    vd <- cp_split_data(pima_val(), prep$spec, "type")$data
    yv <- as.integer(vd$type == "Yes")
    lp <- cp_lp(beta, cp_design(prep$spec, vd))
    m <- cp_binary_metrics(yv, lp)
    vp <- rms::val.prob(plogis(lp), yv, pl = FALSE)
    expect_equal(m[["slope"]], unname(vp["Slope"]), tolerance = 1e-6)
    expect_equal(m[["brier"]], unname(vp["Brier"]), tolerance = 1e-7)
    oracle_citl <- coef(glm(yv ~ offset(lp), family = binomial))[[1]]
    expect_equal(m[["citl"]], oracle_citl, tolerance = 1e-6)
    sm <- predict(loess(yv ~ plogis(lp), degree = 2, span = 0.75))
    expect_equal(m[["ici"]], mean(abs(pmin(pmax(sm, 0), 1) - plogis(lp))), tolerance = 1e-6)
})

test_that("bootstrap optimism matches rms::validate with the same seed", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data)
    k <- prep$spec$vars$glu$knots
    g <- rms::lrm(type ~ npreg + bp + bmi + ped + rms::rcs(glu, k), data = pima_dev(), x = TRUE, y = TRUE)
    set.seed(11); vv <- rms::validate(g, B = 30)
    iv <- cp_internal(X, prep$y, cp_family("binary"), "boot", B = 30, seed = 11, stability = FALSE)
    expect_equal(2 * iv$optimism[["auc"]], vv["Dxy", "optimism"], tolerance = 1e-6)
    expect_equal(iv$optimism[["slope"]], vv["Slope", "optimism"], tolerance = 1e-6)
    expect_equal(iv$optimism[["brier"]], vv["B", "optimism"], tolerance = 1e-6)
    expect_equal(iv$shrinkage, vv["Slope", "index.corrected"], tolerance = 1e-6)
})

test_that("cross-validation equals a hand-made out-of-fold computation", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data); y <- prep$y
    cv <- cp_internal(X, y, cp_family("binary"), "cv", k = 5, R = 1, seed = 3)
    set.seed(3); fold <- .folds(y, 5)
    oof <- numeric(length(y))
    for (f in 1:5) {
        fit <- glm.fit(cbind(1, X[fold != f, ]), y[fold != f], family = binomial())
        oof[fold == f] <- cbind(1, X[fold == f, ]) %*% fit$coefficients
    }
    expect_equal(cv$corrected[["auc"]], cp_binary_metrics(y, oof)[["auc"]], tolerance = 1e-10)
})

test_that("Riley sample size reproduces the published example", {
    ss <- cp_binary_samplesize(24, 0.174, 0.288)
    expect_equal(unname(ss$n), c(623, 662, 221))
})

test_that("net benefit and classification follow their definitions", {
    y <- c(1, 1, 0, 0, 1, 0, 0, 0); p <- c(.9, .6, .55, .2, .3, .1, .4, .05)
    nb <- cp_net_benefit(y, p, 0.5)
    expect_equal(nb$model, 2 / 8 - 1 / 8 * 1)
    expect_equal(nb$all, 3 / 8 - 5 / 8)
    cl <- cp_binary_classify(y, p, 0.5)
    expect_equal(c(cl$tp, cl$fp, cl$fn, cl$tn), c(2, 1, 1, 4))
    expect_equal(cl$sens[1], 2 / 3)
    expect_equal(cp_binary_threshold(y, p, "prevalence", 0.5), 3 / 8)
})
