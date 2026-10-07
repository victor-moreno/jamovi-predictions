pima_age <- function() {
    d <- MASS::Pima.tr
    d$agegrp <- cut(d$age, c(0, 30, 45, 100), labels = c("<30", "30-45", "45+"))
    d
}

prep_int <- function(ints, refs = list(agegrp = "30-45"), scales = list(bmi = "ten")) {
    prep <- cp_prepare(pima_age(), "type", covs = c("bmi", "ped"), splineVars = "glu", factors = "agegrp",
                       refLevels = refs, scales = scales, interactions = ints)
    cp_binary_outcome(prep, "Yes")
}

test_that("spline x factor and covariate interactions reproduce glm", {
    skip_if_not_installed("rms")
    prep <- prep_int(list(c("glu", "agegrp"), c("bmi", "ped")))
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    b <- fam$fit(X, prep$y)
    d <- pima_age(); d$ag2 <- relevel(d$agegrp, "30-45"); d$bmi10 <- d$bmi / 10
    k <- prep$spec$vars$glu$knots
    g <- glm(type ~ bmi10 + ped + rms::rcs(glu, k) * ag2 + bmi10:ped, binomial, d)
    expect_equal(-fam$loss(prep$y, fam$lp(b, X)), as.numeric(logLik(g)), tolerance = 1e-7)
})

test_that("screening chi-square equals the sequential deviance test", {
    prep <- prep_int(list(c("bmi", "agegrp")))
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    sel <- cp_select_interactions(X, prep$y, prep$spec, fam, "lrt", 0)
    d <- pima_age(); d$ag2 <- relevel(d$agegrp, "30-45")
    k <- prep$spec$vars$glu$knots
    a <- glm(type ~ bmi + ped + glu + ag2, binomial, d)
    m0 <- cp_with_interactions(prep$spec, character(0))
    X0 <- cp_design(m0, prep$data)
    l0 <- fam$loss(prep$y, fam$lp(fam$fit(X0, prep$y), X0))
    l1 <- fam$loss(prep$y, fam$lp(fam$fit(X, prep$y), X))
    expect_equal(attr(sel, "screen")$chisq, 2 * (l0 - l1), tolerance = 1e-8)
    expect_equal(attr(sel, "screen")$df, 2)
    expect_length(attr(sel, "selected"), 0)
})

test_that("forward selection keeps a true interaction and drops noise", {
    set.seed(7); n <- 600
    d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)))
    d$y <- factor(rbinom(n, 1, plogis(-0.5 + 0.8 * d$x + 0.3 * d$z + 1.2 * d$x * (d$f == "b"))))
    prep <- cp_binary_outcome(cp_prepare(d, "y", covs = c("x", "z", "w"), factors = "f",
                                         interactions = utils::combn(c("x", "z", "w", "f"), 2, simplify = FALSE)), "1")
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    for (m in c("aic", "bic", "lrt")) {
        sel <- cp_select_interactions(X, prep$y, prep$spec, fam, m, 0.01)
        expect_true("x:f" %in% attr(sel, "selected"))
        expect_true(all(attr(sel, "selected") %in% c("x:f", "x:z", "z:f", "x:w", "w:f", "z:w")))
    }
    expect_equal(attr(cp_select_interactions(X, prep$y, prep$spec, fam, "bic"), "selected"), "x:f")
    fsel <- fam; fsel$fit <- function(X, y) cp_select_interactions(X, y, prep$spec, fam, "bic")
    iv <- cp_internal(X, prep$y, fsel, "boot", B = 15, seed = 2, stability = FALSE)
    expect_gt(iv$selection[["x:f"]], 0.9)
})

test_that("SHAP with interactions equals brute-force interventional Shapley values", {
    prep <- prep_int(list(c("bmi", "agegrp"), c("bmi", "ped"), c("glu", "agegrp")))
    spec <- prep$spec; fam <- cp_family("binary")
    X <- cp_design(spec, prep$data); b <- fam$fit(X, prep$y)
    phi <- cp_shap(b, spec, X)
    vars <- names(spec$vars); m <- length(vars)
    value <- function(i, S) {
        dd <- prep$data
        for (v in vars[S]) dd[[v]] <- prep$data[[v]][i]
        mean(fam$lp(b, cp_design(spec, dd)))
    }
    for (i in c(1, 50)) {
        brute <- vapply(seq_len(m), function(j) {
            others <- setdiff(seq_len(m), j); tot <- 0
            for (s in 0:(m - 1)) for (S in if (s == 0) list(integer(0)) else utils::combn(others, s, simplify = FALSE))
                tot <- tot + factorial(s) * factorial(m - s - 1) / factorial(m) * (value(i, c(S, j)) - value(i, S))
            tot
        }, 0)
        expect_equal(unname(phi[i, ]), brute, tolerance = 1e-10)
    }
    expect_equal(unname(rowSums(phi)) + attr(phi, "base"), fam$lp(b, X), tolerance = 1e-10)
})

test_that("reference levels, scaling and contrasts give the glm odds ratios", {
    skip_if_not_installed("rms")
    prep <- prep_int(list(), scales = list(bmi = "ten", ped = "sd"))
    spec <- prep$spec; fam <- cp_family("binary")
    X <- cp_design(spec, prep$data); b <- fam$fit(X, prep$y)
    eff <- cp_effects(spec, b, fam$vcov(b))
    d <- pima_age(); d$ag2 <- relevel(d$agegrp, "30-45")
    g <- glm(type ~ bmi + ped + rms::rcs(glu, spec$vars$glu$knots) + ag2, binomial, d)
    expect_equal(eff$est[eff$key == "bmi"], exp(10 * coef(g)[["bmi"]]), tolerance = 1e-7)
    expect_equal(eff$est[eff$key == "ped"], exp(sd(d$ped) * coef(g)[["ped"]]), tolerance = 1e-7)
    expect_equal(eff$est[eff$key == "agegrp: <30"], exp(coef(g)[["ag2<30"]]), tolerance = 1e-7)
    q <- quantile(d$glu, c(.25, .75))
    nd <- d[c(1, 1), ]; nd$glu <- q[c(2, 1)]
    expect_equal(eff$est[eff$key == "glu"], unname(exp(diff(rev(predict(g, nd))))), tolerance = 1e-7)
    cr <- cp_crude(spec, prep$data, prep$y, fam)
    expect_equal(cr$est[cr$key == "bmi"], exp(10 * coef(glm(type ~ bmi, binomial, d))[2]), tolerance = 1e-7, ignore_attr = TRUE)
    expect_equal(attr(cr, "pOverall")[["agegrp"]],
                 anova(glm(type ~ agegrp, binomial, d), test = "LRT")[2, "Pr(>Chi)"], tolerance = 1e-7)
})

test_that("effects within subgroups add the interaction coefficient", {
    skip_if_not_installed("rms")
    prep <- prep_int(list(c("bmi", "agegrp")), scales = list())
    spec <- prep$spec; fam <- cp_family("binary")
    X <- cp_design(spec, prep$data); b <- fam$fit(X, prep$y)
    sg <- cp_subgroups(spec, b, fam$vcov(b))
    row <- sg[sg$effect == "bmi: per unit" & sg$within == "agegrp = 45+", ]
    expect_equal(row$est, exp(b[["bmi"]] + b[["bmi × agegrp: 45+"]]), tolerance = 1e-10)
    lt <- cp_lr_tests(spec, X, prep$y, fam)
    d <- pima_age(); d$ag2 <- relevel(d$agegrp, "30-45")
    g1 <- glm(type ~ bmi * ag2 + ped + rms::rcs(glu, spec$vars$glu$knots), binomial, d)
    g0 <- update(g1, . ~ . - bmi:ag2)
    expect_equal(lt$chisq[lt$key == "bmi:agegrp"], deviance(g0) - deviance(g1), tolerance = 1e-7)
})

test_that("the analysis reports screening, selection, subgroups and descriptives", {
    d <- pima_age()
    r <- cpLogistic(data = d, dep = "type", covs = c("bmi", "ped"), splineVars = "glu", factors = "agegrp",
                    refLevels = list(list(var = "agegrp", ref = "30-45")),
                    covScales = list(list(var = "bmi", scale = "ten")),
                    interactions = list(list("bmi", "agegrp"), list("ped", "glu")), intMode = "lrt",
                    descTable = TRUE, predTests = TRUE)
    it <- r$intTable$asDF
    expect_equal(nrow(it), 2)
    or <- r$orTable$asDF
    expect_equal(or$level[1], "per 10 units")
    expect_equal(or$or[or$level == "30-45"], 1)
    expect_equal(r$descTable$asDF$total[1], "200")
    expect_true(nrow(r$testTable$asDF) >= 5)
})
