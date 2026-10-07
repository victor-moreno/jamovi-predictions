test_that("exact Shapley R2 shares add up to the full-model R2", {
    prep <- pima_prep(factors = "agegrp", data = pima_stacked()[1:200, ])
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    sh <- cp_shapley_r2(X, prep$y, prep$spec, fam, "r2mf")
    expect_true(attr(sh, "exact"))
    expect_equal(sum(sh$shapley), attr(sh, "total"), tolerance = 1e-10)
    expect_equal(attr(sh, "nFits"), 2^6 - 1)
})

test_that("with two predictors Shapley is the average of the two entry orders", {
    prep <- cp_binary_outcome(cp_prepare(pima_dev(), "type", covs = c("glu", "bmi")), "Yes")
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary"); y <- prep$y
    r2 <- function(cols) { Xs <- X[, cols, drop = FALSE]; fam$r2(y, fam$lp(fam$fit(Xs, y), Xs), "r2mf") }
    full <- r2(1:2); a <- r2(1); b <- r2(2)
    sh <- cp_shapley_r2(X, y, prep$spec, fam, "r2mf")
    expect_equal(sh$shapley, c((a + full - b) / 2, (b + full - a) / 2), tolerance = 1e-10)
})

test_that("Monte Carlo Shapley approximates the exact values", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    ex <- cp_shapley_r2(X, prep$y, prep$spec, fam, "r2mf")
    mc <- cp_shapley_r2(X, prep$y, prep$spec, fam, "r2mf", maxExact = 2, nPerm = 300)
    expect_false(attr(mc, "exact"))
    expect_equal(sum(mc$shapley), attr(ex, "total"), tolerance = 1e-10)
    expect_equal(mc$shapley, ex$shapley, tolerance = 0.15)
})

test_that("SHAP values add up to the linear predictor minus its mean", {
    prep <- pima_prep(factors = "agegrp", data = pima_stacked()[1:200, ])
    X <- cp_design(prep$spec, prep$data)
    beta <- cp_binary_fit(X, prep$y)
    phi <- cp_shap(beta, prep$spec, X)
    lp <- cp_lp(beta, X)
    expect_equal(unname(rowSums(phi)) + attr(phi, "base"), lp, tolerance = 1e-10)
    expect_equal(attr(phi, "base"), mean(lp), tolerance = 1e-10)
    expect_equal(unname(colMeans(phi)), rep(0, ncol(phi)), tolerance = 1e-10)
})

test_that("permutation importance is positive for a strong predictor and near 0 for noise", {
    set.seed(5)
    d <- data.frame(x = rnorm(400), noise = rnorm(400))
    d$y <- factor(rbinom(400, 1, plogis(2 * d$x)))
    prep <- cp_binary_outcome(cp_prepare(d, "y", covs = c("x", "noise")), "1")
    beta <- cp_binary_fit(cp_design(prep$spec, prep$data), prep$y)
    pi <- cp_perm_importance(beta, prep$spec, prep$data, prep$y, cp_family("binary"), nPerm = 10)
    expect_gt(pi$d_auc[1], 0.2)
    expect_gt(pi$d_brier[1], 0.05)
    expect_lt(abs(pi$d_auc[2]), 0.02)
})
