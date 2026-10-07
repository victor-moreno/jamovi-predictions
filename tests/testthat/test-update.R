test_that("updating methods and tests reproduce glm in external data", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data); fam <- cp_family("binary")
    beta <- fam$fit(X, prep$y)
    val <- cp_split_data(pima_val(), prep$spec, "type")$data
    yv <- as.integer(val$type == "Yes"); Xv <- cp_design(prep$spec, val)
    lp <- fam$lp(beta, Xv)
    u <- cp_update(fam, Xv, yv, lp)
    g1 <- glm(yv ~ 1, offset = lp, family = binomial)
    g2 <- glm(yv ~ lp, family = binomial)
    g3 <- glm(yv ~ Xv, family = binomial)
    expect_equal(u$methods$citl$a, unname(coef(g1)), tolerance = 1e-6)
    expect_equal(c(u$methods$recal$a, u$methods$recal$b), unname(coef(g2)), tolerance = 1e-6)
    l0 <- sum(dbinom(yv, 1, plogis(lp), log = TRUE))
    ll <- c(original = l0, citl = as.numeric(logLik(g1)), recal = as.numeric(logLik(g2)))
    for (m in names(ll)) {
        expect_equal(u$methods[[m]]$test[["stat"]], 2 * (as.numeric(logLik(g3)) - ll[[m]]), tolerance = 1e-6)
        expect_equal(u$methods[[m]]$test[["df"]], length(coef(g3)) - c(original = 0, citl = 1, recal = 2)[[m]])
    }
    expect_true(u$chosen %in% c("original", "citl", "recal", "revision"))
})

test_that("closed testing keeps the original model when nothing changed and updates a shifted intercept", {
    set.seed(5)
    n <- 4000
    x <- matrix(rnorm(n * 2), n, dimnames = list(NULL, c("a", "b")))
    lp <- -1 + x %*% c(0.8, -0.5)
    y <- rbinom(n, 1, plogis(lp))
    fam <- cp_family("binary")
    expect_equal(cp_update(fam, x, y, drop(lp))$chosen, "original")
    y2 <- rbinom(n, 1, plogis(lp + 1))
    expect_equal(cp_update(fam, x, y2, drop(lp))$chosen, "citl")
})

test_that("the analysis fills the updating table for each validation set", {
    d <- pima_stacked()
    r <- cpLogistic(data = d, dep = "type", eventLevel = "Yes", covs = c("glu", "bmi", "ped"),
                    splitVar = "sample", splitLevel = "test", updating = TRUE)
    u <- r$updTable$asDF
    expect_equal(nrow(u), 4)
    expect_equal(sum(u$chosen == "✓", na.rm = TRUE), 1)
    expect_true(all(is.finite(u$auc)))
})
