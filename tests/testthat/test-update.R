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

test_that("Cox updating: slope and likelihood ratios as coxph, baseline re-estimated at the horizon", {
    set.seed(8)
    n <- 3000
    x <- matrix(rnorm(n * 2), n, dimnames = list(NULL, c("a", "b")))
    eta <- drop(x %*% c(0.7, -0.4))
    tm <- rexp(n, 0.1 * exp(eta)); cens <- runif(n, 0, 15)
    y <- cbind(time = pmin(tm, cens), status = as.integer(tm <= cens))
    fam <- cp_family("cox", horizon = 5)
    beta <- fam$fit(x, y); lp <- fam$lp(beta, x)
    u <- cp_update(fam, x, y, lp)
    expect_equal(u$chosen, "original")
    f1 <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ lp)
    expect_equal(u$methods$recal$b, unname(coef(f1)), tolerance = 1e-6)
    f0 <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ offset(lp))
    f2 <- survival::coxph(survival::Surv(y[, 1], y[, 2]) ~ x)
    expect_equal(u$methods$recal$test[["stat"]], 2 * (f1$loglik[2] - f0$loglik), tolerance = 1e-6)
    expect_equal(u$methods$revision$test[["stat"]], 2 * (f2$loglik[2] - f1$loglik[2]), tolerance = 1e-6)
    # a doubled hazard: the baseline needs updating, not the slope
    t2 <- rexp(n, 0.2 * exp(eta)); y2 <- cbind(time = pmin(t2, cens), status = as.integer(t2 <= cens))
    u2 <- cp_update(fam, x, y2, lp)
    expect_equal(u2$chosen, "citl")
    expect_equal(u2$methods$citl$a, log(2), tolerance = 0.1)
})

test_that("the Cox analysis fills the updating table", {
    d <- survival::rotterdam; d$status <- factor(ifelse(d$recur == 1 | d$death == 1, "event", "censored"))
    d$years <- pmin(d$rtime, d$dtime) / 365.25; d$set <- factor(ifelse(d$year < 1990, "early", "late"))
    r <- cpCox(data = d, time = "years", status = "status", eventLevel = "event", horizon = 5,
               covs = c("age", "nodes"), splitVar = "set", splitLevel = "late", updating = TRUE)
    expect_null(r$get("problem"))
    u <- r$updTable$asDF
    expect_equal(nrow(u), 4)
    expect_equal(sum(u$chosen == "✓", na.rm = TRUE), 1)
    expect_true(all(is.finite(u$c)))
})
