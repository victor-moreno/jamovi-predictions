tmpfile <- function(name) file.path(tempdir(), name)

test_that("a logistic model exported as JSON gives back the same predictions", {
    d <- pima_stacked()
    prep <- cp_prepare(d[d$sample == "train", ], "type", covs = c("bp", "bmi"), splineVars = "glu",
                       factors = "agegrp", refLevels = list(agegrp = "30-45"), scales = list(bmi = "ten"),
                       interactions = list(c("bmi", "agegrp")))
    prep <- cp_binary_outcome(prep, "Yes")
    fam <- cp_family("binary")
    X <- cp_design(prep$spec, prep$data); b <- fam$fit(X, prep$y)
    f <- tmpfile("m.json")
    writeLines(cp_model_json(prep$spec, b, "binary", se = sqrt(diag(fam$vcov(b)))), f)
    m <- cp_model_read(f, "m.json")
    expect_equal(m$kind, "spec"); expect_equal(m$family, "binary")
    te <- d[d$sample == "test", ]
    des <- cp_model_design(m, te)
    lp1 <- fam$lp(cp_model_beta(m), des$X)
    lp0 <- fam$lp(b, cp_design(prep$spec, cp_split_data(te, prep$spec, "type")$data))
    expect_equal(unname(lp1), unname(lp0), tolerance = 1e-12)
    expect_equal(m$se, unname(sqrt(diag(fam$vcov(b))))[-1], tolerance = 1e-12)
})

test_that("a Cox model exported as JSON keeps its baseline survival at the horizon", {
    d <- survival::gbsg; d$status <- factor(d$status)
    prep <- cp_cox_outcome(cp_prepare(d, c("rfstime", "status"), covs = c("age", "nodes"), splineVars = "pgr"),
                           "rfstime", "status", "1", 1826)
    fam <- cp_family("cox", horizon = 1826)
    X <- cp_design(prep$spec, prep$data); b <- fam$fit(X, prep$y)
    f <- tmpfile("cox.json")
    writeLines(cp_model_json(prep$spec, b, "cox"), f)
    m <- cp_model_read(f, "cox.json")
    expect_equal(m$horizon, 1826)
    lp <- fam$lp(cp_model_beta(m), cp_model_design(m, d)$X)
    expect_equal(fam$risk(lp), fam$risk(fam$lp(b, X)), tolerance = 1e-10)
})

test_that("a coefficient table with transforms, products and levels reproduces glm", {
    d <- pima_stacked()[1:200, ]
    g <- glm(type ~ glu + log(bmi) + I(age^2) + npreg * ped + agegrp, binomial, d)
    cf <- coef(g); se <- sqrt(diag(vcov(g)))
    terms <- c("(Intercept)", "glu", "log(bmi)", "age^2", "npreg", "ped", "agegrp = 30-45", "agegrp45+", "npreg × ped")
    tab <- data.frame(Term = terms, Estimate = unname(cf), EE = unname(se), check.names = FALSE)
    f <- tmpfile("coefs.csv"); utils::write.csv(tab, f, row.names = FALSE)
    m <- cp_model_read(f, "coefs.csv"); m$family <- "binary"
    expect_equal(m$intercept, unname(cf[1])); expect_equal(m$se, unname(se[-1]))
    des <- cp_model_design(m, d)
    lp <- cp_lp(cp_model_beta(m), des$X)
    expect_equal(unname(lp), unname(predict(g)), tolerance = 1e-10)
    expect_true(all(des$check$status == "OK"))
})

test_that("semicolon tables with decimal commas and unknown terms are handled", {
    f <- tmpfile("semi.csv")
    writeLines(c("term;coef;se", "Intercept;-5,2;0,9", "glu;0,03;0,01", "chol;0,01;0,002"), f)
    m <- cp_model_read(f, "semi.csv"); m$family <- "binary"
    expect_equal(m$estimate, c(0.03, 0.01))
    expect_error(cp_model_design(m, MASS::Pima.tr), class = "cp_user_error", regexp = "'chol'")
    expect_error(cp_model_design(list(terms = "type", kind = "terms"), MASS::Pima.tr),
                 class = "cp_user_error", regexp = "categorical")
})

test_that("a published Cox table with a centred index and baseline survival reproduces survfit", {
    d <- survival::gbsg
    f0 <- survival::coxph(survival::Surv(rfstime, status) ~ age + nodes + hormon, d)
    mns <- colMeans(d[, c("age", "nodes", "hormon")])
    s0 <- summary(survival::survfit(f0, newdata = as.data.frame(t(mns))), times = 1826)$surv
    tab <- data.frame(term = c("age", "nodes", "hormon", "baseline survival", "horizon", "centre"),
                      estimate = c(coef(f0), s0, 1826, sum(coef(f0) * mns)))
    f <- tmpfile("cox.csv"); utils::write.csv(tab, f, row.names = FALSE)
    m <- cp_model_read(f, "cox.csv"); m$family <- "cox"
    lp <- cp_lp(cp_model_beta(m), cp_model_design(m, d)$X)
    pr <- 1 - summary(survival::survfit(f0, newdata = d[1:5, ]), times = 1826)$surv
    expect_equal(.cox_risk(lp[1:5]), as.numeric(pr), tolerance = 1e-8)
})

ext_perf <- function(r) {
    pv <- r$perfVal$asDF
    stats::setNames(pv$est, sub("^ext ", "", unlist(r$perfVal$rowKeys)))
}

test_that("a coefficient table validated with cpPublished equals cpLogistic's external validation", {
    tr <- MASS::Pima.tr; te <- MASS::Pima.te
    g <- glm(type ~ glu + bmi + ped, binomial, tr)
    f <- tmpfile("pima.csv")
    utils::write.csv(data.frame(term = c("Intercept", "glu", "bmi", "ped"), estimate = coef(g),
                                se = sqrt(diag(vcov(g)))), f, row.names = FALSE)
    rp <- cpPublished(data = te, modelFile = list(path = f, filename = "pima.csv"), family = "binary",
                      dep = "type", eventLevel = "Yes", covs = c("glu", "bmi", "ped"), updating = TRUE)
    fe <- tmpfile("pima_te.csv"); utils::write.csv(te, fe, row.names = FALSE)
    rl <- cpLogistic(data = tr, dep = "type", eventLevel = "Yes", covs = c("glu", "bmi", "ped"),
                     extFile = list(path = fe, filename = "pima_te.csv"), updating = TRUE)
    pp <- rp$perf$asDF
    keys <- unlist(rp$perf$rowKeys)
    expect_equal(stats::setNames(pp$est, keys), ext_perf(rl)[keys], tolerance = 1e-8)
    expect_equal(rp$updTable$asDF$stat, rl$updTable$asDF$stat, tolerance = 1e-6)
    expect_equal(rp$modelTable$asDF$ratio[2], unname(exp(coef(g)[2])), tolerance = 1e-12)
})

test_that("a clinpred JSON model (spline, factor) validated with cpPublished equals the external validation", {
    d <- pima_stacked(); tr <- d[d$sample == "train", ]; te <- d[d$sample == "test", ]
    prep <- cp_binary_outcome(cp_prepare(tr, "type", covs = "bmi", splineVars = "glu", factors = "agegrp"), "Yes")
    fam <- cp_family("binary"); X <- cp_design(prep$spec, prep$data); b <- fam$fit(X, prep$y)
    f <- tmpfile("pima.json"); writeLines(cp_model_json(prep$spec, b, "binary"), f)
    rp <- cpPublished(data = te, modelFile = list(path = f, filename = "pima.json"), family = "binary",
                      dep = "type", eventLevel = "Yes", covs = c("bmi", "glu", "agegrp"))
    rl <- cpLogistic(data = d, dep = "type", eventLevel = "Yes", covs = "bmi", splineVars = "glu", factors = "agegrp",
                     splitVar = "sample", splitLevel = "test")
    pp <- rp$perf$asDF
    pv <- rl$perfVal$asDF
    expect_equal(pp$est, pv$est, tolerance = 1e-8)
    r2 <- cpPublished(data = te, modelFile = list(path = f, filename = "pima.json"), family = "cox",
                      dep = "type", eventLevel = "Yes", covs = "bmi", time = "glu")
    expect_match(r2$get("problem")$content, "of type Binary \\(logistic\\)")
    r3 <- cpPublished(data = te, modelFile = list(path = f, filename = "pima.json"), family = "binary",
                      dep = "type", eventLevel = "Yes", covs = "bmi")
    expect_match(r3$get("problem")$content, "add the variables to Predictor variables")
})

test_that("a Cox model from JSON keeps its horizon; without a baseline it is estimated from the data", {
    d <- survival::gbsg; d$status <- factor(d$status, 0:1, c("censored", "event"))
    tr <- d[1:400, ]; te <- d[401:686, ]
    prep <- cp_cox_outcome(cp_prepare(tr, c("rfstime", "status"), covs = c("age", "nodes")), "rfstime", "status", "event", 1500)
    fam <- cp_family("cox", horizon = 1500); X <- cp_design(prep$spec, prep$data); b <- fam$fit(X, prep$y)
    f <- tmpfile("gbsg.json"); writeLines(cp_model_json(prep$spec, b, "cox"), f)
    rp <- cpPublished(data = te, modelFile = list(path = f, filename = "gbsg.json"), family = "cox",
                      dep = "status", eventLevel = "event", time = "rfstime", covs = c("age", "nodes"),
                      dca = TRUE, kmPlot = TRUE, calPlot = TRUE)
    expect_null(rp$get("problem"))
    rc <- cpCox(data = rbind(cbind(tr, s = "a"), cbind(te, s = "b")) |> transform(s = factor(s)),
                time = "rfstime", status = "status", eventLevel = "event", horizon = 1500, covs = c("age", "nodes"),
                splitVar = "s", splitLevel = "b")
    expect_equal(rp$perf$asDF$est, rc$perfVal$asDF$est, tolerance = 1e-8)
    fc <- tmpfile("gbsg.csv")
    utils::write.csv(data.frame(term = c("age", "nodes"), estimate = b[-1]), fc, row.names = FALSE)
    rn <- cpPublished(data = te, modelFile = list(path = fc, filename = "gbsg.csv"), family = "cox",
                      dep = "status", eventLevel = "event", time = "rfstime", covs = c("age", "nodes"), horizon = 1500)
    pn <- rn$perf$asDF; kn <- unlist(rn$perf$rowKeys)
    expect_equal(pn$est[kn == "oe"], 1, tolerance = 0.05)
    expect_match(rn$modelTable$notes$s0$note, "estimated from these data")
    expect_equal(pn$est[kn == "c"], rp$perf$asDF$est[1], tolerance = 1e-10)
})

test_that("a linear coefficient table reproduces lm predictions, and the tables render cleanly", {
    tr <- MASS::Pima.tr; te <- MASS::Pima.te
    m <- lm(glu ~ bmi + age + type, tr)
    f <- tmpfile("lin.csv")
    utils::write.csv(data.frame(term = c("(Intercept)", "bmi", "age", "type = Yes"), estimate = coef(m)), f, row.names = FALSE)
    opts <- cpPublishedOptions$new(modelFile = list(path = f, filename = "lin.csv"), family = "gaussian",
                                   dep = "glu", covs = c("bmi", "age", "type"), updating = TRUE, calPlot = TRUE)
    opts <- enable_output(opts, "pred")
    an <- cpPublishedClass$new(options = opts, data = te, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    v <- unlist(an$results$pred$.__enclos_env__$private$.values[[1]])
    expect_equal(unname(v), unname(predict(m, te)), tolerance = 1e-10)
    pr <- predict(m, te)
    expect_equal(an$results$perf$asDF$est[1], 1 - sum((te$glu - pr)^2) / sum((te$glu - mean(te$glu))^2), tolerance = 1e-10)
    src <- readLines(test_path("test-render.R"))
    eval(parse(text = src[seq_len(grep("^test_that", src)[1] - 1)]))
    expect_equal(heading_problems(an$results), character(0))
    jmvcore:::initProtoBuf()
    expect_gt(length(RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)), 1000)
    grDevices::png(file.path(tempdir(), "pubcal.png"))
    expect_true(an$.__enclos_env__$private$.plotCal(an$results$plots$calibration, ggtheme = ggplot2::theme_grey(), theme = NULL))
    grDevices::dev.off()
})

test_that("a missing-indicator model exported as JSON rebuilds its indicators from other data", {
    d <- pima_stacked()
    d$bmi[c(2, 5, 9, 300, 310)] <- NA
    d$agegrp[c(4, 320)] <- NA
    dev <- d[d$sample == "train", ]
    h <- cp_miss_category(dev, c("bmi", "glu"), "agegrp")
    prep <- cp_binary_outcome(cp_prepare(cp_miss_category_apply(h, dev), "type",
                                         covs = c("glu", "bmi", "bmi (missing)"), factors = "agegrp"), "Yes")
    fam <- cp_family("binary")
    X <- cp_design(prep$spec, prep$data); b <- fam$fit(X, prep$y)
    f <- tmpfile("miss.json")
    writeLines(cp_model_json(prep$spec, b, "binary", missing = h), f)
    m <- cp_model_read(f, "miss.json")
    te <- d[d$sample == "test", c("type", "glu", "bmi", "agegrp")]
    des <- cp_model_design(m, te)
    expect_false(anyNA(des$X))
    lp0 <- fam$lp(b, cp_design(prep$spec, cp_miss_category_apply(h, te)))
    expect_equal(unname(fam$lp(cp_model_beta(m), des$X)), unname(lp0), tolerance = 1e-12)
})

test_that("several data levels falling into the reference of a coefficient table are a warning", {
    f <- tmpfile("ref.csv")
    writeLines(c("term,estimate", "Intercept,-1", "agegrp = 45+,0.5"), f)
    m <- cp_model_read(f, "ref.csv"); m$family <- "binary"
    d <- pima_stacked()
    ck <- cp_model_design(m, d)$check
    expect_equal(ck$status[ck$var == "agegrp"], "Warning")
    expect_match(ck$details[ck$var == "agegrp"], "as reference: 30-45, <30", fixed = TRUE)
})
