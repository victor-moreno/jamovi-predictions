dir <- normalizePath(tempdir())

test_that("csv with semicolons and decimal commas is read", {
    f <- file.path(dir, "semi.csv")
    writeLines(c("glu;bmi;type", "120;33,5;Yes", "90;25,1;No"), f)
    d <- cp_read_external(f, "semi.csv")
    expect_equal(d$bmi, c(33.5, 25.1))
})

test_that("harmonise matches names and levels ignoring case, and reports problems", {
    prep <- pima_prep(factors = "agegrp", data = pima_stacked()[1:200, ])
    ext <- pima_stacked()[201:532, ]
    names(ext)[names(ext) == "glu"] <- "GLU"
    ext$type <- tolower(as.character(ext$type))
    ext$agegrp <- as.character(ext$agegrp)
    ext$agegrp[1:3] <- "unknown"
    ext$bmi[4:5] <- NA
    h <- cp_harmonise(ext, prep$spec, "type", "External")
    ch <- setNames(h$check$status, h$check$var)
    expect_equal(ch[["glu"]], "OK")
    expect_match(h$check$details[h$check$var == "glu"], "GLU")
    expect_equal(ch[["agegrp"]], "Warning")
    expect_equal(h$nMissing, 5)
    expect_equal(levels(h$data$type), c("No", "Yes"))
    expect_equal(sum(h$data$type == "Yes"), sum(ext$type[-(1:5)] == "yes"))
})

test_that("a missing predictor stops with the module's own message", {
    prep <- pima_prep()
    ext <- pima_val(); ext$ped <- NULL
    expect_error(cp_harmonise(ext, prep$spec, "type", "External (x.csv)"),
                 class = "cp_user_error", regexp = "'ped'")
})

test_that("omv and xlsx files are read", {
    skip_if_not_installed("openxlsx")
    f <- file.path(dir, "val.xlsx")
    openxlsx::write.xlsx(pima_val(), f)
    expect_equal(nrow(cp_read_external(f, "val.xlsx")), 332)
    skip_if_not_installed("jmvReadWrite")
    g <- file.path(dir, "val.omv")
    jmvReadWrite::write_omv(pima_val(), g)
    expect_equal(nrow(cp_read_external(g, "val.omv")), 332)
})

test_that("membership c is near 0.5 for random halves and high for shifted data", {
    prep <- pima_prep()
    X <- cp_design(prep$spec, prep$data); y <- prep$y
    set.seed(2); i <- sample(200, 100)
    expect_lt(cp_membership_c(X[i, ], y[i], X[-i, ], y[-i]), 0.7)
    Xs <- X[-i, ]; Xs[, "bmi"] <- Xs[, "bmi"] + 8
    expect_gt(cp_membership_c(X[i, ], y[i], Xs, y[-i]), 0.8)
})
