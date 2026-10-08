# Example data for every clinpred feature (run from the repository root):
#   Rscript data-raw/make_examples.R
# A simulated cardiovascular cohort with three outcomes (10-year event, blood
# pressure at one year, time to event), a non-linear age effect, an age x sex
# interaction, three-level factors, missing values and a period for split
# validation; an external population with a different case-mix and baseline
# risk, written with the problems the variable check reports; and published
# models (coefficient tables, clinpred JSON) for Validate a Published Model.
# Needs the clinpred package (pkgload), survival, jmvReadWrite, openxlsx.

pkgload::load_all(".", quiet = TRUE)
set.seed(2026)

# ---- simulation ------------------------------------------------------------

sim <- function(n, ageMean = 60, pCurrent = 0.2, shift = 0, effect = 1) {
    trunc <- function(x, lo, hi) pmin(pmax(x, lo), hi)
    age <- round(trunc(rnorm(n, ageMean, 10), 35, 85))
    sex <- factor(ifelse(runif(n) < 0.5, "male", "female"), c("female", "male"))
    smoking <- factor(sample(c("never", "former", "current"), n, TRUE, c(0.7 - pCurrent, 0.3, pCurrent)),
                      c("never", "former", "current"))
    bmi <- round(trunc(rnorm(n, 27, 4.5), 17, 50), 1)
    sbp <- round(trunc(120 + 0.5 * (age - 60) + 0.4 * (bmi - 27) + rnorm(n, 0, 15), 90, 220))
    chol <- round(exp(rnorm(n, log(5.4), 0.18)), 2)
    diabetes <- factor(ifelse(runif(n) < plogis(-2.6 + 0.04 * (age - 60) + 0.1 * (bmi - 27)), "yes", "no"), c("no", "yes"))
    male <- sex == "male"
    a <- age - 60
    lp <- effect * (0.055 * a + 0.0012 * a^2 + 0.55 * male + 0.25 * (smoking == "former") +
                    0.75 * (smoking == "current") + 0.18 * (sbp - 130) / 10 + 0.9 * log(chol / 5.4) +
                    0.6 * (diabetes == "yes") + 0.025 * (bmi - 27) - 0.06 * a * male)
    # 10-year event (binary)
    event10 <- factor(ifelse(runif(n) < plogis(-2.3 + shift + lp), "yes", "no"), c("no", "yes"))
    # time to event: Weibull proportional hazards, censoring uniform 3-12 years, administrative at 10
    shape <- 1.3; lambda <- 0.012 * exp(shift)
    tev <- (-log(runif(n)) / (lambda * exp(lp)))^(1 / shape)
    cens <- pmin(runif(n, 3, 12), 10)
    time <- round(pmin(tev, cens), 3)
    status <- factor(ifelse(tev <= cens, "event", "censored"), c("censored", "event"))
    # blood pressure at one year (continuous)
    sbp1y <- round(40 + 0.68 * sbp + 0.12 * (age - 60) + 0.5 * (bmi - 27) + 3 * (smoking == "current") +
                   4 * (diabetes == "yes") + rnorm(n, 0, 9))
    data.frame(age, sex, smoking, sbp, chol, bmi, diabetes, event10, sbp1y, time, status)
}

# the later period has a somewhat lower baseline risk (split validation shows it)
dev <- rbind(sim(1600), sim(800, shift = -0.25))
dev$period <- factor(rep(c("2012-2016", "2017-2019"), c(1600, 800)), c("2012-2016", "2017-2019"))
dev$centre <- factor(sample(c("A", "B", "C"), nrow(dev), TRUE, c(0.5, 0.3, 0.2)))
# missing values (complete-case analysis)
dev$bmi[sample(nrow(dev), 70)] <- NA
dev$chol[sample(nrow(dev), 40)] <- NA

# external population: older, more smokers, higher baseline risk, weaker effects
ext <- sim(1000, ageMean = 66, pCurrent = 0.3, shift = 0.5, effect = 0.85)
ext$bmi[sample(nrow(ext), 80)] <- NA
# what the variable check reports: names in another case, levels in another
# case, an unknown level, values outside the development range
names(ext)[names(ext) == "age"] <- "Age"
names(ext)[names(ext) == "sex"] <- "SEX"
ext$SEX <- ifelse(ext$SEX == "male", "Male", "Female")
ext$smoking <- as.character(ext$smoking); ext$smoking[sample(nrow(ext), 12)] <- "unknown"
ext$sbp[sample(which(ext$sbp > 150), 25)] <- sample(225:240, 25, TRUE)

dir.create("examples", showWarnings = FALSE)
jmvReadWrite::write_omv(dev, "data/heart_cohort.omv", frcWrt = TRUE)
utils::write.csv(ext, "examples/heart_external.csv", row.names = FALSE, na = "")
openxlsx::write.xlsx(ext, "examples/heart_external.xlsx")

# ---- published models (fitted on the 2012-2016 patients) ---------------------

d <- dev[dev$period == "2012-2016", ]
d <- d[stats::complete.cases(d), ]

wtab <- function(tab, file, ...) utils::write.csv(tab, file, row.names = FALSE, ...)

# logistic: transforms, levels, a product
g <- glm(event10 ~ age + I(age^2) + sex + smoking + I(sbp / 10) + log(chol) + diabetes + age:sex, binomial, d)
logTerms <- c("Intercept", "age", "age^2", "sex = male", "smoking = former", "smoking = current",
              "sbp/10", "log(chol)", "diabetes = yes", "age * sex = male")
tab <- data.frame(term = logTerms, estimate = signif(unname(coef(g)), 6), SE = signif(sqrt(diag(vcov(g))), 4))
wtab(tab, "examples/published_logistic.csv")
openxlsx::write.xlsx(tab, "examples/published_logistic.xlsx")
# the same in a European spreadsheet export: semicolons, decimal commas, R's
# own term names (sexmale, smokingcurrent) and the Spanish EE header
rTerms <- c("(Intercept)", "age", "age^2", "sexmale", "smokingformer", "smokingcurrent", "sbp/10", "log(chol)",
            "diabetesyes", "age * sexmale")
eu <- data.frame(term = rTerms, coef = format(tab$estimate, decimal.mark = ",", trim = TRUE),
                 EE = format(tab$SE, decimal.mark = ",", trim = TRUE))
utils::write.table(eu, "examples/published_logistic_semicolon.csv", sep = ";", row.names = FALSE, quote = FALSE)

# linear
m <- lm(sbp1y ~ sbp + age + bmi + smoking + diabetes, d)
wtab(data.frame(term = c("Intercept", "sbp", "age", "bmi", "smoking = former", "smoking = current", "diabetes = yes"),
                estimate = signif(unname(coef(m)), 6), SE = signif(sqrt(diag(vcov(m))), 4)),
     "examples/published_linear.csv")

# Cox, as usually published: coefficients, the baseline survival at 5 years for
# the mean covariate values, and that centre
S <- survival::Surv(d$time, d$status == "event")
X <- cbind(age = d$age, male = as.numeric(d$sex == "male"), former = as.numeric(d$smoking == "former"),
           current = as.numeric(d$smoking == "current"), sbp10 = d$sbp / 10, logchol = log(d$chol),
           diab = as.numeric(d$diabetes == "yes"))
fc <- survival::coxph(S ~ X)
b <- coef(fc); mns <- colMeans(X)
# S0 at the mean covariates = exp(-H0(5) exp(b'mean)), H0 the uncentred baseline
bh <- survival::basehaz(fc, centered = FALSE)
s0 <- exp(-max(bh$hazard[bh$time <= 5]) * exp(sum(b * mns)))
coxTerms <- c("age", "sex = male", "smoking = former", "smoking = current", "sbp/10", "log(chol)", "diabetes = yes")
wtab(data.frame(term = c(coxTerms, "baseline survival", "horizon", "centre"),
                estimate = signif(c(unname(b), s0, 5, sum(b * mns)), 7),
                SE = c(signif(sqrt(diag(vcov(fc))), 4), NA, NA, NA)),
     "examples/published_cox.csv", na = "")
wtab(data.frame(term = coxTerms, estimate = signif(unname(b), 7), SE = signif(sqrt(diag(vcov(fc))), 4)),
     "examples/published_cox_nobaseline.csv")

# clinpred JSON exports (spline for age, sbp per 10 units), as the Export button writes
ver <- as.character(utils::packageVersion("clinpred"))
prep <- cp_prepare(d, "event10", covs = c("sbp", "chol", "bmi"), splineVars = "age",
                   factors = c("sex", "smoking", "diabetes"), scales = list(sbp = "ten"),
                   interactions = list(c("age", "sex")))
prep <- cp_binary_outcome(prep, "yes")
fam <- cp_family("binary"); Xs <- cp_design(prep$spec, prep$data); bs <- fam$fit(Xs, prep$y)
writeLines(cp_model_json(prep$spec, bs, "binary", se = sqrt(diag(fam$vcov(bs))), version = ver),
           "examples/heart_logistic_model.json")
prep <- cp_prepare(d, c("time", "status"), covs = c("sbp", "chol"), splineVars = "age",
                   factors = c("sex", "smoking", "diabetes"), scales = list(sbp = "ten"))
prep <- cp_cox_outcome(prep, "time", "status", "event", 5)
fam <- cp_family("cox", horizon = 5); Xs <- cp_design(prep$spec, prep$data); bs <- fam$fit(Xs, prep$y)
writeLines(cp_model_json(prep$spec, bs, "cox", se = sqrt(diag(fam$vcov(bs))), version = ver),
           "examples/heart_cox_model.json")

cat("development:", nrow(dev), "rows,", sum(dev$event10 == "yes"), "events (10 y),",
    sum(dev$status == "event"), "events (time);",
    "external:", nrow(ext), "rows,", sum(ext$event10 == "yes"), "events\n")
