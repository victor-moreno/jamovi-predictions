# results in English whatever the machine locale (the catalogues translate them otherwise)
Sys.setenv(LANGUAGE = "en")
pima_dev <- function() MASS::Pima.tr
pima_val <- function() MASS::Pima.te
pima_prep <- function(splines = "glu", factors = character(0), data = pima_dev()) {
    prep <- cp_prepare(data, "type", covs = c("npreg", "bp", "bmi", "ped"), splineVars = splines,
                       factors = factors, nknots = 4)
    cp_binary_outcome(prep, "Yes")
}
pima_stacked <- function() {
    d <- rbind(cbind(MASS::Pima.tr, sample = "train"), cbind(MASS::Pima.te, sample = "test"))
    d$sample <- factor(d$sample, c("train", "test"))
    d$agegrp <- cut(d$age, c(0, 30, 45, 100), labels = c("<30", "30-45", "45+"))
    d
}
enable_output <- function(opts, name) {
    p <- opts$.__enclos_env__$private
    p$.options[[name]]$value <- list(value = TRUE, synced = TRUE)
    p$.env[[name]] <- TRUE
    opts
}
