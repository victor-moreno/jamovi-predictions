# Plots. Each takes the plain data kept in the image state (built in .run) and
# jamovi's ggtheme, so it can be drawn and tested without an analysis object.

.set_colours <- c(dev = "#0072B2", split = "#E69F00", ext = "#009E73")

# jamovi's ggtheme is a list: the theme plus two discrete palette scales, which
# would replace the scales set here, so only the theme is added
.jtheme <- function(ggtheme) if (inherits(ggtheme, "theme")) ggtheme else ggtheme[[1]]

.base_theme <- function(ggtheme) {
    list(.jtheme(ggtheme), ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()))
}

# Points of an ROC curve, thinned to at most maxPts
cp_roc_points <- function(y, p, maxPts = 300) {
    o <- order(p, decreasing = TRUE)
    ys <- y[o]; ps <- p[o]
    tpr <- cumsum(ys) / sum(ys); fpr <- cumsum(1 - ys) / sum(1 - ys)
    last <- c(diff(ps) != 0, TRUE)          # one point per distinct threshold
    tpr <- c(0, tpr[last]); fpr <- c(0, fpr[last])
    if (length(tpr) > maxPts) {
        keep <- unique(round(seq(1, length(tpr), length.out = maxPts)))
        tpr <- tpr[keep]; fpr <- fpr[keep]
    }
    list(fpr = fpr, tpr = tpr)
}

cp_plot_roc <- function(sets, ggtheme) {
    df <- do.call(rbind, lapply(sets, function(s)
        data.frame(fpr = s$fpr, tpr = s$tpr, key = s$key,
                   label = sprintf("%s (AUC %.3f)", s$label, s$auc))))
    df$label <- factor(df$label, levels = unique(df$label))
    cols <- stats::setNames(.set_colours[vapply(sets, `[[`, "", "key")], levels(df$label))
    ggplot2::ggplot(df, ggplot2::aes(.data$fpr, .data$tpr, colour = .data$label)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_step(direction = "vh", linewidth = 0.9) +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "1 − Specificity", y = "Sensitivity") +
        # one legend entry per line: data-set labels are long
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1)) +
        .base_theme(ggtheme)
}

cp_plot_cal <- function(sets, ggtheme, xlab = "Predicted risk", ylab = "Observed proportion") {
    curve <- do.call(rbind, lapply(sets, function(s) data.frame(s$curve, label = s$label)))
    grp <- do.call(rbind, lapply(sets, function(s) data.frame(s$groups, label = s$label)))
    hist <- do.call(rbind, lapply(sets, function(s) data.frame(s$hist, label = s$label)))
    lev <- vapply(sets, `[[`, "", "label")
    curve$label <- factor(curve$label, lev); grp$label <- factor(grp$label, lev)
    hist$label <- factor(hist$label, lev)
    lim <- min(1, max(c(grp$obs_hi, grp$pred, curve$x), na.rm = TRUE) * 1.05)
    hist$h <- hist$count / max(hist$count) * 0.12 * lim
    ggplot2::ggplot() +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_rect(data = hist, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                                     ymin = 0, ymax = .data$h),
                           fill = "grey60", alpha = 0.6) +
        ggplot2::geom_ribbon(data = curve, ggplot2::aes(.data$x, ymin = .data$lo, ymax = .data$hi),
                             fill = "#0072B2", alpha = 0.18) +
        ggplot2::geom_line(data = curve, ggplot2::aes(.data$x, .data$y), colour = "#0072B2",
                           linewidth = 0.9) +
        ggplot2::geom_pointrange(data = grp, ggplot2::aes(.data$pred, .data$obs, ymin = .data$obs_lo,
                                                          ymax = .data$obs_hi), size = 0.3) +
        ggplot2::facet_wrap(~label) +
        ggplot2::coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
        ggplot2::labs(x = xlab, y = ylab) +
        .base_theme(ggtheme)
}

cp_plot_risk <- function(sets, ggtheme) {
    df <- do.call(rbind, lapply(sets, function(s)
        data.frame(p = s$p, outcome = s$outcome, label = s$label)))
    df$label <- factor(df$label, vapply(sets, `[[`, "", "label"))
    ggplot2::ggplot(df, ggplot2::aes(.data$p, fill = .data$outcome)) +
        ggplot2::geom_histogram(bins = 30, position = "identity", alpha = 0.55, boundary = 0) +
        ggplot2::scale_fill_manual(values = c("#999999", "#D55E00")) +
        ggplot2::facet_wrap(~label, scales = "free_y") +
        ggplot2::labs(x = "Predicted risk", y = "Count") +
        .base_theme(ggtheme)
}

cp_plot_dca <- function(sets, ggtheme) {
    df <- do.call(rbind, lapply(sets, function(s) {
        nb <- s$nb
        rbind(data.frame(threshold = nb$threshold, nb = nb$model, strategy = if (is.null(nb$base)) "Model" else "Full model", label = s$label),
              # the baseline model, when comparing models
              if (!is.null(nb$base)) data.frame(threshold = nb$threshold, nb = nb$base, strategy = "Baseline model", label = s$label),
              data.frame(threshold = nb$threshold, nb = nb$all, strategy = "Treat all", label = s$label),
              data.frame(threshold = nb$threshold, nb = nb$none, strategy = "Treat none", label = s$label))
    }))
    df$label <- factor(df$label, vapply(sets, `[[`, "", "label"))
    strat <- intersect(c("Model", "Full model", "Baseline model", "Treat all", "Treat none"), unique(df$strategy))
    df$strategy <- factor(df$strategy, strat)
    top <- max(df$nb[df$strategy != "Treat none"], na.rm = TRUE)
    ggplot2::ggplot(df, ggplot2::aes(.data$threshold, .data$nb, colour = .data$strategy,
                                     linetype = .data$strategy)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::scale_colour_manual(values = c(Model = "#0072B2", `Full model` = "#0072B2", `Baseline model` = "#009E73",
                                                `Treat all` = "#D55E00", `Treat none` = "grey40")) +
        ggplot2::scale_linetype_manual(values = c(Model = "solid", `Full model` = "solid", `Baseline model` = "twodash",
                                                  `Treat all` = "longdash", `Treat none` = "dotted")) +
        ggplot2::facet_wrap(~label) +
        ggplot2::coord_cartesian(ylim = c(-0.1 * top, top * 1.05)) +
        ggplot2::labs(x = "Threshold probability", y = "Net benefit") +
        .base_theme(ggtheme)
}

.bar_importance <- function(df, value, lo = NULL, hi = NULL, xlab, ggtheme) {
    df$var <- factor(df$var, df$var[order(df[[value]])])
    p <- ggplot2::ggplot(df, ggplot2::aes(.data[[value]], .data$var)) +
        ggplot2::geom_col(fill = "#0072B2", width = 0.7)
    if (!is.null(lo))
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data[[lo]], xmax = .data[[hi]]),
                                        width = 0.25, orientation = "y")
    .baseline_panel(p + ggplot2::labs(x = xlab, y = NULL) + .base_theme(ggtheme), df)
}

# baseline covariates (column forced) in a panel of their own, below the
# other predictors, as in the forest plot
.baseline_panel <- function(p, df) {
    if (is.null(df$forced) || all(df$forced) || !any(df$forced)) return(p)
    p$data$group <- factor(ifelse(df$forced, "Baseline covariates", "Predictors"), c("Predictors", "Baseline covariates"))
    p + ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y")
}

cp_plot_shapley <- function(df, measureLabel, ggtheme)
    .bar_importance(df, "shapley", xlab = paste0("Share of ", measureLabel), ggtheme = ggtheme)

cp_plot_perm <- function(df, ggtheme, key = "auc", xlab = "Loss in AUC when permuted (mean ± SD)") {
    d <- paste0("d_", key); s <- paste0("sd_", key)
    df$lo <- df[[d]] - df[[s]]; df$hi <- df[[d]] + df[[s]]
    .bar_importance(df, d, "lo", "hi", xlab = xlab, ggtheme = ggtheme)
}

# Beeswarm: one row per predictor, points spread vertically by local density,
# coloured by the predictor's value (low to high)
cp_plot_shap_summary <- function(d, ggtheme, xlab = "SHAP value (contribution to the log-odds)") {
    imp <- tapply(abs(d$phi), d$var, mean)
    # baseline covariates at the bottom, in their own panel
    fv <- if (is.null(d$forced)) character(0) else unique(d$var[d$forced])
    ord <- names(imp)[order(!names(imp) %in% fv, imp)]
    d$var <- factor(d$var, ord)
    d$offset <- 0
    for (v in ord) {
        i <- which(d$var == v)
        if (length(i) > 1 && stats::sd(d$phi[i]) > 0) {
            dens <- stats::density(d$phi[i], n = 256)
            h <- stats::approx(dens$x, dens$y, d$phi[i])$y
            # deterministic spread (golden-ratio sequence): the plot does not change between renders
            jit <- ((seq_along(i) * 0.618034) %% 1) * 2 - 1
            d$offset[i] <- jit * h / max(dens$y) * 0.4
        }
    }
    d$y <- as.integer(d$var) + d$offset
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$phi, .data$y, colour = .data$value)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
        ggplot2::geom_point(size = 1.2, alpha = 0.8) +
        ggplot2::scale_colour_gradient(low = "#0072B2", high = "#D55E00", limits = c(0, 1),
                                       breaks = c(0, 1), labels = c("Low", "High"),
                                       name = "Predictor value") +
        ggplot2::scale_y_continuous(breaks = seq_along(ord), labels = ord) +
        ggplot2::labs(x = xlab, y = NULL) +
        ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0.5, barheight = ggplot2::unit(0.5, "lines"),
                                                          barwidth = ggplot2::unit(10, "lines"))) +
        .jtheme(ggtheme) + ggplot2::theme(legend.position = "bottom")
    .baseline_panel(p, d)
}

# Waterfall: from the mean log-odds of the development data to this person's
cp_plot_shap_patient <- function(w, ggtheme) {
    o <- order(abs(w$phi), decreasing = TRUE)
    lab <- w$label[o]; phi <- w$phi[o]
    end <- w$base + cumsum(phi); start <- c(w$base, utils::head(end, -1))
    # legend labels: c(lowers, raises), for risks unless the state says otherwise
    sg <- w$signs %||% c("Lowers risk", "Raises risk")
    df <- data.frame(label = factor(lab, rev(lab)), start = start, end = end,
                     sign = ifelse(phi > 0, sg[2], sg[1]))
    df$y <- as.integer(df$label)
    ggplot2::ggplot(df) +
        ggplot2::geom_vline(xintercept = w$base, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_rect(ggplot2::aes(xmin = pmin(.data$start, .data$end), xmax = pmax(.data$start, .data$end),
                                        ymin = .data$y - 0.35, ymax = .data$y + 0.35, fill = .data$sign)) +
        ggplot2::scale_fill_manual(values = stats::setNames(c("#0072B2", "#D55E00"), sg)) +
        ggplot2::scale_y_continuous(breaks = df$y, labels = as.character(df$label)) +
        # axis: the scale; subtitle: this person (too long for the axis at jamovi's text size)
        ggplot2::labs(x = (w$axis %||% "Log-odds")[1], y = NULL,
                      subtitle = paste0(w$title, ": ", if (length(w$axis) > 1) w$axis[2] else
                          sprintf("mean %.2f → this person %.2f, risk %.1f%%", w$base, w$lp, 100 * w$risk))) +
        .base_theme(ggtheme)
}

# predictions are stored on a 1/1000 scale of s$lim (see .stability)
cp_plot_instability <- function(s, ggtheme, xlab = "Risk from the developed model",
                                ylab = "Risk from bootstrap models") {
    lim <- s$lim %||% c(0, 1)
    dec <- function(v) lim[1] + v / 1000 * (lim[2] - lim[1])
    pts <- data.frame(p0 = dec(rep(s$p0, ncol(s$pts))), pb = dec(as.vector(s$pts)))
    band <- data.frame(p0 = dec(s$p0), lo = dec(s$lo), hi = dec(s$hi))
    band <- band[order(band$p0), ]
    sm <- function(v) stats::lowess(band$p0, v, f = 0.3)$y
    band$lo_s <- sm(band$lo); band$hi_s <- sm(band$hi)
    ggplot2::ggplot() +
        ggplot2::geom_point(data = pts, ggplot2::aes(.data$p0, .data$pb), alpha = 0.15, size = 0.6,
                            colour = "grey30") +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "#0072B2") +
        ggplot2::geom_line(data = band, ggplot2::aes(.data$p0, .data$lo_s), linetype = "dashed",
                           colour = "#D55E00") +
        ggplot2::geom_line(data = band, ggplot2::aes(.data$p0, .data$hi_s), linetype = "dashed",
                           colour = "#D55E00") +
        ggplot2::coord_equal(xlim = lim, ylim = lim) +
        ggplot2::labs(x = xlab, y = ylab,
                      subtitle = sprintf("MAPE %.3f; dashed: 2.5th and 97.5th percentiles", s$mape)) +
        .base_theme(ggtheme)
}

cp_plot_cal_instability <- function(s, ggtheme, xlab = "Predicted risk",
                                    ylab = "Observed proportion (smoothed)") {
    df <- data.frame(x = rep(s$grid, ncol(s$curves)), y = as.vector(s$curves),
                     b = rep(seq_len(ncol(s$curves)), each = length(s$grid)))
    lim <- if (identical(s$lim %||% c(0, 1), c(0, 1))) c(0, min(1, max(c(s$grid, df$y)) * 1.05))
        else range(s$grid, df$y)
    ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, group = .data$b)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_line(alpha = 0.25, colour = "#0072B2") +
        ggplot2::coord_equal(xlim = lim, ylim = lim) +
        ggplot2::labs(x = xlab, y = ylab,
                      subtitle = "Each line: a bootstrap model applied to the development data") +
        .base_theme(ggtheme)
}

# Continuous outcome: observed against predicted, loess with its CI, identity line
cp_plot_cal_cont <- function(sets, ggtheme) {
    curve <- do.call(rbind, lapply(sets, function(s) data.frame(s$curve, label = s$label)))
    pts <- do.call(rbind, lapply(sets, function(s) data.frame(s$points, label = s$label)))
    lev <- vapply(sets, `[[`, "", "label")
    curve$label <- factor(curve$label, lev); pts$label <- factor(pts$label, lev)
    lim <- range(pts$pred, pts$obs)
    ggplot2::ggplot() +
        ggplot2::geom_point(data = pts, ggplot2::aes(.data$pred, .data$obs), alpha = 0.3, size = 0.9,
                            colour = "grey40") +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_ribbon(data = curve, ggplot2::aes(.data$x, ymin = .data$lo, ymax = .data$hi),
                             fill = "#0072B2", alpha = 0.18) +
        ggplot2::geom_line(data = curve, ggplot2::aes(.data$x, .data$y), colour = "#0072B2",
                           linewidth = 0.9) +
        ggplot2::facet_wrap(~label) +
        ggplot2::coord_cartesian(xlim = lim, ylim = lim) +
        ggplot2::labs(x = "Predicted value", y = "Observed value") +
        .base_theme(ggtheme)
}

# Kaplan–Meier curves by risk group, one facet per data set; the dots at the
# horizon are each group's mean predicted event-free probability
cp_plot_km <- function(st, ggtheme) {
    curves <- do.call(rbind, lapply(st$sets, function(s) data.frame(s$curves, label = s$label)))
    pred <- do.call(rbind, lapply(st$sets, function(s) data.frame(s$pred, label = s$label)))
    lev <- vapply(st$sets, `[[`, "", "label")
    curves$label <- factor(curves$label, lev); pred$label <- factor(pred$label, lev)
    ng <- length(unique(curves$group))
    cols <- grDevices::colorRampPalette(c("#0072B2", "#E69F00", "#D55E00"))(ng)
    curves$group <- factor(curves$group, unique(curves$group)); pred$group <- factor(pred$group, levels(curves$group))
    ggplot2::ggplot(curves, ggplot2::aes(.data$time, .data$surv, colour = .data$group)) +
        ggplot2::geom_vline(xintercept = st$h, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_step(linewidth = 0.8) +
        ggplot2::geom_point(data = pred, size = 2.5, shape = 18) +
        ggplot2::scale_colour_manual(values = cols, name = "Risk group") +
        ggplot2::facet_wrap(~label) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = "Time", y = "Event-free (Kaplan–Meier)",
                      subtitle = "Diamonds: mean predicted at the horizon") +
        .base_theme(ggtheme) +
        # jamovi's 16 pt theme: the group labels need two legend rows
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
        ggplot2::theme(legend.title = ggplot2::element_text())
}

# Cross-validated deviance along the lambda path, with 1 SE bars; dashed lines
# at the minimum and at the largest lambda within 1 SE, the solid one at the
# lambda chosen; numbers of non-zero coefficients on top
cp_plot_enet_cv <- function(st, ggtheme) {
    df <- data.frame(ll = log(st$lambda), cvm = st$cvm, lo = st$cvm - st$cvsd, hi = st$cvm + st$cvsd)
    lab <- df[unique(round(seq(1, nrow(df), length.out = 8))), ]
    lab$nz <- st$nz[as.integer(rownames(lab))]
    top <- max(df$hi, na.rm = TRUE)
    ggplot2::ggplot(df, ggplot2::aes(.data$ll, .data$cvm)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0, colour = "grey60") +
        ggplot2::geom_point(colour = "#D55E00", size = 1.4) +
        ggplot2::geom_vline(xintercept = log(st$lambda[c(st$imin, st$i1se)]), linetype = "dashed", colour = "grey50") +
        ggplot2::geom_vline(xintercept = log(st$chosen), colour = "#0072B2") +
        ggplot2::geom_text(data = lab, ggplot2::aes(.data$ll, top, label = .data$nz), vjust = -0.3,
                           # fallback on the size itself: %||% binds tighter than * and /, so it never applied
                           size = 0.8 * (ggplot2::calc_element("axis.text", .jtheme(ggtheme))$size %||% 11) / ggplot2::.pt) +
        ggplot2::labs(x = "log(λ)", y = st$ylab, subtitle = "Numbers: non-zero coefficients") +
        .jtheme(ggtheme)
}

# ROC curves of the baseline and the full model, one panel per data set
cp_plot_inc_roc <- function(sets, ggtheme) {
    df <- do.call(rbind, lapply(sets, function(s) rbind(
        data.frame(fpr = s$base$fpr, tpr = s$base$tpr, model = "Baseline model", label = s$label),
        data.frame(fpr = s$full$fpr, tpr = s$full$tpr, model = "Full model", label = s$label))))
    df$label <- factor(df$label, vapply(sets, `[[`, "", "label"))
    df$model <- factor(df$model, c("Full model", "Baseline model"))
    aucs <- do.call(rbind, lapply(sets, function(s) data.frame(label = s$label,
        text = sprintf("AUC full %.3f\nbaseline %.3f", s$full$auc, s$base$auc))))
    aucs$label <- factor(aucs$label, levels(df$label))
    th <- .jtheme(ggtheme)
    ggplot2::ggplot(df, ggplot2::aes(.data$fpr, .data$tpr, colour = .data$model, linetype = .data$model)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_step(direction = "vh", linewidth = 0.9) +
        ggplot2::geom_text(data = aucs, ggplot2::aes(1, 0, label = .data$text), inherit.aes = FALSE,
                           hjust = 1, vjust = 0, size = 0.8 * (ggplot2::calc_element("axis.text", th)$size %||% 11) / ggplot2::.pt,
                           colour = ggplot2::calc_element("axis.text", th)$colour %||% "grey30") +
        ggplot2::scale_colour_manual(values = c(`Full model` = "#0072B2", `Baseline model` = "#009E73")) +
        ggplot2::scale_linetype_manual(values = c(`Full model` = "solid", `Baseline model` = "twodash")) +
        ggplot2::facet_wrap(~label) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "1 − Specificity", y = "Sensitivity") +
        .base_theme(ggtheme)
}

# Forest plot of the effects table: estimates with CIs, one colour per model
# (crude, adjusted, selected), baseline covariates in their own panel below
cp_plot_forest <- function(st, blocks, ggtheme) {
    df <- st$rows[st$rows$block %in% blocks, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    lev <- intersect(c("Crude", "Adjusted", "Adjusted (all predictors)", "Selected model"), unique(df$block))
    df$block <- factor(df$block, lev)
    k <- length(lev); off <- if (k > 1) (seq_len(k) - (k + 1) / 2) * 0.22 else 0
    # one row per label, top to bottom in table order; a gap between the panels
    df$y <- match(df$label, rev(st$order)) + off[as.integer(df$block)]
    df$group <- factor(ifelse(df$forced, "Baseline covariates", "Predictors"), c("Predictors", "Baseline covariates"))
    br <- unique(match(df$label, rev(st$order)))
    cols <- c(Crude = "grey45", Adjusted = "#0072B2", `Adjusted (all predictors)` = "#0072B2", `Selected model` = "#D55E00")
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$est, .data$y, colour = .data$block)) +
        ggplot2::geom_vline(xintercept = if (st$exp) 1 else 0, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), width = 0.15, orientation = "y", na.rm = TRUE) +
        ggplot2::geom_point(size = 2.2) +
        ggplot2::scale_colour_manual(values = cols[lev]) +
        ggplot2::scale_y_continuous(breaks = br, labels = rev(st$order)[br], expand = ggplot2::expansion(add = 0.6)) +
        ggplot2::labs(x = st$xlab, y = NULL) +
        .base_theme(ggtheme)
    if (st$exp) p <- p + ggplot2::scale_x_log10()
    if (any(df$forced) && !all(df$forced))
        p <- p + ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y")
    p
}
