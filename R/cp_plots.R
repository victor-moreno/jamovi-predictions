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

cp_plot_cal <- function(sets, ggtheme) {
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
        ggplot2::labs(x = "Predicted risk", y = "Observed proportion") +
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
        rbind(data.frame(threshold = nb$threshold, nb = nb$model, strategy = "Model", label = s$label),
              data.frame(threshold = nb$threshold, nb = nb$all, strategy = "Treat all", label = s$label),
              data.frame(threshold = nb$threshold, nb = nb$none, strategy = "Treat none", label = s$label))
    }))
    df$label <- factor(df$label, vapply(sets, `[[`, "", "label"))
    df$strategy <- factor(df$strategy, c("Model", "Treat all", "Treat none"))
    top <- max(df$nb[df$strategy != "Treat none"], na.rm = TRUE)
    ggplot2::ggplot(df, ggplot2::aes(.data$threshold, .data$nb, colour = .data$strategy,
                                     linetype = .data$strategy)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::scale_colour_manual(values = c("#0072B2", "#D55E00", "grey40")) +
        ggplot2::scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
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
    p + ggplot2::labs(x = xlab, y = NULL) + .base_theme(ggtheme)
}

cp_plot_shapley <- function(df, measureLabel, ggtheme)
    .bar_importance(df, "shapley", xlab = paste0("Share of ", measureLabel), ggtheme = ggtheme)

cp_plot_perm <- function(df, ggtheme) {
    df$lo <- df$d_auc - df$sd_auc; df$hi <- df$d_auc + df$sd_auc
    .bar_importance(df, "d_auc", "lo", "hi", xlab = "Loss in AUC when permuted (mean ± SD)",
                    ggtheme = ggtheme)
}

# Beeswarm: one row per predictor, points spread vertically by local density,
# coloured by the predictor's value (low to high)
cp_plot_shap_summary <- function(d, ggtheme) {
    ord <- names(sort(tapply(abs(d$phi), d$var, mean)))
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
    ggplot2::ggplot(d, ggplot2::aes(.data$phi, .data$y, colour = .data$value)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
        ggplot2::geom_point(size = 1.2, alpha = 0.8) +
        ggplot2::scale_colour_gradient(low = "#0072B2", high = "#D55E00", limits = c(0, 1),
                                       breaks = c(0, 1), labels = c("Low", "High"),
                                       name = "Predictor value") +
        ggplot2::scale_y_continuous(breaks = seq_along(ord), labels = ord) +
        ggplot2::labs(x = "SHAP value (contribution to the log-odds)", y = NULL) +
        .jtheme(ggtheme) + ggplot2::theme(legend.position = "right")
}

# Waterfall: from the mean log-odds of the development data to this person's
cp_plot_shap_patient <- function(w, ggtheme) {
    o <- order(abs(w$phi), decreasing = TRUE)
    lab <- w$label[o]; phi <- w$phi[o]
    end <- w$base + cumsum(phi); start <- c(w$base, utils::head(end, -1))
    df <- data.frame(label = factor(lab, rev(lab)), start = start, end = end,
                     sign = ifelse(phi > 0, "Raises risk", "Lowers risk"))
    df$y <- as.integer(df$label)
    ggplot2::ggplot(df) +
        ggplot2::geom_vline(xintercept = w$base, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_rect(ggplot2::aes(xmin = pmin(.data$start, .data$end), xmax = pmax(.data$start, .data$end),
                                        ymin = .data$y - 0.35, ymax = .data$y + 0.35, fill = .data$sign)) +
        ggplot2::scale_fill_manual(values = c(`Raises risk` = "#D55E00", `Lowers risk` = "#0072B2")) +
        ggplot2::scale_y_continuous(breaks = df$y, labels = as.character(df$label)) +
        ggplot2::labs(x = sprintf("Log-odds  (mean %.2f → this person %.2f, risk %.1f%%)",
                                  w$base, w$lp, 100 * w$risk), y = NULL,
                      subtitle = w$title) +
        .base_theme(ggtheme)
}

cp_plot_instability <- function(s, ggtheme) {
    pts <- data.frame(p0 = rep(s$p0, ncol(s$pts)) / 1000, pb = as.vector(s$pts) / 1000)
    band <- data.frame(p0 = s$p0 / 1000, lo = s$lo / 1000, hi = s$hi / 1000)
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
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "Risk from the developed model", y = "Risk from bootstrap models",
                      subtitle = sprintf("MAPE %.3f; dashed: 2.5th and 97.5th percentiles", s$mape)) +
        .base_theme(ggtheme)
}

cp_plot_cal_instability <- function(s, ggtheme) {
    df <- data.frame(x = rep(s$grid, ncol(s$curves)), y = as.vector(s$curves),
                     b = rep(seq_len(ncol(s$curves)), each = length(s$grid)))
    lim <- min(1, max(c(s$grid, df$y)) * 1.05)
    ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, group = .data$b)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
        ggplot2::geom_line(alpha = 0.25, colour = "#0072B2") +
        ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
        ggplot2::labs(x = "Predicted risk", y = "Observed proportion (smoothed)",
                      subtitle = "Each line: a bootstrap model applied to the development data") +
        .base_theme(ggtheme)
}
