##############################################################################
## Manuscript figures (phontrast 2.1.0)
##
## Reads the pairwise metric tables written by rerun_paper_metrics_phontrast210.R
## and renders the manuscript figures. Run from the repository root AFTER the
## rerun script:  Rscript manuscript/make_figures.R
##
## Requires: ggplot2, scales. Reads from manuscript/outputs/, writes PNG (300 dpi)
## and PDF into manuscript/outputs/figures/.
##
## Figures:
##   fig1_dimensionality_slopegraph  -- headline: correspondence of each metric
##       with KDE overlap in 2-D (PB52) vs 13-D (SBCSAE), showing JSD holding
##       while Pillai and Bhattacharyya degrade with dimensionality.
##   fig2_mfcc_correspondence_panels -- per-metric scatter vs (1 - KDE overlap)
##       in the 13-D space, the evidence behind fig1's right-hand column.
##   fig3_estimator_robustness       -- mc vs legacy JSD across the 13-D pairs,
##       showing the 2.1.0 correction preserves the ranking.
##############################################################################
suppressPackageStartupMessages({library(ggplot2); library(scales)})

BASE <- if (dir.exists("manuscript")) "manuscript" else "."
OUT  <- file.path(BASE, "outputs")
FIG  <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## Okabe-Ito, CVD-safe; fixed order, color follows the metric (never its rank).
PAL <- c(JSD = "#0072B2", Pillai = "#E69F00", Bhattacharyya = "#009E73")
INK <- "#333333"; INK2 <- "#6b6b6b"; GRID <- "#e9e9e6"

theme_pub <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = GRID, linewidth = 0.4),
      axis.title = element_text(color = INK), axis.text = element_text(color = INK2),
      plot.title = element_text(color = INK, face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = INK2, size = rel(0.9)),
      plot.caption = element_text(color = INK2, size = rel(0.75), hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      strip.text = element_text(color = INK, face = "bold"),
      legend.position = "none", plot.margin = margin(12, 16, 10, 12)
    )
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ok <- tryCatch({ ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h,
                          device = cairo_pdf, bg = "white"); TRUE }, error = function(e) FALSE)
  if (!ok) ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  message("wrote ", name, ".{png,pdf}")
}
sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

read_tab <- function(f) {
  p <- file.path(OUT, f)
  if (!file.exists(p)) { message("missing ", p, " -- run rerun_paper_metrics first"); return(NULL) }
  read.csv(p)
}
pb52 <- read_tab("pb52_pairwise_phontrast210.csv")
mfcc <- read_tab("mfcc13_pairwise_phontrast210.csv")

## ---- rho summary (correspondence of each metric with 1 - KDE overlap, mc) ---
rho_row <- function(tab, space, dim_lab) {
  if (is.null(tab)) return(NULL)
  sep <- 1 - tab$ovl_mc
  data.frame(space = space, dim_lab = dim_lab,
             metric = names(PAL),
             rho = c(sp(tab$jsd_mc, sep), sp(tab$pillai, sep), sp(tab$bhatt_dist, sep)))
}
rho <- rbind(rho_row(pb52, "PB52 F1/F2", "2-D"), rho_row(mfcc, "SBCSAE MFCC", "13-D"))

## ---- fig 1: dimensionality slopegraph --------------------------------------
if (!is.null(rho) && length(unique(rho$dim_lab)) == 2) {
  rho$dim_lab <- factor(rho$dim_lab, levels = c("2-D", "13-D"))
  rho$metric  <- factor(rho$metric, levels = names(PAL))
  lab_r <- subset(rho, dim_lab == "13-D"); lab_l <- subset(rho, dim_lab == "2-D")
  # The three 2-D values are near-coincident (all metrics tie in 2-D), so label
  # the cluster once instead of overprinting three values.
  pb_lbl <- sprintf("'all metrics' %%~~%% %.2f", mean(lab_l$rho))
  p1 <- ggplot(rho, aes(dim_lab, rho, color = metric, group = metric)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2.8) +
    geom_text(data = lab_r, aes(label = sprintf("%s  %.3f", metric, rho)),
              hjust = 0, nudge_x = 0.04, size = 3.5, fontface = "bold") +
    annotate("text", x = 1, y = 1.02, label = pb_lbl, parse = TRUE, hjust = 0.5,
             vjust = 0, size = 3.3, color = INK2) +
    scale_color_manual(values = PAL) +
    scale_x_discrete(expand = expansion(mult = c(0.18, 0.42)),
                     labels = c("2-D" = "2-D\nPB52 F1/F2", "13-D" = "13-D\nSBCSAE MFCC")) +
    scale_y_continuous(limits = c(min(rho$rho) - 0.03, 1.04), breaks = seq(0.75, 1, 0.05)) +
    labs(title = "Correspondence with KDE overlap holds for JSD as dimensionality grows",
         subtitle = "Spearman rho between each separation metric and (1 - KDE overlap); Monte-Carlo estimator",
         x = "Acoustic space", y = expression(Spearman~rho),
         caption = "Okabe-Ito CVD-safe palette. Higher = tighter correspondence with distributional overlap.") +
    theme_pub()
  save_fig(p1, "fig1_dimensionality_slopegraph", 7.4, 5.0)
}

## ---- fig 2: MFCC correspondence panels -------------------------------------
if (!is.null(mfcc)) {
  sep <- 1 - mfcc$ovl_mc
  long <- rbind(
    data.frame(metric = "JSD",           x = sep, y = mfcc$jsd_mc),
    data.frame(metric = "Pillai",        x = sep, y = mfcc$pillai),
    data.frame(metric = "Bhattacharyya", x = sep, y = mfcc$bhatt_dist))
  long$metric <- factor(long$metric, levels = names(PAL))
  ann <- do.call(rbind, lapply(names(PAL), function(m) {
    d <- long[long$metric == m, ]; data.frame(metric = factor(m, levels = names(PAL)),
      x = min(d$x, na.rm = TRUE), y = max(d$y, na.rm = TRUE),
      lab = sprintf("rho == %.3f", sp(d$y, d$x))) }))
  p2 <- ggplot(long, aes(x, y, color = metric)) +
    geom_point(alpha = 0.45, size = 0.9) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, color = INK) +
    geom_text(data = ann, aes(label = lab), parse = TRUE, hjust = 0, vjust = 1,
              size = 3.6, fontface = "bold", color = INK) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = PAL) +
    labs(title = "Separation-metric vs. distributional overlap in 13-D MFCC space",
         subtitle = "Each point is one of 351 vowel pairs; JSD tracks overlap most tightly",
         x = "1 - KDE overlap  (greater separation)", y = "Metric value",
         caption = "SBCSAE 13-MFCC. Monte-Carlo JSD and percent-overlap; Pillai and Bhattacharyya are parametric.") +
    theme_pub() + theme(strip.text = element_text(size = rel(1.05)))
  save_fig(p2, "fig2_mfcc_correspondence_panels", 9.2, 3.6)
}

## ---- fig 3: estimator robustness (mc vs legacy) ----------------------------
if (!is.null(mfcc)) {
  r <- sp(mfcc$jsd_mc, mfcc$jsd_legacy)
  lim <- range(c(mfcc$jsd_mc, mfcc$jsd_legacy), na.rm = TRUE) + c(-0.02, 0.02)
  p3 <- ggplot(mfcc, aes(jsd_legacy, jsd_mc)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", color = INK2) +
    geom_point(alpha = 0.5, size = 1.0, color = PAL[["JSD"]]) +
    annotate("text", x = lim[1], y = lim[2],
             label = sprintf("Spearman~rho == %.2f", r), parse = TRUE,
             hjust = 0, vjust = 1, size = 3.8, fontface = "bold", color = INK) +
    coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
    labs(title = "The 2.1.0 correction preserves the JSD ranking",
         subtitle = "Monte-Carlo vs. legacy JSD, 351 SBCSAE vowel pairs",
         x = "JSD, legacy estimator", y = "JSD, Monte-Carlo estimator (2.1.0)",
         caption = "Dashed line: y = x. Ranks agree; the Monte-Carlo estimator compresses the high end.") +
    theme_pub()
  save_fig(p3, "fig3_estimator_robustness", 6.0, 6.2)
}

## ---- fig 4: estimator-influence audit (fairness centerpiece) ---------------
corr <- read_tab("estimator_influence_correspondence.csv")
if (!is.null(corr)) {
  d <- corr[corr$metric %in% c("JSD (mc)", "Pillai", "Bhattacharyya"), ]
  d$metric <- factor(ifelse(d$metric == "JSD (mc)", "JSD", d$metric), levels = names(PAL))
  grp <- c("KDE (JSD)"          = "KDE\n(JSD's\nestimator)",
           "neutral"            = "neutral\n(no shared\nassumptions)",
           "MVN (Pillai/Bhatt)" = "MVN\n(Pillai /\nBhatt)")
  d$grp <- factor(grp[d$shares_machinery], levels = grp)
  ylev <- c("KDE density", "kNN density k=10", "kNN density k=25", "kNN density k=50",
            "PCA-2D grid", "PCA-2D convex hull", "MVN Monte-Carlo")
  d$yardstick <- factor(d$yardstick, levels = rev(ylev))
  dodge <- position_dodge(width = 0.62)
  p4 <- ggplot(d, aes(abs_spearman, yardstick, color = metric)) +
    geom_errorbarh(aes(xmin = spearman_lo, xmax = spearman_hi), height = 0.28,
                   position = dodge, linewidth = 0.5, alpha = 0.55) +
    geom_point(size = 2.7, position = dodge) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = PAL) +
    scale_x_continuous(limits = c(0.70, 1.0), breaks = seq(0.75, 1, 0.05)) +
    coord_cartesian(clip = "off") +
    labs(title = "The best-matching metric depends on the overlap yardstick",
         subtitle = "|Spearman rho| with each overlap estimator (13-D MFCC, 351 pairs, Monte-Carlo JSD)",
         x = "|Spearman rho|  (separation metric vs. overlap)", y = NULL, color = NULL,
         caption = paste("95% bootstrap intervals over vowel pairs. Each metric peaks on its own-assumption yardstick;",
                         "on neutral references they are comparable.")) +
    theme_pub() +
    theme(legend.position = "top", legend.text = element_text(color = INK),
          strip.placement = "outside", panel.spacing.y = unit(6, "pt"),
          strip.text.y.left = element_text(angle = 0, face = "bold", size = rel(0.8), color = INK))
  save_fig(p4, "fig4_estimator_influence_audit", 8.8, 5.8)
}

message("Figures written to ", FIG)
