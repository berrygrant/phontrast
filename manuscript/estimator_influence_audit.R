##############################################################################
## Estimator-influence audit (phontrast 2.1.0), self-contained in R
##
## Fairness check for the metric comparison: JSD is estimated with KDE, and so
## is phontrast's percent-overlap, so correlating JSD against *KDE* overlap is
## partly circular. This audit correlates each separation metric against a panel
## of overlap yardsticks that do NOT share JSD's machinery -- kNN density, MVN
## Monte-Carlo, PCA-2D grid / convex hull -- for BOTH acoustic spaces (2-D PB52
## F1/F2 and 13-D SBCSAE MFCC), recomputed with the corrected 2.1.0 Monte-Carlo
## JSD. Mirrors the LabPhon poster's overlap-estimator robustness analysis; the
## neutral estimators are ported from its reproduce_poster_figures.py
## (overlap_estimators.R) and validated to reproduce its values (rho > 0.98).
##
## Prerequisite: rerun_paper_metrics_phontrast210.R (writes the metric tables).
## Inputs (manuscript/data/): pb52.csv, rerun_mfcc_sampled_tokens_20260624.csv.
## Output: manuscript/outputs/estimator_influence_correspondence.csv
##############################################################################
BASE <- if (dir.exists("manuscript")) "manuscript" else "."
OUT  <- file.path(BASE, "outputs"); DATA <- file.path(BASE, "data")
source(file.path(BASE, "overlap_estimators.R"))

canon <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "")
datasets <- list(
  list(name = "PB52 F1/F2", dim = "2-D", tokens = "pb52.csv",
       feats = c("f1", "f2"), metrics = "pb52_pairwise_phontrast210.csv"),
  list(name = "SBCSAE MFCC", dim = "13-D", tokens = "rerun_mfcc_sampled_tokens_20260624.csv",
       feats = paste0("mfcc", 1:13), metrics = "mfcc13_pairwise_phontrast210.csv"))

## Per-pair metric + overlap table for one dataset.
pair_table <- function(ds) {
  mp <- file.path(OUT, ds$metrics); tp <- file.path(DATA, ds$tokens)
  if (!file.exists(mp)) { message("skip ", ds$name, ": missing ", mp, " (run rerun first)"); return(NULL) }
  if (!file.exists(tp)) { message("skip ", ds$name, ": missing ", tp); return(NULL) }
  met <- read.csv(mp); tok <- read.csv(tp); tok$vowel <- as.character(tok$vowel)
  met$.k <- canon(met$v1, met$v2)
  vs <- sort(unique(tok$vowel)); pairs <- t(combn(vs, 2))
  rows <- vector("list", nrow(pairs)); t0 <- Sys.time()
  for (i in seq_len(nrow(pairs))) {
    v1 <- pairs[i, 1]; v2 <- pairs[i, 2]
    d0 <- as.matrix(tok[tok$vowel == v1, ds$feats]); d1 <- as.matrix(tok[tok$vowel == v2, ds$feats])
    ov <- neutral_overlaps(d0, d1, seed = sum(utf8ToInt(paste0(v1, v2))) + 12345L)
    ov$.k <- canon(v1, v2); rows[[i]] <- ov
    if (i %% 25 == 0) message(ds$name, ": ", i, "/", nrow(pairs), " (",
                              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "m)")
  }
  neu <- do.call(rbind, rows)
  m <- merge(met[, c(".k", "jsd_mc", "jsd_legacy", "pillai", "bhatt_dist", "ovl_mc")], neu, by = ".k")
  if (nrow(m) < nrow(met) * 0.99) stop(ds$name, ": metric/overlap merge lost pairs")
  m$dataset <- ds$name; m$dim <- ds$dim; m
}

## Yardsticks (overlap columns) and whose machinery they share.
yard <- data.frame(
  label = c("KDE density", "kNN density k=10", "kNN density k=25", "kNN density k=50",
            "MVN Monte-Carlo", "PCA-2D grid", "PCA-2D convex hull"),
  col   = c("ovl_mc", "overlap_knn_density_k10", "overlap_knn_density_k25",
            "overlap_knn_density_k50", "overlap_mvn_mc", "overlap_pca2_grid", "overlap_pca2_hull"),
  shares = c("KDE (JSD)", "neutral", "neutral", "neutral", "MVN (Pillai/Bhatt)", "neutral", "neutral"),
  stringsAsFactors = FALSE)
mets <- c(jsd_mc = "JSD (mc)", jsd_legacy = "JSD (legacy)", pillai = "Pillai", bhatt_dist = "Bhattacharyya")

set.seed(2026); B <- 2000
boot_ci <- function(x, o) { n <- length(x); v <- numeric(B)
  for (b in seq_len(B)) { i <- sample.int(n, n, replace = TRUE)
    v[b] <- suppressWarnings(abs(cor(x[i], o[i], method = "spearman"))) }
  quantile(v, c(0.025, 0.975), na.rm = TRUE) }

correspondence <- function(m) {
  rows <- list()
  for (j in seq_len(nrow(yard))) {
    o <- m[[yard$col[j]]]
    for (mc_col in names(mets)) {
      x <- m[[mc_col]]; ok <- is.finite(x) & is.finite(o); if (sum(ok) < 5) next
      ci <- boot_ci(x[ok], o[ok])
      rows[[length(rows) + 1]] <- data.frame(
        dataset = m$dataset[1], dim = m$dim[1], yardstick = yard$label[j],
        shares_machinery = yard$shares[j], metric = mets[[mc_col]], n_pairs = sum(ok),
        abs_spearman = round(abs(cor(x[ok], o[ok], method = "spearman")), 4),
        spearman_lo = round(ci[1], 4), spearman_hi = round(ci[2], 4),
        abs_pearson = round(abs(cor(x[ok], o[ok], method = "pearson")), 4))
    }
  }
  do.call(rbind, rows)
}

tabs <- lapply(datasets, function(ds) { m <- pair_table(ds); if (is.null(m)) NULL else correspondence(m) })
tab <- do.call(rbind, Filter(Negate(is.null), tabs))
if (is.null(tab)) stop("no datasets available")
write.csv(tab, file.path(OUT, "estimator_influence_correspondence.csv"), row.names = FALSE)

for (dn in unique(tab$dataset)) {
  cat("\n== ", dn, " (|Spearman| of metric with overlap yardstick, mc JSD) ==\n", sep = "")
  s <- tab[tab$dataset == dn & tab$metric %in% c("JSD (mc)", "Pillai", "Bhattacharyya"), ]
  for (yl in unique(s$yardstick)) {
    r <- s[s$yardstick == yl, ]; win <- r$metric[which.max(r$abs_spearman)]
    cat(sprintf("  %-20s [%-18s]  JSD %.3f  Pillai %.3f  Bhatt %.3f  -> %s\n",
        yl, r$shares_machinery[1],
        r$abs_spearman[r$metric == "JSD (mc)"], r$abs_spearman[r$metric == "Pillai"],
        r$abs_spearman[r$metric == "Bhattacharyya"], win))
  }
  neu <- s[s$shares_machinery == "neutral", ]
  am <- tapply(neu$abs_spearman, neu$metric, mean)
  cat(sprintf("  neutral avg: JSD %.3f  Pillai %.3f  Bhatt %.3f  -> %s\n",
      am[["JSD (mc)"]], am[["Pillai"]], am[["Bhattacharyya"]], names(which.max(am))))
}
message("\nwrote estimator_influence_correspondence.csv")
