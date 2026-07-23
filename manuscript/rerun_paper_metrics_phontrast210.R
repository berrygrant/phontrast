##############################################################################
## Manuscript metric rerun with phontrast 2.1.0
##
## Regenerates the paper's metric tables with the corrected Monte-Carlo
## estimator (method = "mc", the phontrast 2.1.0 default) alongside
## method = "legacy" (the estimator the 2026-06 LabPhon poster tables were
## built with), so every table carries both columns and the text can cite
## either. Run from the repository root:  Rscript manuscript/rerun_paper_metrics_phontrast210.R
##
## Requires: phontrast (>= 2.1.0), ks, and (for the PB52 block) phonTools.
## Inputs (place under manuscript/data/, which is git-ignored -- see README):
##   - rerun_mfcc_sampled_tokens_20260624.csv  (SBCSAE 13-MFCC sampled tokens;
##       a copy is committed on the labphon_2026 branch under
##       analysis/labphon_2026/data/)
## PB52 F1/F2 is loaded directly from the phonTools package, no file needed.
## Outputs are written to manuscript/outputs/ (also git-ignored).
##############################################################################
suppressPackageStartupMessages(library(phontrast))
stopifnot(packageVersion("phontrast") >= "2.1.0")

BASE    <- if (dir.exists("manuscript")) "manuscript" else "."
DATA    <- file.path(BASE, "data")
OUT     <- file.path(BASE, "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Parametric metrics (Pillai, Bhattacharyya) do not depend on the KDE JSD
## estimator, so they are identical across method = "mc"/"legacy".
pairwise_metrics <- function(tokens, feats, vcol, bw, engine) {
  vs <- sort(unique(tokens[[vcol]]))
  pairs <- t(combn(vs, 2))
  res <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    d <- tokens[tokens[[vcol]] %in% pairs[i, ], c(vcol, feats)]
    res[[i]] <- data.frame(
      v1 = pairs[i, 1], v2 = pairs[i, 2], n_tokens = nrow(d),
      jsd_mc     = jsd_kde_nd(d, feats, vcol, bw = bw, engine = engine, method = "mc"),
      jsd_legacy = jsd_kde_nd(d, feats, vcol, bw = bw, engine = engine, method = "legacy"),
      ovl_mc     = percent_overlap_kde(d, feats, vcol, bw = bw, engine = engine, method = "mc"),
      ovl_legacy = percent_overlap_kde(d, feats, vcol, bw = bw, engine = engine, method = "legacy"),
      pillai     = pillai_overlap(d, feats, vcol)$pillai,
      bhatt_dist = bhattacharyya_mvnorm(d, feats, vcol)$distance
    )
    if (i %% 25 == 0) message(i, "/", nrow(pairs), " pairs")
  }
  do.call(rbind, res)
}

corr_report <- function(tab, label) {
  sp <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
  message(sprintf(
    "[%s] Spearman with (1 - KDE overlap):  mc: JSD %.3f Pillai %.3f Bhatt %.3f  | legacy: JSD %.3f",
    label,
    sp(tab$jsd_mc, 1 - tab$ovl_mc), sp(tab$pillai, 1 - tab$ovl_mc),
    sp(tab$bhatt_dist, 1 - tab$ovl_mc), sp(tab$jsd_legacy, 1 - tab$ovl_legacy)))
}

## ---- 1. PB52 F1/F2 (10 monophthongs, pooled speakers) ---------------------
## Prefer a committed data/pb52.csv (columns vowel, f1, f2); fall back to the
## phonTools package copy.
pb52_path <- file.path(DATA, "pb52.csv")
pb <- if (file.exists(pb52_path)) {
  read.csv(pb52_path)
} else if (requireNamespace("phonTools", quietly = TRUE)) {
  utils::data("pb52", package = "phonTools"); get("pb52")
} else NULL
if (!is.null(pb)) {
  pb$vowel <- as.character(pb$vowel)
  pb52_tab <- pairwise_metrics(pb, c("f1", "f2"), "vowel", bw = "Hpi", engine = "ks")
  write.csv(pb52_tab, file.path(OUT, "pb52_pairwise_phontrast210.csv"), row.names = FALSE)
  corr_report(pb52_tab, "PB52 F1/F2")
} else {
  message("no data/pb52.csv and phonTools not installed; skipping PB52 block.")
}

## ---- 2. SBCSAE 13-MFCC (sampled tokens) -----------------------------------
mfcc_path <- file.path(DATA, "rerun_mfcc_sampled_tokens_20260624.csv")
if (file.exists(mfcc_path)) {
  tok <- read.csv(mfcc_path)
  mfcc_tab <- pairwise_metrics(tok, paste0("mfcc", 1:13), "vowel",
                               bw = "scott.diag", engine = "fast_diag")
  write.csv(mfcc_tab, file.path(OUT, "mfcc13_pairwise_phontrast210.csv"), row.names = FALSE)
  corr_report(mfcc_tab, "SBCSAE 13-MFCC")
} else {
  message("missing ", mfcc_path, "; skipping MFCC block (see README for source).")
}

message("Done. Tables written to ", OUT, " with phontrast ", packageVersion("phontrast"))
