##############################################################################
## Neutral overlap estimators (R port of the poster's reproduce_poster_figures.py)
##
## These overlap yardsticks do NOT use phontrast: kNN density, MVN Monte-Carlo,
## and two PCA-2D support measures. Porting them to R makes the manuscript's
## estimator-influence audit self-contained (no Python) and lets it run on any
## dataset (PB52 F1/F2 as well as SBCSAE 13-MFCC). Definitions follow the poster
## exactly; Monte-Carlo/subsample seeds differ, so values match the poster to
## Monte-Carlo error rather than bit-for-bit.
##############################################################################
suppressPackageStartupMessages({library(FNN); library(mvtnorm)})

.standardize_pair <- function(X0, X1) {
  X <- rbind(X0, X1); mu <- colMeans(X); sdv <- apply(X, 2, stats::sd)
  sdv[sdv == 0] <- 1; Xz <- sweep(sweep(X, 2, mu), 2, sdv, "/")
  list(X0 = Xz[seq_len(nrow(X0)), , drop = FALSE],
       X1 = Xz[nrow(X0) + seq_len(nrow(X1)), , drop = FALSE])
}

## Balance the two categories to equal n (subsample larger, seeded), standardize,
## and stack into a pooled matrix X with 0/1 labels y.
balanced_pair <- function(d0, d1, seed) {
  set.seed(seed); n <- min(nrow(d0), nrow(d1))
  if (nrow(d0) > n) d0 <- d0[sample.int(nrow(d0), n), , drop = FALSE]
  if (nrow(d1) > n) d1 <- d1[sample.int(nrow(d1), n), , drop = FALSE]
  z <- .standardize_pair(d0, d1)
  list(X = rbind(z$X0, z$X1), y = c(rep(0L, nrow(z$X0)), rep(1L, nrow(z$X1))))
}

knn_density_overlap <- function(X, y, k = 25) {
  X0 <- X[y == 0, , drop = FALSE]; X1 <- X[y == 1, , drop = FALSE]
  if (nrow(X0) <= k + 2 || nrow(X1) <= k + 2) return(NA_real_)
  d <- ncol(X)
  log_density <- function(train, eval) {
    dist <- FNN::knnx.dist(data = train, query = eval, k = k + 1)[, k + 1]
    dist <- pmax(dist, .Machine$double.eps)
    log(k) - log(nrow(train)) - d * log(dist)   # unit-ball const cancels on normalize
  }
  norm_ld <- function(lp) { m <- max(lp); w <- exp(lp - m); w / sum(w) }
  p <- norm_ld(log_density(X0, X)); q <- norm_ld(log_density(X1, X))
  sum(pmin(p, q))
}

mvn_monte_carlo_overlap <- function(X, y, n_mc = 10000, ridge = 1e-3, seed = 1) {
  X0 <- X[y == 0, , drop = FALSE]; X1 <- X[y == 1, , drop = FALSE]; d <- ncol(X)
  if (nrow(X0) <= d + 2 || nrow(X1) <= d + 2) return(NA_real_)
  mu0 <- colMeans(X0); mu1 <- colMeans(X1)
  c0 <- stats::cov(X0) + diag(ridge, d); c1 <- stats::cov(X1) + diag(ridge, d)
  set.seed(seed); n0 <- n_mc %/% 2
  s <- rbind(mvtnorm::rmvnorm(n0, mu0, c0), mvtnorm::rmvnorm(n_mc - n0, mu1, c1))
  f0 <- mvtnorm::dmvnorm(s, mu0, c0); f1 <- mvtnorm::dmvnorm(s, mu1, c1)
  den <- f0 + f1; ok <- is.finite(den) & den > 0 & is.finite(f0) & is.finite(f1)
  if (sum(ok) < 10) return(NA_real_)
  post1 <- f1[ok] / den[ok]
  mean(2 * pmin(post1, 1 - post1))
}

pca2_grid_overlap <- function(X, y, bins = 35) {
  if (length(y) < 10) return(NA_real_)
  Z <- stats::prcomp(X, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
  edges <- lapply(1:2, function(j) {
    lo <- min(Z[, j]); hi <- max(Z[, j]); w <- hi - lo
    seq(lo - 0.02 * w - 1e-6, hi + 0.02 * w + 1e-6, length.out = bins + 1)
  })
  binize <- function(zz) {
    i <- findInterval(zz[, 1], edges[[1]], rightmost.closed = TRUE)
    j <- findInterval(zz[, 2], edges[[2]], rightmost.closed = TRUE)
    ok <- i >= 1 & i <= bins & j >= 1 & j <= bins
    tabulate((j[ok] - 1) * bins + i[ok], nbins = bins * bins)
  }
  p <- binize(Z[y == 0, , drop = FALSE]); q <- binize(Z[y == 1, , drop = FALSE])
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  sum(pmin(p / sum(p), q / sum(q)))
}

## convex-polygon area (shoelace) and convex-convex intersection (Sutherland-Hodgman)
.poly_area <- function(P) { n <- nrow(P); if (n < 3) return(0)
  0.5 * abs(sum(P[, 1] * P[c(2:n, 1), 2] - P[c(2:n, 1), 1] * P[, 2])) }
.ccw <- function(P) if (sum(P[, 1] * P[c(2:nrow(P), 1), 2] - P[c(2:nrow(P), 1), 1] * P[, 2]) < 0) P[nrow(P):1, ] else P
.clip <- function(sub, clip) {
  clip <- .ccw(clip); out <- sub
  for (e in seq_len(nrow(clip))) {
    A <- clip[e, ]; B <- clip[if (e < nrow(clip)) e + 1 else 1, ]
    inside <- function(p) (B[1]-A[1])*(p[2]-A[2]) - (B[2]-A[2])*(p[1]-A[1]) >= -1e-12
    inp <- out; out <- matrix(numeric(0), 0, 2); m <- nrow(inp); if (m == 0) break
    for (i in seq_len(m)) {
      cur <- inp[i, ]; prv <- inp[if (i > 1) i - 1 else m, ]
      ci <- inside(cur); pi <- inside(prv)
      isect <- function() { d1 <- cur - prv
        t <- ((A[1]-prv[1])*(A[2]-B[2]) - (A[2]-prv[2])*(A[1]-B[1])) /
             ((d1[1])*(A[2]-B[2]) - (d1[2])*(A[1]-B[1])); prv + t * d1 }
      if (ci) { if (!pi) out <- rbind(out, isect()); out <- rbind(out, cur) }
      else if (pi) out <- rbind(out, isect())
    }
  }
  out
}
pca2_convex_hull_overlap <- function(X, y) {
  if (length(y) < 10) return(NA_real_)
  Z <- stats::prcomp(X, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
  z0 <- Z[y == 0, , drop = FALSE]; z1 <- Z[y == 1, , drop = FALSE]
  h0 <- tryCatch(z0[grDevices::chull(z0), , drop = FALSE], error = function(e) NULL)
  h1 <- tryCatch(z1[grDevices::chull(z1), , drop = FALSE], error = function(e) NULL)
  if (is.null(h0) || is.null(h1) || nrow(h0) < 3 || nrow(h1) < 3) return(NA_real_)
  denom <- min(.poly_area(h0), .poly_area(h1)); if (denom <= 0) return(NA_real_)
  inter <- .clip(.ccw(h0), h1)
  if (nrow(inter) < 3) return(0)
  .poly_area(inter) / denom
}

## Compute all neutral overlaps for one category pair (raw feature matrices).
neutral_overlaps <- function(d0, d1, seed) {
  bp <- balanced_pair(d0, d1, seed); X <- bp$X; y <- bp$y
  data.frame(
    overlap_knn_density_k10 = knn_density_overlap(X, y, 10),
    overlap_knn_density_k25 = knn_density_overlap(X, y, 25),
    overlap_knn_density_k50 = knn_density_overlap(X, y, 50),
    overlap_mvn_mc          = mvn_monte_carlo_overlap(X, y, seed = seed),
    overlap_pca2_grid       = pca2_grid_overlap(X, y),
    overlap_pca2_hull       = pca2_convex_hull_overlap(X, y))
}
