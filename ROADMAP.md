# phontrast roadmap

phontrast grew out of `phonJSD`, a package focused on Jensen-Shannon
divergence, and is being reoriented into a general toolkit for computing and
comparing multiple phonological category **contrast and separation metrics**.
This document records where it is headed.

## Now — 2.0.0 (shipped)

- Renamed `phonJSD` → `phontrast` and reframed the package around multi-metric
  contrast rather than JSD alone.
- Added `phontrast()`, a unified entry point that computes and compares any
  subset of the supported metrics (Jensen-Shannon divergence and distance,
  Pillai-Bartlett trace, Bhattacharyya distance and affinity, Mahalanobis
  distance, proportional overlap) in one call, globally or by group, with
  optional bootstrap intervals.
- Deprecated `compare_overlap_metrics()` in favour of `phontrast()`.

## P1 — Full architectural redesign (next)

Turn the working-but-ad-hoc multi-metric engine into a real, extensible
architecture. Target: additive across the 2.x series (2.1.0 shipped the
corrected Monte-Carlo estimator); the shim removals in step 6 are the trigger
for **3.0.0**.

1. **Metric registry.** Register each metric with metadata: id, label,
   orientation (overlap vs separation), theoretical range, whether it supports
   bootstrap, and modelling assumptions (KDE vs multivariate-normal).
   `phontrast()` dispatches through the registry.
2. **Uniform per-metric contract.** One internal interface,
   `estimate(data, features, group, ...) -> scalar`, that every metric
   implements, so adding a metric means registering a single function.
3. **Pluggable density backends for the distributional metrics.** Decouple the
   density *estimator* from the *metric*: add a `density` argument
   (`"kde"`, the current default; `"mvnorm"`; later `"gmm"`) to the
   Jensen-Shannon and proportional-overlap estimators, so the distribution
   behind each metric is a controlled choice rather than hard-wired to KDE.
   - *Motivation.* The phontrast metric-comparison analyses find that **each
     metric aligns best with the overlap estimator whose structural assumptions
     it shares** — Jensen-Shannon divergence with KDE overlap; Bhattacharyya and
     Pillai with a multivariate-normal reference. Welding JSD to KDE confounds
     the metric with its estimator. Making the estimator pluggable turns it into
     a controlled variable: users can match the estimator to their distributional
     assumptions, and the package can study the metric × estimator interaction
     directly (e.g. Jensen-Shannon divergence computed under a Gaussian fit
     against a Gaussian overlap reference).
   - *Design.* One internal density-model interface — `fit(X)` and
     `logdens(model, points)` — with KDE as the current implementation and a
     multivariate-normal backend (`mvtnorm::dmvnorm`). Jensen-Shannon divergence
     between two Gaussians has no closed form (the mixture is a Gaussian
     mixture), so the Gaussian backend estimates it by Monte-Carlo, reusing the
     existing plug-in averaging without the KDE-specific leave-one-out
     self-kernel term. The manuscript's MVN Monte-Carlo overlap then becomes a
     native backend rather than external analysis code. KDE-only arguments
     (`bw`, `engine`, `eval_on`, `loo`) are ignored with a warning under a
     parametric backend. Estimated effort ~1–1.5 days for the Gaussian backend;
     a Gaussian-mixture backend is a larger, separate step (adds an EM /
     `mclust` dependency and component selection).
4. **Generalised bootstrap.** Lift the currently JSD-centric bootstrap machinery
   to resample *any* registered metric with uniform confidence-interval columns.
5. **Orientation as a first-class concept.** Promote the existing
   `orientation` / `separation_value` / `separation_rank` idea into the type
   system so cross-metric comparison and ranking are coherent by construction.
6. **Consolidate wrappers.** Refactor `speaker_*`, `estimate_*`, `jsd_summary()`,
   and `hier_boot_jsd_model()` into thin shims over the unified core, deprecating
   gradually (removal triggers 3.0.0).
7. **Extensibility.** Document how users register custom metrics.
8. **CRAN.** Submit once the redesigned API is stable.

## Later

- Additional contrast metrics (e.g. energy distance, classifier-based
  separability).
- Richer visualisation of multi-metric comparisons.

Contributions and suggestions are welcome via
<https://github.com/berrygrant/phontrast/issues>.
