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
architecture. Target: additive as **2.1.0**; the shim removals in step 5 are the
trigger for **3.0.0**.

1. **Metric registry.** Register each metric with metadata: id, label,
   orientation (overlap vs separation), theoretical range, whether it supports
   bootstrap, and modelling assumptions (KDE vs multivariate-normal).
   `phontrast()` dispatches through the registry.
2. **Uniform per-metric contract.** One internal interface,
   `estimate(data, features, group, ...) -> scalar`, that every metric
   implements, so adding a metric means registering a single function.
3. **Generalised bootstrap.** Lift the currently JSD-centric bootstrap machinery
   to resample *any* registered metric with uniform confidence-interval columns.
4. **Orientation as a first-class concept.** Promote the existing
   `orientation` / `separation_value` / `separation_rank` idea into the type
   system so cross-metric comparison and ranking are coherent by construction.
5. **Consolidate wrappers.** Refactor `speaker_*`, `estimate_*`, `jsd_summary()`,
   and `hier_boot_jsd_model()` into thin shims over the unified core, deprecating
   gradually (removal triggers 3.0.0).
6. **Extensibility.** Document how users register custom metrics.
7. **CRAN.** Submit once the redesigned API is stable.

## Later

- Additional contrast metrics (e.g. energy distance, classifier-based
  separability).
- Richer visualisation of multi-metric comparisons.

Contributions and suggestions are welcome via
<https://github.com/berrygrant/phontrast/issues>.
