# Issue: SummarizedExperiment input path was fully broken

## Summary
`signature_projection_contributors()` and `spc_heatmap_all()`/`spc_heatmap_sig()` (via `.spc_prepare_heatmap_data()`) all documented support for `SummarizedExperiment` input, alongside `ExpressionSet`. In practice, calling any of them with a real `SummarizedExperiment` errored out.

## Root cause
Two independent bugs stacked in the same code path:
1. `signature_projection_contributors()`'s `sig_score`/`col_ha` validation `stopifnot`s called `Biobase::sampleNames(eset)` *before* the `SummarizedExperiment -> ExpressionSet` conversion ran. `sampleNames()` has no method for `SummarizedExperiment`, so this raised `unable to find an inherited method for function 'sampleNames'...` immediately.
2. Once that was fixed by reordering, a second bug surfaced: the conversion itself passed `colData(eset)`/`rowData(eset)` (an S4 `DFrame`) straight into `Biobase::AnnotatedDataFrame()`, which has no method for `DFrame` either — `unable to find an inherited method for function 'AnnotatedDataFrame'...`. This bug was **latent** (never reached until bug #1 was fixed) and existed in duplicate, unsynchronized form in `spc_heatmaps.R`'s own copy of the same coercion logic.

## Fix implemented
- Reordered `signature_projection_contributors()` so the `SummarizedExperiment -> ExpressionSet` coercion runs before any `sampleNames()`-dependent checks.
- Wrapped `colData()`/`rowData()` in `as.data.frame()` before handing them to `AnnotatedDataFrame()`.
- Deduplicated the coercion into a single shared internal helper, `.as_expressionset()`, in a new `R/utils.R`, used by both `signature_projection_contributors.R` and `spc_heatmaps.R`. This fixed the latent duplicate copy of bug #2 in `spc_heatmaps.R` as a side effect, and removed now-unused `SummarizedExperiment` `importFrom`s from `NAMESPACE`.

## Files changed
- `R/signature_projection_contributors.R`
- `R/spc_heatmaps.R`
- `R/utils.R` (new)
- `NAMESPACE`

## Verification
- Reproduced both bugs directly against a real `SummarizedExperiment` before fixing, confirmed both errors, then confirmed a full `signature_projection_contributors(se, ...)` call succeeds end-to-end after the fix.
- Full `testthat` suite (50/50) passing throughout.

## Notes
Commits: `9702c6b` (bug #1 + #2, first file), `89c020e` (dedup + fix in second file).
