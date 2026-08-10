# Issue: API-clarity, safety, and cleanup pass (post code review)

## Summary
A batch of smaller, independent improvements identified during a full code review of the package, addressed together after the correctness-bug fixes and CI work documented in the sibling `docs/issues/2026-07-14-*` files.

## Items addressed

### 1. Misleading single-signature API
`signature` was documented everywhere as "a list of signatures (at least 1)", but `omics_signature_score()` (and transitively `signature_projection_contributors()`) hard-required `length(signature) == 1` via a bare `stopifnot`. Docs now state the real one-signature-only contract, and the opaque assertion was replaced with `stop("'signature' must be a list containing exactly one signature; multi-signature scoring is not currently supported")`.

### 2. Silent overwriting of reserved eset columns
`signature_projection_contributors()` writes fixed-name columns (`sig_score`, `score_cor`, `pval_cor`, `insig`, `leading_edge`) into the working `eset`'s pData/fData, and `extract_eigengenes()` does the same for a `key` fData column — with no check for a pre-existing column of the same name. Added `.warn_if_overwriting()` to `R/utils.R`, called before each such assignment, so callers whose own data happens to collide with a reserved name get a `warning()` instead of silent data loss.

### 3. Dead code in `.kstest()`
Removed `.ks_anno_lines()` (defined, never called except from a commented-out line, and dragging in an unused `ComplexHeatmap::anno_lines` import) and the `absolute` parameter (documented in-source as "not really implemented, should be removed", with no caller ever passing it). `.kstest()` now always computes the signed enrichment score.

### 4. Unused correlation confidence intervals
`psych::corr.test(...)` was called with its default `ci = TRUE`, computing per-gene confidence intervals that are never read (only `$r`/`$p` are used downstream). Passed `ci = FALSE` instead. Verified byte-identical `$r`/`$p` output with and without `ci`, so this is a pure efficiency win with no behavior change — meaningful at the gene counts (tens of thousands) this function is meant to run on.

## Files changed
- `R/signature_projection_contributors.R`
- `R/omics_signature_score.R`
- `R/extract_eigengenes.R`
- `R/ks_enrichment.R`
- `R/utils.R`
- `NAMESPACE`, `man/*.Rd` (regenerated)

## Verification
- Full `testthat` suite (50/50) passing at every intermediate commit.
- Local `R CMD check`: 0 errors, 0 warnings (one incidental doc regression from regenerating `man/` with a different local roxygen2 version — a missing `fontsize` `@param` on `spc_heatmap_all()`/`spc_heatmap_sig()` — was caught by this check and fixed in the same pass, see commit `8ddc82d`).
- Manually verified the overwrite-warning fires correctly and the `ci = FALSE` correlation values match `ci = TRUE`.

## Notes
Commits (in order): `1e54692` (dead code), `89c020e` (dedup, documented separately), `f7d5e7a` (overwrite warnings), `7d9218b` (API clarity + `ci = FALSE`), `8ddc82d` (fontsize doc regression fix).
