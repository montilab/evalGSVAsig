# Issue: extract_eigengenes() false-rejects non-overlapping signatures

## Summary
`extract_eigengenes(..., method = "multiple")` refuses to run when it detects "overlapping signatures", but the overlap check could trigger even when the signatures had no actual overlap in the data being analyzed.

## Root cause
The overlap check ran on the *raw* `signatures` argument as passed in by the caller, before restricting to features present in `eset`:

```r
if (method == "multiple" && length(signatures) > 1) {
  tmp <- utils::combn(signatures, 2, function(mod) length(intersect(mod[[1]], mod[[2]])), simplify = FALSE)
  if (any(tmp > 0)) stop("method=='multiple' can't be used with overlapping signatures")
}
```

If two signatures happened to share a gene symbol that wasn't even present in `eset` (e.g. a stale/renamed identifier), the check would still fire and reject the call — even though the actual downstream computation, restricted to `eset`'s features, would have been perfectly well-defined and non-overlapping.

## Fix implemented
Moved the overlap check to run after `signatures1` (signatures filtered by `min_size` overlap with `eset`) is computed, and intersect each signature with `fData(eset)$key` before comparing pairwise overlaps:

```r
if (method == "multiple" && length(signatures1) > 1) {
  signatures1_eff <- lapply(signatures1, function(Z) intersect(Z, Biobase::fData(eset)$key))
  tmp <- utils::combn(signatures1_eff, 2, function(mod) length(intersect(mod[[1]], mod[[2]])), simplify = FALSE)
  if (any(tmp > 0)) stop("method=='multiple' can't be used with overlapping signatures")
}
```

## Files changed
- `R/extract_eigengenes.R`

## Verification
- Reproduced the false-rejection with two signatures sharing only an out-of-data gene; confirmed no error after the fix.
- Confirmed a *genuine* in-data overlap still correctly raises the error (no regression on the intended guard).
- Full `testthat` suite (50/50) passing.

## Notes
Commit: `5c0b013`.
